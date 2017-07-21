***, 
!job computes, at a DFT-SAPT level of theory (pbe0/AVTZ+midbond) the
!interaction energy between monomers. eint is calculated as the
!difference between the dimer and monomers' energy.
!
!This template does NOT utilize density fitting, and should probably only be
!necessary for systems that involve lithium.
!
!Need to change the following before using template:
!   *Change dummy atom labels on monomers a and b for both dhf and sapt calc
!   *Update HOMO and IP for each monomer
!   *(Monomer names)
!May need to also change the following:
!   *Memory limits
!   *Monomer/dimer charges

memory,200,m

basis={
set,orbital
default=avtz
! ---------------- even tempered mid bond basis ------------------------
S,Be,EVEN,5,2.5,0.5
P,Be,EVEN,3,2.5,0.5
D,Be,EVEN,1,2.5,0.3
F,Be,EVEN,1,2.5,0.3
}

nosym;noorient;angstrom;
GEOMETRY_BLOCK_GOES_HERE

!sapt files
cd=2000.2
ca=2101.2
cb=2102.2

!=========delta(HF) contribution for higher order interaction terms====
    !dimer
    dummy,Be
    charge=FILL_Q_DIMER
    {hf; save,$cd}
    edm=energy

    !monomer A, FILL_MONA
    dummy, Be, FILL_I_MONB
	charge=FILL_Q_MONA
    {hf; save,$ca}
    ema=energy
    sapt;monomerA

	!monomer B, FILL_MONB
    dummy, Be, FILL_I_MONA
	charge=FILL_Q_MONB
    {hf; save,$cb}
    emb=energy
    sapt;monomerB

    !interaction contributions
    {sapt,SAPT_LEVEL=2,THRSW=1e-12,THRPROD=0d0;intermol,ca=$ca,cb=$cb,icpks=1}
    !calculate high-order terms by subtracting 1st+2nd order energies
    eint_hf=(edm-ema-emb)*1000 mH
    delta_hf=eint_hf-e1pol-e1ex-e2ind-e2exind
    dhf(i)=delta_hf

    data,truncate,$cd,$ca,$cb


!=========DFT-SAPT at second order intermol. perturbation theory====
!sapt files
ca=2103.2
cb=2104.2

!asymptotic correction data calcluated in subdirectory monomer_calcs/
	homo_FILL_MONA =  FILL_HOMO_MONA !HOMO/PBE0 AVTZ
	ip_FILL_MONA =  FILL_IP_MONA              !from experiment
	shift_FILL_MONA =ip_FILL_MONA + homo_FILL_MONA	!shift for bulk xc potential

	homo_FILL_MONB =  FILL_HOMO_MONB !HOMO/PBE0 AVTZ
	ip_FILL_MONB =  FILL_IP_MONB              !from experiment
	shift_FILL_MONB =ip_FILL_MONB + homo_FILL_MONB	!shift for bulk xc potential


    !monomer A, FILL_MONA
    dummy, Be, FILL_I_MONB
	charge=FILL_Q_MONA
	{ks,pbe0; start,atdens; asymp,shift_FILL_MONA; save,$ca}
	sapt;monomerA

	!monomer B, FILL_MONB
    dummy, Be, FILL_I_MONA
	charge=FILL_Q_MONB
	{ks,pbe0; start,atdens; asymp,shift_FILL_MONB; save,$cb}
	sapt;monomerB

	!interaction contributions
	{sapt,SAPT_LEVEL=3,THRSW=1e-12;intermol,ca=$ca,cb=$cb,icpks=0}

	!add high-order approximation to obtain the total interaction energy
	eint_dftsapt=e12tot+dhf(i)

	elst(i)=E1pol;  exch(i)=E1ex
	ind(i)=E2ind;   exind(i)=E2exind
	disp(i)=E2disp; exdisp(i)=E2exdisp
	etot(i)=eint_dftsapt

	data,truncate,$ca
