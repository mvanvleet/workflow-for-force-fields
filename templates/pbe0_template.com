***, 
!job computes, at a DFT-SAPT level of theory (pbe0/AVTZ+midbond) the
!interaction energy between monomers. eint is calculated as the
!difference between the dimer and monomers' energy.
!
!This template utilizes density fitting.
!
!Need to change the following before using template:
!   *Change dummy atom labels on monomers a and b for both dhf and sapt calc
!   *Update HOMO and IP for each monomer
!   *(Monomer names)
!   *Monomer/dimer charges
!May need to also change the following:
!   *Memory limits


memory,200,m
gdirect

include, AVTZ.mbas

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
    {df-hf,basis=jkfit,locorb=0; save,$cd}
    edm=energy

    !monomer A, FILL_MONA
    dummy, Be, FILL_I_MONB
	charge=FILL_Q_MONA
    {df-hf,basis=jkfit,locorb=0; save,$ca}
    ema=energy
    sapt;monomerA

	!monomer B, FILL_MONB
    dummy, Be, FILL_I_MONA
	charge=FILL_Q_MONB
    {df-hf,basis=jkfit,locorb=0; save,$cb}
    emb=energy
    sapt;monomerB

    !interaction contributions
    {sapt,SAPT_LEVEL=2,THRSW=1e-12,THRPROD=0d0;intermol,ca=$ca,cb=$cb,icpks=1,fitlevel=3
    dfit,basis_coul=jkfit,basis_exch=jkfit,cfit_scf=3}

    !calculate high-order terms by subtracting 1st+2nd order energies
    eint_hf=(edm-ema-emb)*1000 mH
    delta_hf=eint_hf-e1pol-e1ex-e2ind-e2exind
    dhf(i)=delta_hf

    data,truncate,$cd,$ca,$cb



!=========DFT-SAPT at second order intermol. perturbation theory====
!sapt files
ca=2103.2
cb=2104.2

!asymptotic correction data calcluated in /home/mvanvleet/data_reference/sapt_monomer_data/ips/
	!shifts for asymptotic correction to xc potential (in Ha)
	eps_homo_pbe0_FILL_MONA =  FILL_HOMO_MONA !HOMO/PBE0 AVTZ
	ip_FILL_MONA =  FILL_IP_MONA              !from experiment
	shift_FILL_MONA =ip_FILL_MONA + eps_homo_pbe0_FILL_MONA	!shift for bulk xc potential

	eps_homo_pbe0_FILL_MONB =  FILL_HOMO_MONB !HOMO/PBE0 AVTZ
	ip_FILL_MONB =  FILL_IP_MONB              !from experiment
	shift_FILL_MONB =ip_FILL_MONB + eps_homo_pbe0_FILL_MONB	!shift for bulk xc potential

	dfit,basis_coul=jkfit,basis_exch=jkfit

    !monomer A, FILL_MONA
    dummy, Be, FILL_I_MONB
	charge=FILL_Q_MONA
	{df-ks,pbex,pw91c,lhf; dftfac,0.75,1.0,0.25; start,atdens; asymp,shift_FILL_MONA; save,$ca}
	sapt;monomerA

	!monomer B, FILL_MONB
    dummy, Be, FILL_I_MONA
	charge=FILL_Q_MONB
	{df-ks,pbex,pw91c,lhf; dftfac,0.75,1.0,0.25; start,atdens; asymp,shift_FILL_MONB; save,$cb}
	sapt;monomerB

	!interaction contributions
	{sapt,SAPT_LEVEL=3,THRSW=1e-12;intermol,ca=$ca,cb=$cb,icpks=0,fitlevel=3,nlexfac=0.0
	dfit,basis_coul=jkfit,basis_exch=jkfit,cfit_scf=3}

	!add high-order approximation to obtain the total interaction energy
	eint_dftsapt=e12tot+dhf(i)

	elst(i)=E1pol;  exch(i)=E1ex
	ind(i)=E2ind;   exind(i)=E2exind
	disp(i)=E2disp; exdisp(i)=E2exdisp
	etot(i)=eint_dftsapt

	data,truncate,$ca
