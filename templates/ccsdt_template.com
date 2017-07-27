***, 
!job computes, at a CCSD(T)-f12 level of theory, the
!interaction energy between monomers. eint is calculated as the
!difference between the dimer and monomers' energy.
!
!This template utilizes density fitting.
!
!Need to change the following before using template:
!   *Change dummy atom labels on monomers a and b for both dhf and sapt calc
!   *(Monomer names)
!   *Monomer/dimer charges
!May need to also change the following:
!   *Memory limits
memory,300,m

include, AVTZ.mbas

nosym;noorient;angstrom
GEOMETRY_BLOCK_GOES_HERE

        !dimer
        dummy,Be
        charge=FILL_Q_DIMER
        {rhf}
        {CCSD(T)-F12,SCALE_TRIP=1,df_basis=jkfit,df_basis_exch=jkfit,ri_basis=ri}
        edm=energy

        !monomer A, FILL_MONA
        dummy, Be, FILL_I_MONB
        charge=FILL_Q_MONA
        {rhf}
        {CCSD(T)-F12,SCALE_TRIP=1,df_basis=jkfit,df_basis_exch=jkfit,ri_basis=ri}
        ema=energy

        !monomer B, FILL_MONB
        dummy, Be, FILL_I_MONA
        charge=FILL_Q_MONB
        {rhf}
        {CCSD(T)-F12,SCALE_TRIP=1,df_basis=jkfit,df_basis_exch=jkfit,ri_basis=ri}
        emb=energy

 eint_ccsdt=(edm-ema-emb)*1000 mH
