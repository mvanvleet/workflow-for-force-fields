***, co2 dimer
memory,300,m

include, AVTZ.mbas

nosym;noorient;angstrom
GEOMETRY_BLOCK_GOES_HERE

        !dimer
        dummy,Be
        charge=0
        {rhf}
        {CCSD(T)-F12,SCALE_TRIP=1,df_basis=jkfit,df_basis_exch=jkfit,ri_basis=ri}
        edm=energy

        !monomer A
        dummy,1,2,3,Be
        charge=0
        {rhf}
        {CCSD(T)-F12,SCALE_TRIP=1,df_basis=jkfit,df_basis_exch=jkfit,ri_basis=ri}
        ema=energy

        !monomer B
        dummy,4,5,6,Be
        charge=0
        {rhf}
        {CCSD(T)-F12,SCALE_TRIP=1,df_basis=jkfit,df_basis_exch=jkfit,ri_basis=ri}
        emb=energy

 eint_ccsdt=(edm-ema-emb)*1000 mH
