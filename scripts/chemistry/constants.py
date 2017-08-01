from math import pi, sqrt

# Constants from Konrad Hinsen's PhysicalQuantities module (1986 CODATA):
clight = 299792458.              # speed of light, m/s
hplanck = 6.6260755e-34          # Planck's constant, J s
hbar = hplanck / (2 * pi)       # Planck's constant / 2pi, J s
_e = 1.60217733e-19              # elementary charge
_me = 9.1093897e-31              # electron mass (kg)
_mp = 1.6726231e-27              # proton mass (kg)
avagadro = 6.0221367e23          # Avogadro's number (particles/mol)
kboltz = 1.380658e-23            # Boltzmann constant, J/K
kb=kboltz
amu = 1.6605402e-27              # atomic mass unit, kg

# Distance unit conversions
bohr2ang = 0.529177249  # Conversion of length from bohr to angstrom
ang2bohr = 1/bohr2ang

# Energy unit conversions
hartree2kcal = 627.5095 # Hartree to kcal/mol conversion
kcal2hartree = 1/hartree2kcal

ev2kcal = 23.061        # Conversion of energy in eV to energy in kcal/mol
kcal2ev = 1/ev2kcal

hartree2joule = 4.3597482e-18   # Hatree to Joule conversion factor
joule2hartree = 1/hartree2joule

hartree2ev = 27.21138505  # Hartree to eV conversion Factor
ev2hartree = 1/hartree2ev

hartree2cm1 = 219474.63     #Hartree to cm^-1 conversion factor
cm12hartree = 1/hartree2cm1

# Mass unit conversions
amu2me = 1822.882       # Conversion from mass in amu to mass in au (m_e)
me2amu = 1/amu2me       # Conversion from mass in au (m_e) to mass in amu


# Time unit conversions
tau2ps = 41341.447      # Conversion from time in au to time in ps
ps2tau = 1/tau2ps       # inverse

# Derived quantities
rgas = kboltz*hartree2kcal*1000.0 # gas constant R = 1.98722 cal/mole/K

del pi, sqrt
