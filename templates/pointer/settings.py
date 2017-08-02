##########################################################################
######################### General Settings ###############################
##########################################################################
# Monomer Names (should match ordering in .sapt file)
mon1                       =   'chloromethane'
mon2                       =   'chloromethane'
                            
# Atomtype settings         
constrained_atomtypes      =    []
                            
# Output options            
file_prefix                =   'fit_exp_'
file_suffix                =   '_unconstrained'

##########################################################################
##########################################################################



##########################################################################
##################### Component-Specific Settings ########################
##########################################################################
# Electrostatic Settings
multipoles_suffix          =   '_ISA-GRID_L2.mom'

# Exchange Settings
fit_bii                    =    True

# Induction Parameters
thole_param                =    0.33
thole_damping_type         =   'thole_tinker'
induction_damping_type     =   'Thole'

# Dispersion parameters
fit_anisotropic_dispersion =    True
fit_isotropic_dispersion   =    False


# If set to true, fits a final A parameter to errors in the total
# energy, in an effort to reduce systematic errors in the total energy
fit_residuals              =    False

##########################################################################
##########################################################################



##########################################################################
##################### Functional Form Settings ###########################
##########################################################################
# radial potentials; see Stone's book for more details.
# Options are 'slater', 'stone', 'born-mayer', or 'lennard-jones'
functional_form           =    'slater'

#   aij: 'saptff', 'waldman-hagler5', 'geometric' 
#   bij: 'saptff', 'waldman-hagler5', 'geometric_mean', 'arithmetic_mean'
#   cij: 'geometric'
aij_combination_rule      =    'geometric'
bij_combination_rule      =    'geometric_mean'
cij_combination_rule      =    'geometric'

##########################################################################
##########################################################################



