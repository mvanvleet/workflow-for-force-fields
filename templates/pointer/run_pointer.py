#!/usr/bin/env python
"""

Last Updated:
"""

# Standard modules
import numpy as np
import sys
import os
import subprocess
import time

# mvanvleet specific modules
from chemistry import io
from chemistry.elementdata import Exponent
from pointer.fit_ff_parameters import FitFFParameters

# Settings files
from input import settings
from input import defaults

###########################################################################
####################### Global Variables ##################################
# List value by which to scale exponents; should be 1.00 for Slater-ISA, and
# 0.84 for Born-Mayer-sISA force fields
functional_form = settings.functional_form.lower()
exp_scale = ( 0.84 if functional_form == 'born-mayer-sisa' else 1.00)

# List which atomtypes should be treated as hard constraints
#atomtype_constraints = ['H','O']
atomtype_constraints = settings.constrained_atomtypes

input_dir = 'input/'

sapt_dir = input_dir
atomtypes_dir = input_dir
exponents_dir = input_dir
drudes_dir = input_dir
dispersion_dir = input_dir
multipoles_dir = input_dir
anisotropy_dir = input_dir
constraints_dir = input_dir
induction_exponents_dir = input_dir

multipoles_suffix = settings.multipoles_suffix

use_slater_ff = (settings.functional_form.lower() == 'slater') 

fit_dispersion = False
fit_isotropic_dispersion = False
if settings.fit_dispersion.lower() == 'anisotropic':
    fit_dispersion = True
elif settings.fit_dispersion.lower() == 'all':
    fit_dispersion = True
    fit_isotropic_dispersion = True
elif settings.fit_dispersion.lower() != 'none':
    sys.exit('Dispersion fit type not recognized')

fit_residuals = settings.fit_residuals
separate_induction_exponents = defaults.separate_induction_exponents

n_a_params = 6 # Aexch, Aelst, Aind, Adhf, Adisp, Aresid by default

# Give parameters for lone pairs or other offsite atoms
lp_constraints = [[0],[0],[0],[0],[0],[0],[0,0,0,0],[1.0],[1.0]]

# Weighting function is defined by the Fermi-Dirac distribution and an
# effective fermi level and temperature. Here we keep the fermi level fixed at
# 0mH and modify the temperature to be half the absolute energy of the well-depth.
scale_weighting_temperature = defaults.scale_weighting_temperature

thole_param = settings.thole_param


mon1 = settings.mon1
mon2 = settings.mon2
###########################################################################
###########################################################################


###########################################################################
###########################################################################

###########################################################################
def create_param_file(mon1,mon2,input_sapt, constraints, cn_coeffs, 
                        drude_charges, outputfile, constrain_Aparams=[],
                        springcon=0.1):
    '''Create the parameter file that gets directly read into POInter
    '''

    # Read in data from .sapt file. Geometry information has lines which each
    # contain exactly 4 blocks of information, allowing us to seperate out that
    # information based on that condition.
    with open(input_sapt,'r') as f:
        # First line contains information about the number of lines in monomer a
        num_atoms_mona=int(f.readline())
        count = 0
        geometry_block = []
        while count < 2:
            line = f.readline()
            if len(line.split()) != 4: 
                count +=1
            else:
                geometry_block.append(line.split())

    # Use SAPT file to get an effective kt for fitting the data (1mH or less)
    with open(input_sapt,'r') as f:
        data = [line.split() for line in f.readlines()]

    etot_min = 0.0
    for i,line in enumerate(data):
        if len(line) > 0 and line[0] == 'E1tot+E2tot':
            etot = float(line[1]) + float(data[i+5][1])
            if etot < etot_min:
                etot_min = etot
    eff_kt = -etot_min/1000*scale_weighting_temperature



    # Use information from .sapt file to extract atom lists
    atom_list_mona=[line[0] for line in geometry_block[0:num_atoms_mona]]
    atom_list_monb=[line[0] for line in geometry_block[num_atoms_mona:]]

    # Write output param file
    with open(outputfile,'w') as f:

        # Settings
        f.write('SETTINGS\n')
        if fit_dispersion:
            f.write('fit_dispersion    True\n')
        if fit_isotropic_dispersion:
            f.write('fit_isotropic_dispersion    True\n')
        if fit_residuals:
            f.write('fit_residuals               True \n')
        if separate_induction_exponents:
            f.write('separate_induction_exponents  True\n')
        f.write('thole_damping_type {}\n'.format(settings.thole_damping_type))
        f.write('induction_damping_type {}\n'.format(defaults.induction_damping_type))
        f.write('thole_param            {:8.6f}\n'.format(thole_param))
        f.write('\n')

        # A coefficients
        f.write('LIST OF PARAMETER CONSTRAINTS\n')
        unique_atoms = set(atom_list_mona + atom_list_monb)

        constrained_Aparams = set(constrain_Aparams)
        unique_atoms = set(atom_list_mona + atom_list_monb)
        constrained_atoms = unique_atoms - (unique_atoms - constrained_Aparams)
        constrained_atoms = list(constrained_atoms)


        f.write('EXCHANGE\n')
        for atom in constrained_atoms:
            template='{:5s}'+'{:12.6f}'*len(constraints[atom][0]) + '\n'
            f.write(template.format(atom, *constraints[atom][0]))

        f.write('ELECTROSTATICS\n')
        for atom in constrained_atoms:
            template='{:5s}'+'{:12.6f}'*len(constraints[atom][1]) + '\n'
            f.write(template.format(atom, *constraints[atom][1]))

        f.write('INDUCTION\n')
        for atom in constrained_atoms:
            template='{:5s}'+'{:12.6f}'*len(constraints[atom][2]) + '\n'
            f.write(template.format(atom, *constraints[atom][2]))

        f.write('DHF\n')
        for atom in constrained_atoms:
            template='{:5s}'+'{:12.6f}'*len(constraints[atom][3]) + '\n'
            f.write(template.format(atom, *constraints[atom][3]))

        if fit_dispersion:
            f.write('DISPERSION\n')
            for atom in constrained_atoms:
                template='{:5s}'+'{:12.6f}'*len(constraints[atom][4]) + '\n'
                f.write(template.format(atom, *constraints[atom][4]))

        if fit_residuals:
            f.write('RESIDUALS\n')
            for atom in constrained_atoms:
                template='{:5s}'+'{:12.6f}'*len(constraints[atom][5]) + '\n'
                f.write(template.format(atom, *constraints[atom][5]))

        # Anisotropic Atoms
        anisotropy_terms1, anisotropy_axes1 = get_anisotropy_terms(mon1)
        anisotropy_terms2, anisotropy_axes2 = get_anisotropy_terms(mon2)
        anisotropy_terms1.update(anisotropy_terms2)
        f.write('\nLIST OF ANISOTROPIC ATOMTYPES\n')
        for k, v in anisotropy_terms1.iteritems():
            template = '{:5s}'*(len(v)+1) + '\n'
            f.write(template.format(k,*v))
        f.write('\n')

        f.write('DEFINE COORDINATE AXES FOR EACH ANISOTROPIC ATOM IN EACH MONOMER\n')
        f.write('ATOM#   AXIS (z or x)   Atomic Indices defining vector (either 2 or 3 integers)\n')
        f.write('monomer 1\n')
        for line in anisotropy_axes1:
            template = '{:3s}'*(len(line)) + '\n'
            f.write(template.format(*line))
        f.write('monomer 2\n')
        for line in anisotropy_axes2:
            template = '{:3s}'*(len(line)) + '\n'
            f.write(template.format(*line))
        f.write('\n')

        # Exponents (B coefficients)
        f.write('EXPONENTS: One for each atom\n')
        for atom in unique_atoms:
            template='{:<5s}'+'{:10.6f}\n'
            f.write(template.format(atom, *constraints[atom][7]))
        f.write('\n')

        # Indcution Exponents
        if separate_induction_exponents:
            f.write('INDUCTION EXPONENTS: One for each atom\n')
            #f.write('monomer 1\n')
            for atom in unique_atoms:
                template='{:<5s}'+'{:10.6f}\n'
                print atom
                print constraints[atom]
                print constraints[atom][-1]
                f.write(template.format(atom, *constraints[atom][-1]))
            f.write('\n')

        # Dispersion (Cn coefficients)
        f.write(' Cn COEFFICIENTS (C6, C8, C10, C12)\n')
        for atom in unique_atoms:
            cn = constraints[atom][6]
            template='{:<5s}'+'{:16.6f}'*len(cn) + '\n'
            f.write(template.format(atom,*cn))
        f.write('\n')

        # Multipole file names
        f.write('MULTIPOLE FILES\n')
        f.write('monomer 1\n')
        multipole_file1 = mon1 + multipoles_suffix + '\n'
        f.write(multipole_file1)
        f.write('monomer 2\n')
        multipole_file2 = mon2 + multipoles_suffix + '\n'
        f.write(multipole_file2)
        f.write('\n')


        # Drude charges and spring coefficients
        f.write('DRUDE PARAMETERS (q   kx    ky    kz) \n')
        template='{:7s} {:10.6f} {:10.6f}{:10.6f}{:10.6f}\n'
        f.write('monomer 1\n')
        for i,line in enumerate(atom_list_mona):
            f.write(template.format(*drude_charges[0][i]))
        f.write('monomer 2\n')
        for i,line in enumerate(atom_list_monb):
            f.write(template.format(*drude_charges[1][i]))
        f.write('\n')

        # Fit Weights
        f.write('FITTING PARAMETERS\n')
        f.write('eff_mu(Ha)  0.000\n')
        f.write('eff_kT(Ha)  {:6.4f}\n'.format(eff_kt))


    return
###########################################################################


###########################################################################
def get_atomtypes(sapt_file):
    '''Read in atomtypes from .sapt file'''
        
    atomtypes = []

    with open(sapt_file,'r') as f:
        # First line contains information about the number of lines in monomer a
        num_atoms_mona=int(f.readline())
        for line in range(num_atoms_mona):
            atom = f.readline().split()[0]
            atomtypes.append(atom)
        num_atoms_monb=int(f.readline())
        for line in range(num_atoms_monb):
            atom = f.readline().split()[0]
            atomtypes.append(atom)

    return atomtypes
###########################################################################
def get_exponents(molecule,exponent_type='density-cutoff'):
    ''' Assuming exponent files of the following form (for each element in the
    file)

Atomtype:  C0
Fitting tail-corrected region using points over the range  1.84  to  6.34
B (two-point fit): 2.3864036044
B (average over range): 2.38639372601
B (r_abs_cutoff minimize RMSE): 2.38125128048
Fitting effective exponents using points over the range  1.14  to  8.72
B (density_abs_cutoff minimize RMSE): 2.38221095365

    extracts exponents for each atom in the molecule.

    Allowed exponent_type values are 'two-point', 'tail-average', and
    'density-cutoff'.
    '''

    # 7 lines per atom
    natomlines = 7

    if exponent_type == 'two-point':
        n_preskip_lines = 2
    elif exponent_type == 'tail-average':
        n_preskip_lines = 3
    elif exponent_type == 'density-cutoff':
        n_preskip_lines = 6
    else:
        sys.exit('Unrecognized exponent type.')

    with open(exponents_dir + molecule + '.exp','r') as f:
        data = [line.split() for line in f.readlines()]
        exponents = [float(i[-1]) for i in data[n_preskip_lines::natomlines] ]

    return exponents
###########################################################################


###########################################################################
def get_induction_exponents(molecule):
    '''
    '''

    n_preskip_lines = 1 #First line is header for each atom

    with open(induction_exponents_dir + molecule + '.indexp','r') as f:
        data = [line.split() for line in f.readlines()]
        for line in data:
            print line
        induction_exponents = [float(i[1]) for i in data[n_preskip_lines:] if i] 

    return induction_exponents
###########################################################################


###########################################################################
def get_anisotropy_terms(molecule):
    with open(anisotropy_dir + molecule + '.axes','r') as f:
        data = [line.split() for line in f.readlines()]

    anisotropy_terms = {}
    anisotropy_axes = []
    count = 1
    for line in data[count:]: # first line is title
        count += 1
        if not line:
            break
        # each anisotropy line is atomtype followed by symmetry group and then
        # spherical harmonic terms
        anisotropy_terms[line[0]] = line[2:]

    count += 2
    for line in data[count:]: # first line is title
        count += 1
        if not line:
            break
        # each anisotropy line is atomtype followed by symmetry group and then
        # spherical harmonic terms
        anisotropy_axes.append(line)

    return anisotropy_terms, anisotropy_axes

###########################################################################


###########################################################################
def get_drude_charges(molecule):
    '''
    '''

    n_preskip_lines = 1 #First line is header for each atom

    with open(drudes_dir + molecule + '.drude','r') as f:
        data = [line.split() for line in f.readlines()]
        drude_charges = [i[0:1] + [float(j) for j in i[1:]] for i in data[n_preskip_lines:] ]

    return drude_charges
###########################################################################


###########################################################################
def get_charges(molecule):
    '''
    '''

    n_preskip_lines = 2 #First two lines are header for each atom

    with open(charges_dir + molecule + '.charges','r') as f:
        data = [line.split() for line in f.readlines()]
        charges = [float(i[-1]) for i in data[n_preskip_lines:] ]

    return charges
###########################################################################


###########################################################################
def read_constrained_params(molecule,param_constraints):
    '''
    '''


    with open(constraints_dir + molecule + '.constraints','r') as f:

        data = [line.split() for line in f.readlines()]
        natomtypes = int(data[0][0])

    # get a params
    iline = 1
    n_header_lines = 1 #First line of each section is a header
    for i_component in range(n_a_params) + [7]:
        iline += n_header_lines
        if i_component == -1:
            iline += 1
        for i_atom in range(natomtypes):
            atom = data[iline][0]
            params = [float(i) for i in data[iline][1:]]
            if atom in param_constraints.keys() and atom in atomtype_constraints:
                param_constraints[atom][i_component] = params
            iline += 1

    return param_constraints
###########################################################################


###########################################################################
def get_dispersion_coeffs(molecule):
    '''
    '''

    # Read in available force field parameters from zif_ff_parameters
    with open(dispersion_dir + molecule + '.disp','r') as f:
        data = f.readlines()
        lines = [line.split() for line in data]

    dispersion_coeffs = {}
    for line in lines:
        atom = line[0] 
        if atom == 'C6' or atom[0:2] == '--':
            # Ignore any header lines
            continue
        cn_data = [float(i) for i in line[1:]]
        dispersion_coeffs[atom] = cn_data

    return dispersion_coeffs
###########################################################################


###########################################################################
###########################################################################

###########################################################################
########################## Main Code ######################################

pwd = os.getcwd()

# Copy relevant files for the dimer pair to the dimer directory
sapt_file = sapt_dir + mon1 + '_' + mon2 + '.sapt'
# subprocess.call(['cp', sapt_file, '.'])

unconstrained_param_file = mon1 + '_' + mon2 + '_unconstrained.param'
# sapt_file = mon1 + '_' + mon2 + '.sapt'


if defaults.exponent_source.upper() == 'IP':
    # Get exponents from atomic ionization potentials (for Born-Mayer-IP force
    # fields only)
    exponents_mon1 = [ Exponent(i) for i in atoms_mon1]
    exponents_mon2 = [ Exponent(i) for i in atoms_mon2]
else:
    # Get exponents from the results of the ISA calculations
    exponents_mon1 = get_exponents(mon1)
    exponents_mon2 = get_exponents(mon2)

if separate_induction_exponents:
    induction_exponents_mon1 = get_induction_exponents(mon1)
    induction_exponents_mon2 = get_induction_exponents(mon2)
else:
    induction_exponents_mon1 = exponents_mon1
    induction_exponents_mon2 = exponents_mon2


# Aparams read in later, so temporariliy set to 0
Aparams_mon1 = [ [0,i,j] for (i,j) in zip(exponents_mon1,induction_exponents_mon1)] 
Aparams_mon2 = [ [0,i,j] for (i,j) in zip(exponents_mon2,induction_exponents_mon2)] 

# Get atomtypes from .sapt file
all_atoms = get_atomtypes(sapt_file)

param_constraints = {}
exponents = exponents_mon1 + exponents_mon2
induction_exponents = induction_exponents_mon1 + induction_exponents_mon2
unique_atoms = set(all_atoms)

# Get dispersion from .disp files
dispersion = get_dispersion_coeffs(mon1)
dispersion_mon2 = get_dispersion_coeffs(mon2)
dispersion.update(dispersion_mon2)

# Input default param values for each atomtype; leave A parameters blank for
# now
for atom in unique_atoms:
    default_a_params = [ [0] for a in range(n_a_params) ]
    cn = dispersion[atom]
    param_constraints[atom] = default_a_params + [cn,
                        [exp_scale*np.mean([exponents[i] 
                        for i in range(len(all_atoms)) 
                        if all_atoms[i] == atom ])],
                        [exp_scale*np.mean([induction_exponents[i] 
                        for i in range(len(all_atoms)) 
                        if all_atoms[i] == atom ])],
                                ]

# Update param_constraints with any hard constraints
for atom in defaults.lone_pair_flags:
    param_constraints[atom] = lp_constraints
param_constraints = read_constrained_params(mon1,param_constraints)

# Get drude charges from the results of the polarization calculations
drude_charges_mon1 = get_drude_charges(mon1)
drude_charges_mon2 = get_drude_charges(mon2)
drude_charges = [drude_charges_mon1, drude_charges_mon2]

# Create the parameter file
create_param_file(mon1,mon2,sapt_file, param_constraints, dispersion, drude_charges, unconstrained_param_file,
        constrain_Aparams=defaults.lone_pair_flags+atomtype_constraints, springcon=0.1)
subprocess.call(['cp',multipoles_dir + mon1 + multipoles_suffix,'.'])
subprocess.call(['cp',multipoles_dir + mon2 + multipoles_suffix,'.'])

# Run the Force Field Fitting Program for each of the parameter files listed
fitfiles = ['edrudes','exchange','electrostatics','induction','dhf','dispersion','total_energy']
if settings.fit_residuals:
    fitfiles += ['residual_energy']
coeffs_fitfile = settings.file_prefix + 'coeffs' + settings.file_suffix + '.out'

# Perform Fitting
FitFFParameters(sapt_file, unconstrained_param_file, coeffs_fitfile,
        slater_correction=use_slater_ff, fit_bii=settings.fit_bii,
        aij_combination_rule=settings.aij_combination_rule,
        bij_combination_rule=settings.bij_combination_rule,
        cij_combination_rule=settings.cij_combination_rule,
        )
for fitfile in fitfiles:
    subprocess.call(['mv', fitfile + '.dat', settings.file_prefix + fitfile + settings.file_suffix + '.dat'])

# Change back to the home directory
os.chdir(pwd)





###########################################################################
###########################################################################
