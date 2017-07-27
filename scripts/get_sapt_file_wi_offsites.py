#!/usr/bin/env python
"""

Last Updated:
"""

# Standard modules
import numpy as np
import sys
import os
import subprocess

# mvanvleet specific modules
from chemistry import io
from chemistry.constants import ang2bohr

###########################################################################
####################### Global Variables ##################################
error_message='''
---------------------------------------------------------------------------
Improperly formatted arguments. Proper usage is as follows:

$ get_sapt_fitting_data.py <dirname>

(<...> indicates required argument, [...] indicates optional argument)
---------------------------------------------------------------------------
    '''

# Required scripts
sapt_script='/home/mvanvleet/scripts/force_field_other/get-sapt-fitting-data'
convert_to_orient_script='convert_xyz_to_orient.py'
convert_to_xyz_script='convert_orient_to_xyz.py'

# File Directories
maindir = os.getcwd().replace("/scripts",'')
geodir = maindir + '/geometries/xyz/'
newgeodir = maindir + '/geometries/xyz_wi_offsites/'
templatesdir = maindir + '/templates/'
inputdir = maindir + '/input/'
saptdir = maindir + '/sapt/'
saptoutdir = saptdir + 'output/'
saptbaddir = saptdir + 'bad_calcs/'

# File names and prefixes
sapt_prefix = '_pbe0_'
grid_settings_file = 'generate_grid_settings.inp'
orient_suffix = '.orient'
offsite_suffix = '_offsites.xyz'
sapt_ofile_prefix = 'offsites_'

# Monomer names
dimer_info_file = 'dimer_info.dat'
with open (inputdir + dimer_info_file,'r') as f:
    data = [ line.split() for line in f.readlines()]
itag = [ i[0] if i else [] for i in data ].index('MonA_Name')
mon1 = data[itag][1]
itag = [ i[0] if i else [] for i in data ].index('MonB_Name')
mon2 = data[itag][1]

dimer = mon1 + '_' + mon2

###########################################################################
###########################################################################

# Some command-line scripts (cp, mv, ls) have argument list limits; max batch
# and chunk prevent us from going over this limit
max_batch = 500
def chunk(seq, size=max_batch):
    return (seq[pos:pos + size] for pos in xrange(0, len(seq), size))


###########################################################################
########################## Main Code ######################################

# Make SAPT fitting data

# Copy .xyz files (*only from good SAPT calcs*) to new directory.
sapt_prefix = mon1 + '_' + mon2 + sapt_prefix
ofiles = subprocess.check_output('echo ' + saptoutdir + sapt_prefix + '*.out',shell=True)
ofiles = ofiles.split()

idnos = [ line.split(sapt_prefix)[-1].replace('.out','') for line in ofiles]
new_xyz = [ geodir + mon1 + '_' + mon2 + '_' + id + '.xyz' for id in idnos]

subprocess.call(['mkdir','-p',newgeodir])
for xyz_batch in chunk(new_xyz):
    subprocess.call(['cp'] + xyz_batch + [newgeodir])

# Copy generate grid settings file from old geometry directory to the new one.
subprocess.call(['cp', geodir +  grid_settings_file, newgeodir])

# Use oriented coordinates to generate new .inp file
xyz1 = io.ReadCoordinates(templatesdir + mon1 + '.xyz')[0]
xyz2 = io.ReadCoordinates(templatesdir + mon2 + '.xyz')[0]
template = '{:5s}' + '{:16.8f}'*3 + '\n'
with open(newgeodir + dimer + '.inp', 'w') as f:
    f.write('Dimer\n')
    f.write(str(len(xyz1))+ '\n')
    for line in xyz1:
        f.write(template.format(*line))
    f.write(str(len(xyz2))+ '\n')
    for line in xyz2:
        f.write(template.format(*line))

# Run convert to orient script
pwd = os.getcwd()
os.chdir(newgeodir)
orient_file = dimer + orient_suffix
subprocess.call(['rm',orient_file])
with open(orient_file,'w') as f:
    #subprocess.call([convert_to_orient_script,dimer + '.inp',dimer])
    subprocess.call([convert_to_orient_script,dimer + '.inp',dimer], stdout=f)


# Move xyz files to new directory
subprocess.call(['mv',dimer + '.inp','no_offsites.inp'])
subprocess.call(['mkdir','-p','old'])
subprocess.call('mv *xyz old/',shell=True)

# Generate new .inp file containing offsite positions
offsites1 = io.ReadCoordinates(templatesdir + mon1 + offsite_suffix )[0]
offsites2 = io.ReadCoordinates(templatesdir + mon2 + offsite_suffix )[0]
template = '{:5s}' + '{:16.8f}'*3 + '\n'
with open(newgeodir + dimer + '.inp', 'w') as f:
    f.write('Dimer\n')
    f.write(str(len(xyz1)+len(offsites1))+ '\n')
    for line in xyz1 + offsites1:
        f.write(template.format(*line))
    f.write(str(len(xyz2)+len(offsites2))+ '\n')
    for line in xyz2 + offsites2:
        f.write(template.format(*line))

subprocess.call([convert_to_xyz_script,orient_file])

os.chdir(pwd)

# Change to the relevant output directory and run the get-sapt-fitting-data
# script
os.chdir(saptoutdir)
sapt_file = mon1 + '_' + mon2 + '.sapt'
print 'Reading in SAPT information.'
with open(saptoutdir + sapt_file,'r') as f:
    sapt_lines = [line.split() for line in f.readlines()]
sapt_lines.append([])
## sapt_lines = subprocess.check_output([sapt_script])
## sapt_lines = [line.split() for line in sapt_lines.split('\n')]
print 'Writing SAPT information to file.'

# Check that the SAPT file has the correct number of atoms for each monomer
nmon1 = len(xyz1)
nmon2 = len(xyz2)
if len(xyz1) != int(sapt_lines[0][0]) or len(xyz2) != int(sapt_lines[nmon1+1][0]):
    print 'nmon1 = ', len(xyz1)
    print 'nmon2 = ', len(xyz2)
    sys.exit('Incorrect number of atoms for one of the monomers.')

# nfiles and nlines are the number of dimer geometries and number of lines per
# information block, respectively
nfiles = sum([1 for i in sapt_lines if i == []]) - 1 #-1 to take care of blank line at end of file
nlines = len(sapt_lines)/nfiles

shift_mon1 = 1
shift_mon2 = len(xyz1) + 2

sapt_energies = []
for i in range(nfiles):
    start = i*nlines + nmon1 + nmon2 + 2
    end = (i+1)*nlines
    sapt_energies.append(sapt_lines[start:end])

# Write SAPT file with new atomtpyes to file.
geofiles = subprocess.check_output('echo ' + newgeodir + '*.xyz',shell=True)
#geofiles = geofiles.split('\n')
geofiles = geofiles.split()
with open(sapt_ofile_prefix + sapt_file,'w') as f:
    for i in range(nfiles):
        # Write geometry
        xyz = io.ReadCoordinates(geofiles[i])[0]
        f.write(str(nmon1 + len(offsites1)) + '\n')
        for j, line in enumerate(xyz):
            if j == len(offsites1) + nmon1:
                f.write(str(nmon2 + len(offsites2)) + '\n')
            template = '{:5s}' + '{:16.8f}'*3 + '\n'
            # Convert to bohr
            line = [line[0]] + [k*ang2bohr for k in line[1:]]
            f.write(template.format(*line))
        # Write energies
        for line in sapt_energies[i]:
            template = '{:20s}'*len(line) + '\n'
            f.write(template.format(*line))
            
print 'The SAPT output file can be found in '
print saptoutdir

###########################################################################
###########################################################################
