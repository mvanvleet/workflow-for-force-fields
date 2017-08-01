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

###########################################################################
####################### Global Variables ##################################
error_message='''
---------------------------------------------------------------------------
Improperly formatted arguments. Proper usage is as follows:

$ get_sapt_fitting_data.py <dirname>

(<...> indicates required argument, [...] indicates optional argument)
---------------------------------------------------------------------------
    '''

maindir = os.getcwd().replace("/scripts",'')
sapt_script= maindir + '/scripts/get-sapt-fitting-data'
#'/home/mvanvleet/scripts/force_field_other/get-sapt-fitting-data'
templatesdir = maindir + '/templates/'
inputdir = maindir + '/input/'
#geometriesdir = maindir + '/geometries/'
saptdir = maindir + '/sapt/'
saptoutdir = saptdir + 'output/'
saptbaddir = saptdir + 'bad_calcs/'


# Slater test removal criterion
rtol = 2.00 #percentage
atol = 1.0 # mH




dimer_info_file = 'dimer_info.dat'
with open (inputdir + dimer_info_file,'r') as f:
    data = [ line.split() for line in f.readlines()]
itag = [ i[0] if i else [] for i in data ].index('MonA_Name')
mon1 = data[itag][1]
itag = [ i[0] if i else [] for i in data ].index('MonB_Name')
mon2 = data[itag][1]

###########################################################################
###########################################################################




###########################################################################
########################## Main Code ######################################


###########################################################################
# Remove slater test fails
names = subprocess.check_output('echo ' + saptoutdir + '*.out', shell=True)
names = names.split()

pwd = os.getcwd()
os.chdir(saptoutdir)

slater_data = subprocess.check_output('grep "RKS STATE 1.1 Energy" -B 9 ' +
'*out | grep "Slater Test"', shell=True)

slater_data = slater_data.split('\n')
slater_data = [ line.split() for line in slater_data ]

eint_data = subprocess.check_output('grep "EINT_DFT" ' + '*out', shell=True)
eint_data = eint_data.split('\n')
eint_data = [ line.split() for line in eint_data ]

os.chdir(pwd)

subprocess.call(['mkdir','-p',saptbaddir])
print 'Removing the following unconverged calculations:'
flag=False
for i,line in enumerate(names):
    # print every third line, which should contain the file name
    ## for eint_line in eint_data[count:]:
    ##     # Slater file should be ordered, but may have more entries than slater
    ##     # test data
    ##     count += 1
    ##     if eint_line[0] == line[0]:
    ##         break
    ## else:
    ##     print 'Couldn\'t find filename:'
    ##     print line
    ##     sys.exit()

    ## print eint_line
    ## print line
    ## print slater_data[i+1]
    ## print slater_data[i+2]

    slater_error1 = abs(float(slater_data[2*i][-1])) # units in Ha
    slater_error2 = abs(float(slater_data[2*i+1][-1]))
    slater_error = 1000*max(slater_error1, slater_error2)
    eint = abs(float(eint_data[i][-2])) # units in mH

    if slater_error/eint > rtol or slater_error > atol:
        print line, slater_error, eint, slater_error/eint*100
        subprocess.call(['mv',line,saptbaddir])
        flag = True

if not flag:
    print 'None'
print 'The above calculations have been moved to '
print saptbaddir
        
###########################################################################

###########################################################################
# Make SAPT fitting data

#For each dimer in the directory, read in atomtypes from the relevant .xyz file 

mon1_xyz, tmp = io.ReadCoordinates(templatesdir + mon1 + '.xyz')
mon2_xyz, tmp = io.ReadCoordinates(templatesdir + mon2 + '.xyz')

mon1_atoms = [i[0] for i in mon1_xyz]
mon2_atoms = [i[0] for i in mon2_xyz]

nmon1 = len(mon1_atoms)
nmon2 = len(mon2_atoms)

# Change to the relevant output directory and run the get-sapt-fitting-data
# script
os.chdir(saptoutdir)
sapt_file = mon1 + '_' + mon2 + '.sapt'
print 'Compiling SAPT information (this may take a minute).'
sapt_lines = subprocess.check_output([sapt_script])
sapt_lines = [line.split() for line in sapt_lines.split('\n')]
print 'Writing SAPT information to file.'

# Check that the SAPT file has the correct number of atoms for each monomer
if nmon1 != int(sapt_lines[0][0]) or nmon2 != int(sapt_lines[nmon1+1][0]):
    print 'nmon1 = ', nmon1
    print 'nmon2 = ', nmon2
    sys.exit('Incorrect number of atoms for one of the monomers.')

# nfiles and nlines are the number of dimer geometries and number of lines per
# information block, respectively
nfiles = sum([1 for i in sapt_lines if i == []]) - 1 #-1 to take care of blank line at end of file
nlines = len(sapt_lines)/nfiles

shift_mon1 = 1
shift_mon2 = nmon1 + 2
for i in xrange(nfiles):
    for iatom,atom in enumerate(mon1_atoms):
        sapt_lines[i*nlines + iatom + shift_mon1][0] = atom
    for iatom,atom in enumerate(mon2_atoms):
        sapt_lines[i*nlines + iatom + shift_mon2][0] = atom

# Write SAPT file with new atomtpyes to file.
with open(sapt_file,'w') as f:
    for line in sapt_lines[:-1]:
        template = '{:20s}'*len(line) + '\n'
        f.write(template.format(*line))

print 'The SAPT output file can be found in '
print saptoutdir

#os.chdir(homedir)

###########################################################################
###########################################################################
