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
#from chemistry import io

###########################################################################
####################### Global Variables ##################################
error_message='''
---------------------------------------------------------------------------
Improperly formatted arguments. Proper usage is as follows:

$ 

(<...> indicates required argument, [...] indicates optional argument)
---------------------------------------------------------------------------
    '''

charges_dir = 'isa/point_charges/'
subprocess.call(['mkdir','-p',charges_dir])
isa_mom_file = 'ISA_L4.mom'
l2_isa_mom_file = 'ISA_L2.mom'
mulfit_out_file = 'ISA_L4_to_L0.mom'
out_dir = '/OUT/'

## mulfit = '/home/mvanvleet/scripts/charge_analyses/mulfit'
## mulfit_template = '/home/mvanvleet/templates/mulfit.inp'

###########################################################################
###########################################################################


###########################################################################
######################## Command Line Arguments ###########################
maindir = os.getcwd().replace("/scripts",'')
templatesdir = maindir + '/templates/'
inputdir = maindir + '/input/'
geometriesdir = maindir + '/geometries/'
isadir = maindir + '/isa/'


###########################################################################
###########################################################################


###########################################################################
########################## Main Code ######################################

# Get mona names
dimer_info_file = 'dimer_info.dat'
with open (inputdir + dimer_info_file,'r') as f:
    data = [ line.split() for line in f.readlines()]
itag = [ i[0] if i else [] for i in data ].index('MonA_Name')
mona = data[itag][1]
itag = [ i[0] if i else [] for i in data ].index('MonB_Name')
monb = data[itag][1]

if mona == monb:
    mons = [mona]
else:
    mons = [mona,monb]

for mon in mons:
    # Get list of atomtypes and elements from .clt file
    clt = isadir + '/' + mon + '.clt'
    with open(clt,'r') as f:
        lines = [line.split() for line in f.readlines()]

    start_keys = [ i[0] if len(i) else '' for i in lines ]
    istart = start_keys.index('I.P.')
    iend = start_keys.index('End',istart)

    atomtypes = set([ i[-1] for i in lines[istart+1:iend] ])


    # Check ISA output for convergence before computing charges
    outfile = isadir + '/' + mon + out_dir + mon + '.out'
    try:
        output = subprocess.check_output(['grep','"NO CONVERGENCE"', outfile])
        print 'WARNING!!!! Stockholder procedure did not converge!'
    except subprocess.CalledProcessError:
        pass

    # Print ISA multipole moments to file
    with open(outfile,'r') as f:
        lines = f.readlines()
    start_flag = '! Multipole moments for'
    start_line = [i for i in xrange(len(lines)) if start_flag in lines[i]][0]
    end_flag = 'Total molecular moments relative to origin'
    end_line = [i for i in xrange(len(lines)) if end_flag in lines[i]][0]

    isa_mom_path = isadir + mon + '_' + isa_mom_file
    print 'Writing ISA charges to:'
    print isa_mom_path

    data = [line.split() for line in lines[start_line:end_line]]
    with open(isa_mom_path,'w') as f:
        for line in lines[start_line:end_line]:
            f.write(line)


    # Print truncated ISA multipole moments to file
    multipoles = [ [ str(i) + str(j) + k for k in ['c','s'] ]
                   if j !=0
                   else
                   [str(i) + str(j)]
                        for i in range(5) 
                        for j in range(i+1)  ]
    multipoles = [ 'Q' + item for sublist in multipoles for item in sublist]

    lines = lines[start_line:end_line]

    l2_isa_mom_path = isadir + mon + '_' + l2_isa_mom_file
    print 'Writing truncated ISA charges to:'
    print l2_isa_mom_path
    with open(l2_isa_mom_path,'w') as f:
        collect_flag = True
        multipole_count = 0
        for i,line in enumerate(data):
            if not line:
                collect_flag = True
                f.write('\n')
            elif line[0] == '!':
                # Skip comments
                f.write(lines[i])
            elif collect_flag == True:
                # Write atom header line
                f.write(lines[i])
                multipole_count = 0
                collect_flag = False
            else:
                # Lines containing multipoles
                for l in line:
                    m = multipoles[multipole_count]
                    if int(m[1]) < 3:
                        f.write('\t\t{:5s} = {} \n'.format(m,l))
                    multipole_count += 1

###########################################################################
###########################################################################
