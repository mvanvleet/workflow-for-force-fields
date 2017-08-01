#!/usr/bin/env python
"""

Last Updated:
"""

# Standard modules
import numpy as np
import sys
import os
# mvanvleet specific modules
from chemistry import io

###########################################################################
####################### Global Variables ##################################
maindir = os.getcwd().replace("/scripts",'')
templatesdir = maindir + '/templates/'
inputdir = maindir + '/input/'
dispersiondir = maindir + '/dispersion/'

dimer_info_file = 'dimer_info.dat'

try:
    dispersion_scale = float(sys.argv[1])
except IndexError:
    dispersion_scale = 1.03

###########################################################################
###########################################################################


###########################################################################
########################## Main Code ######################################

# Read in monomer names from dimer info file
with open (inputdir + dimer_info_file,'r') as f:
    data = [ line.split() for line in f.readlines()]

itag = [ i[0] if i else [] for i in data ].index('MonA_Name')
mon1 = data[itag][1]
itag = [ i[0] if i else [] for i in data ].index('MonB_Name')
mon2 = data[itag][1]

if mon1 == mon2:
    mons = [mon1]
else:
    mons = [mon1, mon2]

# Read in dispersion information, scale, and print to file
for mon in mons:
    with open(dispersiondir + mon + '.cncoeffs','r') as f:
        data = [line.split() for line in f.readlines()]

    with open(dispersiondir + mon + '.disp','w') as f:
        template = '{:5s}' + '{:16.6f}'*4 + '\n'
        flag = False
        for line in data:
            if len(line) > 2 and line[2] == 'C6':
                if line[0] == line[1]:
                    atom = line[0]
                    flag = True
            elif flag:
                c6 = float(line[3])  * dispersion_scale
                c8 = float(line[5])  * dispersion_scale
                c10 = float(line[7]) * dispersion_scale
                c12 = float(line[9]) * dispersion_scale
                f.write(template.format(atom, c6, c8, c10, c12))
                flag = False
        else:
            continue

    sys.exit()




###########################################################################
###########################################################################
