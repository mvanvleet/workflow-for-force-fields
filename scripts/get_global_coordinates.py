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
geometriesdir = maindir + '/geometries/'

#grid_ifile='generate_grid_settings.inp'
dimer_info_file = 'dimer_info.dat'

grid_ofile = 'generate_grid_settings.out'
outline1 = 'Final Oriented coordinates for Monomer A are:'
outline2 = 'Final Oriented coordinates for Monomer B are:'

###########################################################################
###########################################################################


###########################################################################
########################## Main Code ######################################

## # Read in monomer names from dimer info file
with open (templatesdir + dimer_info_file,'r') as f:
    data = [ line.split() for line in f.readlines()]

itag = [ i[0] if i else [] for i in data ].index('MonA_Name')
mon1_name = data[itag][1]
itag = [ i[0] if i else [] for i in data ].index('MonB_Name')
mon2_name = data[itag][1]

# Get geometries from geometries output file
print 'Reading global coordinates from:', geometriesdir + grid_ofile
with open(geometriesdir + grid_ofile,'r') as f:
    lines = f.readlines()

flag1=False
flag2=False
mon1_coords = []
mon2_coords = []
for line in lines:
    line = line.replace('\n','')
    if not line: 
        continue
    elif line[0] == '#':
        flag1 = False
        flag2 = False
    elif line == outline1:
        flag1 = True
    elif line == outline2:
        flag2 = True
    elif flag1:
        line = line.split()
        mon1_coords += [[line[0]] + [float(i) for i in line[1:]]]
    elif flag2:
        line = line.split()
        mon2_coords += [[line[0]] + [float(i) for i in line[1:]]]
    else:
        continue

# Write geometries to global coordinates file
title_text = mon1_name + '; global coordinates'
ofile = templatesdir + mon1_name + '.xyz'
print 'Writing Monomer 1 global coordinates to:',ofile
io.WriteCoordinates(mon1_coords,ofile,title_text=title_text)

title_text = mon2_name + '; global coordinates'
ofile = templatesdir + mon2_name + '.xyz'
print 'Writing Monomer 2 global coordinates to:',ofile
io.WriteCoordinates(mon2_coords,ofile,title_text=title_text)

###########################################################################
###########################################################################
