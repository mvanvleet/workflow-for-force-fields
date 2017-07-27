#!/usr/bin/env python
"""

Last Updated:
"""

# Standard modules
import numpy as np
import sys
import subprocess
import os

###########################################################################
####################### Global Variables ##################################
maindir = os.getcwd().replace("/scripts",'')
templatesdir = maindir + '/templates/'
inputdir = maindir + '/input/'
geometriesdir = maindir + '/geometries/'
saptdir = maindir + '/sapt/monomer_calcs/'

dimer_info_file = 'dimer_info.dat'

ip_file = '''
***, Ionization Energy

memory,100,m

basis=avtz

nosym;noorient;angstrom

geometry={{
{:s}
    }}

charge={:d}
{{rks,pbe0}}
e_atom=energy

charge={:d}
{{ks,pbe0}}
e_ion=energy

e_ionization=e_ion-e_atom
'''


###########################################################################
###########################################################################


###########################################################################
######################## Command Line Arguments ###########################


###########################################################################
###########################################################################


###########################################################################
########################## Main Code ######################################

subprocess.call(['mkdir','-p',saptdir])

# Read in monomer names from dimer file
with open(inputdir + dimer_info_file,'r') as f:
    data = [ line.split() for line in f.readlines()]

itag = [ i[0] if i else [] for i in data ].index('MonA_Name')
mona = data[itag][1]
itag = [ i[0] if i else [] for i in data ].index('MonA_Charge')
q_mona = int(data[itag][1])
itag = [ i[0] if i else [] for i in data ].index('MonB_Name')
monb = data[itag][1]
itag = [ i[0] if i else [] for i in data ].index('MonB_Charge')
q_monb = int(data[itag][1])

# Write and submit monomer A ionization potential file
with open(templatesdir + mona + '.xyz','r') as f:
    lines = f.readlines()
lines[0] = '! ' + lines[0]
lines[1] = '! ' + lines[1]

with open(saptdir + mona + '_ip.com','w') as f:
    f.write(ip_file.format(''.join(lines),q_mona,q_mona + 1))
subprocess.call(['qmolpro2012',saptdir+mona+'_ip.com'])

# Write and submit monomer B ionization potential file
if mona == monb:
    # Don't repeat calculation if uneccesary
    sys.exit()

with open(templatesdir + monb + '.xyz','r') as f:
    lines = f.readlines()
lines[0] = '! ' + lines[0]
lines[1] = '! ' + lines[1]

with open(saptdir + monb + '_ip.com','w') as f:
    f.write(ip_file.format(''.join(lines),q_monb,q_monb + 1))
subprocess.call(['qmolpro2012',saptdir+monb+'_ip.com'])



###########################################################################
###########################################################################
