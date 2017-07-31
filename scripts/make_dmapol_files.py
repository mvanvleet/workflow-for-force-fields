#!/usr/bin/env python
import sys
import os
import subprocess

#Local imports
from chemistry import elementdata, io

'''Usage: python make_dispersion_files.py'''

description='''
This program takes, as input, a generic template file, and creates
input files using geometries from all .gxyz files that are
located in the pwd and match a user-identified prefix.

File prefix is user specified. Template  file should always be of
the form 'template_file'_template.inp . Geometry files should be of the form
'geometry_file'_id.gxyz and should follow standard file format for .gxyz
files. id can either be an identification number or simply a string.

Output files will be of the form 'template_file'_id.com, where id
is an identification tag that matches the correspoinding .gxyz file.
'''

epilog='''
Last Updated: 07/31/14 by mvanvleet

'''

submit_template = '''#!/bin/bash

runcamcasp.py {0} --clt {0}.clt --direct -q default --nproc 20 --ifexists delete
'''

maindir = os.getcwd().replace("/scripts",'')
templatesdir = maindir + '/templates/'
inputdir = maindir + '/input/'
geometriesdir = maindir + '/geometries/'
dispersiondir = maindir + '/dispersion/'

ip_file = templatesdir + 'ips.dat'
template_file = templatesdir + 'dispersion_template.clt'


dimer_info_file = 'dimer_info.dat'
with open (inputdir + dimer_info_file,'r') as f:
    data = [ line.split() for line in f.readlines()]

itag = [ i[0] if i else [] for i in data ].index('MonA_Name')
mona = data[itag][1]
itag = [ i[0] if i else [] for i in data ].index('MonB_Name')
monb = data[itag][1]

itag = [ i[0] if i else [] for i in data ].index('MonA_Charge')
q_mona = int(data[itag][1])
itag = [ i[0] if i else [] for i in data ].index('MonB_Charge')
q_monb = int(data[itag][1])

subprocess.call(['mkdir','-p',dispersiondir])

if mona == monb:
    mons = [mona]
else:
    mons = [mona,monb]

for mon in mons:
    atomtype_file = inputdir + mon + '.atomtypes'
    geometry_file = templatesdir + mon + '.xyz'


    ##################### Begin Script #######################################
    #read in data from input template file
    input_file = file(template_file, 'r')
    prelines = []
    postlines = []
    before_geometry_block = True
    for line in input_file:
        if 'GEOMETRY_BLOCK_GOES_HERE' not in line and before_geometry_block:
            prelines.append(line)
            continue
        elif 'GEOMETRY_BLOCK_GOES_HERE' in line:
            before_geometry_block = False
        elif not before_geometry_block:
            postlines.append(line)
            continue
    input_file.close()

    #read in geometry from corresponding .gxyz file and format for gamess 
    print geometry_file
    xyz = io.ReadCoordinates(geometry_file)[0]
    ## with open(geometry_file,'r') as f:
    ##     geometry_block = f.readlines()[2:] #skip title lines

    # Read in atomtypes from atomtypes file
    with open(atomtype_file,'r') as f:
        atomtypes = [ line.split()[0] for line in f.readlines()[2:]] #skip title lines


    # Ammend geometry_block to include atomtypes
    geometry_block = []
    for i,line in enumerate(xyz):
        template = '{:2} {:>3} {:>14.6f} {:14.6f} {:>14.6f} \tType {:5} \n'
        args1 = [line[0] + str(i), elementdata.AtomicNumber(line[0]), line[1], line[2], line[3], atomtypes[i] ]
        geometry_block.append(template.format(*args1))



    #actually write output file
    output_file_name = dispersiondir + mon + '.clt'
    output_file = file(output_file_name, 'w')
    for line in prelines:
        output_file.write(line)
    for line in geometry_block:
        output_file.write(line)
    for line in postlines:
        output_file.write(line)
    output_file.close()


    # Substitute in ionization potential and molecule name
    #molecule = geometry_file.replace('.gxyz','')
    molecule = mon
    charge = q_mona if mon == mona else q_monb
    ip = subprocess.check_output(['grep','-iw',molecule,ip_file])
    ip = float(ip.split()[2]) # 2nd column contains i.p.

    fill_items = ['FILL_IP', 'FILL_CHARGE', 'FILL_MOLECULE']
    items = [ip, charge, molecule]
    for [fill,item] in zip(fill_items,items):
        subprocess.call(['sed','-i',"s/"+fill+'/'+str(item)+'/',output_file_name])

    print 'Successfully wrote input file.'

    # Make submit script
    with open(dispersiondir + 'submit_' + mon + '.sh','w') as f:
        f.write(submit_template.format(mon))

    #################### End Script ###########################################
