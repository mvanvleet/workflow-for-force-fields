#!/usr/bin/env python
import sys
import os
import subprocess
import math

'''
This program takes, as input, a generic template file, and creates
input files using geometries from all .xyz files that are
located in the pwd and match a user-identified prefix.

File prefix is user specified. Template  file should always be of
the form 'template_file_prefix'_template.com . Geometry files should be of the form
'geometry_file_prefix'_id.xyz and should follow standard file format for .xyz
files. id can either be an identification number or simply a string.

Output files will be of the form 'template_file_prefix'_id.com, where id
is an identification tag that matches the correspoinding .xyz file.

Usage: $ create_molpro_files_from_xyz.py <template_file_prefix> <geometry_file_prefix>

Last Updated: 10/04/13 by mvanvleet
'''

##################### Input Parameters Section ###########################
error_message='''
---------------------------------------------------------------------------
Improperly formatted arguments. Proper usage is as follows:

$ create_molpro_files_from_xyz.py <template_file_prefix> <geometry_file_prefix>
---------------------------------------------------------------------------
    '''
# template_file_prefix refers to input template file (see above for proper usage).
try:
    template_file_prefix=sys.argv[1]
except IndexError:
    print(error_message)
    sys.exit()
# geometry_file_prefix refers to input .gxyz files (see above)
try:
    geometry_file_prefix=sys.argv[2]
except IndexError:
    print(error_message)
    sys.exit()
##################### End Input Parameters Section ########################


##################### Begin Script #######################################
#read in data from input template file
input_file = open(template_file_prefix+'_template.com', 'r')
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

num_geometries = subprocess.getoutput("ls " + geometry_file_prefix + '_*.xyz | grep -c '+  geometry_file_prefix) 
print('Writing ' + str(num_geometries) + ' input file(s).')

#write all relevant .com files
xyzlist = subprocess.getoutput("ls " + geometry_file_prefix +'_*.xyz').split()
idlist = [line.replace(geometry_file_prefix+"_",'').replace('.xyz','') for line in xyzlist]

for ifile in idlist:
    #designate output files with sequential file numbers
    output_file_name = template_file_prefix + '_' +ifile + '.com'

    #read in geometry from corresponding .xyz file and format for molpro
    geometry_file = geometry_file_prefix + '_'+ifile+'.xyz'
    with open(geometry_file,'r') as f:
        geometry_block= f.readlines()
    geometry_block[0] = "!"+geometry_block[0] #comment out first two lines of geometry
    geometry_block[1] = "!"+geometry_block[1]
    geometry_block.insert(0,"geometry={\n")
    geometry_block.append("\t}\n")

    #actually write output file
    output_file = open(output_file_name, 'w')
    for line in prelines:
        output_file.write(line)
    for line in geometry_block:
        output_file.write(line)
    for line in postlines:
        output_file.write(line)
    output_file.close()

print('Successfully wrote input files.')

#################### End Script ###########################################
