#!/usr/bin/env python
"""A module controlling input (from either .xyz or .gxyz) and output (either
as list, .xyz, or .gxyz).

Includes the following functions:
    ReadCoordinates(file,bohr=False)
    WriteCoordinates(coordinates,outputfile='coordinates.xyz',title_text='Title Required',bohr=False): 

Last Updated: 04/27/15 by mvanvleet
"""

import elementdata
import constants 

###########################################################################
def ReadCoordinates(inputfile,bohr=False):
    """Given either a .xyz or .gxyz file, reads in coordinates and returns a
    list of coordinates of the following form-
    [[symbol1,x1,y1,z1],[symbol2,x2,y2,z2],...[symboln,xn,yn,zn]]
    - followed by the input title text. Coordinates and title text returned as
      a two item tuple, with coordinate units in Angstroms.

    inputfile should be a string; symbols are returned as strings, x,y,z as
    floats.
    """
    # Read in .xyz file, extract unformatted list of coordinates:
    with open(inputfile,'r') as f:
        lines=f.readlines()
    geometry=[line.split() for line in lines[2:]]

    # Remove second column (atomic number) from .gxyz files, if applicable:
    if '.gxyz' in inputfile:
        title_text=lines[0].strip('\n')
        geometry=[[geo[i] for i in [0,2,3,4]] for geo in geometry]
    else:
        title_text=lines[1].strip('\n')

    # Convert coordinates to floats
    if bohr:
        c = constants.bohr2ang
        geometry = [[geo[0],float(geo[1])*c,float(geo[2])*c,float(geo[3])*c] for geo in geometry]
    else:
        geometry = [[geo[0],float(geo[1]),float(geo[2]),float(geo[3])] for geo in geometry]

    return geometry,title_text
###########################################################################


###########################################################################
def WriteCoordinates(coordinates,outputfile='coordinates.xyz',\
        title_text='Title Required',bohr=False):
    """Given a list of coordinates and a string containing the name of an
    output file (ex. 'coordinates.xyz'), writes out the corresponding
    outputfile. If the file extension given is .gxyz, the output format is
    .gxyz; else, a .xyz file format is used. The format of the coordinates
    list should be as follows:
    [[symbol1,x1,y1,z1],[symbol2,x2,y2,z2],...[symboln,xn,yn,zn]]


    Options:
    title_text - String giving title text for .xyz or .gxyz file. Default is 'Title Required'.
    bohr - Boolean. Outputs units in bohr when set to true. Default is false.
    """
    if bohr==True:
        title_text+='; units in bohr'
        coordinates = [[item if type(item)==str else item*constants.ang2bohr for item in row] for row in coordinates]

    num_atoms = str(len(coordinates))+"\n" 
    title_text = title_text + "\n"

    #Output file (either as .gxyz or .xyz, depending on specified file extension):
    if '.gxyz' in outputfile:
        for line in coordinates:
            element=line[0]
            line.insert(1,elementdata.AtomicNumber(element))
        coordinates = ['{:2} {:>3} {:>16.8f} {:16.8f} {:>16.8f}\n'.format(*line) \
                        for line in coordinates]
        with open(outputfile,'w') as f:
            f.writelines(title_text)
            f.write('C1\n')
            f.writelines(coordinates)
        
    else:
        coordinates = ['{:2} {:>16.8f} {:16.8f} {:>16.8f}\n'.format(*line) \
                        for line in coordinates]
        with open(outputfile,'w') as f:
            f.write(num_atoms)
            f.writelines(title_text)
            f.writelines(coordinates)

    return
###########################################################################

if __name__=='__main__':
    #Testing conditions
    coordinates,title_text = ReadCoordinates('test.xyz')
    #print ReadCoordinates('test.gxyz')
    WriteCoordinates(coordinates,'test.gxyz',title_text=title_text,bohr=True)

