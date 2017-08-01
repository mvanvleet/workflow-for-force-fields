#!/home/ehermes/local/bin/python
"""Module controlling changes in geometry to .(g)xyz files.

Includes the following functions:
    GetDistance(inputfile,iatom1,iatom2,verbose=False):
    GetCOM(coordinates)


Last Updated: 11/04/14 by mvanvleet
"""

import numpy as np

import io
import stoichiometry
import elementdata


###########################################################################
def GetDistance(inputfile,iatom1,iatom2,verbose=False):
    """Given either a .xyz or .gxyz file [inputfile] and the index of 2 atoms
    [iatom1] and [iatom2], returns the interatomic distance between atom1 and
    atom2.
    
    Alternately, if the verbose option is set to True (default is false), the
    GetDistance function will return
    (name_atom1,name_atom2,interatomic_distance) as a tuple.

    Note that indexing starts at zero. Units are not converted within this
    function (i.e., the output distance units for the program are indentical
    to whatever the input units were).
    """

    coordinates = io.ReadCoordinates(inputfile)[0]
    xyz_atom1 = np.array(coordinates[iatom1][1:])
    xyz_atom2 = np.array(coordinates[iatom2][1:])
    interatomic_distance = np.linalg.norm(xyz_atom1 - xyz_atom2)

    if verbose:
        name_atom1 = coordinates[iatom1][0]
        name_atom2 = coordinates[iatom2][0]
        return name_atom1,name_atom2,interatomic_distance
    else:
        return interatomic_distance
###########################################################################


###########################################################################
def GetCOM(coordinates):
    """Given an array 'coordinates' of the form 
    [[symbol1,x1,y1,z1],[symbol2,x2,y2,z2],...[symboln,xn,yn,zn]]
    which describes the coordinates of a molecule, returns the center of mass
    of said molecule. Units (generally Angstroms or Bohr)  are unchanged from
    input.
    """
    # COM formula: xCOM = sum(m_i*x_i)/M (sum over i=1,N); same for y and z

    # Total Mass of molecule:
    Mass = stoichiometry.MolecularWeight(coordinates)

    # Generate list of each atom's atomic number, xyzcoordinate, and atomic
    # mass:
    atomic_numbers = [elementdata.AtomicNumber(atom[0]) for atom in coordinates]
    xcoords = [atom[1] for atom in coordinates]
    ycoords = [atom[2] for atom in coordinates]
    zcoords = [atom[3] for atom in coordinates]
    masses = [elementdata.Weight(element) for element in atomic_numbers]
    
    # Compute COM:
    xCOM = np.dot(masses,xcoords)/Mass
    yCOM = np.dot(masses,ycoords)/Mass
    zCOM = np.dot(masses,zcoords)/Mass

    return [xCOM,yCOM,zCOM]
###########################################################################


###########################################################################
def shift_coordinates(coordinates, vec):
    """Given an array 'coordinates' of the form 
    [[symbol1,x1,y1,z1],[symbol2,x2,y2,z2],...[symboln,xn,yn,zn]]
    which describes the coordinates of a molecule and a vector which describes
    how much to shift the coordinates of the molecule, returns an updated list
    of coordinates for the molecule in which the coordinates have all been
    shifted by an amount vec.
    Units (generally Angstroms or Bohr)  are unchanged from
    input.
    """

    # Convert coordinates to a numpy 
    coordinates = [ [coord[0]] + [coord[i+1] + vec[i] for i in range(3)]\
                        for coord in coordinates ] 

    return coordinates
###########################################################################


if __name__=='__main__':
    import sys
    try:
        for file in sys.argv[3:]:
            print GetDistance(file,int(sys.argv[1]),int(sys.argv[2]),verbose=False)
    except ValueError:
        print 'Usage: $ geometry.py iatom1 iatom2 xyzfile1 [xyzfile2...xyfileN]'
        print 'Indexing starts at zero.'
    if len(sys.argv) < 4:
        print 'Usage: $ geometry.py iatom1 iatom2 xyzfile1 [xyzfile2...xyfileN]'
        print 'Indexing starts at zero.'
