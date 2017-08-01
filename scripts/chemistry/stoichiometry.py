#!/home/ehermes/local/bin/python
"""A module that contains functions for performing basic stoichiometric tasks.
Includes the following functions:
    MolecularWeight(xyz_file)

Last Updated: 08/20/13 by mvanvleet
"""

from chemistry import constants
from chemistry import elementdata

# Function to obtain molecular weight of a compound from a .xyz file
def MolecularWeight(coordinates):
    """Given an input list of coordinates of the following form:
    [[symbol1,x1,y1,z1],[symbol2,x2,y2,z2],...[symboln,xn,yn,zn]]
    returns the molecular weight (in g/mol) of
    the compound specified by the list.
    """
    elements=[atom[0] for atom in coordinates]

    # Calculate molecular weight
    molecular_weight=0.00
    for element in elements:
        try:
           atno=elementdata.AtomicNumber(element) 
        except KeyError:
            raise KeyError, 'Error: Element name "'+str(element)+\
                    '" not recognized. Check your .xyz file for errors.'
        #(column 3 of element_data list contains molecular weight data)
        molecular_weight += elementdata.Weight(atno) 
    return molecular_weight
        


if __name__=='__main__':
    from chemistry import io

    coordinates=io.ReadCoordinates('test.xyz')
    print MolecularWeight(coordinates)
