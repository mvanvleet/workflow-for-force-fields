#!/usr/bin/env python
import numpy as np

def add_midbond(xyz_file,i_atoms):
    """Given a specified .xyz file and atoms to place midbond functions in
    between, calculates the location of said midbond function and appends
    the position of said midbond function to the end of the .xyz file.

    xyz_file should be a file located in the same directory where this
    code is being run. (At some point I may fix this).
    i_atoms should be a list of atom indices between which a midbond function
    will be placed. Atom indexing begins at zero.
    """

    #define parameters
    append_to_new_file=False 

    #read in data from xyz file
    input_file = open(xyz_file, 'r')
    lines = input_file.readlines()
    num_atoms=lines[0] #store number of atoms
    save_text = lines[1] #comment line in .xyz not tab delimited; stored separately
    split_lines = [ lin.split() for lin in lines ]
    input_file.close()
    
    #create coordinates for midbond function
    coords = [ line[1:] for line in split_lines[2:] ]
    coords = np.array(coords,dtype=np.float)
    pos_midbond = np.average(coords[i_atoms],axis=0)
    ## pos_atom1 = [float(split_lines[i_atom1+2][j]) for j in range(1,4)]
    ## pos_atom2 = [float(split_lines[i_atom2+2][j]) for j in range(1,4)]
    ## pos_midbond = [(x1+x2)/2 for x1,x2 in zip(pos_atom1,pos_atom2)]
    pos_midbond = ["%0.15f" %i for i in pos_midbond]
    pos_midbond.insert(0,'Be')
    
    #update number of atoms line
    num_atoms = str(int(num_atoms)+1) + "\n"
    #write new xyz file with midbond functions
    split_lines.append(pos_midbond)
    split_lines=[map(str,ln) for ln in split_lines]
    coordinates = ['{:2} {:>22} {:>22} {:>22}\n'.format(*line) \
                    for line in split_lines[2:]]
    
    if append_to_new_file==True:
        output_file=open(xyz_file.replace('.xyz','_midbond.xyz','w'))
    else:
        output_file=open(xyz_file,'w')
    output_file.writelines(num_atoms)
    output_file.writelines(save_text)
    output_file.writelines(coordinates)
    output_file.close()
    
if __name__=='__main__':
    import sys
    try:
        for file in sys.argv[3:]:
            add_midbond(file,[int(sys.argv[1]),int(sys.argv[2])])
    except (IndexError,TypeError):
        print 'Error:Incorrectly formatted arguments and/or xyz file.'
        print "Usage: $ python add_midbond.py <atom1 index> <atom2 index> <file.xyz> [file2.xyz] [etc.]"
        print
        raise
