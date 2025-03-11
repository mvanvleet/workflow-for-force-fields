#!/usr/bin/env python

# built-in python modules
import numpy as np
import random 
import math
from math import sqrt,cos,sin,acos,pi
import os 

# python modules created by mvanvleet; these need to be imported into your
# PYTHONPATH for this code to run properly.
from chemistry.stoichiometry import MolecularWeight
from chemistry import elementdata
#from chemistry import io

class GenerateGridPoints():
    """Given a set of properly formatted input files, one a parameter file and the
    other a geometry file (see generate_grid_settings.inp and small.inp,
    respectively, for examples), creates a set of dimer geometries corresponding
    to mutual rotations and translations of the two monomers specified in the
    geometry file.

    The parameter file, assumed to be titled generate_grid_settings.inp, contains
    all relevant parameters and constraints regarding how the two monomers are
    to be rotated and translated relative to one another, as well as how many dimer
    geometry files are created. More detailed information on the parameter file
    and the various parameters contained within can be found in the comments
    section of generate_grid_settings.inp.

    Last Updated: 04/18/15 by mvanvleet
    """

####################################################################################################    
    def __init__(self,parameterfile='generate_grid_settings.inp',run=True,verbose=True):
        """Runs the main program of the GenerateGridPoints class.
        """

        if run:

            print('')
            print("##########################################################################")
            print("##########################################################################")
            print("Welcome to the generate_grid_points program!\n")
            print("For more help text about what's going on,",\
                    "refer to the example input file 'generate_grid_settings.inp'\n")
            self.ReadInput(parameterfile)

            # Orient each monomer according to settings defined in parameterfile
            # Orient Monomer A:
            self.mona = self.OrientMonomer(self.mona,self.scan_vector,self.mona_origin)
            print('')
            if verbose:
                print("##########################################################################")
                template='\t{:2} {:>14.6f} {:14.6f} {:>14.6f}'
                print('Final Oriented coordinates for Monomer A are:')
                for i,conf in enumerate(self.mona):
                    print(f"Conformer {i+1}:")
                    for atom,line in zip(self.mona_elements, conf):
                        print(template.format(atom,*line))

            # Orient Monomer B (choice of scanvector here is irrelevant):
            self.monb = self.OrientMonomer(self.monb,np.array([0,0,1]),self.monb_origin)
            if verbose:
                print("##########################################################################")
                print('Input coordinates for Monomer B are:')
                for i,conf in enumerate(self.monb):
                    print(f"Conformer {i+1}:")
                    for atom,line in zip(self.monb_elements, conf):
                        print(template.format(atom,*line))
                print("##########################################################################")

            # Generate dimer configurations: 
            print('Generating '+str(self.n_points)+' configurations.')
            print('Output files will be of the form: ',self.output_name+'_id#.xyz')
            random.seed()
            for nfile in range(self.n_points):
                reject_point = True
                fill = int(math.log10(self.n_points-1))+1
                filename = self.output_name+'_'+str(nfile).zfill(fill)+'.xyz'

                # Select a random conformer of Monomer A and B
                iconf_mona = random.randrange(self.mona.shape[0])
                iconf_monb = random.randrange(self.monb.shape[0])
                confa = self.mona[iconf_mona]
                confb = self.monb[iconf_monb]
                while reject_point == True:
                    # Choose new location for Monomer B's center of mass (COM) 
                    (r,theta,phi) = self.ChooseTranslationPoint(self.min_r,self.max_r,
                            self.min_theta,self.max_theta,self.min_phi,self.max_phi)
                    drCOM = self.TranslateCOM(r,theta,phi)
                    # Choose rotation vector and angle (for orienting Monomer B)
                    (theta,b,c,d) = self.ChooseRotation()
                    # Determine new coordinates for monomer b, first by rotating
                    # and then by translating the COM
                    rotated_confb = np.array([self.RotatePoint(theta,[b,c,d],atom) for atom in confb])
                    new_confb = rotated_confb + drCOM
                    ## [[rotated_confb[i][j] + drCOM[j] for j in range(3)] for i in range(len(self.monb))]
                    #[self.new_monb[i].insert(0,self.monb[i][0]) for i in range(len(self.monb))]
                    # Employ rejection criteria for the chosen configuration. Upon
                    # rejection, cycle back through the loop. Otherwise, write the
                    # configuration to file.
                    reject_point = self.RejectPoint(confa,new_confb,self.cutoff_min, self.cutoff_max)
                self.WriteCoordinates(iconf_mona, iconf_monb, confa,new_confb,filename)

            print('')
            print('Points successfully generated. Exiting now.')
            print("##########################################################################")
            print("##########################################################################")

        return
####################################################################################################    


####################################################################################################    
    def ReadInput(self,parameterfile='generate_grid.inp'):
        """Read input from two files: parameterfile, which contains default
        settings, and geometryfile, which contains monomer information. Note,
        geometryfile is specified within parameterfile.
        """
    
        # Read in parameter file 
        print(parameterfile)
        with open(parameterfile,'r') as f:
            lines = f.readlines()
            data1=[line.split('#') for line in lines] #eliminate comment lines

        # Read in all parameters
        data=[line[0].split() for line in data1]
        for line in data:
            if len(line) == 0: #ignore blank lines
                continue
            # General Scan Parameters:
            elif 'n_points' in line[0]:
                self.n_points = int(line[1])
            elif 'geometry_file_mona' in line[0]:
                self.geometry_file_mona = line[1] 
                print(self.geometry_file_mona)
            elif 'geometry_file_monb' in line[0]:
                self.geometry_file_monb = line[1]
                print(self.geometry_file_monb)
            elif 'output_name' in line[0]:
                self.output_name = line[1] 

            # Hard Sphere cutoff parameters:
            elif 'cutoff_type' in line[0]:
                self.cutoff_type = line[1].lower()
            elif 'cutoff_min' in line[0]:
                self.cutoff_min = float(line[1]) 
            elif 'cutoff_max' in line[0]:
                self.cutoff_max = float(line[1]) 

            # Constraints on translating Monomer B:
            elif 'min_r' in line[0]:
                self.min_r = float(line[1])
            elif 'max_r' in line[0]:
                self.max_r = float(line[1])
            elif 'min_theta' in line[0]:
                self.min_theta = pi*float(line[1])
            elif 'max_theta' in line[0]:
                self.max_theta = pi*float(line[1])
            elif 'min_phi' in line[0]:
                self.min_phi = pi*float(line[1])
            elif 'max_phi' in line[0]:
                self.max_phi = pi*float(line[1])

            # Origin and Scan Vector Parameters:
            elif 'mona_origin_type' in line[0]:
                origin_typea= int(line[1])
            elif 'mona_origin' in line[0]:
                mona_origin = line[1]
            elif 'monb_origin_type' in line[0]:
                origin_typeb= int(line[1])
            elif 'monb_origin' in line[0]:
                monb_origin = line[1]
            elif 'scan_vector_type' in line[0]:
                scan_vector_type = int(line[1])
            elif 'scan_vector' in line[0]:
                scan_vector = line[1] 
            else:
                print('Unrecognized option ',line[0])

        # Read in XYZ data for each monomer
        self.mona, self.mona_elements, self.natoms_mona = self.readGeometryData(self.geometry_file_mona)
        if self.geometry_file_mona == self.geometry_file_monb:
            self.monb, self.monb_elements, self.natoms_monb = self.mona, self.mona_elements, self.natoms_mona
        else:
            self.monb, self.monb_elements, self.natoms_monb = self.readGeometryData(self.geometry_file_monb)

        # Determine the scan vector based on user input
        v = list(scan_vector.split(','))
        if scan_vector_type == 0:
            v = [int(item) for item in v]
            if v[0] > self.natoms_mona or v[1] > self.natoms_mona:
                raise RuntimeError(f"""You have specified a scan vector that
                involves atom indices larger than the number of atoms in
                Monomer A ({self.geometry_file_mona}). Correct your geometry
                settings file (likely generate_grid_settings.inp).""")
            self.scan_vector = self.mona[:,v[1]-1,:] - self.mona[:,v[0]-1,:] 
        if scan_vector_type == 1:
            self.scan_vector = np.array([float(coord) for coord in v])

        # Choose origin points for monomers a and b based on user input.
        if origin_typea == 2:
            self.mona_origin = np.array(mona_origin.split(','),dtype=float)
        elif origin_typea == 1:
            self.mona_origin=self.mona[:,int(mona_origin)-1,:]
        else:
            self.mona_origin=self.GetCOM(self.mona_elements, self.mona)
        if origin_typeb == 2:
            self.monb_origin = np.array(monb_origin.split(','),dtype=float)
        elif origin_typeb == 1:
            self.monb_origin=self.monb[:,int(monb_origin)-1,:]
        else:
            self.monb_origin=self.GetCOM(self.monb_elements, self.monb)

        print("##########################################################################")
        print("The following scan parameters have been selected:")
        print('Number of dimer configurations:',self.n_points)
        print("##########################################################################")
        print('Input geometry file for mona:',self.geometry_file_mona)
        print('Origin point for monomer a (relative to input coordinates):',self.mona_origin)
        print("##########################################################################")
        print('Input geometry file for monb:',self.geometry_file_monb)
        print('Origin point for monomer b (relative to input coordinates):',self.monb_origin)
        print("##########################################################################")
        print('Scan Vector: ',self.scan_vector)
        print('')
        print('Monomer B will be placed relative to Monomer A according to the following constraints:')
        print("Radius (Angstroms): "+str(self.min_r)+" <= r < "+str(self.max_r))
        print("Azimuthal angle: "+str(self.min_theta)+" <= theta < "+str(self.max_theta))
        print("Polar angle: "+str(self.min_phi)+" <= phi < "+str(self.max_phi))
        if self.cutoff_type == 'absolute':
            print("Minimum allowed separation between any intermonomer contacts: ",\
                self.cutoff_min,'Angstroms')
            print("Maximum allowed separation between monomers (as measured by shortest contact length): ",\
                self.cutoff_max,'Angstroms')
        elif self.cutoff_type == 'vdw':
            print("Minimum allowed separation between any intermonomer contacts: ",\
                self.cutoff_min,'of the Van der Waals radii between atoms')
            print("Maximum allowed separation between monomers (as measured by shortest contact length): ",\
                self.cutoff_max,'of the Van der Waals radii between atoms')
        else:
            print()
            sys.exit('Cutoff type not recognized. Please specify either absolute or vdw.')
        print("##########################################################################")

        return(self.mona,self.monb)
####################################################################################################    
    

#################################################################################################### 
    def readGeometryData(self,geometry_file):
        """Returns the xyz coordinates, element names, and number of atoms
        from a specified geometry file."""

        # Open the conformers file
        with open(geometry_file, 'r') as f:
            lines = f.readlines()

        # Calculate the total number of atoms and conformers based on the file content
        natoms = int(lines[0])  # Number of atoms in each conformer
        print("Number of Atoms in Conformer: ", natoms)
        lines_per_conformer = natoms + 2  # Including the first line and the empty line

        nconformers = len(lines) // lines_per_conformer
        print(f"Number of conformers: {nconformers}")
        xyz = np.zeros((nconformers,natoms,3))
        elements = [line.split()[0] for line in lines[2:2+natoms]]

        # Iterate through the number of conformers and extract their coordinates
        for i in range(nconformers):
            start_index = 2 + i*(natoms + 2)  # Skip the header lines and the first two lines per conformer
            end_index = start_index + natoms
            if natoms != int(lines[start_index - 2]):
                raise RuntimeError(f"""Inconsistency in the number of atoms
                within the conformer file {geometry_file}.  Check your .xyz
                file for errors and make sure each structure within the .xyz
                file has {natoms} atoms.""")
            try:
                xyz[i] = np.array([line.split()[1:] for line in lines[start_index:end_index]],
                        dtype=float) #skip element names, get coordinates only
            except ValueError:
                print(f"""!!!!!!!!!! Error !!!!!!!!!!!!!
                Inconsistent data structure detected within the
                conformer file {geometry_file}. Make sure this file
                follows standard .xyz format.\n\n""")
                raise
            if elements != [line.split()[0] for line in lines[start_index:end_index]]:
                raise RuntimeError(f"""Inconsistency in the element ordering
                within the conformer file {geometry_file}.  Check your .xyz
                file for errors and make sure each structure within the .xyz
                file has the elements listed in the order {elements}.""")

        return xyz, elements, natoms
#################################################################################################### 


#################################################################################################### 
    def chooseConformer(self, monomer):
        import random
        conformer_shape = self.conformers.shape
        random_con = ((random.randint(0, conformer_shape[0] - 1)))
        random_con = self.conformers[random_con]
        #monomer = random_con

        return random_con
#################################################################################################### 


####################################################################################################    
    def GetCOM(self, elements, coordinates):
        """Given a list of elements and array of xyz coordintes, returns the center of mass
        of the molecule. Units (generally Angstroms or Bohr) are unchanged from
        input.
        """
    
        # Total Mass of molecule:
        Mass = MolecularWeight(elements)
    
        # Generate list of each atom's atomic number, xyzcoordinate, and atomic
        # mass:
        atomic_numbers = [elementdata.AtomicNumber(atom) for atom in elements]
        masses = np.array([elementdata.Weight(element) for element in atomic_numbers])

        # COM formula: xCOM = sum(m_i*x_i)/M (sum over i=1,N); same for y and z
        COM = np.einsum('ijk,j',coordinates,masses)/Mass
    
        return COM
####################################################################################################    

    
####################################################################################################    
    def ChooseTranslationPoint(self,r_min=1.0,r_max=5.0,theta_min=0.0,theta_max=2*pi,phi_min=0.0,phi_max=pi):
        """Randomly selects a point (r,theta,phi) along the interval [i_min,i_max]
        for i=r,theta,phi. Returns the tuple (r,theta,phi).

        Algorithm for random sphere point picking taken from
        http://mathworld.wolfram.com/SpherePointPicking.html
        """
        #Randomly choose r: uniform distribution on the interval [r_min,r_max]
        dr = (r_max - r_min)
        r = random.random()*dr + r_min
    
        #Randomly choose theta: uniform distribution on the interval [theta_min,theta_max]
        dtheta = (theta_max - theta_min)
        theta = random.random()*dtheta + theta_min
    
        #Randomly choose phi: weighted distribution on the interval [phi_min,phi_max]
        dcosphi = (cos(phi_max) - cos(phi_min))
        phi = acos(random.random()*dcosphi + cos(phi_min))
    
        return (r,theta,phi)
####################################################################################################    
    

####################################################################################################    
    def TranslateCOM(self,radius=0.0,theta=0.0,phi=0.0):
        """Moves the center of mass (COM) of a monomer to a specified
        coordinate (radius,theta,phi).
    
        Input:
        monomer_coordinates: list of the form
        [[symbol1,x1,y1,z1],[symbol2,x2,y2,z2],...[symboln,xn,yn,zn]]
        where 'symboli' is a string denoting an element and xi,yi,zi are floats
        denoting said element's spatial coordinates.
        translation_vector: a 3-element list of floats. Describes the vector
        (assumed to start at the origin) along which the monomer coordinates will
        be translated.
        radius: float describing the distance along the translation_vector
        direction to move the monomer's center of mass.
        theta: float ranging from 0 to 2pi describing the monomer's COM rotation
        (math convention used here)
        phi: float ranging from 0 to pi describing the monomer's COM rotation
    
        Output:
        Updated monomer_coordinates
        """
    
        drCOM = [radius*sin(phi)*cos(theta),radius*sin(theta)*sin(phi),radius*cos(phi)]
        return drCOM
####################################################################################################    
    

####################################################################################################    
    def ChooseRotation(self,a_min=0,a_max=360, b_min=-1,b_max=1,c_min=-1,c_max=1,d_min=-1,d_max=1):
        """Randomly selects a point (a,b,c,d) along the interval [i_min,i_max] for
        i=a,b,c,d. Returns the tuple (a,b,c,d).
        """
    
        #Randomly choose a: uniform distribution on the interval [a_min,a_max]
        da = (a_max - a_min)
        a = random.random()*da + a_min
        #Randomly choose b: uniform distribution on the interval [b_min,b_max]
        db = (b_max - b_min)
        b = random.random()*db + b_min
        #Randomly choose c: uniform distribution on the interval [c_min,c_max]
        dc = (c_max - c_min)
        c = random.random()*dc + c_min
        #Randomly choose d: uniform distribution on the interval [d_min,d_max]
        dd = (d_max - d_min)
        d = random.random()*dd + d_min
    
        return (a,b,c,d)
####################################################################################################    
    

####################################################################################################    
    def RotatePoint(self,theta=0,vector=[0,0,1],point=[1,2,3]):
        """Given an angle of rotation 'theta' and a vector (b,c,d) about which to
        rotate a point, computes the new position of a point 'point' in 3-space (given as a
        3-membered list) after a rotation of 'theta' degrees about the vector [b,c,d].
    
        This method uses quaternions to accomplish the transformation. For more
        information about the mathematics of quaternions, refer to 
        http://graphics.stanford.edu/courses/cs164-09-spring/Handouts/handout12.pdf
        """
    
        #Compute unit quaternion a+bi+cj+dk
        a = cos(math.radians(theta/2.0))
        [b,c,d] = vector
        if b == c == d == 0.00: #Deal with case where vector is ill-defined
            return point
        norm = sqrt(sin(math.radians(theta/2.0))**2/(b**2+c**2+d**2))
        [b,c,d] = [i*norm for i in [b,c,d]]
    
        # Compute quaternion rotation matrix:
        [a2,b2,c2,d2] = [a**2,b**2,c**2,d**2]
        [ab,ac,ad,bc,bd,cd] = [a*b,a*c,a*d,b*c,b*d,c*d]
    
        rotation = np.array([[ a2+b2-c2-d2 ,  2*bc-2*ad  ,  2*bd+2*ac  ],\
                             [  2*bc+2*ad  , a2-b2+c2-d2 ,  2*cd-2*ab  ],\
                             [  2*bd-2*ac  ,  2*cd+2*ab  , a2-b2-c2+d2 ]])
    
        # Compute rotation of point about the axis
        new_point = np.dot(rotation,point)
        return new_point
####################################################################################################    
    
    
####################################################################################################    
    def OrientMonomer(self,monomer_coordinates,
            scan_vector=np.array([0,0,1]),new_origin=np.array([0,0,0])):
        """Orients monomer a so that the monomer is centered according to
        new_origin and aligned such that scan_vector and the z-axis run
        parallel to one another.
    
        Input:
        monomer_coordinates for monomer a according to the form:
        [[symbol1,x1,y1,z1],[symbol2,x2,y2,z2],...[symboln,xn,yn,zn]]
        where 'symboli' is a string denoting an element and xi,yi,zi are floats
        denoting said element's spatial coordinates.
        scan_vector as a 3-element list. scan_vector is chosen by the user.
        new_origin. A point (relative to the input monomer_coordinates)
        that should serve as the origin. By default new_origin should be
        the center of mass of monomer a.
    
        Output:
        Updated monomer_coordinates.
        """

        trans_coords = monomer_coordinates - new_origin[:,np.newaxis,:]
    
        # Rotate coordinates such that scan_vector is aligned with the z-axis
        if len(scan_vector.shape) == 1:
            scan_vector = scan_vector[np.newaxis,:]
        scanvec = scan_vector/np.linalg.norm(scan_vector,axis=1,keepdims=True) #normalize
        z_axis = np.array([0,0,1])
        rotation_vector = np.cross(scanvec,z_axis)
        rotation_angle = np.arccos(np.dot(scanvec,z_axis))*180/np.pi
        ## print(rotation_angle)
        ## exit()
        ## rotation_angle = math.degrees(acos(np.dot(scanvec,z_axis)))
        if rotation_vector.shape[0] != trans_coords.shape[0]:
            if rotation_vector.shape[0] == 1:
                rotation_vector = np.repeat(rotation_vector, trans_coords.shape[0], axis=0)
                rotation_angle = np.repeat(rotation_angle, trans_coords.shape[0], axis=0)
            else:
                raise RuntimeError("""The given scan vector has incompatible
                dimensions with the given monomer geometry:
                Rotation Vector Shape: {rotation_vector.shape}
                Monomer Coordinates Shape: {trans_coords.shape[::2]}
                """)

        rotated_coords = np.zeros_like(trans_coords)
        for i,conformer in enumerate(trans_coords):
            for j,atom in enumerate(conformer):
                rotated_coords[i,j] = self.RotatePoint(rotation_angle[i], rotation_vector[i], atom)

        return rotated_coords
####################################################################################################    


####################################################################################################    
    def WriteCoordinates(self,iconf_mona, iconf_monb, confa,confb,filename):
        """Write a .xyz file corresponding to a chosen dimer configuration."""

        parentdir = os.path.dirname(filename)
        if parentdir != '':
            os.makedirs(os.path.dirname(filename),exist_ok=True)

        num_atoms = self.natoms_mona + self.natoms_monb
        title_text = 'MonA Conformer #{0}, MonB Conformer #{1}'.format(iconf_mona+1, iconf_monb+1)
        template = '{:2} {:>16.8f} {:16.8f} {:>16.8f}\n'
        with open(filename,'w') as f:
            f.write(f"{num_atoms}\n")
            f.write(title_text+"\n")
            for i,atom in enumerate(confa):
                f.write(template.format(self.mona_elements[i],*atom))
            for i,atom in enumerate(confb):
                f.write(template.format(self.monb_elements[i],*atom))

        return
####################################################################################################    
    
    
####################################################################################################    
    def RejectPoint(self,mona_coords,monb_coords,cutoff_min=1.9,cutoff_max=6.0):
        """Given coordinates for monomers a and b (mona_coords and monb_coords),
        returns True if any intermonomer pairs of atoms are within a distance
        'cutoff' from one another and false otherwise.
        """

        min_separation = cutoff_max + 0.1

        if self.cutoff_type == 'absolute':  # use this option if using an absolute cutoff (in A)
            for a in mona_coords:
                for b in monb_coords:
                    rvec = a - b
                    radius = sqrt(np.dot(rvec,rvec))
                    if radius < cutoff_min:
                        return True
                    min_separation = min(min_separation,radius)

        elif self.cutoff_type == 'vdw': # string should be 'vdw' to indicate use of VdW cutoff radius
            for i,a in enumerate(mona_coords):
                vdw_a = elementdata.VdWRadius(self.mona_elements[i])
                for j,b in enumerate(monb_coords):
                    vdw_b = elementdata.VdWRadius(self.monb_elements[j])
                    #r_cutoff = (vdw_a + vdw_b)*self.vdw_cutoff
                    rvec = a - b
                    radius = sqrt(np.dot(rvec,rvec))
                    vdw_radius = radius/(vdw_a + vdw_b)
                    if vdw_radius < cutoff_min:
                        return True
                    min_separation = min(min_separation,vdw_radius)

        else:
            print(self.cutoff_type, ' not a known cutoff type.')
            sys.exit('Exiting.')

        # Ensure monomers are sufficiently close to one another
        if min_separation > cutoff_max:
            return True
        else:
            return False

####################################################################################################    


if __name__=='__main__':

    GenerateGridPoints()
