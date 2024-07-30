#!/usr/bin/env python

# built-in python modules
import numpy as np
from random import random
import math
from math import sqrt,cos,sin,acos,pi

# python modules created by mvanvleet; these need to be imported into your
# PYTHONPATH for this code to run properly.
from chemistry.stoichiometry import MolecularWeight
from chemistry import elementdata
from chemistry import io

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
    def __init__(self,parameterfile='generate_grid_settings.inp',run=True):
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
            print("##########################################################################")
            print('Final Oriented coordinates for Monomer A are:')
            coordinates = ['{:2} {:>14.6f} {:14.6f} {:>14.6f}'.format(*line) \
                            for line in self.mona]
            for line in coordinates:
                print(line)

            # Orient Monomer B (choice of scanvector here is irrelevant):
            self.monb = self.OrientMonomer(self.monb,[0,0,1],self.monb_origin)
            print("##########################################################################")
            print('Final Oriented coordinates for Monomer B are:')
            coordinates = ['{:2} {:>14.6f} {:14.6f} {:>14.6f}'.format(*line) \
                            for line in self.monb]
            for line in coordinates:
                print(line)
            print("##########################################################################")

            # Generate dimer configurations: 
            print('Generating '+str(self.n_points)+' configurations.')
            print('Output files will be of the form: ',self.output_name+'_id#.xyz')
            for nfile in range(self.n_points):
                reject_point = True
                fill = int(math.log10(self.n_points-1))+1
                filename = self.output_name+'_'+str(nfile).zfill(fill)+'.xyz'
                while reject_point == True:
                    # Choose new location for Monomer B's center of mass (COM) 
                    (r,theta,phi) = self.ChooseTranslationPoint(self.min_r,self.max_r,self.min_theta,self.max_theta,self.min_phi,self.max_phi)
                    drCOM = self.TranslateCOM(r,theta,phi)
                    # Choose rotation vector and angle (for orienting Monomer B)
                    (theta,b,c,d) = self.ChooseRotation()
                    # Determine new coordinates for monomer b, first by rotating
                    # and then by translating the COM
                    rotated_monb = [self.RotatePoint(theta,[b,c,d],point[1:]) for point in self.monb]
                    self.new_monb = [[rotated_monb[i][j] + drCOM[j] for j in range(3)] for i in range(len(self.monb))]
                    [self.new_monb[i].insert(0,self.monb[i][0]) for i in range(len(self.monb))]
                    # Employ rejection criteria for the chosen configuration. Upon
                    # rejection, cycle back through the loop. Otherwise, write the
                    # configuration to file.
                    reject_point = self.RejectPoint(self.mona,self.new_monb,self.cutoff_min, self.cutoff_max)
                dimer = self.mona + self.new_monb
                io.WriteCoordinates(dimer,filename)

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
            elif 'geometry_file' in line[0]:
                self.geometry_file = line[1] 
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

        # Read in monomer information from geometry file
        try:
            with open(self.geometry_file,'r') as f:
                lines = f.readlines()
        except IOError:
            message = 'Error! The geometry file you indicated does not exist!\n'+\
            'Check your input file to make sure your geometry input file names match.'
            raise SystemExit(message)
        data=[line.split() for line in lines[1:]]
        self.natoms_mona = data[0][0]
        mona_geo = []
        monb_geo = []
        monb=False
        for line in data[1:]:
            if len(line) != 1 and monb==False:
                mona_geo.append(line)
            elif monb==True:
                monb_geo.append(line)
            elif len(line) == 1:
                self.natoms_monb = line[0]
                monb=True
        self.mona = [[geo[0],float(geo[1]),float(geo[2]),float(geo[3])] for geo in mona_geo]
        self.monb = [[geo[0],float(geo[1]),float(geo[2]),float(geo[3])] for geo in monb_geo]

        # Determine the scan vector based on user input
        v = list(scan_vector.split(','))
        if scan_vector_type == 0:
            v = [int(item) for item in v]
            self.scan_vector = [self.mona[v[1]-1][i] - self.mona[v[0]-1][i] for i in range(1,4)] # -1 to deal with reindexing
        if scan_vector_type == 1:
            self.scan_vector = [float(coord) for coord in v]

        # Choose origin points for monomers a and b based on user input.
        if origin_typea == 2:
            self.mona_origin = [float(coord) for coord in mona_origin.split(',')]
        elif origin_typea == 1:
            self.mona_origin=self.mona[int(mona_origin)-1][1:]
        else:
            self.mona_origin=self.GetCOM(self.mona)
        if origin_typeb == 2:
            self.monb_origin = [float(coord) for coord in monb_origin.split(',')]
        elif origin_typeb == 1:
            self.monb_origin=self.monb[int(monb_origin)-1][1:]
        else:
            self.monb_origin=self.GetCOM(self.monb)

        print("##########################################################################")
        print("The following scan parameters have been selected:")
        print('Number of dimer configurations:',self.n_points)
        print('Input geometry file:',self.geometry_file)

        print('Origin point for monomer a (relative to input coordinates):',self.mona_origin)
        print('Origin point for monomer b (relative to input coordinates):',self.monb_origin)
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
    def GetCOM(self,coordinates):
        """Given an array 'coordinates' of the form 
        [[symbol1,x1,y1,z1],[symbol2,x2,y2,z2],...[symboln,xn,yn,zn]]
        which describes the coordinates of a molecule, returns the center of mass
        of said molecule. Units (generally Angstroms or Bohr)  are unchanged from
        input.
        """
        # COM formula: xCOM = sum(m_i*x_i)/M (sum over i=1,N); same for y and z
    
        # Total Mass of molecule:
        Mass = MolecularWeight(coordinates)
    
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
        r = random()*dr + r_min
    
        #Randomly choose theta: uniform distribution on the interval [theta_min,theta_max]
        dtheta = (theta_max - theta_min)
        theta = random()*dtheta + theta_min
    
        #Randomly choose phi: weighted distribution on the interval [phi_min,phi_max]
        dcosphi = (cos(phi_max) - cos(phi_min))
        phi = acos(random()*dcosphi + cos(phi_min))
    
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
        a = random()*da + a_min
        #Randomly choose b: uniform distribution on the interval [b_min,b_max]
        db = (b_max - b_min)
        b = random()*db + b_min
        #Randomly choose c: uniform distribution on the interval [c_min,c_max]
        dc = (c_max - c_min)
        c = random()*dc + c_min
        #Randomly choose d: uniform distribution on the interval [d_min,d_max]
        dd = (d_max - d_min)
        d = random()*dd + d_min
    
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
        return new_point.tolist()
####################################################################################################    
    
    
####################################################################################################    
    def OrientMonomer(self,monomer_coordinates,scan_vector=[0,0,1],new_origin=[0,0,0]):
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
        # Write monomer coordinates as a list
        coordinates = [[float(coord[i]) for i in range(1,4)] for coord in monomer_coordinates]
        # Translate coordinates such that new_origin is at the origin
        origin=new_origin
        trans_coords = [[coord[i]-origin[i] for i in range(3)] for coord in coordinates]
    
        # Rotate coordinates such that scan_vector is aligned with the z-axis
        scanvec = scan_vector/np.linalg.norm(scan_vector) #normalize
        z_axis = [0,0,1]
        rotation_vector = np.cross(scanvec,z_axis)
        rotation_angle = math.degrees(acos(np.dot(scanvec,z_axis)))
    
    
        rotated_coords = [self.RotatePoint(rotation_angle,rotation_vector,point) for point in trans_coords]
        [rotated_coords[i].insert(0,monomer_coordinates[i][0]) for i in range(len(monomer_coordinates))]
        return rotated_coords
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
                    rvec = [a[i] - b[i] for i in range(1,4)]
                    radius = sqrt(np.dot(rvec,rvec))
                    if radius < cutoff_min:
                        return True
                    min_separation = min(min_separation,radius)

        elif self.cutoff_type == 'vdw': # string should be 'vdw' to indicate use of VdW cutoff radius
            for a in mona_coords:
                vdw_a = elementdata.VdWRadius(a[0])
                for b in monb_coords:
                    vdw_b = elementdata.VdWRadius(b[0])
                    #r_cutoff = (vdw_a + vdw_b)*self.vdw_cutoff
                    rvec = [a[i] - b[i] for i in range(1,4)]
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
