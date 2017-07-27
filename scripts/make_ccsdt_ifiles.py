#!/usr/bin/env python
"""

Usage:
    $ python make_ccsdt_ifles.py

Last Updated: Nov 16, 2015 by mvanvleet
"""

# Standard modules
import numpy as np
import sys
import os
import subprocess

# mvanvleet specific modules
from chemistry import io, geometry, elementdata
#from generate_geometries.generate_grid_points import GenerateGridPoints
from add_midbond import add_midbond

###########################################################################
####################### Global Variables ##################################
error_message='''
---------------------------------------------------------------------------
Improperly formatted arguments. Proper usage is as follows:

$ 

(<...> indicates required argument, [...] indicates optional argument)
---------------------------------------------------------------------------
    '''

maindir = os.getcwd().replace("/scripts",'')
scriptsdir = maindir + '/scripts/'
templatesdir = maindir + '/templates/'
inputdir = maindir + '/input/'
geometriesdir = maindir + '/geometries/'
ccsdtdir = maindir + '/ccsdt/input/'

dimer_info_file = 'dimer_info.dat'


###########################################################################
###########################################################################




###########################################################################
########################## Main Code ######################################
class CreateSAPTInputFile():
    '''
    '''

###########################################################################
    def __init__(self, dimer_file, input_path='input/', templates_path='templates/'):
        '''
        '''
        # Set some default behavior
        self.ignorecase = True

        # Read in monomer names
        with open (input_path + dimer_file,'r') as f:
            data = [ line.split() for line in f.readlines()]

        itag = [ i[0] if i else [] for i in data ].index('MonA_Name')
        self.mona = data[itag][1]
        itag = [ i[0] if i else [] for i in data ].index('MonB_Name')
        self.monb = data[itag][1]

        itag = [ i[0] if i else [] for i in data ].index('MonA_Charge')
        self.q_mona = int(data[itag][1])
        itag = [ i[0] if i else [] for i in data ].index('MonB_Charge')
        self.q_monb = int(data[itag][1])

        # Set input variables as class variables
        self.dimer_file = dimer_file
        self.templates_path = templates_path
        self.input_path = input_path

        # Get geometries and molecular centers for each monomer
        self.xyz_mona,self.title_mona = self.read_xyz_file(self.mona)
        self.xyz_monb,self.title_monb = self.read_xyz_file(self.monb)

###########################################################################


###########################################################################
    def read_xyz_file(self, monomer_name):
        '''Given a string denoting the monomer name and a path to the
        directory containing that monomer's .xyz file, reads in the geometry
        for that monomer.
        '''

        if self.ignorecase:
            monomer_name = monomer_name.lower()

        ext = '.xyz'
        path = self.templates_path

        return io.ReadCoordinates(path + monomer_name + ext)
###########################################################################


###########################################################################
    def add_midbonds(self,input_dir):
        '''
        '''

        # Read in midbond information
        with open (inputdir + self.dimer_file,'r') as f:
            data = [ line.split() for line in f.readlines()]

        tags = [ i[0] if i else [] for i in data ]
        indices = [i for i,x in enumerate(tags) if x == 'midbond']

        # Read in xyz file(s)
        #xyz_files = subprocess.check_output('ls '+input_dir+'/*xyz',shell=True)
        xyz_files = subprocess.check_output('echo '+input_dir+'/*xyz',shell=True)
        xyz_files = xyz_files.split()

        xyz_coords = []
        for xyz in xyz_files:
            coords, title_text = io.ReadCoordinates(xyz)
            xyz_coords.append(coords)

        # Place midbonds for each .xyz file
        for i in indices:
            i_mon1 = data[i][1]
            i_mon2 = data[i][2]

            print 'Inserting midbonds between',
            if i_mon1.upper() != 'COM':
                print 'atom',
            print i_mon1 + ' on mon A and',
            if i_mon2.upper() != 'COM':
                print 'atom',
            print i_mon2 + ' on mon B.'

            for c,coords in enumerate(xyz_coords):
                if i_mon1.upper() == 'COM':
                    pos1 = geometry.GetCOM(coords[0:len(self.xyz_mona)])
                else:
                    # pos1 = coords[int(i_mon1)+1][1:]
                    pos1 = coords[int(i_mon1)-1][1:]

                if i_mon2.upper() == 'COM':
                    istart = len(self.xyz_mona)
                    pos2 = geometry.GetCOM(coords[istart:])
                else:
                    pos2 = coords[int(i_mon2)+len(self.xyz_mona)-1][1:]

                midbond = ['Be'] + [(i[0] + i[1])/2 for i in zip(pos1,pos2) ]

                coords.append(midbond) 

        # Rewrite xyz files
        for ofile,coords in zip(xyz_files,xyz_coords):
            io.WriteCoordinates(coords,ofile,title_text='Dimer with midbonds')

        return
###########################################################################


###########################################################################
    def fill_ccsdt_template(self,dir,\
            template_file='ccsdt_template.com',\
            ):
        '''
        '''

        template_filepath = self.templates_path + template_file
        filepath = dir+'/'+self.mona+'_'+self.monb+'_'+ template_file

        # Default input fill fields should be constant
        fill_mona = 'FILL_MONA'
        fill_monb = 'FILL_MONB'
        fill_i_mona = 'FILL_I_MONA'
        fill_i_monb = 'FILL_I_MONB'
        fill_q_dimer = 'FILL_Q_DIMER'
        fill_q_mona = 'FILL_Q_MONA'
        fill_q_monb = 'FILL_Q_MONB'
        fill_items = [ fill_mona, fill_monb, fill_i_mona, fill_i_monb,
                        fill_q_dimer, fill_q_mona, fill_q_monb,
                     ]

        mona = self.mona
        monb = self.monb
        i_mona = ','.join(str(i) for i in range(1,1+len(self.xyz_mona)))
        i_monb = ','.join(str(i) for i in range(1+len(self.xyz_mona),1+len(self.xyz_monb + self.xyz_mona)))

        q_mona = self.q_mona
        q_monb = self.q_monb
        q_dimer = q_mona + q_monb

        items = [ mona, monb, i_mona, i_monb,
                  q_dimer, q_mona, q_monb,
                ]

        subprocess.call(['cp',template_filepath,filepath])
        for [fill,item] in zip(fill_items,items):
            subprocess.call(['sed','-i',"s/"+fill+'/'+str(item)+'/g',filepath])

        return
###########################################################################


###########################################################################
    def create_ccsdt_input_files(self,ccsdtdir,geometriesdir,\
            ccsdt_template_file='ccsdt_template.com'):
            
        '''
        '''

        print 'Creating CCSD(T) input file.'
        homedir = os.getcwd()
        input_filepath = ccsdtdir+'/'+self.mona+'_'+self.monb+'_'+ ccsdt_template_file
        
        # Copy xyz files to input directory (using xargs in case argument list
        # is long)
        subprocess.call('echo '+ geometriesdir + '/*xyz  | xargs cp -t '+ ccsdtdir, shell=True)
        os.chdir(ccsdtdir)

        self.add_midbonds(ccsdtdir)

        # Create input files corresponding to all .xyz configurations
        runscript = scriptsdir + 'create_molpro_files_from_xyz.py'
        xyz_prefix = self.mona + '_' + self.monb
        ccsdt_prefix = self.mona + '_' + self.monb + '_' + ccsdt_template_file.replace('_template.com','')

        subprocess.call([runscript, ccsdt_prefix, xyz_prefix])
        os.chdir(homedir)

        return 


###########################################################################
###########################################################################
###########################################################################


###########################################################################
######################## Command Line Arguments ###########################

# Make input directory
subprocess.call(['mkdir','-p',ccsdtdir])

# Create SAPT template file
c = CreateSAPTInputFile(dimer_info_file, inputdir, templatesdir)
c.fill_ccsdt_template(ccsdtdir)

# Copy xyz files over to the input directory, add midbond functions as
# necessary, and create input files for each .xyz configuration
c.create_ccsdt_input_files(ccsdtdir,geometriesdir)

# Copy basis set file over to ccsdt directory
basisfile=templatesdir + '/basis_sets/AVTZ.mbas'
subprocess.call(['cp',basisfile,ccsdtdir])

        

###########################################################################
###########################################################################
