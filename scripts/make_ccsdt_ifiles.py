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
from xyz.add_midbond import add_midbond

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
saptdir = maindir + '/sapt/input/'
ipdir = maindir + '/sapt/monomer_calcs/'

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
    def get_ip_data(self):
        '''
        '''
        print 'Fetching IP data.'
        ofile1 = ipdir + self.mona + '_ip.out'
        ip1 = subprocess.check_output([scriptsdir + '/print_homo_ip.sh',ofile1])

        ofile2 = ipdir + self.monb + '_ip.out'
        ip2 = subprocess.check_output([scriptsdir + '/print_homo_ip.sh',ofile2])

        with open(templatesdir + 'ips.dat','w') as f:
            f.write(self.mona + '\t')
            f.write(ip1)
            f.write('\n')
            f.write(self.monb + '\t')
            f.write(ip2)
            f.write('\n')

        return
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
    def fill_sapt_template(self,dir,\
            template_file='pbe0_template.com',\
            ip_file='templates/ips.dat'):
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
        fill_ip_mona = 'FILL_IP_MONA'
        fill_ip_monb = 'FILL_IP_MONB'
        fill_homo_mona = 'FILL_HOMO_MONA'
        fill_homo_monb = 'FILL_HOMO_MONB'
        fill_items = [ fill_mona, fill_monb, fill_i_mona, fill_i_monb,
                        fill_q_dimer, fill_q_mona, fill_q_monb,
                        fill_ip_mona, fill_ip_monb, 
                        fill_homo_mona, fill_homo_monb ]

        mona = self.mona
        monb = self.monb
        i_mona = ','.join(str(i) for i in range(1,1+len(self.xyz_mona)))
        i_monb = ','.join(str(i) for i in range(1+len(self.xyz_mona),1+len(self.xyz_monb + self.xyz_mona)))

        q_mona = self.q_mona
        q_monb = self.q_monb
        q_dimer = q_mona + q_monb

        if self.ignorecase:
            info_mona = subprocess.check_output(['grep','-iw',self.mona,ip_file])
            info_monb = subprocess.check_output(['grep','-iw',self.monb,ip_file])
        else:
            info_mona = subprocess.check_output(['grep','-w',self.mona,ip_file])
            info_monb = subprocess.check_output(['grep','-w',self.monb,ip_file])

        info_mona = info_mona.split()
        info_monb = info_monb.split()

        ip_mona = float(info_mona[2]) #column 2 for dft ip's, column 3 for exp ip's
        ip_monb = float(info_monb[2])
        homo_mona = float(info_mona[1])
        homo_monb = float(info_monb[1])

        items = [ mona, monb, i_mona, i_monb,
                  q_dimer, q_mona, q_monb,
                  ip_mona, ip_monb, 
                  homo_mona, homo_monb ]

        subprocess.call(['cp',template_filepath,filepath])
        for [fill,item] in zip(fill_items,items):
            subprocess.call(['sed','-i',"s/"+fill+'/'+str(item)+'/g',filepath])

        return


###########################################################################


###########################################################################
    def create_sapt_input_files(self,saptdir,geometriesdir,\
            sapt_template_file='pbe0_template.com'):
            
        '''
        '''

        print 'Creating SAPT input file.'
        homedir = os.getcwd()
        input_filepath = saptdir+'/'+self.mona+'_'+self.monb+'_'+ sapt_template_file
        
        # Copy xyz files to input directory (using xargs in case argument list
        # is long)
        subprocess.call('echo '+ geometriesdir + '/*xyz  | xargs cp -t '+ saptdir, shell=True)
        os.chdir(saptdir)

        self.add_midbonds(saptdir)

        # Create input files corresponding to all .xyz configurations
        runscript = scriptsdir + 'create_molpro_files_from_xyz.py'
        xyz_prefix = self.mona + '_' + self.monb
        sapt_prefix = self.mona + '_' + self.monb + '_' + sapt_template_file.rstrip('_template.com')

        subprocess.call([runscript, sapt_prefix, xyz_prefix])
        os.chdir(homedir)

        return 


###########################################################################
###########################################################################
###########################################################################


###########################################################################
######################## Command Line Arguments ###########################

# Make input directory
subprocess.call(['mkdir','-p',saptdir])

# Create SAPT template file
c = CreateSAPTInputFile(dimer_info_file, inputdir, templatesdir)
c.get_ip_data()
c.fill_sapt_template(saptdir)

# Copy xyz files over to the input directory, add midbond functions as
# necessary, and create input files for each .xyz configuration
c.create_sapt_input_files(saptdir,geometriesdir)

# Copy basis set file over to sapt directory
basisfile=templatesdir + '/basis_sets/AVTZ.mbas'
subprocess.call(['cp',basisfile,saptdir])

        

###########################################################################
###########################################################################
