#!/usr/bin/env python
"""

Last Updated:
"""

# Standard modules
import numpy as np
import sys
import os
import subprocess
from pathlib import Path
import glob
import argparse
import json
from _ctypes import PyObj_FromPtr  # see https://stackoverflow.com/a/15012814/355230
import re

# Dependencies
from chemistry import io

###########################################################################
####################### Global Variables ##################################
# Directories
source_path = Path(__file__).resolve()
scriptsdir = str(source_path.parent)
maindir = os.getcwd() + '/' # directory where user is calling function
sapt_script= scriptsdir + '/get-sapt-fitting-data'
templatesdir = maindir + '/templates/'
inputdir = maindir + '/input/'
saptdir = maindir + '/sapt/'
saptoutdir = saptdir + 'output/'
saptbaddir = saptdir + 'bad_calcs/'





###########################################################################
###########################################################################
# Pretty printing for json object:
# https://stackoverflow.com/questions/42710879/write-two-dimensional-list-to-json-file
# Not required but makes output .json more human-readable for xyz data
class NoIndent(object):
    """ Value wrapper. """
    def __init__(self, value):
        if not isinstance(value, (list, tuple)):
            raise TypeError('Only lists and tuples can be wrapped')
        self.value = value


class MyEncoder(json.JSONEncoder):
    FORMAT_SPEC = '@@{}@@'  # Unique string pattern of NoIndent object ids.
    regex = re.compile(FORMAT_SPEC.format(r'(\d+)'))  # compile(r'@@(\d+)@@')

    def __init__(self, **kwargs):
        # Keyword arguments to ignore when encoding NoIndent wrapped values.
        ignore = {'cls', 'indent'}

        # Save copy of any keyword argument values needed for use here.
        self._kwargs = {k: v for k, v in kwargs.items() if k not in ignore}
        super(MyEncoder, self).__init__(**kwargs)

    def default(self, obj):
        return (self.FORMAT_SPEC.format(id(obj)) if isinstance(obj, NoIndent)
                    else super(MyEncoder, self).default(obj))

    def iterencode(self, obj, **kwargs):
        format_spec = self.FORMAT_SPEC  # Local var to expedite access.

        # Replace any marked-up NoIndent wrapped values in the JSON repr
        # with the json.dumps() of the corresponding wrapped Python object.
        for encoded in super(MyEncoder, self).iterencode(obj, **kwargs):
            match = self.regex.search(encoded)
            if match:
                id = int(match.group(1))
                no_indent = PyObj_FromPtr(id)
                json_repr = json.dumps(no_indent.value, **self._kwargs)
                # Replace the matched id string with json formatted representation
                # of the corresponding Python object.
                encoded = encoded.replace(
                            '"{}"'.format(format_spec.format(id)), json_repr)

            yield encoded
###########################################################################

###########################################################################
def get_monomer_data(dimer_info_file):
    with open (dimer_info_file,'r') as f:
        data = [ line.split() for line in f.readlines()]
    itag = [ i[0] if i else [] for i in data ].index('MonA_Name')
    mon1 = data[itag][1]
    itag = [ i[0] if i else [] for i in data ].index('MonB_Name')
    mon2 = data[itag][1]

    return mon1, mon2
###########################################################################


###########################################################################
def flag_unfinished_calculation(ofile_lines):
    flag=['Molpro','calculation','terminated']

    if ofile_lines[-1] and ofile_lines[-1] == flag:
        return False
    else:
        return True
    ## finished=subprocess.run(['grep','-q',"Molpro calculation terminated",ofile])
    ## if finished.returncode: #flag any non-zero return codes
    ##     flag=True
###########################################################################

###########################################################################
def flag_unconverged_localization(ofile_lines,rtol=2.00, atol=1.0):
    '''Flag calculations whose localized-Hartree Fock routine (part of the
    DF-LHF subroutine in Molpro did not converge.

    Parameters
    ----------
    rtol : float, optional.
        Relative tolerance for localization errors, as a percentage of the
        total interaction energy. Default: 2.0%
    atol : float, optional.
        Absolute tolerance for localization errors in units of mH.
        Default: 1.0 mH

    Returns
    -------
    flag : bool
        True if large localization errors are found, else False.

    '''
    flag=False
    slater_error = 0
    for i,line in enumerate(ofile_lines):
        if line and line[0] == '!RKS' and line[-2] == 'Energy':
            # Check Slater Test energy appearing in the lines directly above
            # this flag (usually appears 12 lines or fewer before
            for lineA in ofile_lines[i-12:i]:
                if lineA and lineA[0:2] == ['Slater','Test']:
                    slater_error = max(slater_error,float(lineA[-1])) #max error in Ha
        elif line and line[0:2] == ['SETTING','EINT_DFTSAPT']:
            eint = abs(float(line[-2].replace('D','E'))) # units in mH

    slater_error *=1000 # convert to mH
    if slater_error/eint > rtol or slater_error > atol:
        print(slater_error, eint, slater_error/eint*100)
        #subprocess.call(['mv',ofile,saptbaddir])
        flag = True

    return flag
###########################################################################


###########################################################################
def get_sapt_data(ofile_lines):

    flag = 'SETTING' # Molpro keyword when a variable is set
    data = {}
    # Some variables set by the program, below, are not needed for later
    # fitting and can be ignored
    ignore = ['CD','CA','CB','EDM','EMA','EMB',
                    'EINT_HF','EPS_HOMO','IP_','SHIFT_','INDEX']
    for line in ofile_lines:
        if line and line[0] == flag: # check if line contains a variable stored by Molpro
            for item in ignore:
                if item in line[1]: break # don't store data from variables in the ignore list
            else:
                try:
                    data[line[1].replace('(I)','')] = float(line[3])
                except ValueError: #ignore non-numeric data
                    continue
    return data
###########################################################################


###########################################################################
def get_xyz_data(ofile_lines):

    flag = ['Geometry','for']
    found_mon1, found_mon2 = False, False
    for i,line in enumerate(ofile_lines):
        if line and line[0] == flag[0] and line[1] == flag[1]: # check if geometry block
            natoms = int(ofile_lines[i+1][0])
            if line[3] == 'A:':
                mon1_xyz = [ [float(x) for x in xyz[2:] ] 
                                for xyz in ofile_lines[i+2:i+2+natoms] ]
                mon1_elements = [ xyz[0] for xyz in ofile_lines[i+2:i+2+natoms] ]
                found_mon1 = True
            elif line[3] == 'B:':
                mon2_xyz = [ [float(x) for x in xyz[2:] ]
                                for xyz in ofile_lines[i+2:i+2+natoms] ]
                mon2_elements = [ xyz[0] for xyz in ofile_lines[i+2:i+2+natoms] ]
                found_mon2 = True
            else:
                print('ack!')
                sys.exit()
        if found_mon1 and found_mon2:
            return mon1_xyz, mon2_xyz, mon1_elements, mon2_elements
    else:
        print('Error: Geometries not located for both monomers')
        sys.exit()
###########################################################################


###########################################################################
def write_json(data,jsonfile):

    with open(jsonfile,'w',encoding='utf-8') as f:
        json.dump(data, f, indent=0, cls=MyEncoder, sort_keys=True)

    return
###########################################################################


###########################################################################
def write_sapt(data,saptfile):

    print('Not yet implemented')
    sys.exit()
###########################################################################

###########################################################################
# Make SAPT fitting data

#For each dimer in the directory, read in atomtypes from the relevant .xyz file 
def make_sapt_file():

    mon1_xyz, tmp = io.ReadCoordinates(templatesdir + mon1 + '.xyz')
    mon2_xyz, tmp = io.ReadCoordinates(templatesdir + mon2 + '.xyz')

    mon1_atoms = [i[0] for i in mon1_xyz]
    mon2_atoms = [i[0] for i in mon2_xyz]

    nmon1 = len(mon1_atoms)
    nmon2 = len(mon2_atoms)

    # Change to the relevant output directory and run the get-sapt-fitting-data
    # script
    os.chdir(saptoutdir)
    sapt_file = mon1 + '_' + mon2 + '.sapt'
    print('Compiling SAPT information (this may take a minute).')
    sapt_lines = subprocess.check_output([sapt_script],text=True)
    sapt_lines = [line.split() for line in sapt_lines.split('\n')]
    print('Writing SAPT information to file.')

    # Check that the SAPT file has the correct number of atoms for each monomer
    if nmon1 != int(sapt_lines[0][0]) or nmon2 != int(sapt_lines[nmon1+1][0]):
        print('nmon1 = ', nmon1)
        print('nmon2 = ', nmon2)
        sys.exit('Incorrect number of atoms for one of the monomers.')

    # nfiles and nlines are the number of dimer geometries and number of lines per
    # information block, respectively
    nfiles = sum([1 for i in sapt_lines if i == []]) - 1 #-1 to take care of blank line at end of file
    nlines = int(len(sapt_lines)/nfiles)

    shift_mon1 = 1
    shift_mon2 = nmon1 + 2
    for i in range(nfiles):
        for iatom,atom in enumerate(mon1_atoms):
            sapt_lines[i*nlines + iatom + shift_mon1][0] = atom
        for iatom,atom in enumerate(mon2_atoms):
            sapt_lines[i*nlines + iatom + shift_mon2][0] = atom

    # Write SAPT file with new atomtpyes to file.
    with open(sapt_file,'w') as f:
        for line in sapt_lines[:-1]:
            template = '{:20s}'*len(line) + '\n'
            f.write(template.format(*line))

    print('The SAPT output file can be found in ')
    print(saptoutdir)

#os.chdir(homedir)
###########################################################################
########################## Main Code ######################################
if __name__ == '__main__':
    ###########################################################################
    ######################## Command Line Arguments ###########################
    parser = argparse.ArgumentParser()

    outfileshelp="List of Molpro output files to work-up.\n (Default: All .out files in sapt/ subdirectory.)"
    parser.add_argument("outfiles", help=outfileshelp, nargs="*",default='OUTFILES_DEFAULT')

    jsonfileshelp="Name of .json file for summary data.\n (Default: [mon1]_[mon2]_sapt.json)"
    parser.add_argument("-o","--jsonfile", help=jsonfileshelp, default='JSONFILE_DEFAULT')

    saptfileshelp="Name of .sapt file for summary data. Usually not needed.\n (Default: [mon1]_[mon2].sapt)"
    parser.add_argument("-s","--saptfile", help=saptfileshelp, default='SAPTFILE_DEFAULT')

    formatfileshelp="Write to old .sapt file format?\n (Default: False)"
    parser.add_argument("-f","--write_oldformat", type=bool, help=formatfileshelp, default=False)

    args = parser.parse_args()

    dimer_info_file = inputdir + 'dimer_info.dat'
    mon1, mon2 = get_monomer_data(dimer_info_file)

    if args.outfiles == 'OUTFILES_DEFAULT':
        # Get list of outfiles, removing any PBS files or monomer calculations from outfiles list
        args.outfiles = [ofile for ofile in glob.glob(saptdir+"*/*.out",root_dir=maindir) 
                            if not ('pbs' in ofile or 'monomer' in ofile)]
    if args.jsonfile == 'JSONFILE_DEFAULT':
        args.jsonfile = f"{mon1}_{mon2}_sapt.json"
    if args.saptfile == 'SAPTFILE_DEFAULT':
        # Get list of outfiles, removing any PBS files or monomer calculations from outfiles list
        args.saptfile = f"{mon1}_{mon2}.sapt"
    ###########################################################################
    ###########################################################################


    all_data = {}
    first_ofile = True
    for ofile in args.outfiles:
        # Read in data from files
        with open(ofile,'r') as f:
            ofile_lines = [line.split() for line in f.readlines()]
        # Check calculation terminated successfully
        # Prune output files to exclude calculations with convergence or
        # completion issues
        if flag_unfinished_calculation(ofile_lines):
            print('Calculation for', ofile, 'did not complete successfully.')
            print('Not including this calculation in the summary data file')
            args.outfiles.remove(ofile)
            continue
        elif flag_unconverged_localization(ofile_lines):
            print('LHF(Local Hartree-Fock exchange functional) not converged for calculation:',ofile)
            print('Not including this calculation in the summary data file')
            args.outfiles.remove(ofile)
            continue

        # Scrape output files for SAPT decomposition
        data = get_sapt_data(ofile_lines)
        mon1_xyz, mon2_xyz, mon1_elements, mon2_elements = get_xyz_data(ofile_lines)
        if first_ofile:
            all_data['mon1_elements'] = NoIndent(mon1_elements)
            all_data['mon2_elements'] = NoIndent(mon2_elements)
            save_mon1_elements = mon1_elements
            save_mon2_elements = mon2_elements
            first_ofile = False
        else:
            error_string='Cannot load SAPT calculations on different dimers into the same .json file.'
            assert save_mon1_elements == mon1_elements, error_string
            assert save_mon2_elements == mon2_elements, error_string

        # Save xyz and energy data into nested dictionary
        ofile_id = os.path.basename(ofile).replace('.out','')
        all_data[ofile_id] = {}
        all_data[ofile_id]['energies'] = data
        all_data[ofile_id]['mon1_xyz'] = [NoIndent(elem) for elem in mon1_xyz]
        all_data[ofile_id]['mon2_xyz'] = [NoIndent(elem) for elem in mon2_xyz]

    # Write data to JSON (new file format) or .sapt file (deprecated file
    # format)
    if args.write_oldformat == True:
        write_sapt(all_data,maindir+args.saptfile)
    else:
        write_json(all_data,maindir+args.jsonfile)

    # sys.exit()
    while False:
        print(data)
        print( data['ETOT'] )
        etot1 = data['ELST'] + data['EXCH'] 
        etot2 = data['IND'] + data['EXIND'] + data['DISP'] + data['EXDISP']
        edhf = data['DELTA_HF']
        assert np.isclose(data['ETOT'],etot1+etot2+edhf)

        ect = data['EINT_CT']
        ect1 = data['EIND1'] - data['EIND2']
        assert np.isclose(ect,ect1)


        print(data['IND'] + data['EXIND'])
        print(data['EIND1'])
        assert np.isclose(data['IND'] + data['EXIND'],data['EIND1'])
        assert np.isclose(data['IND'] + data['EXIND'] - data['EINT_CT'],data['EIND2'])

        ofile_id = ofile.replace('.out','')

        sys.exit()




###########################################################################
###########################################################################
