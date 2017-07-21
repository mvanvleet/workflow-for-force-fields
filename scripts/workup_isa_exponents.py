#!/usr/bin/env python
"""

Last Updated:
"""

# Standard modules
import numpy as np
import sys
import os
from scipy.signal import argrelmin, argrelmax
from scipy.optimize import curve_fit
# mvanvleet specific modules
from chemistry import io
from chemistry.elementdata import VdWRadius

###########################################################################
####################### Global Variables ##################################
error_message='''
---------------------------------------------------------------------------
Improperly formatted arguments. Proper usage is as follows:

$ 

(<...> indicates required argument, [...] indicates optional argument)
---------------------------------------------------------------------------
    '''

# Maximum absolute accepted deviation from linearity
d2dens_abs_cutoff = 0.3
# Maximum relative accepted deviation from linearity (as a percentage of twice
# the smallest exponent in the ISA basis set, i.e. the curvature that would be
# expected of a simple gaussian decay
d2dens_rel_cutoff = 0.95

# When fitting a final 'effective exponent' to the linear data, choose which
# range of r values (in bohr) to fit to:
fit_rstart = 1.0
fit_rstart1 = 5.0
fit_rend = 9.0
# or which values of the density to use as cutoffs (on a log scale)
fit_dstart = -2
fit_dend = -20
#fit_dend = -15

## # When fitting a final 'effective exponent' to the linear data, choose which
## # range of r values (in bohr) to fit to according to a relative VdWRadius
## # cutoff:
## rel_vdw_start = 0.8
## rel_vdw_end = 1.2

# Output file suffix
out_suffix = '_new_exponent_fit.dat'

# Plot file template name
plot_file = 'isa_templates/plot_exponents.plt'

# Outdir prefix
outdir_prefix = '/OUT/'


###########################################################################
###########################################################################


###########################################################################
######################## Command Line Arguments ###########################
maindir = os.getcwd().replace("/scripts",'')
templatesdir = maindir + '/templates/'
geometriesdir = maindir + '/geometries/'
indir = maindir + '/isa/'

###########################################################################
###########################################################################


###########################################################################
class Gaussian():
    '''
    '''

    def __init__(self,exponents):
        '''
        '''
        self.exponents = exponents
        return

    def fit_exp(self,r,*prefactors):
        '''
        '''

        prefactors = np.array(prefactors)

        return np.log(sum(prefactors[:,np.newaxis]*np.exp(-self.exponents[:,np.newaxis]*r**2)))
###########################################################################


###########################################################################
def fit_line(x,a,b):
    ''' '''
    return a*x + b
###########################################################################


###########################################################################
########################## Main Code ######################################
# Get mona names
dimer_info_file = 'dimer_info.dat'
with open (templatesdir + dimer_info_file,'r') as f:
    data = [ line.split() for line in f.readlines()]
itag = [ i[0] if i else [] for i in data ].index('MonA_Name')
mona = data[itag][1]
itag = [ i[0] if i else [] for i in data ].index('MonB_Name')
monb = data[itag][1]

if mona == monb:
    mons = [mona]
else:
    mons = [mona,monb]

for mon in mons:

    # Get list of atomtypes and elements from .clt file
    clt = indir + '/' + mon + '.clt'
    with open(clt,'r') as f:
        lines = [line.split() for line in f.readlines()]

    start_keys = [ i[0] if len(i) else '' for i in lines ]
    istart = start_keys.index('I.P.')
    iend = start_keys.index('End',istart)

    atomtypes = [ i for i in start_keys[istart+1:iend] ]
    atno = [float(i[1]) for i in lines[istart+1:iend] ]
    natoms = len(atomtypes)

    # Change to output directory for carrying out the calculations
    pwd = os.getcwd()
    outdir = indir + '/' + mon + outdir_prefix
    plot_file = pwd + '/' + plot_file
    os.chdir(outdir)

    # Get equations for the shape functions for each atom from outfile
    outfile = mon + '.out'
    with open(outfile,'r') as f:
        lines = [line.split() for line in f.readlines() ]
        # Index the file for the shape functions which are started by the phrase
        # 'Indx Symm'
        flag = 'Indx'
        indices = [ i for i,x in enumerate(lines) if len(x) and x[0] == flag ]
        # We only want the second half of the indices, as the first half contains
        # the shape functions pre-tail iteration convergence
        # indices = indices[natoms:]
        if len(indices) == 2*natoms:
            # We only want the second half of the indices, as the first half contains
            # the shape functions pre-tail iteration convergence
            # indices = indices[natoms:]
            indices = indices[natoms:]
        elif len(indices) == natoms:
            pass
        elif len(indices) != natoms:
            print 'Number of shape functions is not the same as the number of atoms!'
            sys.exit()

    shape_functions = [ [] for i in range(natoms) ]
    for i,index in enumerate(indices):
        count = 1
        line = lines[index + count]
        while line[0][0] != '=': # '=' marks end of the shape functions
            # Append prefactor, Gaussian exponent
            prefactor = float(line[2])*float(line[3])
            exponent = float(line[1])
            shape_functions[i].append([prefactor, exponent])
            count += 1
            line = lines[index + count]

    #shape_functions = np.array(shape_functions)
    summary_file1 = indir + mon + '.exp'
    summary_file2 = indir + mon + '.pre'

    f1 = open(summary_file1,'w')
    f2 = open(summary_file2,'w')

    print 'Writing exponents (B params) file to '
    print summary_file1
    print 'and prefactors (A params in this script, D params in JCTC 2016) file to'
    print summary_file2
    print

    for i,atom_shape in enumerate(shape_functions):
        atom_shape = np.array(atom_shape)
        f1.write('Atomtype: '+ atomtypes[i] + '\n')
        f2.write('Atomtype: '+ atomtypes[i] + '\n')
        print 'Atomtype: ', atomtypes[i]
        r = np.linspace(0,10,501)
        density = np.zeros_like(r)
        ddens = np.zeros_like(r)
        d2dens = np.zeros_like(r)
        for gauss in atom_shape:
            a = gauss[0]
            b = gauss[1]
            expbr2 = np.exp(-b*r**2)
            density += a*expbr2
            ddens += -2*b*a*r*expbr2
            d2dens += 4*a*b**2*r**2*expbr2 - 2*b*a*expbr2

        # Above we obtained expressions for the density and 1st and 2nd
        # derivatives. We're really interested in the natural log of the density,
        # though, and so ammend our expressions accordingly.
        d2dens = -(ddens)**2/(density**2) + d2dens/density
        ddens = ddens/density 
        density = np.log(density)


        # We now want to find the longest possible portion of the shape
        # function that is approximately linear, which we will do by finding the
        # longest set of data points whose second derivative is zero to within
        # some tolerance and whose density is strictly positive.
        d2dens_cutoff = min(d2dens_abs_cutoff, 2*atom_shape[-1,1]*d2dens_rel_cutoff)
        #print 'Density cutoff = ', d2dens_cutoff
        linear = np.logical_and((abs(d2dens) < d2dens_cutoff),np.logical_not(np.isnan(density)))
        # make sure all runs of ones are well-bounded
        bounded = np.hstack(([0], linear , [0]))
        # get 1 at run starts and -1 at run ends
        difs = np.diff(bounded)
        run_starts, = np.where(difs > 0)
        run_ends, = np.where(difs < 0)
        runs = run_ends - run_starts
        i1 = run_starts[np.argmax(runs)]
        i2 = run_ends[np.argmax(runs)] -1 #avoid off by one error

        tail_r1 = r1 = r[i1]
        tail_r2 = r2 = r[i2]
        y1 = density[i1]
        y2 = density[i2]
        template = 'Fitting tail-corrected region using points over the range {:16.8f} to {:16.8f}\n'
        f1.write(template.format(r1,r2))
        f2.write(template.format(r1,r2))
        #print 'Fitting tail-corrected region using points over the range ',r1,' to ',r2

        # Now that we have values for r1 and r2, we can generate values of a and b
        # to fit to y = a*exp(-b*r)
        b = (y1-y2)/(r1 - r2)
        a = (r1*y2 - r2*y1)/(r1 - r2)

        b_two_point = -b
        b_avg = -np.average(ddens[i1:i2])
        template = 'B (two-point fit): {:16.8f}\n'
        f1.write(template.format(b_two_point))
        template = 'B (average over range): {:16.8f}\n'
        f1.write(template.format(b_avg))
        # print 'B (two-point fit):', b_two_point
        # print 'B (average over range):', b_avg

        a_two_point = np.exp(a)

        template = 'A (two-point fit): {:16.8f}\n'
        f2.write(template.format(a_two_point))
        template = 'A (average over range): Not shown\n'
        f2.write(template.format())
        ## print 'A (two-point fit):', a_two_point
        ## print 'A (average over range):', 'Not shown'

        linear_density = a + b*r

        # To see how much basis set errors might be affecting the tail region,
        # refit the tail region of our new linear density using the ISA Gaussian
        # basis set.
        try:
            popt,pcov = curve_fit(Gaussian(atom_shape[:,1]).fit_exp, r[i1:i2],\
                    linear_density[i1:i2], p0=atom_shape[:,0])
            gaussian_fit_to_line = Gaussian(atom_shape[:,1]).fit_exp(r,*popt)

            popt,pcov = curve_fit(Gaussian(atom_shape[:,1]).fit_exp, r[i1:],\
                    linear_density[i1:], p0=atom_shape[:,0])
            gaussian_fit_to_line2 = Gaussian(atom_shape[:,1]).fit_exp(r,*popt)
        except (RuntimeError,TypeError):
            gaussian_fit_to_line = np.zeros_like(r)
            gaussian_fit_to_line2= np.zeros_like(r)

        # Create a 'tail-corrected' ISA density that corrects for basis set
        # effects in the tail region of the density. 
        tail_corrected_density = np.append(density[:i2],linear_density[i2:])

        # Finally, create a new linear fit to the density that minimizes RMSE
        # over a specified range of points (and thus provides an 'effective'
        # atomic exponent:
        i1 = np.argmin(abs(r-fit_rstart))
        i2 = np.argmin(abs(r-fit_rend))
        popt, pcov = curve_fit(fit_line,r[i1:i2],tail_corrected_density[i1:i2])

        linear_density2 = fit_line(r,*popt)

        b_min_rmse = -popt[0]
        #print 'B (r_abs_cutoff minimize RMSE):', b_min_rmse


        # Finally, create a new linear fit to the density that minimizes RMSE
        # over a specified range of points (and thus provides an 'effective'
        # atomic exponent:
        dend = max(fit_dend,tail_corrected_density[-1])
        i1 = np.argmin(abs(tail_corrected_density-fit_dstart))
        i2 = np.argmin(abs(tail_corrected_density-dend))
        r1 = r[i1]
        r2 = r[i2]
        template = 'Fitting effective exponents using points over the range {:16.8f} to {:16.8f}\n'
        f1.write(template.format(r1,r2))
        f2.write(template.format(r1,r2))
        #print 'Fitting effective exponents using points over the range ',r1,' to ',r2
        popt, pcov = curve_fit(fit_line,r[i1:i2],tail_corrected_density[i1:i2])

        linear_density2 = fit_line(r,*popt)

        b_min_rmse = -popt[0]
        a_min_rmse = np.exp(popt[1])
        template = 'B (density_abs_cutoff minimize RMSE): {:16.8f}\n'
        f1.write(template.format(b_min_rmse))
        template = 'A (density_abs_cutoff minimize RMSE): {:16.8f}\n'
        f2.write(template.format(a_min_rmse))
        print 'A (density_abs_cutoff minimize RMSE):', a_min_rmse
        print 'B (density_abs_cutoff minimize RMSE):', b_min_rmse


        ## # Finally, create a new linear fit to the density that minimizes RMSE
        ## # over a specified range of points (and thus provides an 'effective'
        ## # atomic exponent:
        ## vdw_rstart = rel_vdw_start*VdWRadius(atno[i])
        ## vdw_rend = rel_vdw_end*VdWRadius(atno[i])

        ## i1 = np.argmin(abs(r-vdw_rstart))
        ## i2 = np.argmin(abs(r-vdw_rend))
        ## popt, pcov = curve_fit(fit_line,r[i1:i2],tail_corrected_density[i1:i2])

        ## linear_density2 = fit_line(r,*popt)

        ## b_min_rmse = -popt[0]
        ## print 'B (vdw minimize RMSE):', b_min_rmse


        outfile = atomtypes[i] + out_suffix
        with open(outfile,'w') as f:
            f.write('Fitting densities and creating tail-weighted regioin over the range {} to {} bohr:\n'.format(tail_r1,tail_r2))
            f.write('B (two-point fit): ' + str(b_two_point) + '\n')
            f.write('B (average over range): ' + str(b_avg) + '\n')
            f.write('Fitting densities over the range {} to {} bohr:\n'.format(r1,r2))
            f.write('B (minimize RMSE): ' + str(b_min_rmse) + '\n')
            f.write('R  Density   d(Density)   d2(Density)   Basis-corrected Density  Linear Fit1   Linear Fit2\n')
            for line in zip(r, density, ddens, d2dens, tail_corrected_density, linear_density,
                    linear_density2, gaussian_fit_to_line, gaussian_fit_to_line2):
                template = '{:16.6f}'*len(line) + '\n'
                f.write(template.format(*line))

    f1.close()
    f2.close()


os.chdir(pwd)

###########################################################################
###########################################################################
