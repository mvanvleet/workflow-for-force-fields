Purpose
=================

Derive a first-principles, SAPT-based force field.

References
=================
In order of recency and relevance:

* VanVleet2016: 10.1021/acs.jctc.6b00209
* VanVleet2017: TBA
* McDaniel2013: 10.1021/jp3108182
* Schmidt2015: 10.1021/ar500272n
* Yu2011: 10.1021/jp204563n

Overview
=================
To generate a SAPT-based force field, the following inputs are required:
  1. Benchmark dimer energies from SAPT, computed for a variety of dimer
        configurations
  2. Long-range multipole moments, induced dipoles, and dispersion
        parameters, computed from monomer properties (and BS-ISA in particular)
  3. Short-range exponents computed from monomer properties (and BS-ISA in
        particular)
  4. Short-range pre-factors fit to dimer energies

The following scripts are designed to simplify (as much as is possible) the
workflow for force field generation. 

Method
=================
1. Generate the necessary input files upon which the scripts in step #2
depend. The following files must be manually created/edited, and can all be found in the
templates subdirectory (with an example set of input files given for the
pyridine dimer):
     1. dimer_info.dat
        -- For each monomer, list the monomer's name and the charge on the
            monomer. The appropriate file format should be clear from the
            pyridine example.
        -- In the manner described in dimer_info.dat, list all midbonds 
            that should be added between monomers. Midbonds are important for
            running accurate SAPT calcuations; see \cite{Yu2011} for details.
    2. generate_grid_settings.inp
        -- This is the input file for GenerateGridPoints, which generates the
            dimer configurations for running SAPT calculations. The input file
            is commented so as to be self-explanatory; you will need to change
            (at the very least) the 1st, 3rd, and 4th input sections based on
            the identities of the two monomers
    3. MONA_MONB.inp (where MONA and MONB are replaced by the monomer names 
        listed in dimer_info.dat)
        -- This file contains a title line (line 1), and (for each monomer)
            the number of atoms followed by a list of coordinates in .xyz 
            format. See pyridine_pyridine.inp for an example.
    4. MONA.atomtypes, MONB.atomtypes
        -- Each .atomtypes file has the format of a .xyz file, where the
            element names have been replaced by atomtypes. This file will be
            used to generate the CamCASP input files needed for ISA calculations, 
            and is also necessary for pre-processing the input files for force
            field fitting.

2. To generate all files necessary to run force field calculations, run the
following pre-processing scripts (from this main directory):

```bash
$ ./scripts/make_geometries.sh

$ ./scripts/get_global_coordinates.py

$ ./scripts/submit_ip_calcs.py
```

(wait until IP calculation is finished)

```bash
$ ./scripts/make_sapt_ifiles.py

$ ./scripts/make_isa_files.py

$ ./scripts/make_dispersion_files.py
```

3. Submit all SAPT and ISA calculations to relevant locations. At the time of
this writing, SAPT calculations should preferably be run on HCTC (Condor). ISA and 
dispersion calculations should be run on Phoenix using Camcasp 5.8. Copy all
output files back to Pople.

4. Workup the results of the SAPT and ISA calculations by running the
following post-processing scripts:

```bash
$ ./scripts/workup_sapt_energies.py

$ ./scripts/workup_dispersion_files.sh
```

(Depending on the force field, dynamic polarizabilities may need to be added
to templates/dispersion_base_constraints.index before running this script. See
Jesse McDaniel's thesis and \cite{McDaniel2013} for a full description of the
paramterization process for dispersion coefficients.)

```bash
$ ./scripts/workup_drude_files.sh
```

(Depending on the force field, static polarizabilities may need to be added
to templates/drude_base_constraints.index before running this script. See
Jesse McDaniel's thesis and \cite{McDaniel2013} for a full description of the
paramterization process for drude oscillator charges.)

```bash
$ ./scripts/workup_isa_charges.py

$ ./scripts/workup_isa_exponents.py
```

After running these scripts, you should have the SAPT energies, long-range
coefficients, and short-range exponents required to run the force fitting code
(which is needed to generate short-range pre-factors, see
\cite{VanVleet2016}). The proper running of this code is described in the
POInter documentation, see

https://git.chem.wisc.edu/schmidt/force_fields/wikis/home


Overview of Important Files
=================

* dimer_info.dat <- monomer names and midbond positions
* dispersion_template.clt
* generate_grid_settings.inp <- geometry configuration settings
* isa_template.clt
* pbe0_template.com
* pyridine.atomtypes <- change elements to atomtypes; only matters for dispersion
* pyridine_pyridine.inp <- monomer geometries

For most systems, only dimer_info.dat, the .inp files, and the .atomtypes file
will need to be changed. The examples provided for these files should hopefully make the format self-explanatory.


System Requirements
======
Python dependencies:
* numpy
* scipy
* chemistry (mvanvleet package; not standard, so this needs to be downloaded and
added to your `$PYTHONPATH`)

