#!/usr/bin/env python
"""A module that contains a table of chemical information with the following
headings:
   AtNo  Symbol     Name         Atomic Wt(g/mol)   Van der Waals     Covalent        1st Ionization 
                                                       radius(A)      Radius (A)       Potential (eV)

The module also includes the following functions to facilitate extracting
information from the chemical information table:
    AtomicNumber(element_symbol)
    Symbol(atomic_number)
    Name(element)
    NameToSymbol(element_name)
    Weight(element)
    VdWRadius(element)
    CovalentRadius(element)
    IonizationPotential(element)
    Exponent(element)
Input arguments should be of type int (for atomic_number), str (for
element_symbol), or either int or str (for element).

Last Updated: 03/23/15 by mvanvleet
"""

import numpy as np
from constants import ev2hartree

#Important Note: Our lab tends to use Be as a hack to indicate dummy/ghost atoms in
# .xyz files (we rarely use Be as an actual element). Setting
# beryllium_dummy_convention=True will turn on this behavior (default). Feel
# free to set this variable to false if you need to use Be as a real element.
beryllium_dummy_convention=True

##################### Element Data ########################################################################
## Sources: (1) http://www.chem.qmul.ac.uk/iupac/AtWt/index.html           (all data unless otherwise noted)
#           (2) J Phys. Chem A 113 (19): 5806-12; doi:10.1021/jp8111556    (for VdW radii)
#           (3) Dalton Trans., 2008, 2832-2838; doi:10.1039/b801115j       (for covalent radii)
#           (4) https://dept.astro.lsa.umich.edu/~cowley/ionen.htm         (for ionization potentials)

# A note on the covalent radii: In cases where more than one number is given for the covalent
# radii, each number is, respectively:
#           Carbon: Csp3, Csp2, Csp
#           Chromium: High spin (h.s.) Cr, low spin (l.s.) Cr
#           Manganese: h.s. Mn, l.s. Mn
#           Iron: h.s. iron, l.s. iron
# See ref. (3) above for details on how these numbers were obtained.



element_data=[\
#       AtNo  Symbol     Name         Atomic Wt(g/mol)   Van der Waals     Covalent        1st Ionization 
#                                                           radius(A)      Radius (A)       Potential (eV)
        #[ 0,   'X',   'Dummy'         ,   0.000         ,      0.00   ,     0.00          ,  0.00        ],\
        [ 0,   'Du',   'Dummy'         ,   0.000         ,      0.00   ,     0.00          ,  0.00        ],\
        [ 1,   'H',   'Hydrogen'      ,   1.008         ,      1.10   ,     0.31          ,  13.5984     ],\
        [ 2,   'He',  'Helium'        ,   4.002602      ,      1.40   ,     0.28          ,  24.5874     ],\
        [ 3,   'Li',  'Lithium'       ,   6.94          ,      1.81   ,     1.28          ,  5.3917      ],\
        [ 4,   'Be',  'Beryllium'     ,   9.012182      ,      1.53   ,     0.96          ,  9.3227      ],\
        [ 5,   'B',   'Boron'         ,  10.81          ,      1.92   ,     0.84          ,  8.298       ],\
        [ 6,   'C',   'Carbon'        ,  12.011         ,      1.70   , [0.76,0.73,0.69]  ,  11.2603     ],\
        [ 7,   'N',   'Nitrogen'      ,  14.007         ,      1.55   ,     0.71          ,  14.5341     ],\
        [ 8,   'O',   'Oxygen'        ,  15.999         ,      1.52   ,     0.66          ,  13.6181     ],\
        [ 9,   'F',   'Fluorine'      ,  18.9984032     ,      1.47   ,     0.57          ,  17.4228     ],\
        [ 10,  'Ne',  'Neon'          ,  20.1797        ,      1.54   ,     0.58          ,  21.5646     ],\
        [ 11,  'Na',  'Sodium'        ,  22.98976928    ,      2.27   ,     1.66          ,  5.1391      ],\
        [ 12,  'Mg',  'Magnesium'     ,  24.305         ,      1.73   ,     1.41          ,  7.6462      ],\
        [ 13,  'Al',  'Aluminium'     ,  26.9815386     ,      1.84   ,     1.21          ,  5.9858      ],\
        [ 14,  'Si',  'Silicon'       ,  28.085         ,      2.10   ,     1.11          ,  8.1517      ],\
        [ 15,  'P',   'Phosphorus'    ,  30.973762      ,      1.80   ,     1.07          ,  10.4867     ],\
        [ 16,  'S',   'Sulfur'        ,  32.06          ,      1.80   ,     1.05          ,  10.36       ],\
        [ 17,  'Cl',  'Chlorine'      ,  35.45          ,      1.75   ,     1.02          ,  12.9676     ],\
        [ 18,  'Ar',  'Argon'         ,  39.948         ,      1.88   ,     1.06          ,  15.7596     ],\
        [ 19,  'K',   'Potassium'     ,  39.0983        ,      2.75   ,     2.03          ,  4.3407      ],\
        [ 20,  'Ca',  'Calcium'       ,  40.078         ,      2.31   ,     1.76          ,  6.1132      ],\
        [ 21,  'Sc',  'Scandium'      ,  44.955912      ,      None   ,     1.70          ,  6.5615      ],\
        [ 22,  'Ti',  'Titanium'      ,  47.867         ,      None   ,     1.60          ,  6.8281      ],\
        [ 23,  'V',   'Vanadium'      ,  50.9415        ,      None   ,     1.53          ,  6.7462      ],\
        [ 24,  'Cr',  'Chromium'      ,  51.9961        ,      None   ,     1.39          ,  6.7665      ],\
        [ 25,  'Mn',  'Manganese'     ,  54.938045      ,      None   ,    [1.61,1.39]    ,  7.434       ],\
        [ 26,  'Fe',  'Iron'          ,  55.845         ,      None   ,    [1.52,1.32]    ,  7.9024      ],\
        [ 27,  'Co',  'Cobalt'        ,  58.933195      ,      None   ,    [1.50,1.26]    ,  7.881       ],\
        [ 28,  'Ni',  'Nickel'        ,  58.6934        ,      None   ,     1.24          ,  7.6398      ],\
        [ 29,  'Cu',  'Copper'        ,  63.546         ,      None   ,     1.32          ,  7.7264      ],\
        [ 30,  'Zn',  'Zinc'          ,  65.38          ,      None   ,     1.22          ,  9.3942      ],\
        [ 31,  'Ga',  'Gallium'       ,  69.723         ,      1.87   ,     1.22          ,  5.9993      ],\
        [ 32,  'Ge',  'Germanium'     ,  72.630         ,      2.11   ,     1.20          ,  7.8994      ],\
        [ 33,  'As',  'Arsenic'       ,  74.92160       ,      1.85   ,     1.19          ,  9.7886      ],\
        [ 34,  'Se',  'Selenium'      ,  78.96          ,      1.90   ,     1.20          ,  9.7524      ],\
        [ 35,  'Br',  'Bromine'       ,  79.904         ,      1.83   ,     1.20          ,  11.8138     ],\
        [ 36,  'Kr',  'Krypton'       ,  83.798         ,      2.02   ,     1.16          ,  13.9996     ],\
        [ 37,  'Rb',  'Rubidium'      ,  85.4678        ,      3.03   ,     2.20          ,  4.1771      ],\
        [ 38,  'Sr',  'Strontium'     ,  87.62          ,      2.49   ,     1.95          ,  5.6949      ],\
        [ 39,  'Y',   'Yttrium'       ,  88.90585       ,      None   ,     1.90          ,  6.2171      ],\
        [ 40,  'Zr',  'Zirconium'     ,  91.224         ,      None   ,     1.75          ,  6.6339      ],\
        [ 41,  'Nb',  'Niobium'       ,  92.90638       ,      None   ,     1.64          ,  6.7589      ],\
        [ 42,  'Mo',  'Molybdenum'    ,  95.96          ,      None   ,     1.54          ,  7.0924      ],\
        [ 43,  'Tc',  'Technetium'    ,  97             ,      None   ,     1.47          ,  7.28        ],\
        [ 44,  'Ru',  'Ruthenium'     , 101.07          ,      None   ,     1.46          ,  7.3605      ],\
        [ 45,  'Rh',  'Rhodium'       , 102.90550       ,      None   ,     1.42          ,  7.4589      ],\
        [ 46,  'Pd',  'Palladium'     , 106.42          ,      None   ,     1.39          ,  8.3369      ],\
        [ 47,  'Ag',  'Silver'        , 107.8682        ,      None   ,     1.45          ,  7.5762      ],\
        [ 48,  'Cd',  'Cadmium'       , 112.411         ,      None   ,     1.44          ,  8.9938      ],\
        [ 49,  'In',  'Indium'        , 114.818         ,      1.93   ,     1.42          ,  5.7864      ],\
        [ 50,  'Sn',  'Tin'           , 118.710         ,      2.17   ,     1.39          ,  7.3439      ],\
        [ 51,  'Sb',  'Antimony'      , 121.760         ,      2.06   ,     1.39          ,  8.6084      ],\
        [ 52,  'Te',  'Tellurium'     , 127.60          ,      2.06   ,     1.38          ,  9.0096      ],\
        [ 53,  'I',   'Iodine'        , 126.90447       ,      1.98   ,     1.39          ,  10.4513     ],\
        [ 54,  'Xe',  'Xenon'         , 131.293         ,      2.16   ,     1.40          ,  12.1298     ],\
        [ 55,  'Cs',  'Cesium'        , 132.9054519     ,      3.43   ,     2.44          ,  3.8939      ],\
        [ 56,  'Ba',  'Barium'        , 137.327         ,      2.68   ,     2.15          ,  5.2117      ],\
        [ 57,  'La',  'Lanthanum'     , 138.90547       ,      None   ,     2.07          ,  5.5769      ],\
        [ 58,  'Ce',  'Cerium'        , 140.116         ,      None   ,     2.04          ,  5.5387      ],\
        [ 59,  'Pr',  'Praseodymium'  , 140.90765       ,      None   ,     2.03          ,  5.473       ],\
        [ 60,  'Nd',  'Neodymium'     , 144.242         ,      None   ,     2.01          ,  5.525       ],\
        [ 61,  'Pm',  'Promethium'    , 145             ,      None   ,     1.99          ,  5.582       ],\
        [ 62,  'Sm',  'Samarium'      , 150.36          ,      None   ,     1.98          ,  5.6436      ],\
        [ 63,  'Eu',  'Europium'      , 151.964         ,      None   ,     1.98          ,  5.6704      ],\
        [ 64,  'Gd',  'Gadolinium'    , 157.25          ,      None   ,     1.96          ,  6.1501      ],\
        [ 65,  'Tb',  'Terbium'       , 158.92535       ,      None   ,     1.94          ,  5.8638      ],\
        [ 66,  'Dy',  'Dysprosium'    , 162.500         ,      None   ,     1.92          ,  5.9389      ],\
        [ 67,  'Ho',  'Holmium'       , 164.93032       ,      None   ,     1.92          ,  6.0215      ],\
        [ 68,  'Er',  'Erbium'        , 167.259         ,      None   ,     1.89          ,  6.1077      ],\
        [ 69,  'Tm',  'Thulium'       , 168.93421       ,      None   ,     1.90          ,  6.1843      ],\
        [ 70,  'Yb',  'Ytterbium'     , 173.054         ,      None   ,     1.87          ,  6.2542      ],\
        [ 71,  'Lu',  'Lutetium'      , 174.9668        ,      None   ,     1.87          ,  5.4259      ],\
        [ 72,  'Hf',  'Hafnium'       , 178.49          ,      None   ,     1.75          ,  6.8251      ],\
        [ 73,  'Ta',  'Tantalum'      , 180.94788       ,      None   ,     1.70          ,  7.5496      ],\
        [ 74,  'W',   'Tungsten'      , 183.84          ,      None   ,     1.62          ,  7.864       ],\
        [ 75,  'Re',  'Rhenium'       , 186.207         ,      None   ,     1.51          ,  7.8335      ],\
        [ 76,  'Os',  'Osmium'        , 190.23          ,      None   ,     1.44          ,  8.4382      ],\
        [ 77,  'Ir',  'Iridium'       , 192.217         ,      None   ,     1.41          ,  8.967       ],\
        [ 78,  'Pt',  'Platinum'      , 195.084         ,      None   ,     1.36          ,  8.9587      ],\
        [ 79,  'Au',  'Gold'          , 196.966569      ,      None   ,     1.36          ,  9.2255      ],\
        [ 80,  'Hg',  'Mercury'       , 200.592         ,      None   ,     1.32          ,  10.4375     ],\
        [ 81,  'Tl',  'Thallium'      , 204.38          ,      1.96   ,     1.45          ,  6.1082      ],\
        [ 82,  'Pb',  'Lead'          , 207.2           ,      2.02   ,     1.46          ,  7.4167      ],\
        [ 83,  'Bi',  'Bismuth'       , 208.98040       ,      2.07   ,     1.48          ,  7.2856      ],\
        [ 84,  'Po',  'Polonium'      , 209             ,      1.97   ,     1.40          ,  8.417       ],\
        [ 85,  'At',  'Astatine'      , 210             ,      2.02   ,     1.50          ,  9.3         ],\
        [ 86,  'Rn',  'Radon'         , 222             ,      2.20   ,     1.50          ,  10.7485     ],\
        [ 87,  'Fr',  'Francium'      , 223             ,      None   ,     2.60          ,  4.0727      ],\
        [ 88,  'Ra',  'Radium'        , 226             ,      None   ,     2.21          ,  5.2784      ],\
        [ 89,  'Ac',  'Actinium'      , 227             ,      None   ,     2.15          ,  5.17        ],\
        [ 90,  'Th',  'Thorium'       , 232.03806       ,      None   ,     2.06          ,  6.3067      ],\
        [ 91,  'Pa',  'Protactinium'  , 231.03588       ,      None   ,     2.00          ,  None        ],\
        [ 92,  'U',   'Uranium'       , 238.02891       ,      None   ,     1.96          ,  None        ] ]      


# I excluded later elements simply because they're radioactive and useless.
# Feel free to include later if necessary.                            
##         93  Np  Neptunium   [237] 
##         94  Pu  Plutonium   [244]
##         95  Am  Americium   [243]
##         96  Cm  Curium  [247] 
##         97  Bk  Berkelium   [247]
##         98  Cf  Californium [251]
##         99  Es  Einsteinium [252]
##         100 Fm  Fermium [257] 
##         101 Md  Mendelevium [258] 
##         102 No  Nobelium    [259]
##         103 Lr  Lawrencium  [262]
##         104 Rf  Rutherfordium   [267]
##         105 Db  Dubnium [270] 
##         106 Sg  Seaborgium  [271] 
##         107 Bh  Bohrium [270]
##         108 Hs  Hassium [277]
##         109 Mt  Meitnerium  [276] 
##         110 Ds  Darmstadtium    [281]
##         111 Rg  Roentgenium     [282]
##         112 Cn  Copernicium [285] 
##         113 Uut Ununtrium   [285]
##         114 Fl  Flerovium   [289]
##         115 Uup Ununpentium [2889]
##         116 Lv  Livermorium [293]
##         117 Uus Ununseptium [294]
##         118 Uuo Ununoctium  [294]

###########################################################################


############# Element Symbol to Atomic Number Dictionary ##################

# Dictionary Mapping between Element Symbol and Atomic Number 'Z'
symbol_to_atno={\
  'X':   0,\
  'H':   1,\
  'He':  2,\
  'Li':  3,\
  'Be':  4,\
  'B':   5,\
  'C':   6,\
  'N':   7,\
  'O':   8,\
  'F':   9,\
  'Ne':  10,\
  'Na':  11,\
  'Mg':  12,\
  'Al':  13,\
  'Si':  14,\
  'P':   15,\
  'S':   16,\
  'Cl':  17,\
  'Ar':  18,\
  'K':   19,\
  'Ca':  20,\
  'Sc':  21,\
  'Ti':  22,\
  'V':   23,\
  'Cr':  24,\
  'Mn':  25,\
  'Fe':  26,\
  'Co':  27,\
  'Ni':  28,\
  'Cu':  29,\
  'Zn':  30,\
  'Ga':  31,\
  'Ge':  32,\
  'As':  33,\
  'Se':  34,\
  'Br':  35,\
  'Kr':  36,\
  'Rb':  37,\
  'Sr':  38,\
  'Y':   39,\
  'Zr':  40,\
  'Nb':  41,\
  'Mo':  42,\
  'Tc':  43,\
  'Ru':  44,\
  'Rh':  45,\
  'Pd':  46,\
  'Ag':  47,\
  'Cd':  48,\
  'In':  49,\
  'Sn':  50,\
  'Sb':  51,\
  'Te':  52,\
  'I':   53,\
  'Xe':  54,\
  'Cs':  55,\
  'Ba':  56,\
  'La':  57,\
  'Ce':  58,\
  'Pr':  59,\
  'Nd':  60,\
  'Pm':  61,\
  'Sm':  62,\
  'Eu':  63,\
  'Gd':  64,\
  'Tb':  65,\
  'Dy':  66,\
  'Ho':  67,\
  'Er':  68,\
  'Tm':  69,\
  'Yb':  70,\
  'Lu':  71,\
  'Hf':  72,\
  'Ta':  73,\
  'W':   74,\
  'Re':  75,\
  'Os':  76,\
  'Ir':  77,\
  'Pt':  78,\
  'Au':  79,\
  'Hg':  80,\
  'Tl':  81,\
  'Pb':  82,\
  'Bi':  83,\
  'Po':  84,\
  'At':  85,\
  'Rn':  86,\
  'Fr':  87,\
  'Ra':  88,\
  'Ac':  89,\
  'Th':  90,\
  'Pa':  91,\
  'U':   92, }

if beryllium_dummy_convention==True:
    # Reset Be atomic weight to be zero
    element_data[4][3] = 0.0000

###########################################################################


################ Begin Functions ##########################################

def AtomicNumber(element_symbol):
    """Given an element symbol element_symbol as a string, returns the
    corresponding atomic number as an integer"""
    # Properly format element symbol:
    element_symbol=str(element_symbol).translate(None,'0123456789').capitalize()
    return symbol_to_atno[element_symbol]


def Symbol(atomic_number):
    """Given an input integer atomic_number, returns the corresponding element
    symbol as a string.
    """
    return element_data[int(atomic_number)][1]


def Name(element):
    """Given an input elment, which can either be an integer or a string
    (referring to an element's atomic number or symbol, respectively), returns
    the corresponding element name as a string.
    """
    try:
        return element_data[int(element)][2]
    except ValueError:
        return element_data[AtomicNumber(element)][2]
    

def Weight(element):
    """Returns the atomic weight (in g/mol) of a given element. Input can
    either be given as an integer (corresponding to an element's atomic
    number) or as a string (corresponding to the element's symbol).
    """
    try:
        return element_data[int(element)][3]
    except ValueError:
        return element_data[AtomicNumber(element)][3]


def VdWRadius(element):
    """Returns the Van der Waals radius (in Angstroms) of a given element.
    Input can either be given as an integer (corresponding to an element's
    atomic number) or as a string (corresponding to the element's symbol).
    """
    try:
        return element_data[int(element)][4]
    except ValueError:
        return element_data[AtomicNumber(element)][4]


def CovalentRadius(element,upper_bound=True):
    """Returns the covalent radius (in Angstroms) of a given element.
    Input can either be given as an integer (corresponding to an element's
    atomic number) or as a string (corresponding to the element's symbol). In
    cases where there exist multiple covalent radii for an element (true for
    C, Cr, Mn, and Fe), if upped_bound is set to True (default), the upper
    bound for the covalent radius is returned; else a list of all radii is
    given.
    """
    try:
        covradius = element_data[int(element)][5]
    except ValueError:
        covradius = element_data[AtomicNumber(element)][5]
    # For certain elements (C,Cr,Mn,Fe), several covalent radii exist, in which case either
    # return the upper bound for the covalent radius (if upper_bound==True) or the whole list.
    if type(covradius) == list and upper_bound == True:
        return covradius[0]
    else:
        return covradius


def IonizationPotential(element):
    """Returns the 1st ionization potential (in eV) of a given element.
    Input can either be given as an integer (corresponding to an element's
    atomic number) or as a string (corresponding to the element's symbol).
    """
    try:
        return element_data[int(element)][6]
    except ValueError:
        return element_data[AtomicNumber(element)][6]


def Exponent(element):
    """Returns an exponent (based on the 1st ionization potential of the
    element) according to the formula exponent = 2*sqrt(2*IP).  Input can
    either be given as an integer (corresponding to an element's atomic
    number) or as a string (corresponding to the element's symbol).
    """

    ip = IonizationPotential(element)*ev2hartree

    return 2*np.sqrt(2*ip)


def NameToSymbol(element_name):
    """Given an input string of the element's name, returns the corresponding element
    symbol as a string.
    """
    names = [line[2] for line in element_data]
    symbols = [line[1] for line in element_data]
    name_to_symbols = dict(zip(names, symbols))

    return name_to_symbols[element_name]


def ValenceElectrons(element):
    """Returns the number of valence electrons of a given element.
    Input can either be given as an integer (corresponding to an element's
    atomic number) or as a string (corresponding to the element's symbol).
    """
    # Ensure input is given as an integer; convert to integer if not
    if type(element)==str:
        element=AtomicNumber(element)
    if element <= 2: 
        return element
    elif element <= 18:
        return (element - 3)%8 + 1
    elif element <= 54:
        return (element - 37)%18 + 1
    else:
        raise NotImplementedError
    

######################### End Functions ###################################



######################### Testing  ########################################
if __name__=='__main__':
    print 'Atomic Number for Sr is ',AtomicNumber('Sr')
    print 'Symbol for element 90 is ',Symbol(90)
    print 'Element name for Cr is ',Name('Cr')
    print 'Symbol for Lead is ',NameToSymbol('Lead')
    print 'Atomic weight of Mg is ',Weight(12)
    print 'Atomic weight of Li is ',Weight(3)
    print 'Atomic weight of Be is ',Weight(4)
    print 'Van der Waals radius of H is ',VdWRadius(1)
    print 'Van der Waals radius of He is ',VdWRadius('He')
    print 'Covalent radius of C is ',CovalentRadius(6)
    print 'Covalent radius of Si is ',CovalentRadius('si')
    print 'All covalent radii of Fe is ',CovalentRadius('Fe',upper_bound=False)
    print 'Ionization Potential of Mg is ',IonizationPotential('Mg')
    print 'Number of Valence electrons of He is ',ValenceElectrons('He')
    print 'Number of Valence electrons of O is ',ValenceElectrons('O')
    print 'Number of Valence electrons of Zn is ',ValenceElectrons('Zn')

######################### End Testing  #####################################
