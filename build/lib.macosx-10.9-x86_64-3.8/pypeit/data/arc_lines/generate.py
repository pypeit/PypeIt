"""
Module for generating Arc Line lists
  Should be run where it is located (for now)
"""
import numpy as np
import os, imp, glob, pdb, gzip
import subprocess

from astropy import units as u
from astropy.units.quantity import Quantity
from astropy import constants as const
from astropy.io import fits, ascii
from astropy.table import QTable, Column, Table


# def line_data
# def mk_neon

# TODO

#
def line_data(nrows=1):
    ''' Defines the dict for arc line Data

    Parameters:
    ----------
    nrows: int, optional
      Number of rows in Table [default = 1]
    '''
    aldict = {
        'name': ' '*20,       # Name
        'wave': 0.*u.AA,      # Wavelength (Quantity) :: NIST convention (air for >200nm)
        'f':  0.,             # Oscillator strength
#        'gk': 0.,            # Degeneracy of the upper level
#        'Ej': 0./u.cm,       # Energy of lower level (relative to ground state)
#        'Ek': 0./u.cm,       # Energy of upper level (relative to ground state)
#        'Ex': 0./u.cm,       # Excitation energy (cm^-1)
#        'A': 0./u.s,         # Einstein coefficient
#        'gj': 0,             # Lower statistical weight (2J+1)
#        'gk': 0,             # Upper statistical weight (2J+1)
#        'gamma': 0./u.s,     # Sum of A
#        'nj': 0,             # Orbital level of lower state (or vibrational level)
#        'nk': 0,             # Orbital level of upper state (or vibrational level)
#        'Jj': 0.,            # Tot ang mom (z projection) of lower state (or rotation level)
#        'Jk': 0.,            # Tot ang mom (z projection) of upper state (or rotation level)
#        'el': 0,             # Electronic transition (2=Lyman (B-X), 3=Werner (C-X)) 
        'Z': 0,               # Atomic number (for atoms)
        'Am': 0,              # Mass number (often written as "A"; only used for D)
        'ion': 0,             # Ionic state (1=Neutral)
#        'mol': ' '*10,       # Molecular name (H2, HD, CO, C13O)
        'intensity': 0.,      # Intensity -- Usually taken from NIST
        'f_an': 0,            # Flag for analysis (-1=NG, 0=OK, 1=GOOD)
        'f_bl': 0,         # Flag for blending (0=UNKNOWN,1=GOOD,2=SELF,4=OTHER)
        'f_st': 0,      # Flag for brightness (0=UNKNOWN,1=POOR,2=WEAK,3=FAIR,4=GOOD)
        'Ref': ' '*50,        # References
        'Comment': ' '*50,    # Comment
        'group': 0            # Flag for grouping
        }

    # Table
    clms = []
    for key in aldict.keys():
        if type(aldict[key]) is Quantity:
            clm = Column( ([aldict[key].value]*nrows)*aldict[key].unit, name=key)
        else:
            clm = Column( [aldict[key]]*nrows, name=key)
        # Append
        clms.append(clm)
    tbl = Table(clms)
    tbl = tbl[('name','wave','f','Z','Am','ion','intensity','f_an','f_bl',
        'f_st','Ref','Comment','group')]

    return aldict, tbl


#
def mk_neon():
    '''Generate Ne line list from spec2d (DEIMOS) NIST list
        Restricting to >5500A (for now)
    Those values are in Air, converted from Vacuum as given by NIST
        Note that these differ (slightly) from the values in the spec2d NIST blue list
    '''
    outfil = 'Ne_air_linelist.dat'
    # Read spec2d file
    f = open('spec2d_lamp_NIST.dat','r')
    lines = f.readlines()
    f.close()

    # Generate table to fill
    _, ne_table = line_data(500)

    # Loop
    cnt = 0
    for line in lines:
        # Comment
        if line[0] == '#':
            continue
        # Search for Ne
        if 'Ne' not in line[22:28]:
            continue
        wave = float(line[0:11])
        if wave < 5500.:
            continue
        # Fill
        ne_table['wave'][cnt] = wave 
        ne_table['intensity'][cnt] = float(line[11:17])
        # Quality
        qual = line[17:23].strip()
        if qual == 'BLEND':
            ne_table['f_bl'][cnt] = 2
        elif qual == 'GOOD':
            ne_table['f_bl'][cnt] = 1
            ne_table['f_st'][cnt] = 4
            ne_table['f_an'][cnt] = 1
        elif qual == 'FAIR':
            ne_table['f_bl'][cnt] = 1
            ne_table['f_st'][cnt] = 3
        elif qual in ['JUNK','BAD']:
            ne_table['f_an'][cnt] = -1
        else:
            pdb.set_trace()
        # Comment
        if len(line) > 26:
            ne_table['Comment'][cnt] = line[27:].strip()
        else:
            ne_table['Comment'][cnt] = ' '
        # Last bits
        ne_table['Z'][cnt] = 10
        ne_table['Ref'][cnt] = 'NIST,spec2d'
        ne_table['name'][cnt] = 'Ne {:.2f}'.format(wave)

        # Increment
        cnt += 1

    # Cut
    ne_table = ne_table[0:cnt]

    # Write to file
    subt = ne_table[['name','wave','f_an','f_bl','f_st','intensity','Ref','Comment']]
    subt.write(outfil, delimiter='|',format='ascii.fixed_width')
    print('Wrote Ne line list: {:s}'.format(outfil))

#
def roman_to_number(val):
    '''Convert simple Roman numerals to Arabic

    Parameters:
    -------------
    val: str or unicoce
      Roman numeral for conversion
    Returns:
    ------------
    Number
    '''
    r_to_n = dict(I=1, II=2, III=3, IV=4, V=5, VI=6, 
        VII=7, VIII=8, IX=9, X=10)
    try:
        num = r_to_n[val.strip()]
    except KeyError:
        print(val)
        pdb.set_trace()

    return num



###########################################
