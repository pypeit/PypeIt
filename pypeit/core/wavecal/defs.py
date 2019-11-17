""" Module for arcline definitions
"""
from astropy.table import Table

# TODO: This doesn't belong here.
def instruments():
    """ Dict to convert instrument to bitwise flag
    WARNING: Modifying any of the following is a *bad* idea
      Adding is ok

    Returns
    -------
    instr_dict : Table

    """
    instr_dict = {}
    #
    instr_dict['LRISr'] = 2**0
    instr_dict['LRISb'] = 2**1
    instr_dict['Kastb'] = 2**2
    instr_dict['Kastr'] = 2**3
    instr_dict['DEIMOS'] = 2**4
    instr_dict['NIRSPEC'] = 2**5
    instr_dict['GMOS'] = 2**6

    #
    return instr_dict


def lines():
    """ Dict of lines included in this database
    WARNING: Modifying any of the following is a *bad* idea
      Adding is ok

    Returns
    -------
    lamp_dict : dict

    """
    line_dict = {}
    #
    line_dict['ArI'] = 2**0
    line_dict['HgI'] = 2**1
    line_dict['KrI'] = 2**2
    line_dict['NeI'] = 2**3
    line_dict['XeI'] = 2**4
    line_dict['CdI'] = 2**5
    line_dict['ZnI'] = 2**6
    line_dict['HeI'] = 2**7
    line_dict['OH_R24000'] = 2**8
    line_dict['OH_triplespec'] = 2**9
    line_dict['CuI'] = 2**10
    line_dict['ArII'] = 2**11
    line_dict['OH_XSHOOTER'] = 2**12
    line_dict['OH_GNIRS'] = 2**13
    line_dict['OH_NIRES'] = 2**14
    line_dict['ThAr_XSHOOTER_VIS'] = 2**15
    line_dict['OH_GMOS'] = 2**16
    line_dict['OH_MODS'] = 2**17
    line_dict['ThAr_MagE'] = 2**18  # R=4100
    line_dict['OH_FIRE_Echelle'] = 2**19  # R=6000


    #
    return line_dict


def str_len():
    """ Hard-codes length of strings in the database
    WARNING: Modifying any of the following is a *bad* idea

    Returns
    -------
    strlen_dict : dict

    """
    strlen_dict = {}
    # Length of ion name
    strlen_dict['ion'] = 6
    # Length of data file name for line source
    strlen_dict['Source'] = 30
    # Return
    return strlen_dict
