""" Module for arcline definitions
"""
from pypeit.bitmask import BitMask

# TODO: This doesn't belong here.
def instruments():
    """
    Dict to convert instrument to bitwise flag

    .. warning::

        Modifying any of the following is a *bad* idea.  Adding is ok

    Returns
    -------
    instr_dict : dict

    """
    instr_dict = {}
    #
    instr_dict['LRISr'] = 2**0
    instr_dict['LRISb'] = 2**1
    instr_dict['Kastb'] = 2**2
    instr_dict['shane_kast_red'] = 2**3
    instr_dict['shane_kast_red_ret'] = 2**3
    instr_dict['DEIMOS'] = 2**4
    instr_dict['NIRSPEC'] = 2**5
    instr_dict['GMOS'] = 2**6
    instr_dict['DBSP'] = 2**7
    #
    return instr_dict


class LinesBitMask(BitMask):
    """
    Bits for arc lines
    """
    version = '1.0.0'

    def __init__(self):
        mask = dict([
            ('ArI', 'Argon I'),
            ('HgI', 'Comment'),
            ('KrI', 'Comment'),
            ('NeI', 'Comment'),
            ('XeI', 'Comment'),
            ('CdI', 'Comment'),
            ('ZnI', 'Comment'),
            ('HeI', 'Comment'),
            ('FeI', 'Comment'),
            ('FeII', 'Comment'),
            ('ThAr', 'Comment'),
            ('OH', 'Comment'),
            ('UNKNWN', 'Comment'),
            #('OH_R24000', 'Comment'),
            #('OH_triplespec', 'Comment'),
            ('CuI', 'Comment'),
            ('ArII', 'Comment'),
            ('ThI', 'Comment'),
            ('ThII', 'Comment'),
            #('OH_XSHOOTER', 'Comment'),
            #('OH_GNIRS', 'Comment'),
            #('OH_NIRES', 'Comment'),
            #('OH_GMOS', 'Comment'),
            #('OH_MODS', 'Comment'),
            #('OH_FIRE_Echelle', 'Comment'),
            #('Ar_IR_GNIRS', 'Comment'),
            ('Ar', 'This is for Ar_IR_GNIRS which should specify the real ion'),
        ])
        super(LinesBitMask, self).__init__(list(mask.keys()), descr=list(mask.values()))


def lines():
    """
    Dict of line lists included in this database

    .. warning::
        Modifying any of the following is a *bad* idea.  Adding is ok

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
    line_dict['Ar_IR_GNIRS'] = 2**20  # R=6000
    line_dict['FeI'] = 2**21
    line_dict['FeII'] = 2**22
    line_dict['UNKNWN'] = 2**23
    line_dict['Ar_IR_MOSFIRE'] = 2 ** 24
    line_dict['Ne_IR_MOSFIRE'] = 2 ** 25
    line_dict['OH_MOSFIRE_Y'] = 2 ** 26
    line_dict['OH_MOSFIRE_J'] = 2 ** 27
    line_dict['OH_MOSFIRE_H'] = 2 ** 28
    line_dict['OH_MOSFIRE_K'] = 2 ** 29
    line_dict['ThAr_XSHOOTER_UVB'] = 2**30
    #
    return line_dict


def str_len():
    """
    Hard-codes length of strings in the database

    .. warning::

        Modifying any of the following is a *bad* idea.

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
