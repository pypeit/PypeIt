from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
from astropy.table import Table, Column, vstack
import glob, copy
import yaml
from pkg_resources import resource_filename

from pypeit.core.wavecal import waveio

from pypeit import msgs
from pypeit import debugger


def parse_nist(ion):
    """Parse a NIST ASCII table.  Note that the long ---- should have
    been commented out and also the few lines at the start.

    Parameters
    ----------
    ion : str
      Name of ion
    """
    nist_path = resource_filename('pypeit', 'data/arc_lines/NIST/')
    line_file = nist_path+'{:s}_vacuum.ascii'.format(ion)
    nist_tbl = waveio.load_line_list(line_file, NIST=True)

    # Return
    return nist_tbl


def load_arcline_list(lines, disperser, spectrograph, wvmnx=None, modify_parse_dict=None):
    """Loads arc line list from NIST files
    Parses and rejects

    Parameters
    ----------
    lines : list
      List of ions to load
    disperser : str
      Name of the disperser
    wvmnx : list or tuple
      wvmin, wvmax for line list
    modify_parse_dict : dict, optional
      Used to over-ride default settings of parse_dict

    Returns
    -------
    alist : Table
      Table of arc lines
    """
    # Get the parse dict
    parse_dict = load_parse_dict(modify_dict=modify_parse_dict)
    root = resource_filename('pypeit', 'data/arc_lines')
    with open(root+'/rejected_lines.yaml', 'r') as infile:
        rej_dict = yaml.load(infile)

    # TODO JFH This is not elegant, but I cannot figure out how the codes reads in linelists from other directories.
    # It appears these codes were never ported properly
    if 'OH' in lines[0]:
        alist = waveio.load_line_list(lines[0],use_ion = True, NIST = False)
    else:
        # Loop through the NIST Tables
        tbls = []
        for iline in lines:
            # Load
            tbl = parse_nist(iline)
            # Parse
            if iline in parse_dict.keys():
                tbl = parse_nist_tbl(tbl,parse_dict[iline])
            # Reject
            if iline in rej_dict.keys():
                msgs.info("Rejecting select {:s} lines".format(iline))
                tbl = reject_lines(tbl,rej_dict[iline], disperser, spectrograph)
            tbls.append(tbl[['Ion','wave','RelInt']])
        # Stack
        alist = vstack(tbls)

    # wvmnx?
    if wvmnx is not None:
        msgs.info('Cutting down line list by wvmnx: {:g},{:g}'.format(wvmnx[0],wvmnx[1]))
        gdwv = (alist['wave'] >= wvmnx[0]) & (alist['wave'] <= wvmnx[1])
        alist = alist[gdwv]
    # Return
    return alist


def reject_lines(tbl, rej_dict, disperser, spectrograph):
    '''Parses a NIST table using various criteria
    Parameters
    ----------
    tbl : Table
      Read previously from NIST ASCII file
    rej_dict : dict
      Dict of rejected lines
    disperser : str
      Name of the disperser
    spectrograph : str

    Returns
    -------
    tbl : Table
      Rows not rejected
    '''
    msgs.warn("Am not sure this method does anything for real -- JXP 03-Jul-2018")
    msk = tbl['wave'] == tbl['wave']
    # Loop on rejected lines
    for wave in rej_dict.keys():
        close = np.where(np.abs(wave-tbl['wave']) < 0.1)[0]
        if rej_dict[wave] == 'all':
            msk[close] = False
        elif spectrograph in rej_dict[wave].keys():
            if rej_dict[wave][spectrograph] == 'all':
                msk[close] = False
            elif disperser in rej_dict[wave][settings.argflag['run']['spectrograph']]:
                msk[close] = False
    # Return
    return tbl[msk]


def parse_nist_tbl(tbl,parse_dict):
    '''Parses a NIST table using various criteria
    Parameters
    ----------
    tbl : Table
      Read previously from NIST ASCII file
    parse_dict : dict
      Dict of parsing criteria.  Read from load_parse_dict

    Returns
    -------
    tbl : Table
      Rows meeting the criteria
    '''
    # Parse
    gdI = tbl['RelInt'] >= parse_dict['min_intensity']
    try:
        gdA = tbl['Aki'] >= parse_dict['min_Aki']
    except TypeError:
        debugger.set_trace()
    gdw = tbl['wave'] >= parse_dict['min_wave']
    # Combine
    allgd = gdI & gdA & gdw
    # Return
    return tbl[allgd]


def load_parse_dict(modify_dict=None):
    """Dicts for parsing Arc line lists from NIST
    Rejected lines are in the rejected_lines.yaml file

    Parameters:
    modify_dict : dict, optional
      Allows

    """
    dict_parse = dict(min_intensity=0., min_Aki=0., min_wave=0.)
    arcline_parse = {} 
    # ArI
    arcline_parse['ArI'] = copy.deepcopy(dict_parse)
    arcline_parse['ArI']['min_intensity'] = 1000. # NOT PICKING UP REDDEST LINES
    # HgI
    arcline_parse['HgI'] = copy.deepcopy(dict_parse)
    arcline_parse['HgI']['min_intensity'] = 800.
    # HeI
    arcline_parse['HeI'] = copy.deepcopy(dict_parse)
    arcline_parse['HeI']['min_intensity'] = 20.
    # NeI
    arcline_parse['NeI'] = copy.deepcopy(dict_parse)
    arcline_parse['NeI']['min_intensity'] = 500.
    arcline_parse['NeI']['min_Aki']  = 1. # NOT GOOD FOR DEIMOS, DESI, ISIS
    #arcline_parse['NeI']['min_wave'] = 5700. 
    arcline_parse['NeI']['min_wave'] = 5850. # NOT GOOD FOR DEIMOS and others
    # ZnI
    arcline_parse['ZnI'] = copy.deepcopy(dict_parse)
    arcline_parse['ZnI']['min_intensity'] = 50.
    # KrI
    arcline_parse['KrI'] = copy.deepcopy(dict_parse)
    arcline_parse['KrI']['min_Aki'] = 1.  # MAY NOT BE GOOD FOR DEIMOS, DESI, far red of LRISr
    # XeI
    arcline_parse['XeI'] = copy.deepcopy(dict_parse)
    arcline_parse['XeI']['min_intensity'] = 1000.
    #
    if modify_dict is not None:
        for key in modify_dict.keys():  # dict.update doesn't nest
            for ikey in modify_dict[key].keys():
                arcline_parse[key].update(modify_dict[key])
    return arcline_parse


