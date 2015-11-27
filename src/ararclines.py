import numpy as np
import armsgs as msgs
from astropy.table import Table, Column, vstack
import os, glob, copy
import yaml
import pdb
import time

try:
    from xastropy.xutils import xdebug as xdb
except:
    pass

def parse_nist(slf,ion):
    '''Parse a NIST ASCII table.  Note that the long ---- should have
    been commented out and also the few lines at the start.

    Parameters
    ----------
    ion : str
      Name of ion
    '''
    # Root (for development only)
    if slf is None:
        root = "/".join(os.path.dirname(msgs.__file__).split("/")[:-1])
    else:
        root = slf._argflag['run']['pypitdir'] 
    # Find file
    srch_file = root + '/data/arc_lines/NIST/'+ion+'*'
    nist_file = glob.glob(srch_file)
    if len(nist_file) != 1:
        msgs.error("Cannot find NIST file {:s}".format(srch_file))
    # Read
    nist_tbl = Table.read(nist_file[0], format='ascii.fixed_width')
    gdrow = nist_tbl['Observed'] > 0.  # Eliminate dummy lines
    nist_tbl = nist_tbl[gdrow]
    # Now unique values only (no duplicates)
    uniq, indices = np.unique(nist_tbl['Observed'],return_index=True)
    nist_tbl = nist_tbl[indices]
    # Deal with Rel
    agdrel = []
    for row in nist_tbl:
        try:
            gdrel = int(row['Rel.'])
        except:
            try:
                gdrel = int(row['Rel.'][:-1])
            except: 
                gdrel = 0
        agdrel.append(gdrel) 
    agdrel = np.array(agdrel)
    # Remove and add
    nist_tbl.remove_column('Rel.')
    nist_tbl.remove_column('Ritz')
    nist_tbl.add_column(Column(agdrel,name='RelInt'))
    nist_tbl.add_column(Column([ion]*len(nist_tbl), name='Ion', dtype='S5'))
    nist_tbl.rename_column('Observed','wave')
    # Return
    return nist_tbl

def load_arcline_list(slf, idx, lines, wvmnx=None):
    '''Loads arc line list from NIST files
    Parses and rejects

    Parameters
    ----------
    idx : list  
      indices of the arc
    lines : list
      List of ions to load
    wvmnx : list or tuple
      wvmin, wvmax for line list

    Returns
    -------
    alist : Table
      Table of arc lines
    '''
    # Get the parse dict
    parse_dict = load_parse_dict()
    # Read rejection file
    if slf is None:
        root = "/".join(os.path.dirname(msgs.__file__).split("/")[:-1])
    else:
        root = slf._argflag['run']['pypitdir'] 
    with open(root+'/data/arc_lines/rejected_lines.yaml', 'r') as infile:
        rej_dict = yaml.load(infile)
    # Loop through the NIST Tables
    tbls = []
    for iline in lines:
        # Load
        tbl = parse_nist(slf,iline)
        # Parse
        if iline in parse_dict.keys():
            tbl = parse_nist_tbl(tbl,parse_dict[iline])
        # Reject
        if iline in rej_dict.keys():
            msgs.info("Rejecting select {:s} lines".format(iline))
            tbl = reject_lines(slf,tbl,idx,rej_dict[iline])
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


def reject_lines(slf,tbl,idx,rej_dict):
    '''Parses a NIST table using various criteria
    Parameters
    ----------
    tbl : Table
      Read previously from NIST ASCII file
    idx : list  
      indices of the arc
    rej_dict : dict
      Dict of rejected lines

    Returns
    -------
    tbl : Table
      Rows not rejected
    '''
    msk = tbl['wave'] == tbl['wave']
    # Loop on rejected lines
    for wave in rej_dict.keys():
        close = np.where(np.abs(wave-tbl['wave']) < 0.1)[0]
        if rej_dict[wave] == 'all':
            msk[close] = False
        elif slf == None:
            continue
        elif slf._argflag['run']['spectrograph'] in rej_dict[wave].keys():
            if rej_dict[wave][slf._argflag['run']['spectrograph']] == 'all':
                msk[close] = False
            elif slf._fitsdict["disperser"][idx[0]] in rej_dict[wave][slf._argflag['run']['spectrograph']]:
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
    gdA = tbl['Aki'] >= parse_dict['min_Aki']
    gdw = tbl['wave'] >= parse_dict['min_wave']
    # Combine
    allgd = gdI & gdA & gdw
    # Return
    return tbl[allgd]

def load_parse_dict():
    '''Dicts for parsing Arc line lists from NIST
    Rejected lines are in the rejected_lines.yaml file
    '''
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
    arcline_parse['NeI']['min_Aki']  = 1. # NOT GOOD FOR DEIMOS, DESI
    #arcline_parse['NeI']['min_wave'] = 5700. 
    arcline_parse['NeI']['min_wave'] = 5850. # NOT GOOD FOR DEIMOS?
    # ZnI
    arcline_parse['ZnI'] = copy.deepcopy(dict_parse)
    arcline_parse['ZnI']['min_intensity'] = 50.
    #
    return arcline_parse


