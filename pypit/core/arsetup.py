""" Routines for sorting data to be reduced by PYPIT"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import os
import glob
import string

import numpy as np
import yaml

import linetools.utils

from pypit import msgs
#from pypit import arparse as settings
from pypit import arparse
from pypit import arutils
from pypit import arparse
from pypit.core import arsort
from pypit import ardebug as debugger


def dummy_setup_dict(fitstbl, setup):
    """ Generates a dummy setup_dict

    Parameters
    ----------
    fitstbl : Table
    setup : str

    Returns
    -------
    setup_dict : dict
    """
    # setup_dict
    setup_dict = {}
    setup_dict[setup[0]] = {}
    # Fill with dummy dicts
    for ii in range(1,20): # Dummy detectors
        setup_dict[setup[0]][arparse.get_dnum(ii)] = dict(binning='1x1')
    setup_dict[setup[0]][setup[-2:]] = {}
    # Fill up filenames
    setup_dict[setup[0]][setup[-2:]]['sci'] = fitstbl['filename'][fitstbl['science']].tolist()
    # Write
    _ = write_calib(setup_dict)
    return setup_dict

'''
def dummy_setup_dict(filesort, fitsdict):
    """ Generates a dummy setup_dict

    Parameters
    ----------
    filesort : dict
    fitsdict : dict

    Returns
    -------
    setup_dict : dict
    """
    # setup_dict
    setup = settings.argflag['reduce']['masters']['setup']
    setup_dict = {}
    setup_dict[setup[0]] = {}
    # Fill with dummy dicts
    for ii in range(1,20): # Dummy detectors
        setup_dict[setup[0]][arparse.get_dnum(ii)] = dict(binning='1x1')
    setup_dict[setup[0]][setup[-2:]] = {}
    iSCI = filesort['science']
    # Fill up filenames
    setup_dict[setup[0]][setup[-2:]]['sci'] = [fitsdict['filename'][i] for i in iSCI]
    # Write
    calib_file = write_calib(setup_dict)
    return setup_dict
'''


def calib_set(isetup_dict, fitstbl, sci_ID):
    """ Generate a calibration dataset string

    Parameters
    ----------
    isetup_dict : dict
    fitstbl : Table
    sci_ID : int
      ID of the science frame

    Returns
    -------
    cb_str : str
      Ordering is 'aa', 'ab', ...

    """
    cb_strs = []
    default = 'aa'
    for ll in ['a', 'b', 'c', 'd']:
        for lower in string.ascii_lowercase:
            cb_strs.append(ll+lower)
    # Build cbset from sciexp
    new_cbset = {}
    #cbkeys = ['arcs', 'bias', 'trace', 'flat', 'sci', 'cent']
    cbkeys = ['arc', 'bias', 'trace', 'pixelflat', 'science']#, 'cent']
    for cbkey in cbkeys:
        #if len(getattr(sciexp, '_idx_'+cbkey)) != 0:
        #    nms = list(fitsdict['filename'][getattr(sciexp, '_idx_'+cbkey)])
        #    nms.sort()
        #else:
        #    nms = []
        # Grab the names
        idx = arsort.ftype_indices(fitstbl, cbkey, sci_ID)
        names = fitstbl['filename'][idx].tolist()
        # Save
        new_cbset[cbkey] = names

    # Uninitialized?
    if default not in isetup_dict.keys():
        isetup_dict[default] = new_cbset
        return default

    # Not in default
    for cbkey in cbkeys:
        def_names = np.array(isetup_dict[default][cbkey])
        if np.array_equal(def_names, new_cbset[cbkey]):
            _ = new_cbset.pop(cbkey)
    # Science only or exactly the same?
    if len(new_cbset) == 0:
        return default
    elif len(new_cbset) == 1:
        assert list(new_cbset.keys()) == ['science']
        if new_cbset['science'][0] not in isetup_dict[default]['science']:
            isetup_dict[default]['science'] += new_cbset['science']
        return default

    # New calib set?
    for cb_str in cb_strs[1:]:
        mtch = True
        if cb_str not in isetup_dict.keys():
            isetup_dict[cb_str] = new_cbset
            break
        for key in new_cbset:
            if key not in isetup_dict[cb_str]:
                mtch = False
                break
            if key in ['science']:
                continue
            ar_names = np.array(isetup_dict[default][cbkey])
            mtch &= np.array_equal(ar_names, new_cbset[key])
        if mtch:
            # Add sci frames
            for sciframe in new_cbset['science']:
                break
    # Return
    return cb_str

'''
def calib_set(isetup_dict, sciexp, fitsdict):
    """ Generate a calibration dataset string

    Parameters
    ----------
    isetup_dict : dict
    sciexp
    fitsdict : dict

    Returns
    -------
    cb_str : str
      Ordering is 'aa', 'ab', ...

    """
    cb_strs = []
    default = 'aa'
    for ll in ['a', 'b', 'c', 'd']:
        for lower in string.ascii_lowercase:
            cb_strs.append(ll+lower)
    # Build cbset from sciexp
    new_cbset = {}
    cbkeys = ['arcs', 'bias', 'trace', 'flat', 'sci', 'cent']
    for cbkey in cbkeys:
        if len(getattr(sciexp, '_idx_'+cbkey)) != 0:
            nms = list(fitsdict['filename'][getattr(sciexp, '_idx_'+cbkey)])
            nms.sort()
        else:
            nms = []
        new_cbset[cbkey] = nms

    # Uninitialized?
    if default not in isetup_dict.keys():
        isetup_dict[default] = new_cbset
        return default

    # Not in default
    for cbkey in cbkeys:
        def_names = np.array(isetup_dict[default][cbkey])
        if np.array_equal(def_names, new_cbset[cbkey]):
            _ = new_cbset.pop(cbkey)
    # Science only or exactly the same?
    if len(new_cbset) == 0:
        return default
    elif len(new_cbset) == 1:
        assert list(new_cbset.keys()) == ['sci']
        if new_cbset['sci'][0] not in isetup_dict[default]['sci']:
            isetup_dict[default]['sci'] += new_cbset['sci']
        return default

    # New calib set?
    for cb_str in cb_strs[1:]:
        mtch = True
        if cb_str not in isetup_dict.keys():
            isetup_dict[cb_str] = new_cbset
            break
        for key in new_cbset:
            if key not in isetup_dict[cb_str]:
                mtch = False
                break
            if key in ['sci']:
                continue
            ar_names = np.array(isetup_dict[default][cbkey])
            mtch &= np.array_equal(ar_names, new_cbset[key])
        if mtch:
            # Add sci frames
            for sciframe in new_cbset['sci']:
                break
    # Return
    return cb_str
'''


def det_setup(isetup_dict, ddict):
    """ Return detector setup on config dict or add to it if new

    Parameters
    ----------
    isetup_dict : dict
      Selected setup_dict
    ddict : dict
      detector dict

    Returns
    -------
    dkey : str
      Name of the detector string associated to the input ddict
      May be new or previously used

    """
    det_str = [arparse.get_dnum(i+1, prefix=False) for i in range(99)]
    # Init
    for dkey in det_str:
        mtch = True
        if dkey not in isetup_dict.keys():
            isetup_dict[dkey] = ddict
            break
        for key in isetup_dict[dkey].keys():
            mtch &= isetup_dict[dkey][key] == ddict[key]
        if mtch:
            break
    # Return
    return dkey


def instr_setup(sci_ID, det, fitstbl, setup_dict, numamplifiers,
                    must_exist=False, skip_cset=False, config_name=None):
    """ Define instrument config
    Make calls to detector and calib set

    config: A, B, C
    detector setup: 01, 02, 03, ..
    calib set: aa, ab, ac, ..


    Parameters
    ----------
    sci_ID : int
      science frame identifier (binary)
    det : int
      detector identifier
    fitstbl : Table
      contains header info
    setup_dict : dict
    numamplifiers : int
      Number of amplifiers for this detector
    skip_cset : bool, optional
      Skip calib_set;  only used when first generating instrument .setup file
    config_name : str, optional
      Can be used to choose the config value

    Returns
    -------
    setup : str
      Full setup ID, e.g. A_01_aa
    Also fills entries in setup_dict
    """
    # MasterFrame force?  If so, use user input setup updating detector
    #if settings_argflag['reduce']['masters']['force']:
    #    input_setup = settings_argflag['reduce']['masters']['setup']
    #    sdet = settings.get_dnum(det, prefix=None)
    #    setup = '{:s}_{:s}_{:s}'.format(input_setup[0],sdet,input_setup[-2:])
    #    return setup
    # Labels
    cfig_str = string.ascii_uppercase
    cstr = '--'
    # Arc index
    #idx = np.where(fitstbl['arc'] & (fitstbl['sci_idx'] & sci_idx))[0]
    idx = arsort.ftype_indices(fitstbl, 'arc', sci_ID)
    try:
        disp_name = fitstbl["dispname"][idx[0]]
    except:
        debugger.set_trace()
    try:
        disp_angle = fitstbl["dispangle"][idx[0]]
    except KeyError:
        disp_angle = 'none'
    # Common
    dichroic = fitstbl["dichroic"][idx[0]]
    decker = fitstbl["decker"][idx[0]]
    try:
        slitwid = fitstbl["slitwid"][idx[0]]
    except:
        slitwid = 'none'
    try:
        slitlen = fitstbl["slitlen"][idx[0]]
    except:
        slitlen = 'none'

    # Detector -- These may not be set properly from the header alone, e.g. LRIS
    binning = fitstbl["binning"][idx[0]]
    naxis0 = fitstbl["naxis0"][idx[0]]
    naxis1 = fitstbl["naxis1"][idx[0]]
    namp = numamplifiers

    # Generate
    # Don't nest deeper than 1
    cdict = dict(disperser={'name': disp_name,
                            'angle': disp_angle},
                 dichroic=dichroic,
                 slit={'decker': decker,
                       'slitwid': slitwid,
                       'slitlen': slitlen},
                 )
    ddict = {'binning': binning,
              'det': det,
              'namp': namp}
              #'naxis0': naxis0,
              #'naxis1': naxis1}

    def chk_key(val1, val2, tol=1e-3):
        if isinstance(val1,float):
            if np.isclose(val1,val2,rtol=tol):
                return True
            else:
                return False
        else:
            return val1 == val2

    # Configuration
    setup = None
    if len(setup_dict) == 0: #  New one, generate
        if config_name is None:
            setup = 'A'
        else:
            setup = config_name
        # Finish
        setup_dict[setup] = {}
        setup_dict[setup][cstr] = cdict
    else:  # Is it new?
        for ckey in setup_dict.keys():
            mtch = True
            for key in setup_dict[ckey][cstr].keys():
                # Dict?
                if isinstance(setup_dict[ckey][cstr][key], dict):
                    for ikey in setup_dict[ckey][cstr][key].keys():
                        mtch &= chk_key(setup_dict[ckey][cstr][key][ikey],cdict[key][ikey])
                else:
                    mtch &= chk_key(setup_dict[ckey][cstr][key], cdict[key])
            if mtch:
                setup = ckey
                break
        # Augment setup_dict?
        if setup is None:
            if must_exist:
                msgs.error("This setup is not present in the setup_dict.  Something went wrong..")
            maxs = max(setup_dict.keys())
            setup = cfig_str[cfig_str.index(maxs)+1]
            setup_dict[setup] = {}
            setup_dict[setup][cstr] = cdict

    # Detector
    dkey = det_setup(setup_dict[setup], ddict)
    # Calib set
    if not skip_cset:
        calib_key = calib_set(setup_dict[setup], fitstbl, sci_ID)
    else:
        calib_key = '--'

    # Finish and return
    setup = '{:s}_{:s}_{:s}'.format(setup, dkey, calib_key)
    return setup


'''
def instr_setup(sciexp, det, fitsdict, setup_dict, must_exist=False,
                skip_cset=False, config_name=None):
    """ Define instrument config
    Make calls to detector and calib set

    config: A, B, C
    detector setup: 01, 02, 03, ..
    calib set: aa, ab, ac, ..


    Parameters
    ----------
    sciexp : ScienceExposure
    det : int
      detector identifier
    fitsdict : dict
      contains header info
    setup_dict : dict
    skip_cset : bool, optional
      Skip calib_set;  only used when first generating instrument .setup file
    config_name : str, optional
      Can be used to choose the config value

    Returns
    -------
    setup : str
      Full setup ID, e.g. A_01_aa
    Also fills entries in setup_dict
    """
    dnum = settings.get_dnum(det)
    # MasterFrame force?  If so, use user input setup updating detector
    if settings.argflag['reduce']['masters']['force']:
        input_setup = settings.argflag['reduce']['masters']['setup']
        sdet = settings.get_dnum(det, prefix=None)
        setup = '{:s}_{:s}_{:s}'.format(input_setup[0],sdet,input_setup[-2:])
        return setup
    # Labels
    cfig_str = string.ascii_uppercase
    cstr = '--'
    # Arc
    idx = sciexp._idx_arcs
    disp_name = fitsdict["dispname"][idx[0]]
    disp_angle = fitsdict["dispangle"][idx[0]]
    # Common
    dichroic = fitsdict["dichroic"][idx[0]]
    decker = fitsdict["decker"][idx[0]]
    slitwid = fitsdict["slitwid"][idx[0]]
    slitlen = fitsdict["slitlen"][idx[0]]

    # Detector -- These may not be set properly from the header alone, e.g. LRIS
    binning = fitsdict["binning"][idx[0]]
    naxis0 = fitsdict["naxis0"][idx[0]]
    naxis1 = fitsdict["naxis1"][idx[0]]
    namp = settings.spect[dnum]["numamplifiers"]

    # Generate
    # Don't nest deeper than 1
    cdict = dict(disperser={'name': disp_name,
                            'angle': disp_angle},
                 dichroic=dichroic,
                 slit={'decker': decker,
                       'slitwid': slitwid,
                       'slitlen': slitlen},
                 )
    ddict = {'binning': binning,
              'det': det,
              'namp': namp}
              #'naxis0': naxis0,
              #'naxis1': naxis1}

    def chk_key(val1, val2, tol=1e-3):
        if isinstance(val1,float):
            if np.isclose(val1,val2,rtol=tol):
                return True
            else:
                return False
        else:
            return val1 == val2

    # Configuration
    setup = None
    if len(setup_dict) == 0: #  New one, generate
        if config_name is None:
            setup = 'A'
        else:
            setup = config_name
        # Finish
        setup_dict[setup] = {}
        setup_dict[setup][cstr] = cdict
    else:  # Is it new?
        for ckey in setup_dict.keys():
            mtch = True
            for key in setup_dict[ckey][cstr].keys():
                # Dict?
                if isinstance(setup_dict[ckey][cstr][key], dict):
                    for ikey in setup_dict[ckey][cstr][key].keys():
                        mtch &= chk_key(setup_dict[ckey][cstr][key][ikey],cdict[key][ikey])
                else:
                    mtch &= chk_key(setup_dict[ckey][cstr][key], cdict[key])
            if mtch:
                setup = ckey
                break
        # Augment setup_dict?
        if setup is None:
            if must_exist:
                msgs.error("This setup is not present in the setup_dict.  Something went wrong..")
            maxs = max(setup_dict.keys())
            setup = cfig_str[cfig_str.index(maxs)+1]
            setup_dict[setup] = {}
            setup_dict[setup][cstr] = cdict

    # Detector
    dkey = det_setup(setup_dict[setup], ddict)
    # Calib set
    if not skip_cset:
        calib_key = calib_set(setup_dict[setup], sciexp, fitsdict)
    else:
        calib_key = '--'

    # Finish and return
    setup = '{:s}_{:s}_{:s}'.format(setup, dkey, calib_key)
    return setup
'''


def get_setup_file(settings_argflag, spectrograph=None):
    """ Passes back name of setup file
    Also checks for existing setup files

    Parameters
    ----------
    spectrograph : str, optional

    Returns
    -------
    setup_file : str
      Name for the setup file
    nexist : int
      Number of existing setup files (0 or 1)

    """
    if spectrograph is None:
        spectrograph = settings_argflag['run']['spectrograph']
    setup_files = glob.glob('./{:s}*.setups'.format(spectrograph))
    nexist = len(setup_files)
    # Require 1 or 0
    if nexist == 1:
        return setup_files[0], nexist
    elif nexist == 0:
        setup_file = settings_argflag['run']['redname'].replace('.pypit', '.setups')
        if os.path.isfile(setup_file):
            nexist = 1
        #date = str(datetime.date.today().strftime('%Y-%b-%d'))
        #return '{:s}_{:s}.setups'.format(spectrograph,date), nexist
        return setup_file, nexist
    else:
        msgs.error("Found more than one .setup file in the working directory.  Limit to one.")


'''
def load_setup(**kwargs):
    """ Load setup from the disk

    Returns
    -------
    setup_dict : dict
    setup_file : str

    """
    setup_file, nexist = get_setup_file(**kwargs)
    if nexist == 0:
        debugger.set_trace()
        msgs.error("No existing setup file.  Generate one first (e.g. pypit_setup)!")
    # YAML
    with open(setup_file, 'r') as infile:
        setup_dict = yaml.load(infile)
    # Return
    return setup_dict, setup_file
'''


def write_calib(calib_file, setup_dict):
    """ Output setup_dict with full calibrations to hard drive

    Parameters
    ----------
    setup_dict : dict
      setup dict
    """
    # Write
    ydict = arutils.yamlify(setup_dict)
    with open(calib_file, 'w') as yamlf:
        yamlf.write(yaml.dump(ydict))


def write_setup(setup_dict, setup_file=None, use_json=False):
    """ Output setup_dict to hard drive

    Parameters
    ----------
    setup_dict : dict
      setup dict
    use_json : bool, optional
      Output with JSON instead of YAML (not recommended)

    Returns
    -------

    """
    # Write
    if setup_file is None:
        setup_file, nexist = get_setup_file()
    if use_json:
        gddict = linetools.utils.jsonify(setup_dict)
        linetools.utils.savejson(setup_file, gddict, easy_to_read=True)
    else: # YAML
        ydict = arutils.yamlify(setup_dict)
        with open(setup_file, 'w') as yamlf:
            yamlf.write(yaml.dump(ydict))


def load_sorted(sorted_file):
    """ Load a .sorted file (mainly to generate a .pypit file)

    Parameters
    ----------
    sorted_file : str

    Returns
    -------
    all_setups : list
     list of all setups eg. ['A','B']
    all_setuplines : list
     list of lists of all setup lines which describe the setup
    all_setupfiles : list
     list of lists of all setup files including the header
    """
    all_setups, all_setuplines, all_setupfiles = [], [], []
    try:
        with open(sorted_file,'r') as ff:
            # Should begin with ####
            fline = ff.readline()
            if fline[0:4] != '####':
                msgs.error('Bad .sorted fomatting')
            # Loop on setups
            while fline[0:5] != '##end':
                # Setup lines
                setuplines = []
                while fline[0:2] != '#-':
                    fline = ff.readline()
                    # Setup name
                    if 'Setup' in fline:
                        all_setups.append(fline[6:].strip())
                    #
                    setuplines.append(fline)
                all_setuplines.append(setuplines[:-1])
                # Data files
                datafiles = []
                while fline[0:2] != '##':
                    fline = ff.readline()
                    datafiles.append(fline)
                all_setupfiles.append(datafiles[:-1])
    except:
        debugger.set_trace()
    # Return
    return all_setups, all_setuplines, all_setupfiles


def write_sorted(group_file, fitstbl, group_dict, setup_dict):
    """ Write the .sorted file

    Parameters
    ----------
    group_file : str
    fitstbl : Table
    group_dict : dict
    setup_dict : dict

    Returns
    -------

    """
    # Setup
    srt_tbl = fitstbl.copy()
    srt_tbl['frametype'] = arsort.build_frametype_list(fitstbl)
    # Output file
    ff = open(group_file, 'w')
    # Keys
    setups = list(group_dict.keys())
    setups.sort()
    ftypes = list(group_dict[setups[0]].keys())
    ftypes.sort()
    # Loop on Setup
    asciiord = np.array(['filename', 'date', 'frameno', 'frametype',
                         'target', 'exptime', 'dispname', 'decker', 'AB_frame'])
    for setup in setups:
        ff.write('##########################################################\n')
        in_setup = []
        ff.write('Setup {:s}\n'.format(setup))
        ydict = arutils.yamlify(setup_dict[setup])
        ff.write(yaml.dump(ydict))
        ff.write('#---------------------------------------------------------\n')
        # ID files
        for key in ftypes:
            # Sort
            gfiles = group_dict[setup][key]
            gfiles.sort()
            for ifile in gfiles:
                mt = np.where(srt_tbl['filename'] == ifile)[0]
                if (len(mt) > 0) and (mt not in in_setup):
                    in_setup.append(mt[0])
        # Write overlapping keys
        gdkeys = np.in1d(asciiord, np.array(srt_tbl.keys()))
        subtbl = srt_tbl[asciiord[gdkeys].tolist()][np.array(in_setup)]
        subtbl.write(ff, format='ascii.fixed_width')
    ff.write('##end\n')
    ff.close()

def build_group_dict(fitstbl, setupIDs, all_sci_idx, all_sci_ID):
    """ Generate a group dict
    Only used for generating the .sorted output file

    Parameters
    ---------
    filetbl : Table
    setupIDs : list
    all_sci_idx : ndarray
      all the science frame indices in the fitstbl

    Returns
    -------
    group_dict : dict
    """
    ftype_keys = arsort.ftype_list + ['failures']

    group_dict = {}
    for sc,setupID in enumerate(setupIDs):
        #scidx = sciexp[sc]._idx_sci[0]
        scidx = all_sci_idx[sc]
        sci_ID = all_sci_ID[sc]
        # Set group_key
        config_key = setupID[0]
        # Plan init
        if config_key not in group_dict.keys():
            group_dict[config_key] = {}
            #for key in filesort.keys():
            for key in ftype_keys:
                if key not in ['unknown', 'dark']:
                    group_dict[config_key][key] = []
                group_dict[config_key]['sciobj'] = []
                group_dict[config_key]['stdobj'] = []
        # Fill group_dict too
        #for key in filesort.keys():
        for key in ftype_keys:
            if key in ['unknown', 'dark', 'failures']:
                continue
            #for idx in settings_spect[key]['index'][sc]:
            indices = arsort.ftype_indices(fitstbl, key, sci_ID)
            #for idx in settings_spect[key]['index'][sc]:
            for idx in indices:
                # Only add if new
                if fitstbl['filename'].data[idx] not in group_dict[config_key][key]:
                    group_dict[config_key][key].append(fitstbl['filename'].data[idx])
                    if key == 'standard':  # Add target name
                        group_dict[config_key]['stdobj'].append(fitstbl['target'].data[idx])
                if key == 'science':  # Add target name
                    group_dict[config_key]['sciobj'].append(fitstbl['target'].data[scidx])

    return group_dict

'''
def build_group_dict(filesort, setupIDs, sciexp, fitsdict):
    """ Generate a group dict
    Only used for generating the .sorted output file

    Parameters
    ---------
    filesort : dict
    setupIDs : list
    sciexp : list
    fitsdict : dict

    Returns
    -------
    group_dict : dict
    """

    group_dict = {}
    for sc,setupID in enumerate(setupIDs):
        scidx = sciexp[sc]._idx_sci[0]
        # Set group_key
        config_key = setupID[0]
        # Plan init
        if config_key not in group_dict.keys():
            group_dict[config_key] = {}
            for key in filesort.keys():
                if key not in ['unknown', 'dark']:
                    group_dict[config_key][key] = []
                group_dict[config_key]['sciobj'] = []
                group_dict[config_key]['stdobj'] = []
        # Fill group_dict too
        for key in filesort.keys():
            if key in ['unknown', 'dark', 'failures']:
                continue
            for idx in settings.spect[key]['index'][sc]:
                # Only add if new
                if fitsdict['filename'][idx] not in group_dict[config_key][key]:
                    group_dict[config_key][key].append(fitsdict['filename'][idx])
                    if key == 'standard':  # Add target name
                        group_dict[config_key]['stdobj'].append(fitsdict['target'][idx])
                if key == 'science':  # Add target name
                    group_dict[config_key]['sciobj'].append(fitsdict['target'][scidx])

    return group_dict
'''
