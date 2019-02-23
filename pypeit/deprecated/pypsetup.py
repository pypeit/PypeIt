"""
Core IO methods for the PypeIt setup.

THIS MODULE IS NOW DEPRECATED AND LIKELY TO BE REMOVED
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import os
import glob
import string

import numpy as np
import yaml

import linetools.utils

from pypeit import msgs
from pypeit import utils
from pypeit.core import parse
from pypeit.core import framematch
from pypeit import debugger


def dummy_setup_dict(file_list, setup):
    """
    Generates a dummy setup_dict.

    .. todo::
        - Describe how this is used.
        - Improve the description of setup and/or point to the function
          used to generate it.

    Args:
        file_list (:obj:`list`):
            List of files to set as science exposures.
        setup (:obj:`str`):
            String representation of the instrument setup.
    
    Returns:
        dict: Dictionary with the instrument setup.
    """
    # setup_dict
    setup_dict = {}
    setup_dict[setup[0]] = {}
    # Fill with dummy dicts
    for ii in range(1,20): # Dummy detectors
        setup_dict[setup[0]][parse.get_dnum(ii)] = dict(binning='1x1')
    setup_dict[setup[0]][setup[-2:]] = {}
    # Fill up filenames
    setup_dict[setup[0]][setup[-2:]]['sci'] = file_list
    # Write
    # TODO: Why is this called here?
    _ = write_calib(setup_dict)
    return setup_dict


def calib_set(isetup_dict, fitstbl, sci_ID):
    """ Generate a calibration dataset string

    Parameters
    ----------
    isetup_dict : dict
    fitstbl : PypeItMetaData
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
    cbkeys = ['arc', 'bias', 'trace', 'pixelflat', 'science']
    for cbkey in cbkeys:
        new_cbset[cbkey], _ = fitstbl.find_frame_files(cbkey, sci_ID=sci_ID)
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
    det_str = [parse.get_dnum(i+1, prefix=False) for i in range(99)]
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


def is_equal(val1, val2, tol=1e-3):
    if isinstance(val1,float):
        return np.isclose(val1,val2,rtol=tol)
    return val1 == val2


def instr_setup(sci_ID, det, fitstbl, setup_dict=None, must_exist=False, skip_cset=False,
                config_name=None, copy=False):
    """
    DEPRECATED

    Define the instrument configuration.

    .. todo::
        - This needs to be made a class object and pulled out of core.

    configuration ID: A, B, C
    detector number: 01, 02, 03, ..
    calibration ID: aa, ab, ac, ..

    Args:
        sci_ID (:obj:`int`):
            The selected science frame (binary)
        det (:obj:`int`):
            The 1-indexed detector
        fitstbl (:class:`pypeit.metadata.PypeItMetaData`):
            The fits file metadata used by PypeIt
        setup_dict (:obj:`dict`, optional):
            The dictionary with the instrument configurations that have
            already been parsed for this execution of PypeIt.  If None,
            the dictionary is instantiated from scratch; otherwise, any
            new setup information is added to the output dictionary.
        must_exist (:obj:`bool`, optional);
            The setup must exist in the provided `setup_dict`.  If not,
            the function will raise an error.
        skip_cset (:obj:`bool`, optional):
            Skip setting the calibration identifier.  This should only
            be True when first generating instrument .setup file.
        config_name (:obj:`str`, optional):
            Can be used to choose a specific configuration ID.
        copy (:obj:`bool`, optional):
            Do not perform any in-place additions to the setup
            dictionary.  Instead copy any input dictionary and return
            the modified dictionary.  By default, modifications are made
            to the input dictionary, which is returned instead of a
            modified copy.

    Returns:
        str, dict: Returns the string identifier for the instrument
        configuration and a dictionary with the configuration data.  The
        latter is either a new object, a modified copy of the input
        dictionary, or the input dictionary after it has been modified
        in place.
    """
    debugger.set_trace()

    # Labels
    cfig_str = string.ascii_uppercase
    cstr = '--'

    # Find the first arc exposure tied to this science exposure
    idx = np.where(fitstbl.find_frames('arc', sci_ID=sci_ID))[0][0]

    # Use this exposure to pull the relevant information for the instrument setup
    dispname = fitstbl['dispname'][idx] if 'dispname' in fitstbl.keys() else 'none'
    dispangle = fitstbl['dispangle'][idx] if 'dispangle' in fitstbl.keys() else 'none'
    dichroic = fitstbl['dichroic'][idx] if 'dichroic' in fitstbl.keys() else 'none'
    decker = fitstbl['decker'][idx] if 'decker' in fitstbl.keys() else 'none'
    slitwid = fitstbl['slitwid'][idx] if 'slitwid' in fitstbl.keys() else 'none'
    slitlen = fitstbl['slitlen'][idx] if 'slitlen' in fitstbl.keys() else 'none'
    binning = fitstbl['binning'][idx] if 'binning' in fitstbl.keys() else 'none'

    # Generate the relevant dictionaries
    cdict = dict(disperser={'name':dispname, 'angle':dispangle},
                 dichroic=dichroic,
                 slit={'decker':decker, 'slitwid':slitwid, 'slitlen':slitlen})
    ddict = {'binning':binning,
             'det':det,
             'namp':fitstbl.spectrograph.detector[det-1]['numamplifiers']}

    # Configuration
    setupID = None
    if setup_dict is None:
        # Generate new configuration dictionary
        setupID = 'A' if config_name is None else config_name
        _setup_dict = dict()
        _setup_dict[setupID] = {}
        _setup_dict[setupID][cstr] = cdict
    else:
        # Try to find the setup in the existing configuration dictionary
        for ckey in setup_dict.keys():
            mtch = True
            for key in setup_dict[ckey][cstr].keys():
                # Dict?
                if isinstance(setup_dict[ckey][cstr][key], dict):
                    for ikey in setup_dict[ckey][cstr][key].keys():
                        mtch &= is_equal(setup_dict[ckey][cstr][key][ikey],cdict[key][ikey])
                else:
                    mtch &= is_equal(setup_dict[ckey][cstr][key], cdict[key])
            if mtch:
                setupID = ckey
                break

        # Augment setup_dict?
        _setup_dict = setup_dict.copy() if copy else setup_dict

        if setupID is None:
            if must_exist:
                msgs.error('This setup ID is not present in the setup_dict.')
            maxs = max(_setup_dict.keys())
            setupID = cfig_str[cfig_str.index(maxs)+1]
            _setup_dict[setupID] = {}
            _setup_dict[setupID][cstr] = cdict

    # Detector
    dkey = det_setup(_setup_dict[setupID], ddict)
    # Calib set
    if not skip_cset:
        calib_key = calib_set(_setup_dict[setupID], fitstbl, sci_ID)
    else:
        calib_key = '--'

    # Finish and return
    return '{:s}_{:s}_{:s}'.format(setupID, dkey, calib_key), _setup_dict


# TODO: This is out of date!
#def get_setup_file(settings_argflag, spectrograph=None):
#    """ Passes back name of setup file
#    Also checks for existing setup files
#
#    Parameters
#    ----------
#    spectrograph : str, optional
#
#    Returns
#    -------
#    setup_file : str
#      Name for the setup file
#    nexist : int
#      Number of existing setup files (0 or 1)
#
#    """
#    if spectrograph is None:
#        spectrograph = settings_argflag['run']['spectrograph']
#    setup_files = glob.glob('./{:s}*.setups'.format(spectrograph))
#    nexist = len(setup_files)
#    # Require 1 or 0
#    if nexist == 1:
#        return setup_files[0], nexist
#    elif nexist == 0:
#        setup_file = settings_argflag['run']['redname'].replace('.pypeit', '.setups')
#        if os.path.isfile(setup_file):
#            nexist = 1
#        #date = str(datetime.date.today().strftime('%Y-%b-%d'))
#        #return '{:s}_{:s}.setups'.format(spectrograph,date), nexist
#        return setup_file, nexist
#    else:
#        msgs.error("Found more than one .setup file in the working directory.  Limit to one.")


def write_calib(calib_file, setup_dict):
    """ Output setup_dict with full calibrations to hard drive

    Parameters
    ----------
    setup_dict : dict
      setup dict
    """
    # Write
    ydict = utils.yamlify(setup_dict)
    with open(calib_file, 'w') as yamlf:
        yamlf.write(yaml.dump(ydict))


def write_setup(setup_dict, ofile, use_json=False):
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
#    if setup_file is None:
#        setup_file, nexist = get_setup_file()
    if use_json:
        gddict = linetools.utils.jsonify(setup_dict)
        linetools.utils.savejson(setup_file, gddict, easy_to_read=True)
        return

    ydict = utils.yamlify(setup_dict)
    with open(ofile, 'w') as yamlf:
        yamlf.write(yaml.dump(ydict))


def load_sorted(sorted_file):
    """ Load a .sorted file (mainly to generate a .pypeit file)

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
        with open(sorted_file, 'r') as ff:
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
    fitstbl : PypeItMetaData
    group_dict : dict
    setup_dict : dict

    Returns
    -------

    """
    # TODO: Not sure we should do this.  Instead just check that
    # frametype is in the relevant keys
    if 'frametype' not in fitstbl.keys():
        fitstbl.get_frame_types()
    # Output file
    ff = open(group_file, 'w')
    # Keys
    setups = list(group_dict.keys())
    setups.sort()
    ftypes = list(group_dict[setups[0]].keys())
    ftypes.sort()
    # Loop on Setup
    asciiord = np.array(fitstbl.spectrograph.metadata_keys())
    for setup in setups:
        ff.write('##########################################################\n')
        in_setup = []
        ff.write('Setup {:s}\n'.format(setup))
        ydict = utils.yamlify(setup_dict[setup])
        ff.write(yaml.dump(ydict))
        ff.write('#---------------------------------------------------------\n')
        # ID files
        for key in ftypes:
            # Sort
            gfiles = group_dict[setup][key]
            gfiles.sort()
            for ifile in gfiles:
                mt = np.where(fitstbl['filename'] == ifile)[0]
                if (len(mt) > 0) and (mt not in in_setup):
                    in_setup.append(mt[0])
        # Write overlapping keys
        gdkeys = np.in1d(asciiord, np.array(fitstbl.keys()))
        subtbl = fitstbl[asciiord[gdkeys].tolist()][np.array(in_setup)]
        subtbl.write(ff, format='ascii.fixed_width')
    ff.write('##end\n')
    ff.close()


def build_group_dict(fitstbl, setupIDs, all_sci_idx, all_sci_ID):
    """
    Build a dictionary with the list of exposures of each type for a
    given instrumental setup.

    Args:
        fitstbl (:class:`pypeit.metadata.PypeItMetaData`):
            The metadata table for the fits files to reduce.
        setupIDs (:obj:`list`):
            The list of setups.
        all_sci_idx (:obj:`list`):
            The indices of the science frames in the data table.
            TODO: Why is this needed?
        all_sci_ID (:obj:`list`):
            The ID number assigned to each science frame.

    Returns:
        dict: A dictionary with the list of file names associated with
        each instrument configuration.  It also provides the list of
        science and standard star targets. TODO: How are the latter
        used, if at all?
    """
    type_keys = framematch.FrameTypeBitMask().keys() + ['None']

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
            for key in type_keys:
#                if key not in ['None', 'dark']:
#                    group_dict[config_key][key] = []
                group_dict[config_key][key] = []
                group_dict[config_key]['sciobj'] = []
                group_dict[config_key]['stdobj'] = []
        # Fill group_dict too
        for key in type_keys:
#            if key in ['None', 'dark']:
#                continue
            indices = np.where(fitstbl.find_frames(key))[0] if key in ['None', 'dark'] \
                        else np.where(fitstbl.find_frames(key, sci_ID=sci_ID))[0]
            for idx in indices:
                # Only add if new
                if fitstbl['filename'][idx] not in group_dict[config_key][key]:
                    group_dict[config_key][key].append(fitstbl['filename'][idx])
                    # TODO: How is this used?
                    if key == 'standard' and 'target' in fitstbl.keys():  # Add target name
                        group_dict[config_key]['stdobj'].append(fitstbl['target'][idx])
                # TODO: How is this used?
                if key == 'science' and 'target' in fitstbl.keys():  # Add target name
                    group_dict[config_key]['sciobj'].append(fitstbl['target'][scidx])

    return group_dict

