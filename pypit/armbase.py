from __future__ import (print_function, absolute_import, division, unicode_literals)

import sys
import os
import numpy as np

from pypit import arparse as settings
from pypit import armsgs
from pypit import arsetup
from pypit import arsort
from pypit import arsciexp
from pypit import arparse

# Logging
msgs = armsgs.get_logger()

from pypit import ardebug as debugger


def setup_science(fitsdict):
    """ Create an exposure class for every science frame
    Also links to standard star frames and calibrations

    Parameters
    ----------
    fitsdict : dict
      Contains relevant information from fits header files

    Returns
    -------
    sciexp : list
      A list containing all science exposure classes
    """
    # Init
    if settings.argflag['run']['calcheck'] or settings.argflag['run']['setup']:
        do_qa = False
        bad_to_unknown = True
    else:
        do_qa = True
        bad_to_unknown = False
    if settings.argflag['run']['setup']:
        skip_cset = True
    else:
        skip_cset = False

    # Sort the data
    msgs.bug("Files and folders should not be deleted -- there should be an option to overwrite files automatically if they already exist, or choose to rename them if necessary")
    filesort = arsort.sort_data(fitsdict, flag_unknown=bad_to_unknown)
    # Write out the details of the sorted files into a .lst file
    if settings.argflag['output']['sorted'] is not None:
        srt_tbl = arsort.sort_write(fitsdict, filesort)

    # Match calibration frames to science frames
    _ = arsort.match_science(fitsdict, filesort)
    # Make directory structure for different objects
    if do_qa:
        sci_targs = arsort.make_dirs(fitsdict, filesort)

    # Create the list of science exposures
    numsci = np.size(settings.spect['science']['index'])
    sciexp = []
    for i in range(numsci):
        sciexp.append(arsciexp.ScienceExposure(i, fitsdict, do_qa=do_qa))

    # Generate setup dict
    setup_dict = {}
    if settings.argflag['reduce']['masters']['force']:
        # Check that setup was input
        if len(settings.argflag['reduce']['masters']['setup']) == 0:
            msgs.error("Need to specify --  reduce masters setup  -- in your PYPIT file!")
        # setup_dict
        setup = settings.argflag['reduce']['masters']['setup']
        setup_dict = {}
        setup_dict[setup[0]] = {}
        for ii in range(1,20): # Dummy detectors
            setup_dict[setup[0]][arparse.get_dnum(ii)] = dict(binning='1x1')
        setup_dict[setup[0]][setup[-2:]] = {}
        iSCI = filesort['science']
        setup_dict[setup[0]][setup[-2:]]['sci'] = [fitsdict['filename'][i] for i in iSCI]
        # Write
        calib_file = arsetup.write_calib(setup_dict)
        return sciexp, setup_dict

    # Run through the setups to fill setup_dict
    setupIDs = []
    for sc in range(numsci):
        for kk in range(settings.spect['mosaic']['ndet']):
            # Use user-supplied setup name (useful for some book-keeping)?
            try:
                cname = settings.argflag['setup']['name']
            except KeyError:
                cname = None
            setupID = arsetup.instr_setup(sciexp[sc], kk+1, fitsdict, setup_dict, skip_cset=skip_cset, config_name=cname)
            if kk == 0: # Only save the first detector for run setup
                setupIDs.append(setupID)

    # Calib IDs
    group_dict = {}
    if settings.argflag['run']['setup']: # Collate all matching files
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
            #debugger.set_trace()
        # Write .sorted file
        if len(group_dict) > 0:
            arsetup.write_sorted(srt_tbl, group_dict, setup_dict)
        else:
            msgs.warn("No group dict entries and therefore no .sorted file")

    # Write setup -- only if not present
    setup_file, nexist = arsetup.get_setup_file()
    # Write calib file (if not in setup mode)
    if not settings.argflag['run']['setup']:
        calib_file = arsetup.write_calib(setup_dict)
    else:
        arsetup.write_setup(setup_dict)
    # Finish calcheck or setup
    if settings.argflag['run']['calcheck'] or settings.argflag['run']['setup']:
        if settings.argflag['run']['calcheck']:
            msgs.info("Inspect the .calib file: {:s}".format(calib_file))
            msgs.info("*********************************************************")
            msgs.info("Calibration check complete and successful!")
            msgs.info("Set 'run calcheck False' to continue with data reduction")
            msgs.info("*********************************************************")
            # Instrument specific (might push into a separate file)
            if settings.argflag['run']['spectrograph'] in ['keck_lris_blue']:
                if settings.argflag['reduce']['flatfield']['useframe'] in ['pixelflat']:
                    msgs.warn("We recommend a slitless flat for your instrument.")
            return 'calcheck', None
        elif settings.argflag['run']['setup']:
            for idx in filesort['failures']:
                msgs.warn("No Arc found: Skipping object {:s} with file {:s}".format(fitsdict['target'][idx],fitsdict['filename'][idx]))
            msgs.info("Setup is complete.")
            msgs.info("Inspect the .setups file: {:s}".format(setup_file))
            return 'setup', None
        else:
            msgs.error("Should not get here")
    return sciexp, setup_dict


def UpdateMasters(sciexp, sc, det, ftype=None, chktype=None):
    """ Update the master calibrations for other science targets

    If they will use an identical master frame

    Parameters
    ----------
    sciexp : list
      A list containing all science exposure classes
    sc : int
      Index of sciexp for the science exposure currently being reduced
    det : int
      detector index (starting from 1)
    ftype : str
      Describes the type of Master frame being udpated
    chktype : str
      Describes the subtype of Master frame being updated
    """
    numsci = len(sciexp)
    if ftype == "arc":
        chkarr = sciexp[sc]._idx_arcs
    elif ftype == "bias": chkarr = sciexp[sc]._idx_bias
    elif ftype == "readnoise": chkarr = sciexp[sc]._idx_rn
    elif ftype == "flat":
        if chktype == "trace": chkarr = sciexp[sc]._idx_trace
        elif chktype == "pixelflat": chkarr = sciexp[sc]._idx_flat
        elif chktype == "pinhole": chkarr = sciexp[sc]._idx_cent
        else:
            msgs.bug("I could not update frame of type {0:s} and subtype {1:s}".format(ftype, chktype))
            return
    elif ftype == "standard": chkarr = sciexp[sc]._idx_std
    else:
        msgs.bug("I could not update frame of type: {0:s}".format(ftype))
        return
    if ftype == "flat":
        # First check flats of the same type
        for i in range(sc+1, numsci):
            # Check if an *identical* master frame has already been produced
            if chktype == "trace": chkfarr = sciexp[i]._idx_trace
            elif chktype == "pixelflat": chkfarr = sciexp[i]._idx_flat
            elif chktype == "pinhole": chkfarr = sciexp[i]._idx_cent
            else:
                msgs.bug("I could not update frame of type {0:s} and subtype {1:s}".format(ftype, chktype))
                return
            if np.array_equal(chkarr, chkfarr) and sciexp[i].GetMasterFrame(chktype, det, mkcopy=False) is None:
                msgs.info("Updating master {0:s} frame for science target {1:d}/{2:d}".format(chktype, i+1, numsci))
                sciexp[i].SetMasterFrame(sciexp[sc].GetMasterFrame(chktype, det), chktype, det)
        # Now check flats of a different type
        origtype = chktype
        if chktype == "trace": chktype = "pixelflat"
        elif chktype == "pixelflat": chktype = "trace"
        for i in range(sc, numsci):
            # Check if an *identical* master frame has already been produced
            if chktype == "trace": chkfarr = sciexp[i]._idx_trace
            elif chktype == "pixelflat": chkfarr = sciexp[i]._idx_flat
            elif chktype == "pinhole": chkfarr = sciexp[i]._idx_cent
            else:
                msgs.bug("I could not update frame of type {0:s} and subtype {1:s}".format(ftype, chktype))
                return
            if np.array_equal(chkarr, chkfarr) and sciexp[i].GetMasterFrame(chktype, det, mkcopy=False) is None:
                msgs.info("Updating master {0:s} frame for science target {1:d}/{2:d}".format(chktype, i+1, numsci))
                sciexp[i].SetMasterFrame(sciexp[sc].GetMasterFrame(origtype, det), chktype, det)
    else:
        for i in range(sc+1, numsci):
            # Check if an *identical* master frame has already been produced
            if ftype == "arc":
                chkfarr = sciexp[i]._idx_arcs
            elif ftype == "bias": chkfarr = sciexp[i]._idx_bias
            elif ftype == "standard": chkfarr = sciexp[i]._idx_std
            else:
                msgs.bug("I could not update frame of type: {0:s}".format(ftype))
                return
            if np.array_equal(chkarr, chkfarr) and sciexp[i].GetMasterFrame(ftype, det, mkcopy=False) is None:
                msgs.info("Updating master {0:s} frame for science target {1:d}/{2:d}".format(ftype, i+1, numsci))
                sciexp[i].SetMasterFrame(sciexp[sc].GetMasterFrame(ftype, det), ftype, det)
    return

