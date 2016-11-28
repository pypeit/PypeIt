from __future__ import (print_function, absolute_import, division, unicode_literals)

import sys
import os
import numpy as np
import yaml

from pypit import arparse as settings
from pypit import armsgs
from pypit import arsort
from pypit import arsciexp
from pypit import arutils

# Logging
msgs = armsgs.get_logger()

try:
    from xastropy.xutils import xdebug as debugger
except ImportError:
    import pdb as debugger


def SetupScience(fitsdict):
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
    # Sort the data
    msgs.bug("Files and folders should not be deleted -- there should be an option to overwrite files automatically if they already exist, or choose to rename them if necessary")
    filesort = arsort.sort_data(fitsdict, flag_unknown=bad_to_unknown)
    # Write out the details of the sorted files
    if settings.argflag['output']['sorted'] is not None:
        arsort.sort_write(fitsdict, filesort)
    # Match calibration frames to science frames
    arsort.match_science(fitsdict, filesort)
    # Make directory structure for different objects
    if do_qa:
        sci_targs = arsort.make_dirs(fitsdict, filesort)
    # Create the list of science exposures
    numsci = np.size(filesort['science'])
    sciexp = []
    for i in range(numsci):
        sciexp.append(arsciexp.ScienceExposure(i, fitsdict, do_qa=do_qa))
    # Generate setup and group dicts
    setup_dict = {}
    group_dict = {}
    for sc in range(numsci):
        scidx = sciexp[sc]._idx_sci[0]
        # Run setup
        setup = arsort.instr_setup(sciexp[sc], 1, fitsdict, setup_dict, skip_cset=True)
        # Set group_key
        group_key = setup[0]
        # Plan init
        if group_key not in group_dict.keys():
            group_dict[group_key] = {}
            for key in filesort.keys():
                if key not in ['unknown', 'dark']:
                    group_dict[group_key][key] = []
                group_dict[group_key]['sciobj'] = []
                group_dict[group_key]['stdobj'] = []
        # Run through the setups
        for kk in range(settings.spect['mosaic']['ndet']):
            _ = arsort.instr_setup(sciexp[sc], kk+1, fitsdict, setup_dict, skip_cset=True)
            # Fill group_dict too
            if kk==0:
                for key in filesort.keys():
                    if key in ['unknown', 'dark']:
                        continue
                    for idx in settings.spect[key]['index'][sc]:
                        # Only add if new
                        if fitsdict['filename'][idx] not in group_dict[group_key][key]:
                            group_dict[group_key][key].append(fitsdict['filename'][idx])
                            if key == 'standard':  # Add target name
                                group_dict[group_key]['stdobj'].append(fitsdict['target'][idx])
                        if key == 'science':  # Add target name
                            group_dict[group_key]['sciobj'].append(fitsdict['target'][scidx])
    # Write setup -- only if not present
    setup_file, nexist = arsort.get_setup_file()
    if nexist == 0:
        arsort.write_setup(setup_dict)
    elif nexist == 1: # Compare
        pass
        #prev_setup_dict = arsort.load_setup()
        #if arsort.compare_setup(setup_dict, prev_setup_dict) is False:
        #    msgs.error("Existing setup (from disk) does not match new one.  Regenerate setup file")
    # Write group file
    group_file = settings.argflag['run']['redname'].replace('.pypit', '.group')
    ydict = arutils.yamlify(group_dict)
    with open(group_file, 'w') as yamlf:
        yamlf.write( yaml.dump(ydict))#, default_flow_style=True) )
    # Finish calcheck or setup
    if settings.argflag['run']['calcheck'] or settings.argflag['run']['setup']:
        msgs.info("Inspect the setup file: {:s}".format(setup_file))
        msgs.info("Inspect the group file: {:s}".format(group_file))
        if settings.argflag['run']['calcheck']:
            msgs.info("Calibration check complete. Change 'run calcheck' flag to False to continue with data reduction")
        elif settings.argflag['run']['setup']:
            msgs.info("Setup is complete. Change 'run setup' to False to continue with data reduction")
            return 'setup'
        else:
            sys.exit()
    return sciexp


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

