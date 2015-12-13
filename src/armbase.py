import sys
import numpy as np
from collections import OrderedDict
import armsgs
import arsort
import arsciexp

# Logging
msgs = armsgs.get_logger()

try:
    from xastropy.xutils import xdebug as xdb
except:
    pass

def SetupScience(argflag, spect, fitsdict):
    """ Create an exposure class for every science frame
    Also links to standard star frames

    Parameters
    ----------
    argflag : dict
      Arguments and flags used for reduction
    spect : dict
      Properties of the spectrograph.
    fitsdict : dict
      Contains relevant information from fits header files

    Returns
    -------
    sciexp : list
      A list containing all science exposure classes
    """
    # Sort the data
    msgs.bug("Files and folders should not be deleted -- there should be an option to overwrite files automatically if they already exist, or choose to rename them if necessary")
    filesort = arsort.sort_data(argflag, spect, fitsdict)
    # Write out the details of the sorted files
    if argflag['out']['sorted'] is not None:
        arsort.sort_write(argflag['out']['sorted'], spect, fitsdict, filesort)
    # Match calibration frames to science frames
    spect = arsort.match_science(argflag, spect, fitsdict, filesort)
    # If the user is only debugging, then exit now
    if argflag['run']['calcheck']:
        msgs.info("Calibration check complete. Change the 'calcheck' flag to continue with data reduction")
        sys.exit()
    # Make directory structure for different objects
    sci_targs = arsort.make_dirs(argflag, fitsdict, filesort)
    # Create the list of science exposures
    numsci = np.size(filesort['science'])
    sciexp = []
    for i in xrange(numsci):
        sciexp.append(arsciexp.ScienceExposure(i, argflag, spect, fitsdict))
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
    elif ftype == "flat":
        if chktype == "trace": chkarr = sciexp[sc]._idx_trace
        elif chktype == "pixflat": chkarr = sciexp[sc]._idx_flat
        else:
            msgs.bug("I could not update frame of type {0:s} and subtype {1:s}".format(ftype, chktype))
            return
    elif ftype == "standard": chkarr = sciexp[sc]._idx_std
    else:
        msgs.bug("I could not update frame of type: {0:s}".format(ftype))
        return
    if ftype == "flat":
        # First check flats of the same type
        for i in xrange(sc+1, numsci):
            # Check if an *identical* master frame has already been produced
            if chktype == "trace": chkfarr = sciexp[i]._idx_trace
            elif chktype == "pixflat": chkfarr = sciexp[i]._idx_flat
            else:
                msgs.bug("I could not update frame of type {0:s} and subtype {1:s}".format(ftype, chktype))
                return
            if np.array_equal(chkarr, chkfarr) and sciexp[i].GetMasterFrame(chktype, det, msgs, copy=False) is None:
                msgs.info("Updating master {0:s} frame for science target {1:d}/{2:d}".format(chktype, i+1, numsci))
                sciexp[i].SetMasterFrame(sciexp[sc].GetMasterFrame(chktype, det), chktype, det)
        # Now check flats of a different type
        origtype = chktype
        if chktype == "trace": chktype = "pixflat"
        elif chktype == "pixflat": chktype = "trace"
        for i in xrange(sc, numsci):
            # Check if an *identical* master frame has already been produced
            if chktype == "trace": chkfarr = sciexp[i]._idx_trace
            elif chktype == "pixflat": chkfarr = sciexp[i]._idx_flat
            else:
                msgs.bug("I could not update frame of type {0:s} and subtype {1:s}".format(ftype, chktype))
                return
            if np.array_equal(chkarr, chkfarr) and sciexp[i].GetMasterFrame(chktype, det, copy=False) is None:
                msgs.info("Updating master {0:s} frame for science target {1:d}/{2:d}".format(chktype, i+1, numsci))
                sciexp[i].SetMasterFrame(sciexp[sc].GetMasterFrame(origtype, det), chktype, det)
    else:
        for i in xrange(sc+1, numsci):
            # Check if an *identical* master frame has already been produced
            if ftype == "arc":
                chkfarr = sciexp[i]._idx_arcs
            elif ftype == "bias": chkfarr = sciexp[i]._idx_bias
            elif ftype == "standard": chkfarr = sciexp[i]._idx_std
            else:
                msgs.bug("I could not update frame of type: {0:s}".format(ftype))
                return
            if np.array_equal(chkarr, chkfarr) and sciexp[i].GetMasterFrame(ftype, det, copy=False) is None:
                msgs.info("Updating master {0:s} frame for science target {1:d}/{2:d}".format(ftype, i+1, numsci))
                sciexp[i].SetMasterFrame(sciexp[sc].GetMasterFrame(ftype, det), ftype, det)
    return


def instconfig(slf, det, scidx, fitsdict):
    """ Returns a unique config string for the current slf

    Parameters
    ----------
    scidx: int
       Exposure index (max=9999)
    """
    config_dict = OrderedDict()
    config_dict['S'] = 'slitwid'
    config_dict['D'] = 'dichroic'
    config_dict['G'] = 'disperser'
    config_dict['T'] = 'cdangle'
    #
    config = ''
    for key in config_dict.keys():
        try:
            comp = str(fitsdict[config_dict[key]][scidx])
        except KeyError:
            comp = '0'
        #
        val = ''
        for s in comp:
            if s.isdigit():
                val = val + s
        config = config + key+'{:s}-'.format(val)
    # Binning
    try:
        binning = slf._spect['det'][det-1]['binning']
    except KeyError:
        binning = '1x1'
    val = ''
    for s in binning:
        if s.isdigit():
            val = val + s
    config = config + 'B{:s}'.format(val)
    # Return
    return config

    """
    msgs.warn("Flat indexing needs to be improved in arsort.setup")
    fidx = slf._name_flat.index(slf._mspixflat_name)
    if fidx > 9:
        msgs.error("Not ready for that many flats!")
    aidx = slf._name_flat.index(slf._mspixflat_name)
    setup = 10*(aidx+1) + fidx
    return setup
    """
