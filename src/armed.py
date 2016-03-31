import numpy as np
import armasters
import armbase
import artrace
import armsgs
import arproc

try:
    from xastropy.xutils.xdebug import set_trace
#    from xastropy.xutils import xdebug as xdb
except ImportError:
    from pdb import set_trace

# Logging
msgs = armsgs.get_logger()

def ARMED(argflag, spect, fitsdict, reuseMaster=False):
    """
    Automatic Reduction and Modeling of Echelle Data

    Parameters
    ----------
    argflag : dict
      Arguments and flags used for reduction
    spect : dict
      Properties of the spectrograph.
    fitsdict : dict
      Contains relevant information from fits header files
    msgs : class
      Messages class used to log data reduction process
    reuseMaster : bool
      If True, a master frame that will be used for another science frame
      will not be regenerated after it is first made.
      This setting comes with a price, and if a large number of science frames are
      being generated, it may be more efficient to simply regenerate the master
      calibrations on the fly.

    Returns
    -------
    status : int
      Status of the reduction procedure
      0 = Successful execution
      1 = ...
    """
    status = 0

    # Create a list of science exposure classes
    sciexp = armbase.SetupScience(argflag, spect, fitsdict)
    numsci = len(sciexp)

    # Create a list of master calibration frames
    masters = armasters.MasterFrames(spect['mosaic']['ndet'])

    # Start reducing the data
    for sc in range(numsci):
        slf = sciexp[sc]
        scidx = slf._idx_sci[0]
        msgs.info("Reducing file {0:s}, target {1:s}".format(fitsdict['filename'][scidx], slf._target_name))
        # Loop on Detectors
        for kk in xrange(slf._spect['mosaic']['ndet']):
            det = kk + 1  # Detectors indexed from 1
            ###############
            # Get amplifier sections
            fitsdict = arproc.get_ampsec_trimmed(slf, fitsdict, det, scidx)
            ###############
            # Generate master bias frame
            update = slf.MasterBias(fitsdict, det)
            if update and reuseMaster:
                armbase.UpdateMasters(sciexp, sc, det, ftype="bias")
            ###############
            # Generate a bad pixel mask (should not repeat)
            update = slf.BadPixelMask(det)
            if update and reuseMaster:
                armbase.UpdateMasters(sciexp, sc, det, ftype="arc")
            ###############
            # Estimate gain and readout noise for the amplifiers
            msgs.work("Estimate Gain and Readout noise from the raw frames...")
            ###############
            # Generate a master arc frame
            update = slf.MasterArc(fitsdict, det)
            if update and reuseMaster:
                armbase.UpdateMasters(sciexp, sc, det, ftype="arc")
            ###############
            # Determine the dispersion direction (and transpose if necessary)
            slf.GetDispersionDirection(fitsdict, det, scidx)
            if slf._bpix[det-1] is None:
                slf.SetFrame(slf._bpix, np.zeros((slf._nspec[det-1], slf._nspat[det-1])), det)
            ###############
            # Generate a master trace frame
            update = slf.MasterTrace(fitsdict, det)
            if update and reuseMaster:
                armbase.UpdateMasters(sciexp, sc, det, ftype="flat", chktype="trace")
            ###############
            # Generate an array that provides the physical pixel locations on the detector
            slf.GetPixelLocations(det)
            ###############
            # Determine the edges of the spectrum (spatial)
            set_trace()
            lordloc, rordloc, extord = artrace.trace_orders(slf, slf._mstrace[det-1], det, pcadesc="PCA trace of the slit edges")
            slf.SetFrame(slf._lordloc, lordloc, det)
            slf.SetFrame(slf._rordloc, rordloc, det)


    return status
