""" Primary module for guiding the reduction of echelle data
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
#from pypit import arflux
#from pypit import arload
from pypit import armasters
from pypit import armbase
from pypit import armsgs
from pypit import arproc
#from pypit import arsave
from pypit import arsort
from pypit import artrace
#from pypit import arqa

from linetools import utils as ltu

try:
    from xastropy.xutils import xdebug as debugger
except ImportError:
    import pdb as debugger

# Logging
msgs = armsgs.get_logger()


def ARMED(argflag, spect, fitsdict, reuseMaster=False, reloadMaster=True):
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

    # Use Masters?  Requires setup file
    setup_file = argflag['out']['sorted']+'.setup'
    try:
        calib_dict = ltu.loadjson(setup_file)
    except:
        msgs.info("No setup file {:s} for MasterFrames".format(setup_file))
        calib_dict = {}
    else:
        argflag['masters']['setup_file'] = setup_file

    # Start reducing the data
    for sc in range(numsci):
        slf = sciexp[sc]
        scidx = slf._idx_sci[0]
        msgs.info("Reducing file {0:s}, target {1:s}".format(fitsdict['filename'][scidx], slf._target_name))
        msgs.sciexp = slf  # For QA writing on exit, if nothing else.  Could write Masters too
        if reloadMaster and (sc > 0):
            slf._argflag['masters']['use'] = True
        # Loop on Detectors
        for kk in range(slf._spect['mosaic']['ndet']):
            det = kk + 1  # Detectors indexed from 1
            slf.det = det
            ###############
            # Get amplifier sections
            arproc.get_ampsec_trimmed(slf, fitsdict, det, scidx)
            # Setup
            setup = arsort.calib_setup(slf, sc, det, fitsdict, calib_dict, write=False)
            slf._argflag['masters']['setup'] = setup
            ###############
            # Generate master bias frame
            update = slf.MasterBias(fitsdict, det)
            if update and reuseMaster:
                armbase.UpdateMasters(sciexp, sc, det, ftype="bias")
            ###############
            # Generate a bad pixel mask (should not repeat)
            update = slf.BadPixelMask(fitsdict, det)
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
            if slf._bpix[det-1] is None:  # Needs to be done here after nspec is set
                slf.SetFrame(slf._bpix, np.zeros((slf._nspec[det-1], slf._nspat[det-1])), det)
            ###############
            # Generate a master trace frame
            update = slf.MasterTrace(fitsdict, det)
            if update and reuseMaster:
                armbase.UpdateMasters(sciexp, sc, det, ftype="flat", chktype="trace")
            ###############
            # Generate an array that provides the physical pixel locations on the detector
            slf.GetPixelLocations(det)
            if ('trace'+slf._argflag['masters']['setup'] not in slf._argflag['masters']['loaded']):
                ###############
                # Determine the edges of the spectrum (spatial)
                lordloc, rordloc, extord = artrace.trace_slits(slf, slf._mstrace[det-1], det, pcadesc="PCA trace of the slit edges")
                slf.SetFrame(slf._lordloc, lordloc, det)
                slf.SetFrame(slf._rordloc, rordloc, det)


    return status
