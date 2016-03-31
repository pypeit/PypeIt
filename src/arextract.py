import numpy as np
from astropy import units as u
import arcyutils
import armsgs

# Logging
msgs = armsgs.get_logger()

try:
    from xastropy.xutils.xdebug import set_trace
except ImportError:
    from pdb import set_trace

# MASK VALUES FROM EXTRACTION
# 0 
# 2**0 = Flagged as bad detector pixel
# 2**1 = Flagged as affected by Cosmic Ray 
# 2**5 = Flagged as NAN (from something gone wrong)

mask_flags = dict(bad_pix=2**0, CR=2**1, NAN=2**5)


def boxcar(slf, det, specobjs, sciframe, varframe, skyframe, crmask, scitrace):
    """ Perform boxcar extraction on the traced objects.

    Parameters
    ----------
    det : int
      Detector index
    specobjs : list
      list of SpecObj objects
    sciframe : ndarray
      science frame
    varframe : ndarray
      variance image
    bgframe : ndarray
      sky background frame
    crmask : int ndarray
        mask of cosmic ray hits
    scitrace : dict
      traces, object and background trace images

    Returns
    -------
    Nothing.  slf._specobjs.boxcar is updated
    """
    bgfitord = 1  # Polynomial order used to fit the background
    nobj = scitrace['traces'].shape[1]
    cr_mask = 1.0-crmask
    bgfit = np.linspace(0.0, 1.0, sciframe.shape[1])
    # Loop on Objects
    for o in range(nobj):
        #set_trace()
        msgs.info("Performing boxcar extraction on object {0:d}/{1:d}".format(o+1,nobj))
        # Fit the background
        msgs.info("   Fitting the background")
        bgframe = arcyutils.func2d_fit_val(bgfit, sciframe, scitrace['background'][:,:,o]*cr_mask, bgfitord)
        # Weights
        weight = scitrace['object'][:,:,o]
        sumweight = np.sum(weight, axis=1)
        # Generate wavelength array (average over the pixels)
        wvsum = np.sum(slf._mswave[det-1]*weight, axis=1)
        wvsum /= sumweight
        # Generate sky spectrum (flux per pixel)
        skysum = np.sum(skyframe*weight, axis=1)
        skysum /= sumweight
        # Total the object flux
        msgs.info("   Summing object counts")
        scisum = np.sum((sciframe-bgframe)*weight, axis=1)
        # Total the variance array
        msgs.info("   Summing variance array")
        varsum = np.sum(varframe*weight, axis=1)
        # Mask 
        boxmask = np.zeros_like(wvsum, dtype=np.int)
        # Bad detector pixels
        BPs = np.sum(weight*slf._bpix[det-1], axis=1)
        bp = BPs > 0.
        boxmask[bp] += mask_flags['bad_pix']
        # CR
        CRs = np.sum(weight*cr_mask, axis=1)
        cr = CRs > 0.
        boxmask[cr] += mask_flags['CR']
        # NAN
        NANs = np.isnan(scisum)
        if np.sum(NANs) > 0:
            msgs.warn("   NANs in the spectrum somehow..")
            boxmask[NANs] += mask_flags['NANs']
            scisum[NANs] = 0.
            varsum[NANs] = 0.
            skysum[NANs] = 0.
        # Check on specobjs
        if not specobjs[o].check_trace(scitrace['traces'][:, o]):
            set_trace()
            msgs.error("Bad match to specobj in boxcar!")
        # Fill
        specobjs[o].boxcar['wave'] = wvsum*u.AA  # Yes, units enter here
        specobjs[o].boxcar['counts'] = scisum
        specobjs[o].boxcar['var'] = varsum
        specobjs[o].boxcar['sky'] = skysum  # per pixel
        specobjs[o].boxcar['mask'] = boxmask
    # Return
    return None
