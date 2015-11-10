import numpy as np
import scipy
import armsgs as msgs
from astropy import units as u
import arcyutils
from matplotlib import pyplot as plt
import pdb

try:
    from xastropy.xutils import xdebug as xdb
except:
    pass

# MASK VALUES FROM EXTRACTION
# 0 
# 2**0 = Flagged as bad detector pixel
# 2**1 = Flagged as affected by Cosmic Ray 
# 2**5 = Flagged as NAN (from something gone wrong)

mask_flags = dict(bad_pix=2**0, CR=2**1, NAN=2**5)

def boxcar(slf, sciframe, varframe, skyframe, crmask, scitrace):
    """
    Perform boxcar extraction on the traced objects. 

    Parameters:
    -----------
    slf
    sciframe: ndarray
      science frame
    varframe: ndarray
      variance image
    bgframe: ndarray
      sky background frame
    crmask: int ndarray
        mask of cosmic ray hits
    scitrace: dict
      traces, object and background trace images

    Returns:
    -----------
    Nothing.  slf._specobjs.boxcar is updated
    
    """
    bgfitord = 1  # Polynomial order used to fit the background
    nobj = scitrace['traces'].shape[1]
    cr_mask = 1.0-crmask
    bgfit = np.linspace(0.0, 1.0, sciframe.shape[1])
    # Loop on Objects
    for o in range(nobj):
        #pdb.set_trace()
        msgs.info("Performing boxcar extraction on object {0:d}/{1:d}".format(o+1,nobj))
        # Fit the background
        msgs.info("   Fitting the background")
        bgframe = arcyutils.func2d_fit_val(bgfit, sciframe, scitrace['background'][:,:,o]*cr_mask, bgfitord)
        # Weights
        weight = scitrace['object'][:,:,o]
        sumweight = np.sum(weight, axis=1)
        # Generate wavelength array (average over the pixels)
        wvsum = np.sum(slf._mswvimg*weight, axis=1)
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
        boxmask = np.zeros_like(wvsum).astype(int)
        # Bad detector pixels
        BPs = np.sum(weight*slf._bpix,axis=1)
        bp = BPs > 0.
        boxmask[bp] += mask_flags['bad_pix']
        # CR
        CRs = np.sum(weight*cr_mask,axis=1)
        cr = CRs > 0.
        boxmask[cr] += mask_flags['CR']
        # NAN
        NANs = np.isnan(scisum)
        if np.sum(NANs) > 0:
            msgs.warn("   NANs in the spectrum somehow..")
            boxmask[NANs] += mask_flags['NANs']
            scisum[NANs] = 0.
            varsum[NANs] = 0.
            sky[NANs] = 0.
        # Check on specobjs
        if not slf._specobjs[o].check_trace(scitrace['traces'][:,o]):
            xdb.set_trace()
            msgs.error("Bad match to specobj in boxcar!")
        # Fill
        slf._specobjs[o].boxcar['wave'] = wvsum*u.AA # Yes, units enter here
        slf._specobjs[o].boxcar['counts'] = scisum
        slf._specobjs[o].boxcar['var'] = varsum
        slf._specobjs[o].boxcar['sky'] = skysum # per pixel
        slf._specobjs[o].boxcar['mask'] = boxmask 
    # Return
    return None
