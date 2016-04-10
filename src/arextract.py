import numpy as np
from astropy import units as u
import arcyutils
import arutils
import armsgs
import pdb

# Logging
msgs = armsgs.get_logger()

try:
    from xastropy.xutils import xdebug as debugger
except:
    import pdb as debugger

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
    bgcorr = np.zeros_like(cr_mask)
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
        # Update background correction image
        tmp = scitrace['background'][:,:,o] + scitrace['object'][:,:,o]
        gdp = np.where(tmp > 0)
        bgcorr[gdp] = bgframe[gdp]
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
            debugger.set_trace()
            msgs.error("Bad match to specobj in boxcar!")
        # Fill
        specobjs[o].boxcar['wave'] = wvsum*u.AA  # Yes, units enter here
        specobjs[o].boxcar['counts'] = scisum
        specobjs[o].boxcar['var'] = varsum
        specobjs[o].boxcar['sky'] = skysum  # per pixel
        specobjs[o].boxcar['mask'] = boxmask
    # Return
    return bgcorr

def obj_profiles(slf, det, specobjs, sciframe, varframe, crmask, scitrace,
                 COUNT_LIM=15., pickle_file=None):
    """ Derive spatial profiles for each object

    Parameters
    ----------
    slf
    det
    specobjs
    sciframe
    varframe
    skyframe
    crmask
    scitrace

    Returns
    -------

    """
    tilts = slf._tilts[det-1]
    if False:
        import pickle
        args = [det, specobjs, sciframe, varframe, crmask, scitrace, tilts]
        msgs.warn("Pickling in the profile code")
        with open("trc_pickle.p",'wb') as f:
            pickle.dump(args,f)
    if pickle_file is not None:
        f = open(pickle_file,'r')
        args = pickle.load(f)
        f.close()
        det, specobjs, sciframe, varframe, crmask, scitrace, tilts = args
    #
    ximg = np.outer(np.ones(sciframe.shape[0]), np.arange(sciframe.shape[1]))
    dypix = 1./sciframe.shape[0]
    # Loop
    nobj = scitrace['traces'].shape[1]
    for o in range(nobj):
        # Calculate tilts for the object trace
        # Using closest pixe for now
        msgs.work("Use 2D spline to evaluate tilts")
        xtrc = np.round(scitrace['traces'][:,o]).astype(int)
        trc_tilt = tilts[np.arange(tilts.shape[0]),xtrc]
        trc_tilt_img = np.outer(trc_tilt, np.ones(sciframe.shape[1]))
        # Slit image  (should worry about changing plate scale)
        dy = (tilts - trc_tilt_img)/dypix  # Pixels
        dx = ximg - np.outer(scitrace['traces'][:,o],np.ones(sciframe.shape[1]))
        slit_img = np.sqrt(dx**2 - dy**2)
        neg = dx < 0.
        slit_img[neg] *= -1
        # Object pixels
        weight = scitrace['object'][:,:,o]
        # Identify good rows
        gdrow = np.where(specobjs[o].boxcar['counts'] > COUNT_LIM)[0]
        # Normalized image
        norm_img = sciframe / np.outer(specobjs[o].boxcar['counts'], np.ones(sciframe.shape[1]))
        # Slit image
        # Eliminate rows with CRs (wipes out boxcar)
        crspec = np.sum(crmask*weight,axis=1)
        cr_rows = np.where(crspec > 0)[0]
        for row in cr_rows:
            weight[row,:] = 0.
        #
        if len(gdrow) > 100:  # Good S/N regime
            msgs.info("Good S/N for profile")
            # Eliminate low count regions
            badrow = np.where(specobjs[o].boxcar['counts'] < COUNT_LIM)[0]
            for row in badrow:
                weight[row,:] = 0.
            # Extract profile
            gdprof = weight > 0
            slit_val = slit_img[gdprof]
            flux_val = norm_img[gdprof]
            # Fit Gaussian
            #gauss, flag = arutils.gauss_fit(slit_val, flux_val, 0.)
            gfit = arutils.robust_polyfit(slit_val, flux_val, 3,
                                          weights=None, function='gaussian')

            debugger.set_trace()
            debugger.xplot(slit_val, flux_val, scatter=True)
        elif len(gdrow) > 10:  #
            debugger.set_trace()


