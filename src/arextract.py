import numpy as np
from astropy import units as u
import arcyutils
import artrace
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

def obj_profiles(slf, det, specobjs, sciframe, varframe, skyframe, crmask, scitrace,
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
    import pickle
    if False:
        tilts = slf._tilts[det-1]
        args = [det, specobjs, sciframe, varframe, skyframe, crmask, scitrace, tilts]
        msgs.warn("Pickling in the profile code")
        with open("trc_pickle.p",'wb') as f:
            pickle.dump(args,f)
        debugger.set_trace()
    if pickle_file is not None:
        f = open(pickle_file,'r')
        args = pickle.load(f)
        f.close()
        det, specobjs, sciframe, varframe, skyframe, crmask, scitrace, tilts = args
        slf = None
    else:
        tilts = slf._tilts[det-1]
    #
    sigframe = np.sqrt(varframe)
    # Loop
    nobj = scitrace['traces'].shape[1]
    scitrace['opt_profile'] = []
    for o in range(nobj):
        # Calculate tilts for the object trace
        # Using closest pixel for now
        slit_img = artrace.slit_image(slf, det, scitrace, o, tilts=tilts)
        # Object pixels
        weight = scitrace['object'][:,:,o]
        # Identify good rows
        gdrow = np.where(specobjs[o].boxcar['counts'] > COUNT_LIM)[0]
        # Normalized image
        norm_img = sciframe / np.outer(specobjs[o].boxcar['counts'], np.ones(sciframe.shape[1]))
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
            weight_val = sciframe[gdprof]/sigframe[gdprof]  # S/N
            # Fit Gaussian
            fdict = dict(func='gaussian', deg=3)
            mask, gfit = arutils.robust_polyfit(slit_val, flux_val, fdict['deg'],
                                                function=fdict['func'],
                                                weights=weight_val, maxone=False)
            msgs.work("Consider flagging CRs here")
            # Record
            fdict['param'] = gfit
            fdict['mask'] = mask
            scitrace['opt_profile'].append(fdict)
            # Plot?
            if False:
                gdp = mask == 0
                mn = np.min(slit_val[gdp])
                mx = np.max(slit_val[gdp])
                xval = np.linspace(mn,mx,1000)
                gauss = arutils.func_val(gfit, xval, fdict['func'])
                debugger.xplot(slit_val[gdp], flux_val[gdp], xtwo=xval, ytwo=gauss, scatter=True)
                debugger.set_trace()
        elif len(gdrow) > 10:  #
            debugger.set_trace()
    return scitrace['opt_profile']


def optimal_extract(slf, det, specobjs, sciframe, varframe, skyframe, crmask, scitrace,
                 pickle_file=None, profiles=None):
    """ Preform optimal extraction
    Standard Horne approach

    Parameters
    ----------
    slf
    det
    specobjs
    sciframe
    varframe
    crmask
    scitrace
    COUNT_LIM
    pickle_file

    Returns
    -------

    """
    import pickle
    if pickle_file is not None:
        f = open(pickle_file,'r')
        args = pickle.load(f)
        f.close()
        det, specobjs, sciframe, varframe, skyframe, crmask, scitrace, tilts = args
        slf = None
        scitrace['opt_profile'] = profiles
    else:
        tilts = None
    # Setup
    ivar = np.zeros_like(skyframe)
    cr_mask = 1.0-crmask
    msgs.warn("Should include RN too?")
    gdv = skyframe > 0.
    ivar[gdv] = 1./skyframe[gdv]
    # Loop
    nobj = scitrace['traces'].shape[1]
    for o in range(nobj):
        # Fit dict
        fit_dict = scitrace['opt_profile'][o]
        # Slit image
        slit_img = artrace.slit_image(slf, det, scitrace, o, tilts=tilts)
        msgs.warn("Turn off tilts")
        # Object pixels (avoiding CRs)
        weight = scitrace['object'][:,:,o] * cr_mask
        gdo = (weight > 0) & (ivar > 0)
        # Profile image
        prof_img = np.zeros_like(weight)
        prof_img[gdo] = arutils.func_val(fit_dict['param'], slit_img[gdo],
                                         fit_dict['func'])
        # Normalize
        norm_prof = np.sum(prof_img, axis=1)
        prof_img /= np.outer(norm_prof, np.ones(prof_img.shape[1]))
        # Mask
        mask = np.zeros_like(prof_img)
        mask[gdo] = 1.
        # Optimal flux
        opt_num = np.sum(mask * sciframe * ivar * prof_img, axis=1)
        opt_den = np.sum(mask * ivar * prof_img**2, axis=1)
        opt_flux = opt_num / opt_den
        # Optimal wave
        opt_num = np.sum(mask * slf._mswave[det-1] * ivar * prof_img**2, axis=1)
        opt_wave = opt_num / opt_den
        # Optimal variance

        # Save
        specobjs[o].optimal['wave'] = opt_wave*u.AA  # Yes, units enter here
        specobjs[o].optimal['counts'] = opt_flux
        #specobjs[o].boxcar['var'] = varsum
        #specobjs[o].boxcar['sky'] = skysum  # per pixel

        debugger.set_trace()
