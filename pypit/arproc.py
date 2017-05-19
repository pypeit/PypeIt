from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
from scipy.signal import savgol_filter
import scipy.signal as signal
import scipy.ndimage as ndimage
import scipy.interpolate as interp
from matplotlib import pyplot as plt
from pypit import arextract
from pypit import arlris
from pypit import armsgs
from pypit import artrace
from pypit import arutils
from pypit import arparse as settings
from pypit import arspecobj
from pypit import arqa
from pypit import arpca
from pypit import arwave

from pypit import ardebug as debugger

# Logging
msgs = armsgs.get_logger()


def background_subtraction(slf, sciframe, varframe, slitn, det, refine=0.0):
    """ Generate a frame containing the background sky spectrum

    Parameters
    ----------
    slf : Class
      Science Exposure Class
    sciframe : ndarray
      science frame
    varframe : ndarray
      variance frame
    slitn : int
      Slit number
    det : int
      Detector index
    refine : float or ndarray
      refine the object traces. This should be a small value around 0.0.
      If a float, a constant offset will be applied.
      Otherwise, an array needs to be specified of the same length as
      sciframe.shape[0] that contains the refinement of each pixel along
      the spectral direction.

    Returns
    -------
    bgframe : ndarray
      An image, the same size as sciframe, that contains
      the background spectrum within the specified slit.
    nl : int
      number of pixels from the left slit edge to use as background pixels
    nr : int
      number of pixels from the right slit edge to use as background pixels
    """
    # Obtain all pixels that are within the slit edges, and are not masked
    word = np.where((slf._slitpix[det - 1] == slitn + 1) & (slf._scimask[det - 1] == 0))
    if word[0].size == 0:
        msgs.warn("There are no pixels in slit {0:d}".format(slitn))
        debugger.set_trace()
        nl, nr = 0, 0
        return np.zeros_like(sciframe), nl, nr
    # Calculate the oversampled object profiles
    xedges, modvals = object_profile(slf, sciframe, slitn, det, refine=refine, factor=3)
    bincent = 0.5*(xedges[1:]+xedges[:-1])
    npix = slf._pixwid[det - 1][slitn]
    tilts = slf._tilts[det - 1].copy()
    lordloc = slf._lordloc[det - 1][:, slitn]
    rordloc = slf._rordloc[det - 1][:, slitn]
    # For each pixel, calculate the fraction along the slit's spatial direction
    spatval = (word[1] - lordloc[word[0]] + refine) / (rordloc[word[0]] - lordloc[word[0]])
    # Cumulative sum and normalize
    csum = np.cumsum(modvals)
    csum -= csum[0]
    csum /= csum[-1]
    # Find a first guess of the edges of the object profile - assume this is the innermost 90 percent of the flux
    argl = np.argmin(np.abs(csum - 0.05))
    argr = np.argmin(np.abs(csum - 0.95))
    # Considering the possible background pixels that are left of the object,
    # find the first time where the object profile no longer decreases as you
    # move toward the edge of the slit. This is the beginning of the noisy
    # object profile, which is where the object can no longer be distinguished
    # from the background.
    wl = np.where((modvals[1:] < modvals[:-1]) & (bincent[1:] < bincent[argl]))
    wr = np.where((modvals[1:] > modvals[:-1]) & (bincent[1:] > bincent[argr]))
    nl, nr = 0, 0
    if wl[0].size != 0:
        # This is the index of the first time where the object profile
        # no longer decreases as you move towards the slit edge
        nl = np.max(wl[0])
    if wr[0].size != 0:
        # This is the index of the first time where the object profile
        # no longer decreases as you move towards the slit edge
        nr = npix - np.min(wr[0])
    if nl+nr < 5:
        msgs.warn("The object profile appears to extrapolate to the edge of the detector")
        msgs.info("A background subtraction will not be performed for slit {0:d}".format(slitn+1))
        nl, nr = 0, 0
        return np.zeros_like(sciframe), nl, nr
    # Find background pixels and fit
    wbgpix = np.where((spatval <= float(nl)/npix) | (spatval >= float(nr)/npix))
    if settings.argflag['reduce']['skysub']['method'].lower() == 'bspline':
        msgs.info("Using bspline sky subtraction")
        srt = np.argsort(tilts[wbgpix])
        ivar = arutils.calc_ivar(varframe)
        # Perform a weighted b-spline fit to the sky background pixels
        mask, bspl = arutils.robust_polyfit(tilts[wbgpix][srt], sciframe[wbgpix][srt], 3, function='bspline',
                                            weights=np.sqrt(ivar)[wbgpix][srt], sigma=5.,
                                            maxone=False, **settings.argflag['reduce']['skysub']['bspline'])
        bgf_flat = arutils.func_val(bspl, tilts.flatten(), 'bspline')
        bgframe = bgf_flat.reshape(tilts.shape)
        if msgs._debug['sky_sub']:
            def plt_bspline_sky(tilts, scifrcp, bgf_flat):
                # Setup
                srt = np.argsort(tilts.flatten())
                # Plot
                plt.close()
                plt.clf()
                ax = plt.gca()
                ax.scatter(tilts[gdp]*tilts.shape[0], scifrcp[gdp], marker='o')
                ax.plot(tilts.flatten()[srt]*tilts.shape[0], bgf_flat[srt], 'r-')
                plt.show()
            plt_bspline_sky(tilts, sciframe, bgf_flat)
            debugger.set_trace()
    else:
        msgs.error('Not ready for this method for skysub {:s}'.format(
                settings.argflag['reduce']['skysub']['method'].lower()))
    if np.any(np.isnan(bgframe)):
        msgs.warn("NAN in bgframe.  Replacing with 0")
        bad = np.isnan(bgframe)
        bgframe[bad] = 0.
    return bgframe, nl, nr


def badpix(det, frame, sigdev=10.0):
    """
    frame is a master bias frame
    sigdev is the number of standard deviations away from the median that a pixel needs to be in order to be classified as a bad pixel
    """
    dnum = settings.get_dnum(det)
    bpix = np.zeros_like(frame, dtype=np.int)
    subfr, tframe, temp = None, None, None
    for i in range(settings.spect[dnum]['numamplifiers']):
        datasec = "datasec{0:02d}".format(i+1)
        x0, x1 = settings.spect[dnum][datasec][0][0], settings.spect[dnum][datasec][0][1]
        y0, y1 = settings.spect[dnum][datasec][1][0], settings.spect[dnum][datasec][1][1]
        xv = np.arange(x0, x1)
        yv = np.arange(y0, y1)
        # Construct an array with the rows and columns to be extracted
        w = np.ix_(xv,yv)
        tframe = frame[w]
        temp = np.abs(np.median(tframe)-tframe)
        sigval = max(np.median(temp)*1.4826, 1.4826)
        ws = np.where(temp > sigdev*sigval)
        subfr = np.zeros(tframe.shape, dtype=np.int)
        subfr[ws] = 1
        bpix[w] = subfr
    del subfr, tframe, temp
    # Finally, trim the bad pixel frame
    bpix = trim(bpix, det)
    msgs.info("Identified {0:d} bad pixels".format(int(np.sum(bpix))))
    return bpix


def bg_subtraction(slf, det, sciframe, varframe, crpix, tracemask=None,
                   rejsigma=3.0, maskval=-999999.9):
    """ Extract a science target and background flux
    :param slf:
    :param sciframe:
    :param varframe:
    :return:
    """
    from pypit import arcyutils
    from pypit import arcyproc
    # Set some starting parameters (maybe make these available to the user)
    msgs.work("Should these parameters be made available to the user?")
    polyorder, repeat = 5, 1
    # Begin the algorithm
    errframe = np.sqrt(varframe)
    # Find which pixels are within the order edges
    msgs.info("Identifying pixels within each order")
    ordpix = arcyutils.order_pixels(slf._pixlocn[det-1],
                                    slf._lordloc[det-1]*0.95+slf._rordloc[det-1]*0.05,
                                    slf._lordloc[det-1]*0.05+slf._rordloc[det-1]*0.95)
    msgs.info("Applying bad pixel mask")
    ordpix *= (1-slf._bpix[det-1].astype(np.int)) * (1-crpix.astype(np.int))
    if tracemask is not None: ordpix *= (1-tracemask.astype(np.int))
    # Construct an array of pixels to be fit with a spline
    msgs.bug("Remember to include the following in a loop over order number")
    #whord = np.where(ordpix != 0)
    o = 0 # order=1
    whord = np.where(ordpix == o+1)
    tilts = slf._tilts[det-1].copy()
    xvpix  = tilts[whord]
    scipix = sciframe[whord]
    varpix = varframe[whord]
    xargsrt = np.argsort(xvpix,kind='mergesort')
    sxvpix  = xvpix[xargsrt]
    sscipix = scipix[xargsrt]
    svarpix = varpix[xargsrt]
    # Reject deviant pixels -- step through every 1.0/sciframe.shape[0] in sxvpix and reject significantly deviant pixels
    edges = np.linspace(min(0.0,np.min(sxvpix)),max(1.0,np.max(sxvpix)),sciframe.shape[0])
    fitcls = np.zeros(sciframe.shape[0])
    #if tracemask is None:
    if True:
        maskpix = np.zeros(sxvpix.size)
        msgs.info("Identifying pixels containing the science target")
        msgs.work("Speed up this step in cython")
        for i in range(sciframe.shape[0]-1):
            wpix = np.where((sxvpix>=edges[i]) & (sxvpix<=edges[i+1]))
            if (wpix[0].size>5):
                txpix = sxvpix[wpix]
                typix = sscipix[wpix]
                msk, cf = arutils.robust_polyfit(txpix, typix, 0, sigma=rejsigma)
                maskpix[wpix] = msk
                #fitcls[i] = cf[0]
                wgd=np.where(msk == 0)
                szt = np.size(wgd[0])
                if szt > 8:
                    fitcls[i] = np.mean(typix[wgd][szt//2-3:szt//2+4]) # Average the 7 middle pixels
                    #fitcls[i] = np.mean(np.random.shuffle(typix[wgd])[:5]) # Average the 5 random pixels
                else:
                    fitcls[i] = cf[0]
    else:
        msgs.work("Speed up this step in cython")
        for i in range(sciframe.shape[0]-1):
            wpix = np.where((sxvpix >= edges[i]) & (sxvpix <= edges[i+1]))
            typix = sscipix[wpix]
            szt = typix.size
            if szt > 8:
                fitcls[i] = np.mean(typix[szt//2-3:szt//2+4])  # Average the 7 middle pixels
            elif szt != 0:
                fitcls[i] = np.mean(typix)
            else:
                fitcls[i] = 0.0
        # Trace the sky lines to get a better estimate of the tilts
        scicopy = sciframe.copy()
        scicopy[np.where(ordpix == 0)] = maskval
        scitilts, _ = artrace.model_tilt(slf, det, scicopy, guesstilts=tilts.copy(), censpec=fitcls, maskval=maskval, plotQA=True)
        xvpix  = scitilts[whord]
        scipix = sciframe[whord]
        varpix = varframe[whord]
        mskpix = tracemask[whord]
        xargsrt = np.argsort(xvpix, kind='mergesort')
        sxvpix  = xvpix[xargsrt]
        sscipix = scipix[xargsrt]
        svarpix = varpix[xargsrt]
        maskpix = mskpix[xargsrt]
    # Check the mask is reasonable
    scimask = sciframe.copy()
    rxargsrt = np.argsort(xargsrt, kind='mergesort')
    scimask[whord] *= (1.0-maskpix)[rxargsrt]
    #arutils.ds9plot(scimask)
    # Now trace the sky lines to get a better estimate of the spectral tilt during the observations
    scifrcp = scimask.copy()
    scifrcp[whord] += (maskval*maskpix)[rxargsrt]
    scifrcp[np.where(ordpix == 0)] = maskval
    # Check tilts? -- Can also be error in flat fielding or slit illumination
    if msgs._debug['sky_sub']:
        gdp = scifrcp != maskval
        #debugger.xplot(tilts[gdp]*tilts.shape[0], scifrcp[gdp], scatter=True)
        idx = 1893
        if True:
            plt.clf()
            ax = plt.gca()
            ax.scatter(tilts[idx-2,:], scifrcp[idx-2,:], color='green')
            ax.scatter(tilts[idx-1,:], scifrcp[idx-1,:], color='blue')
            ax.scatter(tilts[idx,:], scifrcp[idx,:], color='red')
            ax.scatter(tilts[idx+1,:], scifrcp[idx+1,:], color='orange')
            ax.set_ylim(0., 3000)
            plt.show()
            debugger.set_trace()
    #
    msgs.info("Fitting sky background spectrum")
    if settings.argflag['reduce']['skysub']['method'].lower() == 'polyscan':
        polypoints = 5
        nsmth = 15
        bgmodel = arcyproc.polyscan_fitsky(tilts.copy(), scifrcp.copy(), 1.0/errframe, maskval, polyorder, polypoints, nsmth, repeat)
        bgpix = bgmodel[whord]
        sbgpix = bgpix[xargsrt]
        wbg = np.where(sbgpix != maskval)
        # Smooth this spectrum
        polyorder = 1
        xpix = sxvpix[wbg]
        maxdiff = np.sort(xpix[1:]-xpix[:-1])[xpix.size-sciframe.shape[0]-1] # only include the next pixel in the fit if it is less than 10x the median difference between all pixels
        msgs.info("Generating sky background image")
        #if msgs._debug['sky_sub']:
        #    debugger.set_trace()
        #    debugger.xplot(sxvpix[wbg]*tilts.shape[0], sbgpix[wbg], scatter=True)
        bgscan = arcyutils.polyfit_scan_lim(sxvpix[wbg], sbgpix[wbg].copy(), np.ones(wbg[0].size,dtype=np.float), maskval, polyorder, sciframe.shape[1]/3, repeat, maxdiff)
        # Restrict to good values
        gdscan = bgscan != maskval
        if msgs._debug['sky_sub']:
            debugger.set_trace()
            debugger.xplot(sxvpix[wbg[0][gdscan]]*tilts.shape[0], sbgpix[wbg[0][gdscan]], scatter=True)
        if np.sum(~gdscan) > 0:
            msgs.warn("At least one masked value in bgscan")
        # Generate
        bgframe = np.interp(tilts.flatten(), sxvpix[wbg[0][gdscan]], bgscan[gdscan]).reshape(tilts.shape)
    elif settings.argflag['reduce']['skysub']['method'].lower() == 'bspline':
        msgs.info("Using bspline sky subtraction")
        gdp = scifrcp != maskval
        srt = np.argsort(tilts[gdp])
        #bspl = arutils.func_fit(tilts[gdp][srt], scifrcp[gdp][srt], 'bspline', 3,
        #                        **settings.argflag['reduce']['skysub']['bspline'])
        ivar = arutils.calc_ivar(varframe)
        mask, bspl = arutils.robust_polyfit(tilts[gdp][srt], scifrcp[gdp][srt], 3, function='bspline',
                                            weights=np.sqrt(ivar)[gdp][srt], sigma=5.,
                                            maxone=False, **settings.argflag['reduce']['skysub']['bspline'])
        bgf_flat = arutils.func_val(bspl, tilts.flatten(), 'bspline')
        bgframe = bgf_flat.reshape(tilts.shape)
        if msgs._debug['sky_sub']:
            def plt_bspline_sky(tilts, scifrcp, bgf_flat, maskval):
                # Setup
                srt = np.argsort(tilts.flatten())
                # Plot
                plt.close()
                plt.clf()
                ax = plt.gca()
                ax.scatter(tilts[gdp]*tilts.shape[0], scifrcp[gdp], marker='o')
                ax.plot(tilts.flatten()[srt]*tilts.shape[0], bgf_flat[srt], 'r-')
                plt.show()
            plt_bspline_sky(tilts, scifrcp, bgf_flat, maskval)
            debugger.set_trace()
    else:
        msgs.error('Not ready for this method for skysub {:s}'.format(
                settings.argflag['reduce']['skysub']['method'].lower()))
    if np.sum(np.isnan(bgframe)) > 0:
        msgs.warn("NAN in bgframe.  Replacing with 0")
        bad = np.isnan(bgframe)
        bgframe[bad] = 0.
    if msgs._debug['sky_sub']:
        debugger.set_trace()
        #debugger.show_image(sciframe-bgframe)
    # Plot to make sure that the result is good
    #arutils.ds9plot(bgframe)
    #arutils.ds9plot(sciframe-bgframe)
    return bgframe


def error_frame_postext(sciframe, idx, fitsdict):
    # Dark Current noise
    dnoise = settings.spect['det']['darkcurr'] * float(fitsdict["exptime"][idx])/3600.0
    # The effective read noise
    rnoise = settings.spect['det']['ronoise']**2 + (0.5*settings.spect['det']['gain'])**2
    errframe = np.zeros_like(sciframe)
    w = np.where(sciframe != -999999.9)
    errframe[w] = np.sqrt(sciframe[w] + rnoise + dnoise)
    w = np.where(sciframe == -999999.9)
    errframe[w] = 999999.9
    return errframe


def flatfield(slf, sciframe, flatframe, det, snframe=None,
              varframe=None, slitprofile=None):
    """ Flat field the input image
    Parameters
    ----------
    slf
    sciframe : 2d image
    flatframe : 2d image
    snframe : 2d image, optional
    det : int
      Detector index
    varframe : ndarray
      variance image
    slitprofile : ndarray
      slit profile image

    Returns
    -------
    flat-field image
    and updated sigma array if snframe is input
    or updated variance array if varframe is input

    """
    if (varframe is not None) & (snframe is not None):
        msgs.error("Cannot set both varframe and snframe")
    if slitprofile is not None:
        flatframe *= slitprofile
    # New image
    retframe = np.zeros_like(sciframe)
    w = np.where(flatframe > 0.0)
    retframe[w] = sciframe[w]/flatframe[w]
    if w[0].size != flatframe.size:
        ww = np.where(flatframe <= 0.0)
        slf._bpix[det-1][ww] = 1.0
    # Variance?
    if varframe is not None:
        retvar = np.zeros_like(sciframe)
        retvar[w] = varframe[w]/flatframe[w]**2
        return retframe, retvar
    # Error image
    if snframe is None:
        return retframe
    else:
        errframe = np.zeros_like(sciframe)
        wnz = np.where(snframe>0.0)
        errframe[wnz] = retframe[wnz]/snframe[wnz]
        return retframe, errframe


def flatnorm(slf, det, msflat, maskval=-999999.9, overpix=6, plotdesc=""):
    """ Normalize the flat-field frame

    Parameters
    ----------
    slf : class
      An instance of the Science Exposure class
    det : int
      Detector number
    msflat : ndarray
      Flat-field image
    maskval : float
      Global floating point mask value used throughout the code
    overpix : int
      overpix/2 = the number of pixels to extend beyond each side of the order trace
    plotdesc : str
      A title for the plotted QA

    Returns
    -------
    msnormflat : ndarray
      The normalized flat-field frame
    msblaze : ndarray
      A 2d array containing the blaze function for each slit
    """
    from pypit import arcyutils
    from pypit import arcyextract
    from pypit import arcyproc
    dnum = settings.get_dnum(det)

    msgs.info("Normalizing the master flat field frame")
    norders = slf._lordloc[det-1].shape[1]
    # First, determine the relative scale of each amplifier (assume amplifier 1 has a scale of 1.0)
    if (settings.spect[dnum]['numamplifiers'] > 1) & (norders > 1):
        sclframe = get_ampscale(slf, det, msflat)
        # Divide the master flat by the relative scale frame
        msflat /= sclframe
    else:
        sclframe = np.ones(msflat, dtype=np.float)
    # Determine the blaze
    polyord_blz = 2  # This probably doesn't need to be a parameter that can be set by the user
    # Look at the end corners of the detector to get detector size in the dispersion direction
    #xstr = slf._pixlocn[det-1][0,0,0]-slf._pixlocn[det-1][0,0,2]/2.0
    #xfin = slf._pixlocn[det-1][-1,-1,0]+slf._pixlocn[det-1][-1,-1,2]/2.0
    #xint = slf._pixlocn[det-1][:,0,0]
    # Find which pixels are within the order edges
    msgs.info("Identifying pixels within each order")
    ordpix = arcyutils.order_pixels(slf._pixlocn[det-1], slf._lordloc[det-1], slf._rordloc[det-1])
    msgs.info("Applying bad pixel mask")
    ordpix *= (1-slf._bpix[det-1].astype(np.int))
    mskord = np.zeros(msflat.shape)
    msgs.info("Rectifying the orders to estimate the background locations")
    #badorders = np.zeros(norders)
    msnormflat = maskval*np.ones_like(msflat)
    msblaze = maskval*np.ones((msflat.shape[0],norders))
    msgs.work("Must consider different amplifiers when normalizing and determining the blaze function")
    msgs.work("Multiprocess this step to make it faster")
    flat_ext1d = maskval*np.ones((msflat.shape[0],norders))
    for o in range(norders):
        if settings.argflag["reduce"]["flatfield"]["method"].lower() == "bspline":
            msgs.info("Deriving blaze function of slit {0:d} with a bspline".format(o+1))
            tilts = slf._tilts[det - 1].copy()
            gdp = (msflat != maskval) & (ordpix == o + 1)
            srt = np.argsort(tilts[gdp])
            everyn = settings.argflag['reduce']['flatfield']['params'][0]
            if everyn > 0.0 and everyn < 1.0:
                everyn *= msflat.shape[0]
                everyn = int(everyn + 0.5)
            everyn *= slf._pixwid[det - 1][o]
            if np.where(gdp)[0].size < 2*everyn:
                msgs.warn("Not enough pixels in slit {0:d} to fit a bspline")
                continue
            bspl = arutils.func_fit(tilts[gdp][srt], msflat[gdp][srt], 'bspline', 3, everyn=everyn)
            model_flat = arutils.func_val(bspl, tilts.flatten(), 'bspline')
            model = model_flat.reshape(tilts.shape)
            word = np.where(ordpix == o + 1)
            msnormflat[word] = msflat[word] / model[word]
            msblaze[:, o] = arutils.func_val(bspl, np.linspace(0.0, 1.0, msflat.shape[0]), 'bspline')
            mskord[word] = 1.0
            flat_ext1d[:, o] = np.sum(msflat * mskord, axis=1) / np.sum(mskord, axis=1)
            mskord *= 0.0
        elif settings.argflag["reduce"]["flatfield"]["method"].lower() == "polyscan":
            # Rectify this order
            recframe = arcyextract.rectify(msflat, ordpix, slf._pixcen[det - 1][:, o], slf._lordpix[det - 1][:, o],
                                           slf._rordpix[det - 1][:, o], slf._pixwid[det - 1][o] + overpix, maskval)
            polyorder = settings.argflag["reduce"]["flatfield"]["params"][0]
            polypoints = settings.argflag["reduce"]["flatfield"]["params"][1]
            repeat = settings.argflag["reduce"]["flatfield"]["params"][2]
            # Take the median along the spatial dimension
            flatmed = np.median(recframe, axis=1)
            # Perform a polynomial fitting scheme to determine the blaze profile
            xarray = np.arange(flatmed.size, dtype=np.float)
            weight = flatmed.copy()
            msgs.work("Routine doesn't support user parameters yet")
            msgs.bug("Routine doesn't support user parameters yet")
            blazet = arcyutils.polyfit_scan(xarray, flatmed.copy(), weight, maskval, polyorder, polypoints, repeat)
             # Remove the masked endpoints
            outx, outy, outm, lox, hix = arcyproc.remove_maskedends(xarray, flatmed, blazet, maskval)
            # Inspect the end points and extrapolate from the best fitting end pixels
            derv = (outm[1:]-outm[:-1])/(outx[1:]-outx[:-1])
            dervx = 0.5*(outx[1:]+outx[:-1])
            derv2 = (derv[1:]-derv[:-1])/(dervx[1:]-dervx[:-1])
            medv = np.median(derv2)
            madv = 1.4826*np.median(np.abs(derv2-medv))
            blaze = arcyproc.blaze_fitends(outx, outy, outm, derv2-medv, madv, polyord_blz, polypoints)
            #plt.plot(xarray,flatmed,'k-',drawstyle='steps')
            #plt.plot(xarray, blaze, 'r-')
            #plt.show()
            #np.savetxt("check_blaze_ord{0:d}.txt".format(o),np.transpose((xarray,flatmed)))
            # Divide the flat by the fitted flat profile
            finalblaze = np.ones(recframe.shape[0])
            finalblaze[lox:hix] = blaze.copy()
            blazenrm = finalblaze.reshape((finalblaze.size, 1)).repeat(recframe.shape[1], axis=1)
            recframe /= blazenrm
            # Store the blaze for this order
            msblaze[lox:hix,o] = blaze.copy()
            flat_ext1d[:,o] = flatmed.copy()
            # Sort the normalized frames along the dispersion direction
            recsort = np.sort(recframe, axis=0)
            # Find the mean value, but only consider the "innermost" 50 per cent of pixels (i.e. the pixels closest to 1.0)
            recmean = arcyproc.scale_blaze(recsort, maskval)
            #rows = np.arange(recsort.shape[0]/4,(3*recsort.shape[0])/4,dtype=np.int)
            #w = np.ix_(rows,np.arange(recframe.shape[1]))
            #recmean = np.mean(recsort[w],axis=0)
            for i in range(recmean.size):
                recframe[:, i] /= recmean[i]
            # Undo the rectification
            normflat_unrec = arcyextract.rectify_undo(recframe, slf._pixcen[det-1][:,o], slf._lordpix[det-1][:,o],
                                                      slf._rordpix[det-1][:,o], slf._pixwid[det-1][o], maskval,
                                                      msflat.shape[0], msflat.shape[1])
            # Apply the normalized flatfield for this order to the master normalized frame
            msnormflat = arcyproc.combine_nrmflat(msnormflat, normflat_unrec, slf._pixcen[det-1][:,o],
                                                  slf._lordpix[det-1][:,o], slf._rordpix[det-1][:,o],
                                                  slf._pixwid[det-1][o]+overpix, maskval)
        else:
            msgs.error("Flatfield method {0:s} is not supported".format(settings.argflag["reduce"]["flatfield"]["method"]))
    # Send the blaze away to be plotted and saved
    if "2dpca" in settings.argflag["reduce"]["flatfield"].keys():
        if settings.argflag["reduce"]["flatfield"]["2dpca"] >= 1:
            msgs.info("Performing a 2D PCA on the blaze fits")
            msblaze = arpca.pca2d(msblaze, settings.argflag["reduce"]["flatfield"]["2dpca"])
    # Plot the blaze model
    if not msgs._debug['no_qa']:
        msgs.info("Saving blaze fits to QA")
        arqa.plot_orderfits(slf, msblaze, flat_ext1d, desc=plotdesc, textplt="Order")
    # If there is more than 1 amplifier, apply the scale between amplifiers to the normalized flat
    if (settings.spect[dnum]['numamplifiers'] > 1) & (norders > 1):
        msnormflat *= sclframe
    return msnormflat, msblaze


def get_ampscale(slf, det, msflat):
    """ Normalize the flat-field frame

    Parameters
    ----------
    slf : class
      An instance of the Science Exposure class
    det : int
      Detector number
    msflat : ndarray
      Flat-field image

    Returns
    -------
    sclframe : ndarray
      A frame to scale all amplifiers to the same counts at the amplifier borders
    """
    dnum = settings.get_dnum(det)

    sclframe = np.ones_like(msflat)
    ampdone = np.zeros(settings.spect[dnum]['numamplifiers'], dtype=int) # 1 = amplifiers have been assigned a scale
    ampdone[0]=1
    while np.sum(ampdone) != settings.spect[dnum]['numamplifiers']:
        abst, bbst, nbst, n0bst, n1bst = -1, -1, -1, -1, -1 # Reset the values for the most overlapping amplifier
        for a in range(0, settings.spect[dnum]['numamplifiers']): # amplifier 'a' is always the reference amplifier
            if ampdone[a] == 0: continue
            for b in range(0, settings.spect[dnum]['numamplifiers']):
                if ampdone[b] == 1 or a == b: continue
                tstframe = np.zeros_like(msflat)
                tstframe[np.where(slf._datasec[det-1] == a+1)] = 1
                tstframe[np.where(slf._datasec[det-1] == b+1)] = 2
                # Determine the total number of adjacent edges between amplifiers a and b
                n0 = np.sum(tstframe[1:,:]-tstframe[:-1,:])
                n1 = np.sum(tstframe[:,1:]-tstframe[:,:-1])
                if (abs(n0)+abs(n1)) > nbst:
                    n0bst = n0
                    n1bst = n1
                    nbst = abs(n0)+abs(n1)
                    abst = a
                    bbst = b
        # Determine the scaling factor for these two amplifiers
        tstframe = np.zeros_like(msflat)
        tstframe[np.where(slf._datasec[det-1] == abst+1)] = 1
        tstframe[np.where(slf._datasec[det-1] == bbst+1)] = 2
        if abs(n0bst) > abs(n1bst):
            # The amplifiers overlap on the zeroth index
            w = np.where(tstframe[1:,:]-tstframe[:-1,:] != 0)
            sclval = np.median(msflat[w[0][0]+1, w[1]])/np.median(msflat[w[0][0], w[1]])
            # msflat[w[0][0], w[1][0:50]] = 1.0E10
            # msflat[w[0][0]-1, w[1][0:50]] = -1.0E10
            # arutils.ds9plot(msflat)
            if n0bst > 0:
                # Then pixel w[0][0] falls on amplifier a
                sclval = sclframe[w[0][0], w[1]] * sclval
            else:
                # pixel w[0][0] falls on amplifier b
                sclval = sclframe[w[0][0]+1, w[1]] / sclval
        else:
            # The amplifiers overlap on the first index
            w = np.where(tstframe[:,1:]-tstframe[:,:-1] != 0)
            sclval = np.median(msflat[w[0], w[1][0]+1]/msflat[w[0], w[1][0]])
            if n1bst > 0:
                # Then pixel w[1][0] falls on amplifier a
                sclval = sclframe[w[0], w[1][0]] * sclval
            else:
                # pixel w[1][0] falls on amplifier b
                sclval = sclframe[w[0], w[1][0]+1] / sclval
        # Finally, apply the scale factor thwe amplifier b
        w = np.where(slf._datasec[det-1] == bbst+1)
        sclframe[w] = np.median(sclval)
        ampdone[bbst] = 1
    return sclframe


def get_datasec_trimmed(slf, fitsdict, det, scidx):
    """
     Generate a frame that identifies each pixel to an amplifier, and then trim it to the data sections.
     This frame can be used to later identify which trimmed pixels correspond to which amplifier

    Parameters
    ----------
    slf : class
      An instance of the ScienceExposure class
    fitsdict : dict
      Contains relevant information from fits header files
    det : int
      Detector number, starts at 1
    scidx : int
      Index of science frame

    Returns
    -------
    fitsdict : dict
      Updates to the input fitsdict
    """
    dnum = settings.get_dnum(det)

    # Get naxis0, naxis1, datasec, oscansec, ampsec for specific instruments
    if settings.argflag['run']['spectrograph'] in ['lris_blue', 'lris_red']:
        msgs.info("Parsing datasec and oscansec from headers")
        temp, head0, secs = arlris.read_lris(fitsdict['directory'][scidx]+
                                             fitsdict['filename'][scidx],
                                             det)
        # Naxis
        fitsdict['naxis0'][scidx] = temp.shape[0]
        fitsdict['naxis1'][scidx] = temp.shape[1]
        # Loop on amplifiers
        for kk in range(settings.spect[dnum]['numamplifiers']):
            datasec = "datasec{0:02d}".format(kk+1)
            settings.spect[dnum][datasec] = settings.load_sections(secs[0][kk])
            oscansec = "oscansec{0:02d}".format(kk+1)
            settings.spect[dnum][oscansec] = settings.load_sections(secs[1][kk])
    # For convenience
    naxis0, naxis1 = int(fitsdict['naxis0'][scidx]), int(fitsdict['naxis1'][scidx])
    # Initialize the returned array
    retarr = np.zeros((naxis0, naxis1))
    for i in range(settings.spect[dnum]['numamplifiers']):
        datasec = "datasec{0:02d}".format(i+1)
        x0, x1 = settings.spect[dnum][datasec][0][0], settings.spect[dnum][datasec][0][1]
        y0, y1 = settings.spect[dnum][datasec][1][0], settings.spect[dnum][datasec][1][1]
        if x0 < 0: x0 += naxis0
        if x1 <= 0: x1 += naxis0
        if y0 < 0: y0 += naxis1
        if y1 <= 0: y1 += naxis1
        # Fill in the pixels for this amplifier
        xv = np.arange(x0, x1)
        yv = np.arange(y0, y1)
        w = np.ix_(xv, yv)
        try:
            retarr[w] = i+1
        except IndexError:
            debugger.set_trace()
        # Save these locations for trimming
        if i == 0:
            xfin = xv.copy()
            yfin = yv.copy()
        else:
            xfin = np.unique(np.append(xfin, xv.copy()))
            yfin = np.unique(np.append(yfin, yv.copy()))
    # Construct and array with the rows and columns to be extracted
    w = np.ix_(xfin, yfin)
    slf._datasec[det-1] = retarr[w]
    return


def get_wscale(slf):
    """
    This routine calculates the wavelength array based on the sampling size (in km/s) of each pixel.
    It conveniently assumes a standard reference wavelength of 911.75348 A
    """

    lam0 = 911.75348
    step = 1.0 + settings.argflag['reduce']['pixelsize']/299792.458
    # Determine the number of pixels from lam0 that need to be taken to reach the minimum wavelength of the spectrum
    msgs.work("No orders should be masked -- remove this code when the auto wavelength ID routine is fixed, and properly extrapolates.")
    w = np.where(slf._waveids!=-999999.9)
    nmin = int(np.log10(np.min(slf._waveids[w])/lam0)/np.log10(step) )
    nmax = int(1.0 + np.log10(np.max(slf._waveids[w])/lam0)/np.log10(step) ) # 1.0+ is to round up
    wave = np.min(slf._waveids[w]) * (step**np.arange(1+nmax-nmin))
    msgs.info("Extracted wavelength range will be: {0:.5f} - {1:.5f}".format(wave.min(),wave.max()))
    msgs.info("Total number of spectral pixels in the extracted spectrum will be: {0:d}".format(1+nmax-nmin))
    return wave


def object_profile(slf, sciframe, slitn, det, refine=0.0, factor=3):
    """ Generate an array of the object profile

    Parameters
    ----------
    slf : Class
      Science Exposure Class
    sciframe : ndarray
      science frame
    slitn : int
      Slit number
    det : int
      Detector index
    refine : float or ndarray
      refine the object traces. This should be a small value around 0.0.
      If a float, a constant offset will be applied.
      Otherwise, an array needs to be specified of the same length as
      sciframe.shape[0] that contains the refinement of each pixel along
      the spectral direction.
    factor : int, optional
      Sampling factor. factor=1 samples the object profile
      with the number of pixels along the length of the slit.
      factor=2 samples with twice the number of pixels along
      the length of the slit, etc.

    Returns
    -------
    xedges : ndarray
      bin edges
    profile : ndarray
      object profile
    """
    # Obtain the indices of the pixels that are in slit number 'slitn', and are not masked
    word = np.where((slf._slitpix[det - 1] == slitn + 1) & (slf._scimask[det - 1] == 0))
    if word[0].size == 0:
        msgs.warn("There are no pixels in slit {0:d}".format(slitn))
        return None, None
    # Determine the width of the slit in pixels, and calculate the
    # number of bins needed to oversample the object profile.
    npix = slf._pixwid[det-1][slitn]
    nbins = factor*npix
    # Extract the left and right order locations, and estimate the spatial positions
    # of all pixels within the slit.
    lordloc = slf._lordloc[det - 1][:, slitn]
    rordloc = slf._rordloc[det - 1][:, slitn]
    spatval = (word[1] - lordloc[word[0]] + refine) / (rordloc[word[0]] - lordloc[word[0]])
    # Create an array to store the oversampled object profile
    profile = np.zeros(nbins)
    # Determine the bin edges of the oversampled array
    xedges = np.linspace(np.min(spatval), np.max(spatval), nbins+1)
    # Assign each detector pixel within the slit to an oversampled pixel
    groups = np.digitize(spatval, xedges)
    flxfr = sciframe[word]
    # For each oversampled pixel, calculate the median flux
    msgs.work("It might be a good idea to use a weighted mean (where weights=flux), instead of the median here")
    for mm in range(1, xedges.size):
        profile[mm - 1] = np.median(flxfr[groups == mm])
    return xedges, profile


def reduce_prepare(slf, sciframe, scidx, fitsdict, det, standard=False):
    """ Prepare the Run standard extraction steps on a frame

    Parameters
    ----------
    sciframe : image
      Bias subtracted image (using arload.load_frame)
    scidx : int
      Index of the frame
    fitsdict : dict
      Contains relevant information from fits header files
    det : int
      Detector index
    standard : bool, optional
      Standard star frame?
    """
    # Check inputs
    if not isinstance(scidx, (int,np.integer)):
        raise IOError("scidx needs to be an int")
    # Convert ADUs to electrons
    sciframe *= gain_frame(slf, det)
    # Mask
    slf._scimask[det-1] = np.zeros_like(sciframe).astype(int)
    msgs.info("Masking bad pixels")
    slf.update_sci_pixmask(det, slf._bpix[det-1], 'BadPix')
    # Variance
    msgs.info("Generate raw variance frame (from detected counts [flat fielded])")
    rawvarframe = variance_frame(slf, det, sciframe, scidx, fitsdict)
    ###############
    # Subtract off the scattered light from the image
    msgs.work("Scattered light subtraction is not yet implemented...")
    ###############
    # Flat field the science frame (and variance)
    if settings.argflag['reduce']['flatfield']['perform']:
        msgs.info("Flat fielding the science frame")
        sciframe, rawvarframe = flatfield(slf, sciframe, slf._mspixelflatnrm[det-1], det,
                                          varframe=rawvarframe, slitprofile=slf._slitprof[det-1])
    else:
        msgs.info("Not performing a flat field calibration")
    if not standard:
        slf._sciframe[det-1] = sciframe
        slf._rawvarframe[det-1] = rawvarframe
    ###############
    # Identify cosmic rays
    msgs.work("Include L.A.Cosmic arguments in the settings files")
    if True: crmask = lacosmic(slf, fitsdict, det, sciframe, scidx, grow=1.5)
    else: crmask = np.zeros(sciframe.shape)
    # Mask
    slf.update_sci_pixmask(det, crmask, 'CR')
    return sciframe, rawvarframe, crmask


def reduce_echelle(slf, sciframe, scidx, fitsdict, det,
                   standard=False, triml=1, trimr=1):
    """ Run standard extraction steps on an echelle frame

    Parameters
    ----------
    sciframe : image
      Bias subtracted image (using arload.load_frame)
    scidx : int
      Index of the frame
    fitsdict : dict
      Contains relevant information from fits header files
    det : int
      Detector index
    standard : bool, optional
      Standard star frame?
    triml : int (optional)
      Number of pixels to trim from the left slit edge
    trimr : int (optional)
      Number of pixels to trim from the right slit edge
    """
    msgs.work("Multiprocess this algorithm")
    nspec = sciframe.shape[0]
    nord = slf._lordloc[det-1].shape[1]
    # Prepare the frames for tracing and extraction
    sciframe, rawvarframe, crmask = reduce_prepare(slf, sciframe, scidx, fitsdict, det, standard=standard)
    bgframe = np.zeros_like(sciframe)
    bgnl, bgnr = np.zeros(nord, dtype=np.int), np.zeros(nord, dtype=np.int)
    skysub = True
    if settings.argflag['reduce']['skysub']['perform']:
        # Identify background pixels, and generate an image of the sky spectrum in each slit
        for o in range(nord):
            tbgframe, nl, nr = background_subtraction(slf, sciframe, rawvarframe, o, det)
            bgnl[o], bgnr[o] = nl, nr
            bgframe += tbgframe
            if nl == 0 and nr == 0:
                # If just one slit cannot do sky subtraction, don't do sky subtraction
                msgs.warn("A sky subtraction will not be performed")
                skysub = False
                bgframe = np.zeros_like(sciframe)
                modelvarframe = rawvarframe.copy()
                break
        if skysub:
            # Provided the for loop above didn't break early, model the variance frame
            modelvarframe = variance_frame(slf, det, sciframe, scidx, fitsdict, skyframe=bgframe)
    else:
        modelvarframe = rawvarframe.copy()
        bgframe = np.zeros_like(sciframe)
    if not standard:  # Need to save
        slf._modelvarframe[det - 1] = modelvarframe
        slf._bgframe[det - 1] = bgframe
    # Obtain a first estimate of the object trace then
    # fit the traces and perform a PCA for the refinements
    trccoeff = np.zeros((settings.argflag['trace']['object']['order']+1, nord))
    trcxfit = np.arange(nspec)
    for o in range(nord):
        trace, error = artrace.trace_weighted(sciframe-bgframe, slf._lordloc[det-1][:, o], slf._rordloc[det-1][:, o],
                                              mask=slf._scimask[det-1], wght="flux")
        trace -= slf._lordloc[det-1][:, o]
        trace /= (slf._rordloc[det-1][:, o]-slf._lordloc[det-1][:, o])
        msk, trccoeff[:, o] = arutils.robust_polyfit(trcxfit, trace,
                                                     settings.argflag['trace']['object']['order'],
                                                     function=settings.argflag['trace']['object']['function'],
                                                     weights=1.0 / error ** 2, minv=0.0, maxv=nspec-1.0)
    refine = 0.0
    if settings.argflag['trace']['object']['method'] == "pca":
        # Identify the orders to be extrapolated during reconstruction
        orders = 1.0 + np.arange(nord)
        msgs.info("Performing a PCA on the object trace")
        ofit = settings.argflag['trace']['object']['params']
        lnpc = len(ofit) - 1
        msgs.work("May need to do a check here to make sure ofit is reasonable")
        xcen = trcxfit[:, np.newaxis].repeat(nord, axis=1)
        trccen = arutils.func_val(trccoeff, trcxfit, settings.argflag['trace']['object']['function'],
                                  minv=0.0, maxv=nspec-1.0).T
        fitted, outpar = arpca.basis(xcen, trccen, trccoeff, lnpc, ofit, skipx0=False,
                                     function=settings.argflag['trace']['object']['function'])
        if not msgs._debug['no_qa']:
            arpca.pc_plot(slf, outpar, ofit, pcadesc="PCA of object trace")
        # Extrapolate the remaining orders requested
        trccen, outpar = arpca.extrapolate(outpar, orders, function=settings.argflag['trace']['object']['function'])
        #refine = trccen-trccen[nspec//2, :].reshape((1, nord))
    else:
        msgs.error("Not ready for object trace method:" + msgs.newline() +
                   settings.argflag['trace']['object']['method'])
    # Construct the left and right traces of the object profile
    # The following code ensures that the fraction of the slit
    # containing the object remains constant along the spectral
    # direction
    trcmean = np.mean(trccen, axis=0)
    trobjl = (trcmean - (1+bgnl)/slf._pixwid[det - 1].astype(np.float)).reshape((1, nord)).repeat(nspec, axis=0)
    trobjl = trccen - trobjl
    trobjr = (-trcmean + (slf._pixwid[det - 1]-bgnr-1)/slf._pixwid[det - 1].astype(np.float)).reshape((1, nord)).repeat(nspec, axis=0)
    trobjr = trccen + trobjr
    # Convert trccen to the actual trace locations
    trccen *= (slf._rordloc[det - 1] - slf._lordloc[det - 1])
    trccen += slf._lordloc[det - 1]
    trobjl *= (slf._rordloc[det - 1] - slf._lordloc[det - 1])
    trobjl += slf._lordloc[det - 1]
    trobjr *= (slf._rordloc[det - 1] - slf._lordloc[det - 1])
    trobjr += slf._lordloc[det - 1]

    # Generate an image of pixel weights for each object
    msgs.work("Eventually allow ARMED to find multiple objects in the one slit")
    nobj = 1
    rec_obj_img = np.zeros(sciframe.shape+(nobj,))
    rec_bg_img = np.zeros(sciframe.shape+(nobj,))
    for o in range(nord):
        # Prepare object/background regions
        objl = np.array([bgnl[o]])
        objr = np.array([slf._pixwid[det - 1][o]-bgnr[o]-triml-trimr])
        bckl = np.zeros((slf._pixwid[det - 1][o]-triml-trimr, 1))
        bckr = np.zeros((slf._pixwid[det - 1][o]-triml-trimr, 1))
        bckl[:bgnl[o]] = 1
        if bgnr[o] != 0:
            bckr[-bgnr[o]:] = 1
        tobj_img, tbg_img = artrace.trace_objbg_image(slf, det, sciframe-bgframe, o,
                                                      [objl, objr], [bckl, bckr],
                                                      triml=triml, trimr=trimr)
        rec_obj_img += tobj_img
        rec_bg_img += tbg_img

    # Create trace dict
    scitrace = artrace.trace_object_dict(nobj, trccen[:, 0].reshape(trccen.shape[0], 1),
                                         object=rec_obj_img, background=rec_bg_img)
    for o in range(1, nord):
        scitrace = artrace.trace_object_dict(nobj, trccen[:, o].reshape(trccen.shape[0], 1),
                                             tracelist=scitrace)

    # Save the quality control
    if not msgs._debug['no_qa']:
        arqa.obj_trace_qa(slf, sciframe, trobjl, trobjr, None, root="object_trace", normalize=False)

    # Finalize the Sky Background image
    if settings.argflag['reduce']['skysub']['perform'] and (nobj > 0) and skysub:
        # Identify background pixels, and generate an image of the sky spectrum in each slit
        bgframe = np.zeros_like(sciframe)
        for o in range(nord):
            tbgframe, nl, nr = background_subtraction(slf, sciframe, rawvarframe, o, det, refine=refine)
            bgnl[o], bgnr[o] = nl, nr
            bgframe += tbgframe
        modelvarframe = variance_frame(slf, det, sciframe, scidx, fitsdict, skyframe=bgframe)

    # Perform an optimal extraction
    return reduce_frame(slf, sciframe, rawvarframe, modelvarframe, bgframe, scidx, fitsdict, det, crmask,
                        scitrace=scitrace, standard=standard)


def reduce_multislit(slf, sciframe, scidx, fitsdict, det, standard=False):
    """ Run standard extraction steps on an echelle frame

    Parameters
    ----------
    sciframe : image
      Bias subtracted image (using arload.load_frame)
    scidx : int
      Index of the frame
    fitsdict : dict
      Contains relevant information from fits header files
    det : int
      Detector index
    standard : bool, optional
      Standard star frame?
    """
    sciframe, rawvarframe, crmask = reduce_prepare(slf, sciframe, scidx, fitsdict, det, standard=standard)

    ###############
    # Estimate Sky Background
    if settings.argflag['reduce']['skysub']['perform']:
        # Perform an iterative background/science extraction
        if msgs._debug['obj_profile'] and False:
            msgs.warn("Reading background from 2D image on disk")
            from astropy.io import fits
            datfil = settings.argflag['run']['directory']['science']+'/spec2d_{:s}.fits'.format(slf._basename.replace(":","_"))
            hdu = fits.open(datfil)
            bgframe = hdu[1].data - hdu[2].data
        else:
            msgs.info("First estimate of the sky background")
            bgframe = bg_subtraction(slf, det, sciframe, rawvarframe, crmask)
        #bgframe = bg_subtraction(slf, det, sciframe, varframe, crmask)
        modelvarframe = variance_frame(slf, det, sciframe, scidx, fitsdict, skyframe=bgframe)
    else:
        modelvarframe = rawvarframe.copy()
        bgframe = np.zeros_like(sciframe)
    if not standard:  # Need to save
        slf._modelvarframe[det - 1] = modelvarframe
        slf._bgframe[det - 1] = bgframe

    ###############
    # Estimate trace of science objects
    scitrace = artrace.trace_object(slf, det, sciframe-bgframe, modelvarframe, crmask, doqa=False)# (not standard))
    if scitrace is None:
        msgs.info("Not performing extraction for science frame"+msgs.newline()+fitsdict['filename'][scidx[0]])
        debugger.set_trace()
        #continue
    ###############
    # Finalize the Sky Background image
    if settings.argflag['reduce']['skysub']['perform'] & (scitrace['nobj'] > 0):
        # Perform an iterative background/science extraction
        msgs.info("Finalizing the sky background image")
        trcmask = scitrace['object'].sum(axis=2)
        trcmask[np.where(trcmask>0.0)] = 1.0
        bgframe = bg_subtraction(slf, det, sciframe, modelvarframe, crmask, tracemask=trcmask)
        # Redetermine the variance frame based on the new sky model
        modelvarframe = variance_frame(slf, det, sciframe, scidx, fitsdict, skyframe=bgframe)
        # Save
        if not standard:
            slf._modelvarframe[det-1] = modelvarframe
            slf._bgframe[det-1] = bgframe

    ###############
    # Flexure down the slit? -- Not currently recommended
    if settings.argflag['reduce']['flexure']['method'] == 'slitcen':
        flex_dict = arwave.flexure_slit(slf, det)
        if not msgs._debug['no_qa']:
            arqa.flexure(slf, det, flex_dict, slit_cen=True)

    # Perform an optimal extraction
    msgs.work("For now, perform extraction -- really should do this after the flexure+heliocentric correction")
    return reduce_frame(slf, sciframe, rawvarframe, modelvarframe, bgframe, scidx, fitsdict, det, crmask, standard=standard)


def reduce_frame(slf, sciframe, rawvarframe, modelvarframe, bgframe, scidx, fitsdict, det, crmask,
                 scitrace=None, standard=False):
    """ Run standard extraction steps on a frame

    Parameters
    ----------
    sciframe : image
      Bias subtracted, trimmed, and flatfielded image
    rawvarframe : ndarray
      Variance array using the raw detector counts
    modelvarframe : ndarray
      Model variance array using the raw detector counts and an image of the sky background frame.
    bgframe : ndarray
      Sky background image
    scidx : int
      Index of the frame
    fitsdict : dict
      Contains relevant information from fits header files
    det : int
      Detector index
    scitrace : list of dict
      List containing dictionaries of the object trace parameters
    standard : bool, optional
      Standard star frame?
    """

    ###############
    # Determine the final trace of the science objects
    if scitrace is None:
        msgs.info("Performing final object trace")
        scitrace = artrace.trace_object(slf, det, sciframe-bgframe, modelvarframe, crmask, doqa=(not standard))
    if standard:
        slf._msstd[det-1]['trace'] = scitrace
        specobjs = arspecobj.init_exp(slf, scidx, det, fitsdict, scitrace, objtype='standard')
        slf._msstd[det-1]['spobjs'] = specobjs
    else:
        slf._scitrace[det-1] = scitrace
        # Generate SpecObjExp list
        specobjs = arspecobj.init_exp(slf, scidx, det, fitsdict, scitrace, objtype='science')
        slf._specobjs[det-1] = specobjs

    ###############
    # Extract
    noobj = True
    for sl in range(len(scitrace)):
        if scitrace[sl]['nobj'] != 0:
            noobj = False
    if noobj is True:
        msgs.warn("No objects to extract for science frame"+msgs.newline()+fitsdict['filename'][scidx])
        return True

    # Boxcar
    msgs.info("Performing boxcar extraction")
    bgcorr_box = arextract.boxcar(slf, det, specobjs, sciframe-bgframe,
                                  rawvarframe, bgframe, crmask, scitrace)

    # Optimal
    if not standard:
        msgs.info("Attempting optimal extraction with model profile")
        arextract.obj_profiles(slf, det, specobjs, sciframe-bgframe-bgcorr_box,
                               modelvarframe, bgframe+bgcorr_box, crmask, scitrace)
        newvar = arextract.optimal_extract(slf, det, specobjs, sciframe-bgframe-bgcorr_box,
                                           modelvarframe, bgframe+bgcorr_box, crmask, scitrace)
        msgs.work("Should update variance image (and trace?) and repeat")
        #
        arextract.obj_profiles(slf, det, specobjs, sciframe-bgframe-bgcorr_box,
                               newvar, bgframe+bgcorr_box, crmask, scitrace)
        finalvar = arextract.optimal_extract(slf, det, specobjs, sciframe-bgframe-bgcorr_box,
                                             newvar, bgframe+bgcorr_box, crmask, scitrace)
        slf._modelvarframe[det-1] = finalvar.copy()

    # Flexure correction?
    if settings.argflag['reduce']['flexure']['perform'] and (not standard):
        if settings.argflag['reduce']['flexure']['method'] is not None:
            flex_dict = arwave.flexure_obj(slf, det)
            arqa.flexure(slf, det, flex_dict)

    # Correct Earth's motion
    if (settings.argflag['reduce']['calibrate']['refframe'] in ['heliocentric', 'barycentric']) and \
       (settings.argflag['reduce']['calibrate']['wavelength'] != "pixel"):
        if settings.argflag['science']['extraction']['reuse']:
            msgs.warn("{0:s} correction will not be applied if an extracted science frame exists, and is used".format(settings.argflag['reduce']['calibrate']['refframe']))
        msgs.info("Performing a {0:s} correction".format(settings.argflag['reduce']['calibrate']['refframe']))
        # Load the header for the science frame
        arwave.geomotion_correct(slf, det, fitsdict)
    else:
        msgs.info("A heliocentric correction will not be performed")

    # Final
    if not standard:
        slf._bgframe[det-1] += bgcorr_box
    # Return
    return True


def slit_pixels(slf, frameshape, det):
    """ Generate an image indicating the slit associated with each pixel.

    Parameters
    ----------
    slf : class
      Science Exposure Class
    frameshape : tuple
      A two element tuple providing the shape of a trace frame.
    det : int
      Detector index

    Returns
    -------
    msordloc : ndarray
      An image assigning each pixel to a slit number. A zero value indicates
      that this pixel does not belong to any slit.
    """

    from pypit import arcytrace
    nslits = slf._lordloc[det - 1].shape[1]
    msordloc = np.zeros(frameshape)
    for o in range(nslits):
        lordloc = slf._lordloc[det - 1][:, o]
        rordloc = slf._rordloc[det - 1][:, o]
        ordloc = arcytrace.locate_order(lordloc, rordloc, frameshape[0], frameshape[1],
                                        settings.argflag['trace']['slits']['pad'])
        word = np.where(ordloc != 0)
        if word[0].size == 0:
            msgs.warn("There are no pixels in slit {0:d}".format(o + 1))
            continue
        msordloc[word] = o + 1
    return msordloc


def slit_profile(slf, mstrace, det, ntcky=None):
    """ Generate an image of the spatial slit profile.

    Parameters
    ----------
    slf : class
      Science Exposure Class
    mstrace : ndarray
      Master trace frame that is used to trace the slit edges.
    det : int
      Detector index
    ntcky : int
      Number of bspline knots in the spectral direction.

    Returns
    -------
    slit_profile : ndarray
      An image containing the slit profile
    mstracenrm : ndarray
      The input trace frame, normalized by the blaze function (but still contains the slit profile)
    msblaze : ndarray
      A model of the blaze function of each slit
    blazeext : ndarray
      The blaze function extracted down the centre of the slit
    """
    dnum = settings.get_dnum(det)
    nslits = slf._lordloc[det - 1].shape[1]

    # First, determine the relative scale of each amplifier (assume amplifier 1 has a scale of 1.0)
    if (settings.spect[dnum]['numamplifiers'] > 1) & (nslits > 1):
        sclframe = get_ampscale(slf, det, mstrace)
        # Divide the master flat by the relative scale frame
        mstrace /= sclframe

    mstracenrm = mstrace.copy()
    msblaze = np.ones_like(slf._lordloc[det - 1])
    blazeext = np.ones_like(slf._lordloc[det - 1])
    slit_profiles = np.ones_like(mstrace)
    # Set the number of knots in the spectral direction
    if ntcky is None:
        if settings.argflag["reduce"]["flatfield"]["method"] == "bspline":
            ntcky = settings.argflag["reduce"]["flatfield"]["params"][0]
            if settings.argflag["reduce"]["flatfield"]["params"][0] < 1.0:
                ntcky = int(1.0/ntcky)+0.5
        else:
            ntcky = 20
    else:
        if ntcky < 1.0:
            ntcky = int(1.0 / ntcky) + 0.5
    msgs.work("Multiprocess this step")
    for o in range(nslits):
        if settings.argflag["reduce"]["slitprofile"]["perform"]:
            msgs.info("Deriving the spatial profile and blaze function of slit {0:d}".format(o+1))
        else:
            msgs.info("Deriving the blaze function of slit {0:d}".format(o + 1))
        lordloc = slf._lordloc[det - 1][:, o]
        rordloc = slf._rordloc[det - 1][:, o]
        word = np.where(slf._slitpix[det - 1] == o+1)
        if word[0].size <= (ntcky+1)*(2*slf._pixwid[det - 1][o]+1):
            msgs.warn("There are not enough pixels in slit {0:d}".format(o+1))
            continue
        spatval = (word[1] - lordloc[word[0]])/(rordloc[word[0]] - lordloc[word[0]])
        specval = slf._tilts[det-1][word]
        fluxval = mstrace[word]
        ntckx = 2*slf._pixwid[det - 1][o]
        if not settings.argflag["reduce"]["slitprofile"]["perform"]:
            # The slit profile is not needed, so just do the quickest possible fit
            ntckx = 3
        tckx = np.linspace(np.min(spatval), np.max(spatval), ntckx)
        tcky = np.linspace(np.min(specval), np.max(specval), ntcky)

        # Derive the blaze function
        wsp = np.where((spatval > 0.25) & (spatval < 0.75))
        srt = np.argsort(specval[wsp])
        xb, xe = min(specval[wsp][srt][0], tcky[0]), max(specval[wsp][srt][-1], tcky[-1])
        mask, blzspl = arutils.robust_polyfit(specval[wsp][srt], fluxval[wsp][srt], 3, function='bspline',
                                              sigma=5., maxone=False, xmin=xb, xmax=xe, knots=tcky[1:-1])
        blz_flat = arutils.func_val(blzspl, specval, 'bspline')
        msblaze[:, o] = arutils.func_val(blzspl, np.linspace(0.0, 1.0, msblaze.shape[0]), 'bspline')
        blazeext[:, o] = mstrace[(np.arange(mstrace.shape[0]), np.round(0.5*(lordloc+rordloc)).astype(np.int),)]
        # Calculate the slit profile
        sprof_fit = fluxval / (blz_flat + (blz_flat == 0.0))
        srt = np.argsort(spatval)
        xb, xe = min(spatval[srt][0], tckx[0]), max(spatval[srt][-1], tckx[-1])
        mask, sltspl = arutils.robust_polyfit(spatval[srt], sprof_fit[srt], 3, function='bspline',
                                              sigma=5., maxone=False, xmin=xb, xmax=xe, knots=tckx[1:-1])
        slt_flat = arutils.func_val(sltspl, spatval, 'bspline')
        modvals = blz_flat * slt_flat
        # Normalize to the value at the centre of the slit
        nrmvals = blz_flat * arutils.func_val(sltspl, 0.5, 'bspline')
        if settings.argflag["reduce"]["slitprofile"]["perform"]:
            # Leave slit_profiles as ones if the slitprofile is not being determined, otherwise, set the model.
            slit_profiles[word] = modvals/nrmvals
        mstracenrm[word] /= nrmvals
        if msgs._debug['slit_profile'] and o == 30:
            debugger.set_trace()
            model = np.zeros_like(mstrace)
            model[word] = modvals
            diff = mstrace - model
            import astropy.io.fits as pyfits
            hdu = pyfits.PrimaryHDU(mstrace)
            hdu.writeto("mstrace_{0:02d}.fits".format(det), overwrite=True)
            hdu = pyfits.PrimaryHDU(model)
            hdu.writeto("model_{0:02d}.fits".format(det), overwrite=True)
            hdu = pyfits.PrimaryHDU(diff)
            hdu.writeto("diff_{0:02d}.fits".format(det), overwrite=True)
    # Return
    return slit_profiles, mstracenrm, msblaze, blazeext


def sn_frame(slf, sciframe, idx):

    # Dark Current noise
    dnoise = settings.spect['det']['darkcurr'] * float(slf._fitsdict["exptime"][idx])/3600.0
    # The effective read noise
    rnoise = np.sqrt(settings.spect['det']['ronoise']**2 + (0.5*settings.spect['det']['gain'])**2)
    errframe = np.abs(sciframe) + rnoise + dnoise
    # If there are negative pixels, mask them as bad pixels
    w = np.where(errframe <= 0.0)
    if w[0].size != 0:
        msgs.warn("The error frame is negative for {0:d} pixels".format(w[0].size)+msgs.newline()+"Are you sure the bias frame is correct?")
        msgs.info("Masking these {0:d} pixels".format(w[0].size))
        errframe[w]  = 0.0
        slf._bpix[w] = 1.0
    w = np.where(errframe > 0.0)
    snframe = np.zeros_like(sciframe)
    snframe[w] = sciframe[w]/np.sqrt(errframe[w])
    return snframe


def lacosmic(slf, fitsdict, det, sciframe, scidx, maxiter=1, grow=1.5, maskval=-999999.9,
             simple_var=False):
    """
    Identify cosmic rays using the L.A.Cosmic algorithm
    U{http://www.astro.yale.edu/dokkum/lacosmic/}
    (article : U{http://arxiv.org/abs/astro-ph/0108003})
    This routine is mostly courtesy of Malte Tewes

    :param grow: Once CRs are identified, grow each CR detection by all pixels within this radius
    :return: mask of cosmic rays (0=no CR, 1=CR)
    """
    from pypit import arcyutils
    from pypit import arcyproc
    dnum = settings.get_dnum(det)

    msgs.info("Detecting cosmic rays with the L.A.Cosmic algorithm")
    msgs.work("Include these parameters in the settings files to be adjusted by the user")
    sigclip = 5.0
    sigfrac = 0.3
    objlim  = 5.0
    # Set the settings
    scicopy = sciframe.copy()
    crmask = np.cast['bool'](np.zeros(sciframe.shape))
    sigcliplow = sigclip * sigfrac

    # Determine if there are saturated pixels
    satpix = np.zeros_like(sciframe)
    satlev = settings.spect[dnum]['saturation']*settings.spect[dnum]['nonlinear']
    wsat = np.where(sciframe >= satlev)
    if wsat[0].size == 0: satpix = None
    else:
        satpix[wsat] = 1.0
        satpix = np.cast['bool'](satpix)

    # Define the kernels
    laplkernel = np.array([[0.0, -1.0, 0.0], [-1.0, 4.0, -1.0], [0.0, -1.0, 0.0]])  # Laplacian kernal
    growkernel = np.ones((3,3))
    for i in range(1, maxiter+1):
        msgs.info("Convolving image with Laplacian kernel")
        # Subsample, convolve, clip negative values, and rebin to original size
        #set_trace()
        subsam = arutils.subsample(scicopy)
        conved = signal.convolve2d(subsam, laplkernel, mode="same", boundary="symm")
        cliped = conved.clip(min=0.0)
        lplus = arutils.rebin(cliped, np.array(cliped.shape)/2.0)

        msgs.info("Creating noise model")
        # Build a custom noise map, and compare  this to the laplacian
        m5 = ndimage.filters.median_filter(scicopy, size=5, mode='mirror')
        if simple_var:
            noise = np.sqrt(np.abs(m5)) #variance_frame(slf, det, m5, scidx, fitsdict))
        else:
            noise = np.sqrt(variance_frame(slf, det, m5, scidx, fitsdict))
        msgs.info("Calculating Laplacian signal to noise ratio")

        # Laplacian S/N
        s = lplus / (2.0 * noise)  # Note that the 2.0 is from the 2x2 subsampling

        # Remove the large structures
        sp = s - ndimage.filters.median_filter(s, size=5, mode='mirror')

        msgs.info("Selecting candidate cosmic rays")
        # Candidate cosmic rays (this will include HII regions)
        candidates = sp > sigclip
        nbcandidates = np.sum(candidates)

        msgs.info("{0:5d} candidate pixels".format(nbcandidates))

        # At this stage we use the saturated stars to mask the candidates, if available :
        if satpix is not None:
            msgs.info("Masking saturated pixels")
            candidates = np.logical_and(np.logical_not(satpix), candidates)
            nbcandidates = np.sum(candidates)

            msgs.info("{0:5d} candidate pixels not part of saturated stars".format(nbcandidates))

        msgs.info("Building fine structure image")

        # We build the fine structure image :
        m3 = ndimage.filters.median_filter(scicopy, size=3, mode='mirror')
        m37 = ndimage.filters.median_filter(m3, size=7, mode='mirror')
        f = m3 - m37
        f /= noise
        f = f.clip(min=0.01)

        msgs.info("Removing suspected compact bright objects")

        # Now we have our better selection of cosmics :
        cosmics = np.logical_and(candidates, sp/f > objlim)
        nbcosmics = np.sum(cosmics)

        msgs.info("{0:5d} remaining candidate pixels".format(nbcosmics))

        # What follows is a special treatment for neighbors, with more relaxed constains.

        msgs.info("Finding neighboring pixels affected by cosmic rays")

        # We grow these cosmics a first time to determine the immediate neighborhod  :
        growcosmics = np.cast['bool'](signal.convolve2d(np.cast['float32'](cosmics), growkernel, mode="same", boundary="symm"))

        # From this grown set, we keep those that have sp > sigmalim
        # so obviously not requiring sp/f > objlim, otherwise it would be pointless
        growcosmics = np.logical_and(sp > sigclip, growcosmics)

        # Now we repeat this procedure, but lower the detection limit to sigmalimlow :

        finalsel = np.cast['bool'](signal.convolve2d(np.cast['float32'](growcosmics), growkernel, mode="same", boundary="symm"))
        finalsel = np.logical_and(sp > sigcliplow, finalsel)

        # Unmask saturated pixels:
        if satpix != None:
            msgs.info("Masking saturated stars")
            finalsel = np.logical_and(np.logical_not(satpix), finalsel)

        ncrp = np.sum(finalsel)

        msgs.info("{0:5d} pixels detected as cosmics".format(ncrp))

        # We find how many cosmics are not yet known :
        newmask = np.logical_and(np.logical_not(crmask), finalsel)
        nnew = np.sum(newmask)

        # We update the mask with the cosmics we have found :
        crmask = np.logical_or(crmask, finalsel)

        msgs.info("Iteration {0:d} -- {1:d} pixels identified as cosmic rays ({2:d} new)".format(i, ncrp, nnew))
        if ncrp == 0: break
    # Additional algorithms (not traditionally implemented by LA cosmic) to remove some false positives.
    msgs.work("The following algorithm would be better on the rectified, tilts-corrected image")
    filt  = ndimage.sobel(sciframe, axis=1, mode='constant')
    filty = ndimage.sobel(filt/np.sqrt(np.abs(sciframe)), axis=0, mode='constant')
    filty[np.where(np.isnan(filty))]=0.0
    sigimg  = arcyproc.cr_screen(filty,0.0)
    sigsmth = ndimage.filters.gaussian_filter(sigimg,1.5)
    sigsmth[np.where(np.isnan(sigsmth))]=0.0
    sigmask = np.cast['bool'](np.zeros(sciframe.shape))
    sigmask[np.where(sigsmth>sigclip)] = True
    crmask = np.logical_and(crmask, sigmask)
    msgs.info("Growing cosmic ray mask by 1 pixel")
    crmask = arcyutils.grow_masked(crmask.astype(np.float), grow, 1.0)
    return crmask


def gain_frame(slf, det):
    """ Generate a gain image from the spect dict

    Parameters
    ----------
    slf
    det

    Returns
    -------
    gain_img : ndarray

    """
    dnum = settings.get_dnum(det)

    # Loop on amplifiers
    gain_img = np.zeros_like(slf._datasec[det-1])
    for ii in range(settings.spect[dnum]['numamplifiers']):
        amp = ii+1
        try:
            amppix = slf._datasec[det-1] == amp
            gain_img[amppix] = settings.spect[dnum]['gain'][amp - 1]
        except IndexError:
            debugger.set_trace()
    # Return
    return gain_img


def rn_frame(slf, det):
    """ Generate a RN image

    Parameters
    ----------
    slf
    det

    Returns
    -------
    rn_img : ndarray
      Read noise *variance* image (i.e. RN**2)
    """
    dnum = settings.get_dnum(det)

    # Loop on amplifiers
    rnimg = np.zeros_like(slf._datasec[det-1])
    for ii in range(settings.spect[dnum]['numamplifiers']):
        amp = ii+1
        amppix = slf._datasec[det-1] == amp
        rnimg[amppix] = (settings.spect[dnum]['ronoise'][ii]**2 +
                         (0.5*settings.spect[dnum]['gain'][ii])**2)
    # Return
    return rnimg


def sub_overscan(frame, det):
    """
    Subtract overscan

    Parameters
    ----------
    frame : ndarray
      frame which should have the overscan region subtracted
    det : int
      Detector Index

    Returns
    -------
    frame : ndarray
      The input frame with the overscan region subtracted
    """
    dnum = settings.get_dnum(det)

    for i in range(settings.spect[dnum]['numamplifiers']):
        # Determine the section of the chip that contains data
        datasec = "datasec{0:02d}".format(i+1)
        dx0, dx1 = settings.spect[dnum][datasec][0][0], settings.spect[dnum][datasec][0][1]
        dy0, dy1 = settings.spect[dnum][datasec][1][0], settings.spect[dnum][datasec][1][1]
        if dx0 < 0: dx0 += frame.shape[0]
        if dx1 <= 0: dx1 += frame.shape[0]
        if dy0 < 0: dy0 += frame.shape[1]
        if dy1 <= 0: dy1 += frame.shape[1]
        xds = np.arange(dx0, dx1)
        yds = np.arange(dy0, dy1)
        # Determine the section of the chip that contains the overscan region
        oscansec = "oscansec{0:02d}".format(i+1)
        ox0, ox1 = settings.spect[dnum][oscansec][0][0], settings.spect[dnum][oscansec][0][1]
        oy0, oy1 = settings.spect[dnum][oscansec][1][0], settings.spect[dnum][oscansec][1][1]
        if ox0 < 0: ox0 += frame.shape[0]
        if ox1 <= 0: ox1 += min(frame.shape[0], dx1)  # Truncate to datasec
        if oy0 < 0: oy0 += frame.shape[1]
        if oy1 <= 0: oy1 += min(frame.shape[1], dy1)  # Truncate to datasec
        xos = np.arange(ox0, ox1)
        yos = np.arange(oy0, oy1)
        w = np.ix_(xos, yos)
        oscan = frame[w]
        # Make sure the overscan section has at least one side consistent with datasec
        if dx1-dx0 == ox1-ox0:
            osfit = np.median(oscan, axis=1)  # Mean was hit by CRs
        elif dy1-dy0 == oy1-oy0:
            osfit = np.median(oscan, axis=0)
        elif settings.argflag['reduce']['overscan']['method'].lower() == "median":
            osfit = np.median(oscan)
        else:
            msgs.error("Overscan sections do not match amplifier sections for amplifier {0:d}".format(i+1))
        # Fit/Model the overscan region
        if settings.argflag['reduce']['overscan']['method'].lower() == "polynomial":
            c = np.polyfit(np.arange(osfit.size), osfit, settings.argflag['reduce']['overscan']['params'][0])
            ossub = np.polyval(c, np.arange(osfit.size))#.reshape(osfit.size,1)
        elif settings.argflag['reduce']['overscan']['method'].lower() == "savgol":
            ossub = savgol_filter(osfit, settings.argflag['reduce']['overscan']['params'][1], settings.argflag['reduce']['overscan']['params'][0])
        elif settings.argflag['reduce']['overscan']['method'].lower() == "median":  # One simple value
            ossub = osfit * np.ones(1)
        else:
            msgs.warn("Overscan subtraction method {0:s} is not implemented".format(settings.argflag['reduce']['overscan']['method']))
            msgs.info("Using a linear fit to the overscan region")
            c = np.polyfit(np.arange(osfit.size), osfit, 1)
            ossub = np.polyval(c, np.arange(osfit.size))#.reshape(osfit.size,1)
        # Determine the section of the chip that contains data for this amplifier
        wd = np.ix_(xds, yds)
        ossub = ossub.reshape(osfit.size, 1)
        if wd[0].shape[0] == ossub.shape[0]:
            frame[wd] -= ossub
        elif wd[1].shape[1] == ossub.shape[0]:
            frame[wd] -= ossub.T
        elif settings.argflag['reduce']['overscan']['method'].lower() == "median":
            frame[wd] -= osfit
        else:
            msgs.error("Could not subtract bias from overscan region --"+msgs.newline()+"size of extracted regions does not match")
    # Return
    del xds, yds, xos, yos, oscan
    return frame


def trim(frame, det):
    dnum = settings.get_dnum(det)
    for i in range(settings.spect[dnum]['numamplifiers']):
        datasec = "datasec{0:02d}".format(i+1)
        x0, x1 = settings.spect[dnum][datasec][0][0], settings.spect[dnum][datasec][0][1]
        y0, y1 = settings.spect[dnum][datasec][1][0], settings.spect[dnum][datasec][1][1]
        if x0 < 0:
            x0 += frame.shape[0]
        if x1 <= 0:
            x1 += frame.shape[0]
        if y0 < 0:
            y0 += frame.shape[1]
        if y1 <= 0:
            y1 += frame.shape[1]
        if i == 0:
            xv = np.arange(x0, x1)
            yv = np.arange(y0, y1)
        else:
            xv = np.unique(np.append(xv, np.arange(x0, x1)))
            yv = np.unique(np.append(yv, np.arange(y0, y1)))
    # Construct and array with the rows and columns to be extracted
    w = np.ix_(xv, yv)
#	if len(file.shape) == 2:
#		trimfile = file[w]
#	elif len(file.shape) == 3:
#		trimfile = np.zeros((w[0].shape[0],w[1].shape[1],file.shape[2]))
#		for f in range(file.shape[2]):
#			trimfile[:,:,f] = file[:,:,f][w]
#	else:
#		msgs.error("Cannot trim {0:d}D frame".format(int(len(file.shape))))
    try:
        return frame[w]
    except:
        msgs.bug("Odds are datasec is set wrong. Maybe due to transpose")
        debugger.set_trace()
        msgs.error("Cannot trim file")


def variance_frame(slf, det, sciframe, idx, fitsdict=None, skyframe=None, objframe=None):
    """ Calculate the variance image including detector noise
    Parameters
    ----------
    fitsdict : dict, optional
      Contains relevant information from fits header files
    objframe : ndarray, optional
      Model of object counts
    Returns
    -------
    variance image : ndarray
    """
    dnum = settings.get_dnum(det)

    # The effective read noise (variance image)
    rnoise = rn_frame(slf, det)
    if skyframe is not None:
        if objframe is None:
            objframe = np.zeros_like(skyframe)
        varframe = np.abs(skyframe + objframe - np.sqrt(2)*np.sqrt(rnoise)) + rnoise
        return varframe
    else:
        scicopy = sciframe.copy()
        # Dark Current noise
        dnoise = (settings.spect[dnum]['darkcurr'] *
                  float(fitsdict["exptime"][idx])/3600.0)
        # Return
        return np.abs(scicopy) + rnoise + dnoise
