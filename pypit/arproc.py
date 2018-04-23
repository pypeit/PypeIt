from __future__ import (print_function, absolute_import, division, unicode_literals)

import time
import inspect

import numpy as np

from scipy import signal, ndimage, interpolate

from matplotlib import pyplot as plt
from matplotlib import gridspec, font_manager

from astropy import units
from astropy.io import fits

from pypit import msgs

from pypit import arextract
from pypit import arlris
from pypit import ardeimos
from pypit import armsgs
from pypit import artrace
from pypit import arutils
from pypit import arparse as settings
from pypit import arspecobj
from pypit import arqa
from pypit import arpca
from pypit import arwave

from pypit import arcytrace
from pypit import arcyutils
from pypit import arcyproc

from pypit import ardebug as debugger

try:
    basestring
except NameError:
    basestring = str

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
    oversampling_factor = 3 # should be an integer according to the description in object_profile()
    xedges, modvals = object_profile(slf, sciframe, slitn, det, refine=refine, factor=oversampling_factor)
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
        nl_index = np.max(wl[0])
        # Calculate nl, defined as:
        # "number of pixels from the left slit edge to use as background pixels",
        # which is just nl_index with the sampling factor taken out
        nl_index_origscale = int(nl_index/oversampling_factor+0.5)
        nl = nl_index_origscale
    if wr[0].size != 0:
        # This is the index of the first time where the object profile
        # no longer decreases as you move towards the slit edge
        nr_index = np.min(wr[0])
        # Calculate nr, defined as:
        # "number of pixels from the right slit edge to use as background pixels",
        # which is npix minus nr_index with the sampling factor taken out
        nr_index_origscale = int(nr_index/oversampling_factor+0.5)
        nr = npix - nr_index_origscale
    if nl+nr < 5:
        msgs.warn("The object profile appears to extrapolate to the edge of the slit")
        msgs.info("A background subtraction will not be performed for slit {0:d}".format(slitn+1))
        nl, nr = 0, 0
        return np.zeros_like(sciframe), nl, nr
    # Find background pixels and fit
    wbgpix_spatval = np.where((spatval <= float(nl)/npix) | (spatval >= float(npix-nr)/npix)) # this cannot be used to index the 2D array tilts
    wbgpix = (word[0][wbgpix_spatval], word[1][wbgpix_spatval]) # this may be approproate for indexing the 2D array tilts
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


def bias_subtract(rawframe, msbias):
    if type(msbias) is np.ndarray:
        msgs.info("Subtracting bias image from raw frame")
        newframe = rawframe-msbias  # Subtract the master bias frame
    elif isinstance(msbias, basestring):
        if msbias == "overscan":
            #def sub_overscan(frame, numamplifiers, datasec, oscansec, settings=None):
            newframe = sub_overscan(rawframe)
        else:
            msgs.error("Could not subtract bias level with the input bias approach") #when loading {0:s} frames".format(frametype))
    return newframe


def bg_subtraction(slf, det, sciframe, varframe, crpix, tracemask=None,
                   rejsigma=3.0, maskval=-999999.9):
    """ Extract a science target and background flux
    :param slf:
    :param sciframe:
    :param varframe:
    :return:
    """
    # Set some starting parameters (maybe make these available to the user)
    msgs.work("Should these parameters be made available to the user?")
    polyorder, repeat = 5, 1
    # Begin the algorithm
    errframe = np.sqrt(varframe)
    # Find which pixels are within the order edges
    msgs.info("Identifying pixels within each order")
#    print('calling order_pixels')
#    t = time.clock()
#    _ordpix = arcyutils.order_pixels(slf._pixlocn[det-1],
#                                     slf._lordloc[det-1]*0.95+slf._rordloc[det-1]*0.05,
#                                     slf._lordloc[det-1]*0.05+slf._rordloc[det-1]*0.95)
#    print('Old order_pixels: {0} seconds'.format(time.clock() - t))
#    t = time.clock()
    ordpix = new_order_pixels(slf._pixlocn[det-1],
                              slf._lordloc[det-1]*0.95+slf._rordloc[det-1]*0.05,
                              slf._lordloc[det-1]*0.05+slf._rordloc[det-1]*0.95)
#    print('New order_pixels: {0} seconds'.format(time.clock() - t))
#    assert np.sum(_ordpix != ordpix) == 0, \
#                    'Difference between old and new order_pixels'

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
    xargsrt = np.argsort(xvpix, kind='mergesort')
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
    if settings.argflag['reduce']['skysub']['method'].lower() == 'bspline':
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

    *** CAUTION ***  This function might be deprecated.

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

#    print('calling order_pixels')
#    t = time.clock()
#    _ordpix = arcyutils.order_pixels(slf._pixlocn[det-1], slf._lordloc[det-1], slf._rordloc[det-1])
#    print('Old order_pixels: {0} seconds'.format(time.clock() - t))
#    t = time.clock()
    ordpix = new_order_pixels(slf._pixlocn[det-1], slf._lordloc[det-1], slf._rordloc[det-1])
#    print('New order_pixels: {0} seconds'.format(time.clock() - t))
#    assert np.sum(_ordpix != ordpix) == 0, \
#                    'Difference between old and new order_pixels'

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
#        arqa.plot_orderfits(slf, msblaze, flat_ext1d, desc=plotdesc, textplt="Order")
        artrace.plot_orderfits(slf, msblaze, flat_ext1d, desc=plotdesc, textplt="Order")
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


def wrapper_get_datasec_trimmed(slf, fitsdict, det, scidx):
    dnum = settings.get_dnum(det)
    spectrograph = settings.argflag['run']['spectrograph']
    scifile = fitsdict['directory'][scidx]+fitsdict['filename'][scidx]
    numamplifiers = settings.spect[dnum]['numamplifiers']

    if spectrograph in ['keck_lris_blue', 'keck_lris_red', 'keck_deimos']:
        # Grab
        datasec, oscansec, naxis0, naxis1 = get_datasec(spectrograph, scifile, numamplifiers, det=det)
        # Fill (for backwards compatability)
        for kk in range(numamplifiers):
            datasec = "datasec{0:02d}".format(kk+1)
            settings.spect[dnum][datasec] = datasec[kk]
            oscansec = "oscansec{0:02d}".format(kk+1)
            settings.spect[dnum][oscansec] = oscansec[kk]
        fitsdict['naxis0'][scidx] = naxis0
        fitsdict['naxis1'][scidx] = naxis1


def get_datasec(spectrograph, scifile, numamplifiers, det=None):
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
    # Get naxis0, naxis1, datasec, oscansec, ampsec for specific instruments
    datasec, oscansec, naxis0, naxis1 = [], [], 0, 0
    if spectrograph in ['keck_lris_blue', 'keck_lris_red']:
        msgs.info("Parsing datasec and oscansec from headers")
        temp, head0, secs = arlris.read_lris(scifile, det)
    elif spectrograph in ['keck_deimos']:
        msgs.info("Parsing datasec and oscansec from headers")
        # TODO -- This should not be reading in the full DEIMOS detector but currently is
        temp, head0, secs = ardeimos.read_deimos(scifile)
    else:  # Other instruments are set in their settings file
        msgs.warn("Should not have called get_datasec!")
        return datasec, oscansec, naxis0, naxis1

    naxis0 = temp.shape[0]
    naxis1 = temp.shape[1]
    for kk in range(numamplifiers):
        #datasec = "datasec{0:02d}".format(kk+1)
        #settings.spect[dnum][datasec] = settings.load_sections(secs[0][kk], fmt_iraf=False)
        datasec.append(settings.load_sections(secs[0][kk], fmt_iraf=False))
        #oscansec = "oscansec{0:02d}".format(kk+1)
        #settings.spect[dnum][oscansec] = settings.load_sections(secs[1][kk], fmt_iraf=False)
        oscansec.append(settings.load_sections(secs[1][kk], fmt_iraf=False))
    # Return
    return datasec, oscansec, naxis0, naxis1

def pix_to_amp(spectrograph, naxis0, naxis1, datasec=None, scifile=None, det=None):
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
    return retarr[w]


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


def flexure_qa(slf, det, flex_list, slit_cen=False):
    """ QA on flexure measurement

    Parameters
    ----------
    slf
    det
    flex_list : list
      list of dict containing flexure results
    slit_cen : bool, optional
      QA on slit center instead of objects

    Returns
    -------

    """
    plt.rcdefaults()
    plt.rcParams['font.family']= 'times new roman'

    # Grab the named of the method
    method = inspect.stack()[0][3]
    #
    for sl in range(len(slf._specobjs[det-1])):
        # Setup
        if slit_cen:
            nobj = 1
            ncol = 1
        else:
            nobj = len(slf._specobjs[det-1][sl])
            if nobj == 0:
                continue
            ncol = min(3, nobj)
        #
        nrow = nobj // ncol + ((nobj % ncol) > 0)

        # Get the flexure dictionary
        flex_dict = flex_list[sl]

        # Outfile
        outfile = arqa.set_qa_filename(slf._basename, method+'_corr', det=det,
                                       slit=slf._specobjs[det-1][sl][0].slitid)

        plt.figure(figsize=(8, 5.0))
        plt.clf()
        gs = gridspec.GridSpec(nrow, ncol)

        # Correlation QA
        for o in range(nobj):
            ax = plt.subplot(gs[o//ncol, o % ncol])
            # Fit
            fit = flex_dict['polyfit'][o]
            xval = np.linspace(-10., 10, 100) + flex_dict['corr_cen'][o] #+ flex_dict['shift'][o]
            #model = (fit[2]*(xval**2.))+(fit[1]*xval)+fit[0]
            model = arutils.func_val(fit, xval, 'polynomial')
            mxmod = np.max(model)
            ylim = [np.min(model/mxmod), 1.3]
            ax.plot(xval-flex_dict['corr_cen'][o], model/mxmod, 'k-')
            # Measurements
            ax.scatter(flex_dict['subpix'][o]-flex_dict['corr_cen'][o],
                       flex_dict['corr'][o]/mxmod, marker='o')
            # Final shift
            ax.plot([flex_dict['shift'][o]]*2, ylim, 'g:')
            # Label
            if slit_cen:
                ax.text(0.5, 0.25, 'Slit Center', transform=ax.transAxes, size='large', ha='center')
            else:
                ax.text(0.5, 0.25, '{:s}'.format(slf._specobjs[det-1][sl][o].idx), transform=ax.transAxes, size='large', ha='center')
            ax.text(0.5, 0.15, 'flex_shift = {:g}'.format(flex_dict['shift'][o]),
                    transform=ax.transAxes, size='large', ha='center')#, bbox={'facecolor':'white'})
            # Axes
            ax.set_ylim(ylim)
            ax.set_xlabel('Lag')

        # Finish
        plt.tight_layout(pad=0.2, h_pad=0.0, w_pad=0.0)
        plt.savefig(outfile, dpi=600)
        plt.close()

        # Sky line QA (just one object)
        if slit_cen:
            o = 0
        else:
            o = 0
            specobj = slf._specobjs[det-1][sl][o]
        sky_spec = flex_dict['sky_spec'][o]
        arx_spec = flex_dict['arx_spec'][o]

        # Sky lines
        sky_lines = np.array([3370.0, 3914.0, 4046.56, 4358.34, 5577.338, 6300.304,
                  7340.885, 7993.332, 8430.174, 8919.610, 9439.660,
                  10013.99, 10372.88])*units.AA
        dwv = 20.*units.AA
        gdsky = np.where((sky_lines > sky_spec.wvmin) & (sky_lines < sky_spec.wvmax))[0]
        if len(gdsky) == 0:
            msgs.warn("No sky lines for Flexure QA")
            return
        if len(gdsky) > 6:
            idx = np.array([0, 1, len(gdsky)//2, len(gdsky)//2+1, -2, -1])
            gdsky = gdsky[idx]

        # Outfile
        outfile = arqa.set_qa_filename(slf._basename, method+'_sky', det=det,
                                       slit=slf._specobjs[det-1][sl][0].slitid)
        # Figure
        plt.figure(figsize=(8, 5.0))
        plt.clf()
        nrow, ncol = 2, 3
        gs = gridspec.GridSpec(nrow, ncol)
        if slit_cen:
            plt.suptitle('Sky Comparison for Slit Center', y=1.05)
        else:
            plt.suptitle('Sky Comparison for {:s}'.format(specobj.idx), y=1.05)

        for ii, igdsky in enumerate(gdsky):
            skyline = sky_lines[igdsky]
            ax = plt.subplot(gs[ii//ncol, ii % ncol])
            # Norm
            pix = np.where(np.abs(sky_spec.wavelength-skyline) < dwv)[0]
            f1 = np.sum(sky_spec.flux[pix])
            f2 = np.sum(arx_spec.flux[pix])
            norm = f1/f2
            # Plot
            ax.plot(sky_spec.wavelength[pix], sky_spec.flux[pix], 'k-', label='Obj',
                    drawstyle='steps-mid')
            pix2 = np.where(np.abs(arx_spec.wavelength-skyline) < dwv)[0]
            ax.plot(arx_spec.wavelength[pix2], arx_spec.flux[pix2]*norm, 'r-', label='Arx',
                    drawstyle='steps-mid')
            # Axes
            ax.xaxis.set_major_locator(plt.MultipleLocator(dwv.value))
            ax.set_xlabel('Wavelength')
            ax.set_ylabel('Counts')

        # Legend
        plt.legend(loc='upper left', scatterpoints=1, borderpad=0.3,
                   handletextpad=0.3, fontsize='small', numpoints=1)

        # Finish
        plt.savefig(outfile, dpi=800)
        plt.close()
        #plt.close()

    plt.rcdefaults()

    return


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
        medpix = flxfr[groups == mm]
        if medpix.size == 0:
            profile[mm - 1] = 0.0
        else:
            profile[mm - 1] = np.median(medpix)
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
            word = np.where((slf._slitpix[det - 1] == o + 1) & (slf._scimask[det - 1] == 0))
            if word[0].size == 0:
                msgs.warn("There are no pixels in slit {0:d}".format(o+1))
                continue
            tbgframe, nl, nr = background_subtraction(slf, sciframe, rawvarframe, o, det)
            bgnl[o], bgnr[o] = nl, nr
            bgframe += tbgframe
            if nl == 0 and nr == 0:
                pass
                # If just one slit cannot do sky subtraction, don't do sky subtraction
                # msgs.warn("A sky subtraction will not be performed")
                # skysub = False
                # bgframe = np.zeros_like(sciframe)
                # modelvarframe = rawvarframe.copy()
                # break
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
    extrap_slit = np.zeros(nord)
    for o in range(nord):
        trace, error = artrace.trace_weighted(sciframe-bgframe, slf._lordloc[det-1][:, o], slf._rordloc[det-1][:, o],
                                              mask=slf._scimask[det-1], wght="flux")
        if trace is None:
            extrap_slit[o] = 1
            continue
        # Find only the good pixels
        w = np.where((error != 0.0) & (~np.isnan(error)))
        if w[0].size <= 2*settings.argflag['trace']['object']['order']:
            extrap_slit[o] = 1
            continue
        # Convert the trace locations to be a fraction of the slit length,
        # measured from the left slit edge.
        trace -= slf._lordloc[det-1][:, o]
        trace /= (slf._rordloc[det-1][:, o]-slf._lordloc[det-1][:, o])
        try:
            msk, trccoeff[:, o] = arutils.robust_polyfit(trcxfit[w], trace[w],
                                                     settings.argflag['trace']['object']['order'],
                                                     function=settings.argflag['trace']['object']['function'],
                                                     weights=1.0 / error[w] ** 2, minv=0.0, maxv=nspec-1.0)
        except:
            msgs.info("arproc.reduce_echelle")
            debugger.set_trace()
    refine = 0.0
    if settings.argflag['trace']['object']['method'] == "pca":
        # Identify the orders to be extrapolated during reconstruction
        orders = 1.0 + np.arange(nord)
        msgs.info("Performing a PCA on the object trace")
        ofit = settings.argflag['trace']['object']['params']
        lnpc = len(ofit) - 1
        maskord = np.where(extrap_slit == 1)[0]

        xcen = trcxfit[:, np.newaxis].repeat(nord, axis=1)
        trccen = arutils.func_val(trccoeff, trcxfit, settings.argflag['trace']['object']['function'],
                                  minv=0.0, maxv=nspec-1.0).T
        if np.sum(1.0 - extrap_slit) > ofit[0] + 1:
            fitted, outpar = arpca.basis(xcen, trccen, trccoeff, lnpc, ofit, skipx0=False, mask=maskord,
                                         function=settings.argflag['trace']['object']['function'])
            if not msgs._debug['no_qa']:
#                arqa.pca_plot(slf, outpar, ofit, "Object_Trace", pcadesc="PCA of object trace")
                arpca.pca_plot(slf, outpar, ofit, "Object_Trace", pcadesc="PCA of object trace")
            # Extrapolate the remaining orders requested
            trccen, outpar = arpca.extrapolate(outpar, orders, function=settings.argflag['trace']['object']['function'])
            #refine = trccen-trccen[nspec//2, :].reshape((1, nord))
        else:
            msgs.warn("Could not perform a PCA on the object trace" + msgs.newline() +
                      "Not enough well-traced orders")
            msgs.info("Using direct determination of the object trace instead")
            pass
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

    # Generate an image of pixel weights for each object. Each weight can
    # take any floating point value from 0 to 1 (inclusive). For the rec_obj_img,
    # a weight of 1 means that the pixel is fully contained within the object
    # region, and 0 means that the pixel is fully contained within the background
    # region. The opposite is true for the rec_bg_img array. A pixel that is on
    # the border of object/background is assigned a value between 0 and 1.
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
        artrace.obj_trace_qa(slf, sciframe, trobjl, trobjr, None, det,
                             root="object_trace", normalize=False)

    # Finalize the Sky Background image
    if settings.argflag['reduce']['skysub']['perform'] and (nobj > 0) and skysub:
        msgs.info("Finalizing the sky background image")
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
    # Find objects and estimate their traces
    scitrace = artrace.trace_object(slf, det, sciframe-bgframe, modelvarframe, crmask,
                                    bgreg=20, doqa=False, standard=standard)
    if scitrace is None:
        msgs.info("Not performing extraction for science frame"+msgs.newline()+fitsdict['filename'][scidx[0]])
        debugger.set_trace()
        #continue

    # Make sure that there are objects
    noobj = True
    for sl in range(len(scitrace)):
        if scitrace[sl]['nobj'] != 0:
            noobj = False
    if noobj is True:
        msgs.warn("No objects to extract for science frame" + msgs.newline() + fitsdict['filename'][scidx])
        return True

    ###############
    # Finalize the Sky Background image
    if settings.argflag['reduce']['skysub']['perform']:
        # Perform an iterative background/science extraction
        msgs.info("Finalizing the sky background image")
        # Create a trace mask of the object
        trcmask = np.zeros_like(sciframe)
        for sl in range(len(scitrace)):
            trcmask += scitrace[sl]['object'].sum(axis=2)
        trcmask[np.where(trcmask > 0.0)] = 1.0
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
        #if not msgs._debug['no_qa']:
#        arqa.flexure(slf, det, flex_dict, slit_cen=True)
        flexure_qa(slf, det, flex_dict, slit_cen=True)

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
        scitrace = artrace.trace_object(slf, det, sciframe-bgframe, modelvarframe, crmask,
                                        bgreg=20, doqa=(not standard), standard=standard)
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

        # KBW: Using variance_frame() in arextract leads to a circular
        # import.  I've changed the arextract.optimal_extract() function
        # to return the object model, then the last step of generating
        # the new variance image is done here.

        msgs.info("Attempting optimal extraction with model profile")
        arextract.obj_profiles(slf, det, specobjs, sciframe-bgframe-bgcorr_box,
                               modelvarframe, bgframe+bgcorr_box, crmask, scitrace, doqa=False)
#        newvar = arextract.optimal_extract(slf, det, specobjs, sciframe-bgframe-bgcorr_box,
#                                           modelvarframe, bgframe+bgcorr_box, crmask, scitrace)
        obj_model = arextract.optimal_extract(slf, det, specobjs, sciframe-bgframe-bgcorr_box,
                                              modelvarframe, bgframe+bgcorr_box, crmask, scitrace)
        newvar = variance_frame(slf, det, sciframe-bgframe-bgcorr_box, -1,
                                skyframe=bgframe+bgcorr_box, objframe=obj_model)
        msgs.work("Should update variance image (and trace?) and repeat")
        #
        arextract.obj_profiles(slf, det, specobjs, sciframe-bgframe-bgcorr_box,
                               newvar, bgframe+bgcorr_box, crmask, scitrace)
#        finalvar = arextract.optimal_extract(slf, det, specobjs, sciframe-bgframe-bgcorr_box,
#                                             newvar, bgframe+bgcorr_box, crmask, scitrace)
        obj_model = arextract.optimal_extract(slf, det, specobjs, sciframe-bgframe-bgcorr_box,
                                              newvar, bgframe+bgcorr_box, crmask, scitrace)
        finalvar = variance_frame(slf, det, sciframe-bgframe-bgcorr_box, -1,
                                  skyframe=bgframe+bgcorr_box, objframe=obj_model)
        slf._modelvarframe[det-1] = finalvar.copy()

    # Flexure correction?
    if settings.argflag['reduce']['flexure']['perform'] and (not standard):
        if settings.argflag['reduce']['flexure']['method'] is not None:
            flex_dict = arwave.flexure_obj(slf, det)
            #if not msgs._debug['no_qa']:
#            arqa.flexure(slf, det, flex_dict)
            flexure_qa(slf, det, flex_dict)

    # Correct Earth's motion
    if (settings.argflag['reduce']['calibrate']['refframe'] in ['heliocentric', 'barycentric']) and \
       (settings.argflag['reduce']['calibrate']['wavelength'] != "pixel"):
        if settings.argflag['science']['extraction']['reuse']:
            msgs.warn("{0:s} correction will not be applied if an extracted science frame exists, and is used".format(settings.argflag['reduce']['calibrate']['refframe']))
        if slf._specobjs[det-1] is not None:
            msgs.info("Performing a {0:s} correction".format(settings.argflag['reduce']['calibrate']['refframe']))
            arwave.geomotion_correct(slf, det, fitsdict)
        else:
            msgs.info("There are no objects on detector {0:d} to perform a {1:s} correction".format(
                det, settings.argflag['reduce']['calibrate']['refframe']))
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

    nslits = slf._lordloc[det - 1].shape[1]
    msordloc = np.zeros(frameshape)
    for o in range(nslits):
        lordloc = slf._lordloc[det - 1][:, o]
        rordloc = slf._rordloc[det - 1][:, o]
#        print('calling locate_order')
#        t = time.clock()
#        _ordloc = arcytrace.locate_order(lordloc, rordloc, frameshape[0], frameshape[1],
#                                         settings.argflag['trace']['slits']['pad'])
#        print('Old locate_order: {0} seconds'.format(time.clock() - t))
#        t = time.clock()
        ordloc = new_locate_order(lordloc, rordloc, frameshape[0], frameshape[1],
                                  settings.argflag['trace']['slits']['pad'])
#        print('New locate_order: {0} seconds'.format(time.clock() - t))
#        assert np.sum(_ordloc != ordloc) == 0, \
#                    'Difference between old and new locate_order'
        word = np.where(ordloc != 0)
        if word[0].size == 0:
            msgs.warn("There are no pixels in slit {0:d}".format(o + 1))
            continue
        msordloc[word] = o + 1
    return msordloc


def new_locate_order(lordloc, rordloc, sz_x, sz_y, pad):
    """ Generate a boolean image that identifies which pixels
    belong to the slit associated with the supplied left and
    right slit edges.

    Parameters
    ----------
    lordloc : ndarray
      Location of the left slit edges of 1 slit
    rordloc : ndarray
      Location of the right slit edges of 1 slit
    sz_x : int
      The size of an image in the spectral (0th) dimension
    sz_y : int
      The size of an image in the spatial (1st) dimension
    pad : int
      Additional pixels to pad the left and right slit edges

    Returns
    -------
    orderloc : ndarray
      An image the same size as the input frame, containing values from 0-1.
      0 = pixel is not in the specified slit
      1 = pixel is in the specified slit
    """
    ow = (rordloc-lordloc)/2.0
    oc = (rordloc+lordloc)/2.0
    ymin = (oc-ow).astype(int)-pad
    ymax = (oc+ow).astype(int)+1+pad
    indx = np.invert((ymax < 0) | (ymin >= sz_y))
    ymin[ymin < 0] = 0
    ymax[ymax > sz_y-1] = sz_y-1
    indx &= (ymax > ymin)

    orderloc = np.zeros((sz_x,sz_y), dtype=int)
    for x in np.arange(sz_x)[indx]:
        orderloc[x,ymin[x]:ymax[x]] = 1
    return orderloc



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
    extrap_slit : ndarray
      Mask indicating if a slit is well-determined (0) or poor (1). If the latter, the slit profile
      and blaze function for those slits should be extrapolated or determined from another means
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
    ntcky = int(ntcky)
    # Set the number of knots in the spatial direction
    ntckx = 2 * np.max(slf._pixwid[det - 1])
    if not settings.argflag["reduce"]["slitprofile"]["perform"]:
        # The slit profile is not needed, so just do the quickest possible fit
        ntckx = 3

    extrap_slit = np.zeros(nslits, dtype=np.int)

    # Calculate the slit and blaze profiles
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
            extrap_slit[o] = 1.0
            continue
        spatval = (word[1] - lordloc[word[0]])/(rordloc[word[0]] - lordloc[word[0]])
        specval = slf._tilts[det-1][word]
        fluxval = mstrace[word]

        # Only use pixels where at least half the slit is on the chip
        cordloc = 0.5 * (lordloc[word[0]] + rordloc[word[0]])
        wcchip = ((cordloc > 0.0) & (cordloc < mstrace.shape[1]-1.0))

        # Derive the blaze function
        wsp = np.where((spatval > 0.25) & (spatval < 0.75) & wcchip)
        if wsp[0].size <= (ntcky+1)*(2*slf._pixwid[det - 1][o]+1):
            msgs.warn("There are not enough pixels in slit {0:d}".format(o+1))
            extrap_slit[o] = 1.0
            continue
        if (np.min(word[0]) > 0) or (np.max(word[0]) < mstrace.shape[0]-1):
            extrap_slit[o] = 1.0
        tcky = np.linspace(min(0.0, np.min(specval[wsp])), max(1.0, np.max(specval[wsp])), ntcky)
        tcky = tcky[np.where((tcky > np.min(specval[wsp])) & (tcky < np.max(specval[wsp])))]
        srt = np.argsort(specval[wsp])
        # Only perform a bspline if there are enough pixels for the specified knots
        if tcky.size >= 2:
            yb, ye = min(np.min(specval), tcky[0]), max(np.max(specval), tcky[-1])
            mask, blzspl = arutils.robust_polyfit(specval[wsp][srt], fluxval[wsp][srt], 3, function='bspline',
                                                  sigma=5., maxone=False, xmin=yb, xmax=ye,
                                                  everyn=specval[wsp].size//tcky.size)  # knots=tcky)
            blz_flat = arutils.func_val(blzspl, specval, 'bspline')
            msblaze[:, o] = arutils.func_val(blzspl, np.linspace(0.0, 1.0, msblaze.shape[0]), 'bspline')
        else:
            mask, blzspl = arutils.robust_polyfit(specval[wsp][srt], fluxval[wsp][srt], 2, function='polynomial',
                                                  sigma=5., maxone=False)
            blz_flat = arutils.func_val(blzspl, specval, 'polynomial')
            msblaze[:, o] = arutils.func_val(blzspl, np.linspace(0.0, 1.0, msblaze.shape[0]), 'polynomial')
            extrap_slit[o] = 1.0

        # Extract a spectrum of the trace frame
        xext = np.arange(mstrace.shape[0])
        yext = np.round(0.5 * (lordloc + rordloc)).astype(np.int)
        wcc = np.where((yext > 0) & (yext < mstrace.shape[1] - 1.0))
        blazeext[wcc[0], o] = mstrace[(xext[wcc], yext[wcc],)]
        if wcc[0].size != mstrace.shape[0]:
            extrap_slit[o] = 1.0

        # Calculate the slit profile
        sprof_fit = fluxval / (blz_flat + (blz_flat == 0.0))
        wch = np.where(wcchip)
        tckx = np.linspace(min(0.0, np.min(spatval[wch])), max(1.0, np.max(spatval[wch])), ntckx)
        tckx = tckx[np.where((tckx > np.min(spatval[wch])) & (tckx < np.max(spatval[wch])))]
        srt = np.argsort(spatval[wch])
        # Only perform a bspline if there are enough pixels for the specified knots
        if tckx.size >= 1:
            xb, xe = min(np.min(spatval), tckx[0]), max(np.max(spatval), tckx[-1])
            mask, sltspl = arutils.robust_polyfit(spatval[wch][srt], sprof_fit[wch][srt], 3, function='bspline',
                                                  sigma=5., maxone=False, xmin=xb, xmax=xe,
                                                  everyn=spatval[wch].size//tckx.size)  #, knots=tckx)
            slt_flat = arutils.func_val(sltspl, spatval, 'bspline')
            sltnrmval = arutils.func_val(sltspl, 0.5, 'bspline')
        else:
            mask, sltspl = arutils.robust_polyfit(spatval[srt], sprof_fit[srt], 2, function='polynomial',
                                                  sigma=5., maxone=False)
            slt_flat = arutils.func_val(sltspl, spatval, 'polynomial')
            sltnrmval = arutils.func_val(sltspl, 0.5, 'polynomial')
            extrap_slit[o] = 1.0

        modvals = blz_flat * slt_flat
        # Normalize to the value at the centre of the slit
        nrmvals = blz_flat * sltnrmval
        if settings.argflag["reduce"]["slitprofile"]["perform"]:
            # Leave slit_profiles as ones if the slitprofile is not being determined, otherwise, set the model.
            slit_profiles[word] = modvals/nrmvals
        mstracenrm[word] /= nrmvals
        if msgs._debug['slit_profile']:
            debugger.set_trace()
            model = np.zeros_like(mstrace)
            model[word] = modvals
            diff = mstrace - model
            hdu = fits.PrimaryHDU(mstrace)
            hdu.writeto("mstrace_{0:02d}.fits".format(det), overwrite=True)
            hdu = fits.PrimaryHDU(model)
            hdu.writeto("model_{0:02d}.fits".format(det), overwrite=True)
            hdu = fits.PrimaryHDU(diff)
            hdu.writeto("diff_{0:02d}.fits".format(det), overwrite=True)

    # Return
    return slit_profiles, mstracenrm, msblaze, blazeext, extrap_slit


def slit_profile_pca(slf, mstrace, det, msblaze, extrap_slit, slit_profiles):
    """ Perform a PCA analysis on the spatial slit profile and blaze function.

    Parameters
    ----------
    slf : class
      Science Exposure Class
    slit_profile : ndarray
      An image containing the slit profile
    det : int
      Detector index
    msblaze : ndarray
      A model of the blaze function of each slit
    extrap_slit : ndarray
      Mask indicating if a slit is well-determined (0) or poor (1). If the latter, the slit profile
      and blaze function for those slits should be extrapolated or determined from another means
    slit_profiles : ndarray
      An image containing the slit profile

    Returns
    -------
    slit_profiles : ndarray
      An image containing the slit profile
    mstracenrm : ndarray
      The input trace frame, normalized by the blaze function (but still contains the slit profile)
    extrap_blz : ndarray
      A model of the blaze function of each slit
    """
    #################
    # Parameters to include in settings file
    fitfunc = "legendre"
    ordfit = 4
    ofit = [2, 3, 3, 2, 2]
    sordfit = 2
    sofit = [1, 3, 1]
    #################

    nslits = extrap_slit.size
    gds = np.where(extrap_slit == 0)
    maskord = np.where(extrap_slit == 1)[0]
    specfit = np.arange(mstrace.shape[0])
    nspec = np.max(slf._pixwid[det - 1])*10
    spatbins = np.linspace(-0.25, 1.25, nspec + 1)
    # Perform a PCA on the spectral (i.e. blaze) function
    blzmxval = np.ones((1, nslits))
    lorr = 0
    for o in range(0, nslits):
        # if extrap_slit[o] == 1:
        #     continue
        # Find which pixels are on the slit
        wch = np.where((slf._lordloc[det-1][:, o] > 0.0) &
                       (slf._rordloc[det - 1][:, o] < mstrace.shape[1]-1.0))
        cordloc = np.round(0.5 * (slf._lordloc[det - 1][:, o] + slf._rordloc[det - 1][:, o])).astype(np.int)
        if wch[0].size < mstrace.shape[0]:
            # The entire order is not on the chip
            if cordloc[int(0.5*mstrace.shape[0])] < mstrace.shape[1]/2:
                lorr = -1  # Once a full order is found, go left
                continue
            else:
                lorr = +1  # Go right
        else:
            blzmxval[0, o] = np.median(mstrace[wch[0], cordloc[wch]])
        if lorr == -1:
            # A full order has been found, go back and fill in the gaps
            for i in range(1, o+1):
                wch = np.where((slf._lordloc[det-1][:, o-i] > 0.0) &
                               (slf._rordloc[det-1][:, o-i] < mstrace.shape[1] - 1.0))
                # Calculate the previous order flux
                cordloc = np.round(0.5 * (slf._lordloc[det-1][:, o-i+1] + slf._rordloc[det-1][:, o-i+1])).astype(np.int)
                prval = mstrace[wch[0], cordloc[wch]]
                # Calculate the current order flux
                cordloc = np.round(0.5 * (slf._lordloc[det-1][:, o-i] + slf._rordloc[det-1][:, o-i])).astype(np.int)
                mnval = mstrace[wch[0], cordloc[wch]]
                wnz = np.where(prval != 0.0)
                blzmxval[0, o-i] = blzmxval[0, o-i+1] * np.median(mnval[wnz] / prval[wnz])
            lorr = 0
        elif lorr == +1:
            # Calibrate the current order with the previous one
            mnval = mstrace[wch[0], cordloc[wch]]
            cordloc = np.round(0.5 * (slf._lordloc[det - 1][:, o-1] + slf._rordloc[det - 1][:, o-1])).astype(np.int)
            prval = mstrace[wch[0], cordloc[wch]]
            wnz = np.where(prval != 0.0)
            blzmxval[0, o] = blzmxval[0, o-1] * np.median(mnval[wnz] / prval[wnz])
            lorr = 0

    # Check for nan values (i.e. when median is given a zero element array)
    blznan = np.isnan(blzmxval[0, :])
    if np.any(blznan):
        # Find the acceptable values and linearly interpolate
        blzx = np.arange(nslits)
        wnnan = np.where(~blznan)
        fblz = interpolate.interp1d(blzx[wnnan], blzmxval[0, wnnan],
                                    kind="linear", bounds_error=False, fill_value="extrapolate")
        blzmxval = fblz(blzx).reshape(blzmxval.shape)
    elif np.all(blznan):
        msgs.bug("All of the blaze values are NaN... time to debug")
        debugger.set_trace()

    # Calculate the mean blaze function of all good orders
    blzmean = np.mean(msblaze[:, gds[0]], axis=1)
    blzmean /= np.max(blzmean)
    blzmean = blzmean.reshape((blzmean.size, 1))
    msblaze /= blzmean
    msblaze /= blzmxval
    # Fit the blaze functions
    fitcoeff = np.ones((ordfit+1, nslits))
    for o in range(nslits):
        if extrap_slit[o] == 1:
            continue
        wmask = np.where(msblaze[:, o] != 0.0)[0]
        null, bcoeff = arutils.robust_polyfit(specfit[wmask], msblaze[wmask, o],
                                              ordfit, function=fitfunc, sigma=2.0,
                                              minv=0.0, maxv=mstrace.shape[0])
        fitcoeff[:, o] = bcoeff

    lnpc = len(ofit) - 1
    xv = np.arange(mstrace.shape[0])
    blzval = arutils.func_val(fitcoeff, xv, fitfunc,
                              minv=0.0, maxv=mstrace.shape[0] - 1).T
    # Only do a PCA if there are enough good orders
    if np.sum(1.0 - extrap_slit) > ofit[0] + 1:
        # Perform a PCA on the tilts
        msgs.info("Performing a PCA on the spectral blaze function")
        ordsnd = np.arange(nslits) + 1.0
        xcen = xv[:, np.newaxis].repeat(nslits, axis=1)
        fitted, outpar = arpca.basis(xcen, blzval, fitcoeff, lnpc, ofit, x0in=ordsnd, mask=maskord, skipx0=False,
                                     function=fitfunc)
        if not msgs._debug['no_qa']:
#            arqa.pca_plot(slf, outpar, ofit, "Blaze_Profile", pcadesc="PCA of blaze function fits")
            arpca.pca_plot(slf, outpar, ofit, "Blaze_Profile",
                           pcadesc="PCA of blaze function fits")
        # Extrapolate the remaining orders requested
        orders = 1.0 + np.arange(nslits)
        extrap_blz, outpar = arpca.extrapolate(outpar, orders, function=fitfunc)
        extrap_blz *= blzmean
        extrap_blz *= blzmxval
    else:
        msgs.warn("Could not perform a PCA on the order blaze function" + msgs.newline() +
                  "Not enough well-traced orders")
        msgs.info("Using direct determination of the blaze function instead")
        extrap_blz = msblaze*blzmean

    # Normalize the trace frame, but don't remove the slit profile
    mstracenrm = mstrace.copy()
    for o in range(nslits):
        word = np.where(slf._slitpix[det - 1] == o+1)
        specval = slf._tilts[det-1][word]
        blzspl = interpolate.interp1d(np.linspace(0.0, 1.0, mstrace.shape[0]), extrap_blz[:, o],
                                      kind="linear", fill_value="extrapolate")
        mstracenrm[word] /= blzspl(specval)

    # Now perform a PCA on the spatial (i.e. slit) profile
    # First generate the original model of the spatial slit profiles
    msslits = np.zeros((nspec, nslits))
    mskslit = np.ones((nspec, nslits))
    for o in range(nslits):
        if extrap_slit[o] == 1:
            continue
        word = np.where(slf._slitpix[det - 1] == o+1)
        spatval = (word[1] + 0.5 - slf._lordloc[det-1][:, o][word[0]]) /\
                  (slf._rordloc[det-1][:, o][word[0]] - slf._lordloc[det-1][:, o][word[0]])
        groups = np.digitize(spatval, spatbins)
        modelw = slit_profiles[word]
        for mm in range(1, spatbins.size):
            tmp = modelw[groups == mm]
            if tmp.size != 0.0:
                msslits[mm - 1, o] = tmp.mean()
            else:
                mskslit[mm - 1, o] = 0.0

    # Calculate the spatial profile of all good orders
    sltmean = np.mean(msslits[:, gds[0]], axis=1)
    sltmean = sltmean.reshape((sltmean.size, 1))
    msslits /= (sltmean + (sltmean == 0))

    # Fit the spatial profiles
    spatfit = 0.5*(spatbins[1:]+spatbins[:-1])
    fitcoeff = np.ones((sordfit+1, nslits))
    for o in range(nslits):
        if extrap_slit[o] == 1:
            continue
        wmask = np.where(mskslit[:, o] == 1.0)[0]
        null, bcoeff = arutils.robust_polyfit(spatfit[wmask], msslits[wmask, o],
                                              sordfit, function=fitfunc, sigma=2.0,
                                              minv=spatfit[0], maxv=spatfit[-1])
        fitcoeff[:, o] = bcoeff

    lnpc = len(sofit) - 1
    sltval = arutils.func_val(fitcoeff, spatfit, fitfunc,
                              minv=spatfit[0], maxv=spatfit[-1]).T
    # Only do a PCA if there are enough good orders
    if np.sum(1.0 - extrap_slit) > sofit[0] + 1:
        # Perform a PCA on the tilts
        msgs.info("Performing a PCA on the spatial slit profiles")
        ordsnd = np.arange(nslits) + 1.0
        xcen = spatfit[:, np.newaxis].repeat(nslits, axis=1)
        fitted, outpar = arpca.basis(xcen, sltval, fitcoeff, lnpc, sofit, x0in=ordsnd, mask=maskord, skipx0=False,
                                     function=fitfunc)
        if not msgs._debug['no_qa']:
#            arqa.pca_plot(slf, outpar, sofit, "Slit_Profile", pcadesc="PCA of slit profile fits")
            arpca.pca_plot(slf, outpar, sofit, "Slit_Profile", pcadesc="PCA of slit profile fits")
        # Extrapolate the remaining orders requested
        orders = 1.0 + np.arange(nslits)
        extrap_slt, outpar = arpca.extrapolate(outpar, orders, function=fitfunc)
        extrap_slt *= sltmean
        extrap_slt *= mskslit
    else:
        msgs.warn("Could not perform a PCA on the spatial slit profiles" + msgs.newline() +
                  "Not enough well-traced orders")
        msgs.info("Using direct determination of the slit profiles instead")
        extrap_slt = (msslits*mskslit)*sltmean

    # Normalize the trace frame, but don't remove the slit profile
    slit_profiles = np.ones_like(mstrace)
    for o in range(nslits):
        lordloc = slf._lordloc[det - 1][:, o]
        rordloc = slf._rordloc[det - 1][:, o]
        word = np.where(slf._slitpix[det - 1] == o+1)
        spatval = (word[1] - lordloc[word[0]])/(rordloc[word[0]] - lordloc[word[0]])

        sltspl = interpolate.interp1d(spatfit, extrap_slt[:, o],
                                      kind="linear", fill_value="extrapolate")
        slit_profiles[word] = sltspl(spatval)

    return slit_profiles, mstracenrm, extrap_blz



def slit_profile_qa(slf, mstrace, model, lordloc, rordloc, msordloc, textplt="Slit", maxp=16,
                    desc=""):
    """ Generate a QA plot for the slit profile of each slit

    Parameters
    ----------
    slf : class
      Science Exposure class
    mstrace : ndarray
      trace frame
    model : ndarray
      model of slit profiles, same shape as frame.
    lordloc : ndarray
      left edge locations of all slits
    rordloc : ndarray
      right edge locations of all slits
    msordloc : ndarray
      An array the same size as frame that determines which pixels contain a given order.
    textplt : str, optional
      A string printed above each panel
    maxp : int, (optional)
      Maximum number of panels per page
    desc : str, (optional)
      A description added to the top of each page
    """

    plt.rcdefaults()
    plt.rcParams['font.family']= 'times new roman'

    # Outfile
    method = inspect.stack()[0][3]
    outroot = arqa.set_qa_filename(slf.setup, method)

    npix, nord = lordloc.shape
    nbins = 40
    bins = np.linspace(-0.25, 1.25, nbins+1)
    pages, npp = arqa.get_dimen(nord, maxp=maxp)
    # Loop through all pages and plot the results
    ndone = 0
    axesIdx = True
    for i in range(len(pages)):
        f, axes = plt.subplots(pages[i][1], pages[i][0])
        ipx, ipy = 0, 0
        for j in range(npp[i]):
            if pages[i][0] == 1 and pages[i][1] == 1: axesIdx = False
            elif pages[i][1] == 1: ind = (ipx,)
            elif pages[i][0] == 1: ind = (ipy,)
            else: ind = (ipy, ipx)
            # Get data to be plotted
            word = np.where(msordloc == ndone+j+1)
            if word[0].size == 0:
                msgs.warn("There are no pixels in slit {0:d}".format(ndone + j + 1))
                # Delete the axis
                if pages[i][1] == 1: ind = (ipx,)
                elif pages[i][0] == 1: ind = (ipy,)
                else: ind = (ipy, ipx)
                f.delaxes(axes[ind])
                ipx += 1
                if ipx == pages[i][0]:
                    ipx = 0
                    ipy += 1
                continue
            spatval = (word[1] + 0.5 - lordloc[:, ndone+j][word[0]]) / (rordloc[:, ndone+j][word[0]] - lordloc[:, ndone+j][word[0]])
            fluxval = mstrace[word]
            mxval = 1.25
            modvals = np.zeros(nbins)
            if axesIdx:
                cnts, xedges, yedges, null = axes[ind].hist2d(spatval, fluxval, bins=bins, cmap=plt.cm.gist_heat_r)
                groups = np.digitize(spatval, xedges)
                modelw = model[word]
                for mm in range(1, xedges.size):
                    modvals[mm-1] = modelw[groups == mm].mean()
                axes[ind].plot(0.5*(xedges[1:]+xedges[:-1]), modvals, 'b-', linewidth=2.0)
                axes[ind].plot([0.0, 0.0], [0.0, mxval], 'r-')
                axes[ind].plot([1.0, 1.0], [0.0, mxval], 'r-')
            else:
                cnts, xedges, yedges, null = axes.hist2d(spatval, fluxval, bins=bins, cmap=plt.cm.gist_heat_r)
                groups = np.digitize(spatval, xedges)
                modelw = model[word]
                for mm in range(1, xedges.size):
                    modvals[mm-1] = modelw[groups == mm].mean()
                axes.plot(0.5*(xedges[1:]+xedges[:-1]), modvals, 'b-', linewidth=2.0)
                axes.plot([0.0, 0.0], [0.0, mxval], 'r-')
                axes.plot([1.0, 1.0], [0.0, mxval], 'r-')
            if axesIdx:
                axes[ind].axis([xedges[0], xedges[-1], 0.0, 1.1*mxval])
                axes[ind].set_title("{0:s} {1:d}".format(textplt, ndone+j+1))
                axes[ind].tick_params(labelsize=10)
            else:
                axes.axis([xedges[0], xedges[-1], 0.0, 1.1*mxval])
                axes.set_title("{0:s} {1:d}".format(textplt, ndone+j+1))
                axes.tick_params(labelsize=10)
            ipx += 1
            if ipx == pages[i][0]:
                ipx = 0
                ipy += 1
        # Delete the unnecessary axes
        if axesIdx:
            for j in range(npp[i], axes.size):
                if pages[i][1] == 1: ind = (ipx,)
                elif pages[i][0] == 1: ind = (ipy,)
                else: ind = (ipy, ipx)
                f.delaxes(axes[ind])
                ipx += 1
                if ipx == pages[i][0]:
                    ipx = 0
                    ipy += 1
        ndone += npp[i]
        # Save the figure
        if axesIdx: axsz = axes.size
        else: axsz = 1.0
        if pages[i][1] == 1 or pages[i][0] == 1: ypngsiz = 11.0/axsz
        else: ypngsiz = 11.0*axes.shape[0]/axes.shape[1]
        f.set_size_inches(11.0, ypngsiz)
        if desc != "":
            pgtxt = ""
            if len(pages) != 1:
                pgtxt = ", page {0:d}/{1:d}".format(i+1, len(pages))
            f.suptitle(desc + pgtxt, y=1.02, size=16)
        f.tight_layout()
        outfile = outroot+'{:03d}.png'.format(i)
        plt.savefig(outfile, dpi=200)
        plt.close()
        f.clf()
    del f

    plt.rcdefaults()

    return


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
        if satpix is not None:
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

    # Old cr_screen can yield nan pixels, new one does not, meaning that
    # there are differences between the returned arrays.  For all
    # non-nan pixels, the two algorithms are identical.
#    print('calling cr_screen')
#    t = time.clock()
#    _sigimg  = arcyproc.cr_screen(filty,0.0)
#    print('Old cr_screen: {0} seconds'.format(time.clock() - t))
#    print(np.sum(np.invert(np.isfinite(_sigimg))))
#    t = time.clock()
    sigimg  = new_cr_screen(filty)
#    print('New cr_screen: {0} seconds'.format(time.clock() - t))
#    print(np.sum(np.invert(np.isfinite(sigimg))))
#    if np.sum(_sigimg != sigimg) != 0:
#        plt.imshow(_sigimg, origin='lower', interpolation='nearest', aspect='auto')
#        plt.show()
#        plt.imshow(sigimg, origin='lower', interpolation='nearest', aspect='auto')
#        plt.show()
#        plt.imshow(_sigimg-sigimg, origin='lower', interpolation='nearest', aspect='auto')
#        plt.colorbar()
#        plt.show()
#        r = np.ma.divide(_sigimg,sigimg)-1
#        print(np.ma.sum(r))
#        plt.imshow(np.ma.divide(_sigimg,sigimg)-1, origin='lower', interpolation='nearest', aspect='auto')
#        plt.colorbar()
#        plt.show()
#    assert np.sum(_sigimg != sigimg) == 0, 'Difference between old and new cr_screen'

#    print(sigimg.shape)
#    print(new_sigimg.shape)
#    print(np.ma.mean(np.ma.log10(sigimg)))
#    print(np.ma.mean(np.ma.log10(new_sigimg)))
#    print(np.ma.mean(np.ma.log10(sigimg)-np.ma.log10(new_sigimg)))
#
#    plt.imshow(np.ma.log10(sigimg), origin='lower', interpolation='nearest', aspect='auto')
#    plt.colorbar()
#    plt.show()
#
#    plt.imshow(np.ma.log10(new_sigimg), origin='lower', interpolation='nearest', aspect='auto')
#    plt.colorbar()
#    plt.show()
#
#    plt.imshow(np.ma.log10(new_sigimg)-np.ma.log10(sigimg), origin='lower', interpolation='nearest', aspect='auto')
#    plt.colorbar()
#    plt.show()

    sigsmth = ndimage.filters.gaussian_filter(sigimg,1.5)
    sigsmth[np.where(np.isnan(sigsmth))]=0.0
    sigmask = np.cast['bool'](np.zeros(sciframe.shape))
    sigmask[np.where(sigsmth>sigclip)] = True
    crmask = np.logical_and(crmask, sigmask)
    msgs.info("Growing cosmic ray mask by 1 pixel")
#    print('calling grow_masked')
#    t = time.clock()
#    _crmask = arcyutils.grow_masked(crmask.astype(np.float), grow, 1.0)
#    print('Old grow_masked: {0} seconds'.format(time.clock() - t))
#    t = time.clock()
    crmask = new_grow_masked(crmask.astype(np.float), grow, 1.0)
#    print('New grow_masked: {0} seconds'.format(time.clock() - t))
#    assert np.sum(_crmask != crmask) == 0, 'Difference between old and new grow_masked'

    return crmask


def new_cr_screen(a, mask_value=0.0, spatial_axis=1):
    r"""
    Calculate the significance of pixel deviations from the median along
    the spatial direction.

    No type checking is performed of the input array; however, the
    function assumes floating point values.

    Args:
        a (numpy.ndarray): Input 2D array
        mask_value (float): (**Optional**) Values to ignore during the
            calculation of the median.  Default is 0.0.
        spatial_axis (int): (**Optional**) Axis along which to calculate
            the median.  Default is 1.

    Returns:
        numpy.ndarray: Returns a map of :math:`|\Delta_{i,j}|/\sigma_j`,
        where :math:`\Delta_{i,j}` is the difference between the pixel
        value and the median along axis :math:`i` and :math:`\sigma_j`
        is robustly determined using the median absolute deviation,
        :math:`sigma_j = 1.4826 MAD`.
    """
    # Check input
    if len(a.shape) != 2:
        msgs.error('Input array must be two-dimensional.')
    if spatial_axis not in [0,1]:
        msgs.error('Spatial axis must be 0 or 1.')

    # Mask the pixels equal to mask value: should use np.isclose()
    _a = np.ma.MaskedArray(a, mask=(a==mask_value))
    # Get the median along the spatial axis
    meda = np.ma.median(_a, axis=spatial_axis)
    # Get a robust measure of the standard deviation using the median
    # absolute deviation; 1.4826 factor is the ratio of sigma/MAD
    d = np.absolute(_a - meda[:,None])
    mada = 1.4826*np.ma.median(d, axis=spatial_axis)
    # Return the ratio of the difference to the standard deviation
    return np.ma.divide(d, mada[:,None]).filled(mask_value)


def new_grow_masked(img, grow, growval):

    if not np.any(img == growval):
        return img

    _img = img.copy()
    sz_x, sz_y = img.shape
    d = int(1+grow)
    rsqr = grow*grow

    # Grow any masked values by the specified amount
    for x in range(sz_x):
        for y in range(sz_y):
            if img[x,y] != growval:
                continue

            mnx = 0 if x-d < 0 else x-d
            mxx = x+d+1 if x+d+1 < sz_x else sz_x
            mny = 0 if y-d < 0 else y-d
            mxy = y+d+1 if y+d+1 < sz_y else sz_y

            for i in range(mnx,mxx):
                for j in range(mny, mxy):
                    if (i-x)*(i-x)+(j-y)*(j-y) <= rsqr:
                        _img[i,j] = growval
    return _img


def new_order_pixels(pixlocn, lord, rord):
    """
    Based on physical pixel locations, determine which pixels are within the orders
    """

    sz_x, sz_y, _ = pixlocn.shape
    sz_o = lord.shape[1]

    outfr = np.zeros((sz_x, sz_y), dtype=int)

    for y in range(sz_y):
        for o in range(sz_o):
            indx = (lord[:,o] < rord[:,o]) & (pixlocn[:,y,1] > lord[:,o])\
                        & (pixlocn[:,y,1] < rord[:,o])
            indx |= ( (lord[:,o] > rord[:,o]) & (pixlocn[:,y,1] < lord[:,o]) 
                        & (pixlocn[:,y,1] > rord[:,o]) )
            if np.any(indx):
                # Only assign a single order to a given pixel
                outfr[indx,y] = o+1
                break

    return outfr


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


def sub_overscan(frame, numamplifiers, datasec, oscansec, settings=None):
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
    #dnum = settings.get_dnum(det)
    if settings is None:
        settings = dict(recduce={'overscan': {'method': 'savgol', 'params': [5,65]}})


    for i in range(settings.spect[dnum]['numamplifiers']):
        # Determine the section of the chip that contains data
        datasec = "datasec{0:02d}".format(i+1)
        #dx0, dx1 = settings.spect[dnum][datasec][0][0], settings.spect[dnum][datasec][0][1]
        #dy0, dy1 = settings.spect[dnum][datasec][1][0], settings.spect[dnum][datasec][1][1]
        dx0, dx1 = datasec[i][0][0], datasec[i][0][1]
        dy0, dy1 = datasec[i][1][0], datasec[i][1][1]
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
            ossub = signal.savgol_filter(osfit,
                                         settings.argflag['reduce']['overscan']['params'][1],
                                         settings.argflag['reduce']['overscan']['params'][0])
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
