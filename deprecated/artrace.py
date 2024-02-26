from __future__ import (print_function, absolute_import, division, unicode_literals)

import inspect
import copy

import numpy as np

from scipy import interpolate

import matplotlib.pyplot as plt
from matplotlib import cm, font_manager

from astropy.stats import sigma_clip
from astropy.convolution import convolve, Gaussian1DKernel
from pypeit import msgs
from pypeit.core import qa
from pypeit.core import plot
from pypeit.core import arc
from pypeit import utils
from pypeit.core import pca
from pypeit.filter import BoxcarFilter
from pypeit import debugger

from pypeit.par.parset import ParSet

try:
    from pypeit import ginga
except ImportError:
    pass

# Testing
import time

# This code is now deprecated and one should be using detect_lines for peak finding.
###########
def fit_min(xarr, yarr, xguess, width=None):

    errcode = 0
    # Edges
    if width is None:
        xleft, xright = np.min(xarr), np.max(xarr)
    else:
        xleft = xguess - width
        xright = xguess + width
    idx = np.where((xarr >= xleft) & (xarr <= xright))[0]

    # Setup
    thisx = xarr[idx]
    thisy = yarr[idx]

    # Guess for Gaussian
    guess = np.max(thisy), 0., width/2.

    # Fit with Gaussian
    try:
        coeff = func_fit(thisx-xguess, thisy, 'gaussian', 3, guesses=guess)
    except RuntimeError:  # Bad fit
        errcode = -1
        return xguess, 0., errcode
    sigma = coeff[2]
    xbest = xguess + coeff[1]

    # Could/should add a bunch of sanity checks
    # Insist on it being a minimum
    if coeff[0] > 0.:
        errcode = -4
    if (xbest < xleft) or (xbest > xright):
        errcode = -6
    # Return
    return xbest, sigma, errcode

# This code is now deprecated and one should be using detect_lines for peak finding.
def find_nminima(yflux, xvec=None, nfind=10, nsmooth=None, minsep=5, width=5):
    """ Find minima in an input 1D array
    Parameters
    ----------
    yflux : ndarray
    xvec : ndarray, optional
      Assumed to be ascending
    nfind : int, optional
      Number of peaks to find in the input array
    nsmooth : int, optional
      Smooth by a Gaussian with kenrel of nsmooth
    minsep : int, optional
      Minimum separation between peaks
    width : int, optional
      Width around a putative peak to fit a Gaussian

    Returns
    -------
    peaks: ndarray
      x values of the peaks
    sigmas: ndarray
      sigma widths of the Gaussian fits to each peak
    ledges: ndarray
      left edges of each peak;  defined to be at least minsep away
      from the peak and where the slope of the data switches
    redges: ndarray
      right edges of each peak;  defined to be at least minsep away
      from the peak and where the slope of the data switches
    """
    # Init
    if xvec is None:
        xvec = np.arange(len(yflux))
    # Gaussian smooth
    if nsmooth is not None:
        yflux = convolve(yflux, Gaussian1DKernel(nsmooth))#, **kwargs)

    # ycopy, yderiv, ydone
    ycopy = yflux.copy()
    yderiv = np.roll(ycopy,1)-ycopy
    yderiv[0] = 0.
    yderiv[-1] = 0.
    ydone = np.max(ycopy)

    # Find first one
    peaks, sigmas, ledges, redges = [], [], [], []
    npeak = 0
    for kk in range(nfind):
        imin = np.argmin(ycopy)
        xbest, sigma, errcode = fit_min(xvec, ycopy, xvec[imin], width=width)
        #
        noldpeak = npeak
        npeak = len(peaks)
        # Find edges and
        # Block out pixels within minsep and 2*minsep
        x1 = (xvec < xvec[imin]-minsep) & (np.roll(yderiv,1) < 0.)
        if np.any(x1):
            ix1 = np.where(x1)[0][-1]
        else:
            ix1 = 0
        x2 = (xvec > xvec[imin]+minsep) & (yderiv > 0.)  # Scans until increasing
        if np.any(x2):
            ix2 = np.where(x2)[0][0]
        else:
            ix2 = len(xvec)
        ycopy[ix1:ix2] = ydone
        # Save
        if npeak == 0:  # Always grab at least one
            peaks.append(xbest)
            sigmas.append(sigma)
            ledges.append(ix1)
            redges.append(ix2-1)
        else:  # Check it is minsep away (seems like it will always be)
            xmin = np.min(np.abs(np.array(peaks-xbest)))
            if (xmin > minsep) & (errcode >= 0):
                peaks.append(xbest)
                sigmas.append(sigma)
                ledges.append(ix1)
                redges.append(ix2-1)
        # Any more to look for?
        if not np.any(ycopy < ydone):
            npeak = nfind
    return np.array(peaks), np.array(sigmas), np.array(ledges), np.array(redges)


def trace_objbg_image(lordloc, sciframe, slitn, objreg, bgreg, trim=2, triml=None, trimr=None):
    """ Creates an image with weights corresponding to object or background pixels.

    Each weight can take any floating point value from 0 to 1 (inclusive). For the
    rec_obj_img, a weight of 1 means that the pixel is fully contained within the
    object region, and 0 means that the pixel is fully contained within the
    background region. The opposite is true for the rec_bg_img array. A pixel that
    is on the border of object/background is assigned a value between 0 and 1, based
    on the percentage overlap with the object/background regions.

    Parameters
    ----------
    lordloc : ndarray
    sciframe: numpy ndarray
      Science frame
    slitn : int
      Slit (or order) number
    objreg : list
      A two element list containing arrays that indicate the left and right
      pixels of the object locations for each object.
    bgreg : list
      A two element list containing arrays that indicate the background regions
      around each object (i.e. [bgleft, bgright]). Each array is 2D and contains
      a list of pixels for each object. The size of each array is [npix, nobj].
    trim : int (optional)
      Number of pixels to trim from the left and right slit edges.
      To separately specify how many pixels to trim from the left
      and right slit edges, use triml and trimr.
    triml : int (optional)
      Number of pixels to trim from the left slit edge
    trimr : int (optional)
      Number of pixels to trim from the right slit edge

    Returns
    -------
    rec_obj_img : ndarray
      An image containing weights to be used for the object
    rec_bg_img : ndarray
      An image containing weights to be used for the background
    """
    # Get the number of objects in the slit
    nobj = len(objreg[0])
    # Set the number of pixels to trim of the left and right edges
    if triml is None:
        triml = trim
    if trimr is None:
        trimr = trim
    # Generate an array of spatial pixels for each row on the detector
    spatdir = np.arange(sciframe.shape[1])[np.newaxis, :].repeat(sciframe.shape[0], axis=0)
    # Make an image of pixel weights for each object
    msgs.info("Creating an image weighted by object pixels for {:d} objects".format(nobj))
    # This array can be a *big* memory hog..
    rec_obj_img = np.zeros((sciframe.shape[0], sciframe.shape[1], nobj)).astype(np.float32)
    for o in range(nobj):
        lobj = lordloc[:, slitn] + triml + objreg[0][o] - 1.0
        robj = lordloc[:, slitn] + trimr + objreg[1][o]
        rec_obj_img[:, :, o] = np.clip(spatdir - lobj.reshape(sciframe.shape[0], 1), 0.0, 1.0) - \
                               np.clip(spatdir - robj.reshape(sciframe.shape[0], 1), 0.0, 1.0)
    # Make an image of pixel weights for the background region of each object
    msgs.info("Creating an image weighted by background pixels")
    rec_bg_img = np.zeros((sciframe.shape[0], sciframe.shape[1], nobj)).astype(np.float32)
    for o in range(nobj):
        wll = np.where(bgreg[0][1:, o] > bgreg[0][:-1, o])[0]
        wlr = np.where(bgreg[0][1:, o] < bgreg[0][:-1, o])[0]
        if len(wll) < len(wlr): #< len(wlr): # JXP kludge
            wll = np.concatenate([np.array([0]).astype(int), wll])
        # Background regions to the left of object
        for ii in range(wlr.size):
            lobj = lordloc[:, slitn] + triml + wll[ii]
            robj = lordloc[:, slitn] + trimr + wlr[ii]
            rec_bg_img[:, :, o] += np.clip(spatdir - lobj.reshape(sciframe.shape[0], 1), 0.0, 1.0) - \
                                   np.clip(spatdir - robj.reshape(sciframe.shape[0], 1), 0.0, 1.0)
        wrl = np.where(bgreg[1][1:, o] > bgreg[1][:-1, o])[0]
        wrr = np.where(bgreg[1][1:, o] < bgreg[1][:-1, o])[0]
        if len(wrr) < len(wrl): # JXP kludge
            wrr = np.concatenate([wrr, np.array([len(bgreg[1][1:,o])-1]).astype(int)])
        # Background regions to the right of object
        for ii in range(wrl.size):
            lobj = lordloc[:, slitn] + triml + wrl[ii]
            robj = lordloc[:, slitn] + trimr + wrr[ii]
            rec_bg_img[:, :, o] += np.clip(spatdir - lobj.reshape(sciframe.shape[0], 1), 0.0, 1.0) - \
                                   np.clip(spatdir - robj.reshape(sciframe.shape[0], 1), 0.0, 1.0)
    return rec_obj_img, rec_bg_img


def trace_object_dict(nobj, traces, object=None, background=None, params=None, tracelist=None):
    """ Creates a list of dictionaries, which contain the object traces in each slit

    Parameters
    ----------
    nobj : int
      Number of objects in this slit
    traces : numpy ndarray
      the trace of each object in this slit
    object: numpy ndarray (optional)
      An image containing weights to be used for the object.
    background : numpy ndarray (optional)
      An image containing weights to be used for the background
    params : dict
      A dictionary containing some of the important parameters used in
      the object tracing.
    tracelist : list of dict
      A list containing a trace dictionary for each slit

    To save memory, the object and background images for multiple slits
    can be stored in a single object image and single background image.
    In this case, you must only set object and background for the zeroth
    slit (with 'None' set for all other slits). In this case, the object
    and background images must be set according to the slit locations
    assigned in the slf._slitpix variable.

    Returns
    -------
    tracelist : list of dict
      A list containing a trace dictionary for each slit
    """
    # Create a dictionary with all the properties of the object traces in this slit
    newdict = dict({})
    newdict['nobj'] = nobj
    newdict['traces'] = traces
    newdict['object'] = object
    newdict['params'] = params
    newdict['background'] = background
    if tracelist is None:
        tracelist = []
    tracelist.append(newdict)
    return tracelist


# TODO: Is this ever called?
def trace_objects_in_slits(slf, det, sciframe, varframe, crmask, doqa=False, **kwargs):
    """ Wrapper for trace_objects_in_slit

    Parameters
    ----------
    slf : SciExposure class
    det : int
    sciframe : ndarray
    varframe : ndarray
    crmask : ndarray
    doqa : generate PNG
    kwargs : optional
      Passed to trace_objects_in_slit

    Returns
    -------
    tracelist : list
       Contains all the trace information for all slits

    """
    nslit = len(slf._maskslits[det-1])
    gdslits = np.where(~slf._maskslits[det-1])[0]
    tracelist = []

    # Loop on good slits
    for slit in range(nslit):
        if slit not in gdslits:
            tracelist.append({})
            continue
        tlist = trace_objects_in_slit(slf, det, slit, sciframe, varframe, crmask, **kwargs)
        # Append
        tracelist += tlist

    # QA?
    if doqa:
        obj_trace_qa(slf, sciframe, det, tracelist, root="object_trace", normalize=False)

    # Return
    return tracelist


def trace_objects_in_slit(det, slitn, tslits_dict, sciframe, skyframe, varframe, crmask,
                          tracefunc, traceorder, tracefind, nsmooth, trim=2, triml=None,
                          trimr=None, sigmin=2.0, bgreg=None, maskval=-999999.9, xedge=0.03,
                          standard=False, debug=False, manual=None, maxnumber=None):
    """ Finds objects, and traces their location on the detector

    Parameters
    ----------
    det : int
      Index of the detector
    slitn : int
      Slit (or order) number
    tslits_dict : dict
      Slits definition
    sciframe : numpy ndarray
      Science frame
    skyframe : numpy ndarray
      Science frame
    varframe : numpy ndarray
      Variance frame
    settings_trace : dict
      Settings for object finding + tracing
    crmask : numpy ndarray
      Mask or cosmic rays
    trim : int (optional)
      Number of pixels to trim from the left and right slit edges.
      To separately specify how many pixels to trim from the left
      and right slit edges, use triml and trimr.
    triml : int (optional)
      Number of pixels to trim from the left slit edge
    trimr : int (optional)
      Number of pixels to trim from the right slit edge
    sigmin : float
      Significance threshold to eliminate CRs
    bgreg : int
      Number of pixels to use in the background region
    maskval : float
      Placeholder value used to mask pixels
    xedge : float
      Trim objects within xedge % of the slit edge
    doqa : bool
      Should QA be output?

    Returns
    -------
    tracelist : list
      A single item list which is a dictionary containing the object trace information
    """
    # TODO -- Synchronize and avoid duplication in the usage of triml, trimr, trim, and xedge

    # Unpack the trace_slits dict (for now)
    pixwid = tslits_dict['pixwid']
    lordloc = tslits_dict['lcen']
    rordloc = tslits_dict['rcen']

    #
    skysub = sciframe - skyframe

    # Find the trace of each object
    if triml is None:
        triml = trim
    if trimr is None:
        trimr = trim
    npix = int(pixwid[slitn] - triml - trimr)
    if bgreg is None:
        bgreg = npix
    # Store the trace parameters
    tracepar = dict(smthby=7, rejhilo=1, bgreg=bgreg, triml=triml, trimr=trimr,
                    tracefunc=tracefunc, traceorder=traceorder, xedge=xedge)
    # Interpolate the science array onto a new grid (with constant spatial slit length)
    msgs.info("Rectifying science frame")
    xint = np.linspace(0.0, 1.0, skysub.shape[0])
    yint = np.linspace(0.0, 1.0, skysub.shape[1])
    scispl = interpolate.RectBivariateSpline(xint, yint, skysub, bbox=[0.0, 1.0, 0.0, 1.0],
                                             kx=1, ky=1, s=0)
    varspl = interpolate.RectBivariateSpline(xint, yint, varframe, bbox=[0.0, 1.0, 0.0, 1.0],
                                             kx=1, ky=1, s=0)
    crmspl = interpolate.RectBivariateSpline(xint, yint, crmask, bbox=[0.0, 1.0, 0.0, 1.0], kx=1,
                                             ky=1, s=0)
    xx, yy = np.meshgrid(np.linspace(0.0, 1.0, skysub.shape[0]), np.linspace(0.0, 1.0, npix), indexing='ij')
    ro = (rordloc[:, slitn] - trimr).reshape((-1, 1)) / (skysub.shape[1] - 1.0)
    lo = (lordloc[:, slitn] + triml).reshape((-1, 1)) / (skysub.shape[1] - 1.0)
    vv = (lo+(ro-lo)*yy).flatten()
    xx = xx.flatten()
    recsh = (skysub.shape[0], npix)
    rec_skysub = scispl.ev(xx, vv).reshape(recsh)
    rec_varframe = varspl.ev(xx, vv).reshape(recsh)
    rec_crmask   = crmspl.ev(xx, vv).reshape(recsh)
    # Update the CR mask to ensure it only contains 1's and 0's
    rec_crmask[np.where(rec_crmask > 0.2)] = 1.0
    rec_crmask[np.where(rec_crmask <= 0.2)] = 0.0
    msgs.info("Estimating object profiles")
    # Avoid any 0's in varframe
    zero_var = rec_varframe == 0.
    if np.any(zero_var):
        rec_varframe[zero_var] = 1.
        rec_crmask[zero_var] = 1.
    # Smooth the S/N frame
    smthby, rejhilo = tracepar['smthby'], tracepar['rejhilo']
    tmp = np.ma.MaskedArray(rec_skysub/np.sqrt(rec_varframe), mask=rec_crmask.astype(bool))
    # TODO: Add rejection to BoxcarFilter?
    rec_sigframe_bin = BoxcarFilter(smthby).smooth(tmp.T).T

    #rec_varframe_bin = arcyutils.smooth_x(rec_varframe, 1.0-rec_crmask, smthby, rejhilo, maskval)
    #rec_sigframe_bin = np.sqrt(rec_varframe_bin/(smthby-2.0*rejhilo))
    #sigframe = rec_skysub*(1.0-rec_crmask)/rec_sigframe_bin
    ww = np.where(rec_crmask == 0.0)
    med, mad = utils.robust_meanstd(rec_sigframe_bin[ww])
    ww = np.where(rec_crmask == 1.0)
    rec_sigframe_bin[ww] = maskval
    srtsig = np.sort(rec_sigframe_bin,axis=1)
    ww = np.where(srtsig[:, -2] > med + sigmin*mad)
    mask_sigframe = np.ma.array(rec_sigframe_bin, mask=rec_crmask, fill_value=maskval)
    # Clip image along spatial dimension -- May eliminate bright emission lines
    clip_image = sigma_clip(mask_sigframe[ww[0], :], axis=0, sigma=4.)
    # Collapse along the spectral direction to get object profile
    #trcprof = np.ma.mean(mask_sigframe[ww[0],:], axis=0).filled(0.0)
    trcprof = np.ma.mean(clip_image, axis=0).filled(0.0)
    trcxrng = np.arange(npix)/(npix-1.0)
    msgs.info("Identifying objects that are significantly detected")
    # Find significantly detected objects
    mskpix, coeff = utils.robust_polyfit(trcxrng, trcprof, 1+npix//40,
                                           function='legendre', sigma=2.0, minv=0.0, maxv=1.0)
    backg = utils.func_val(coeff, trcxrng, 'legendre', minv=0.0, maxv=1.0)
    trcprof -= backg
    wm = np.where(mskpix == 0)
    if wm[0].size == 0:
        msgs.warn("No objects found")
        return trace_object_dict(0, None)
    med, mad = utils.robust_meanstd(trcprof[wm])
    trcprof -= med
    # Gaussian smooth
    #from scipy.ndimage.filters import gaussian_filter1d
    #smth_prof = gaussian_filter1d(trcprof, fwhm/2.35)
    # Define all 5 sigma deviations as objects (should make the 5 user-defined)
    if tracefind == 'nminima':
        trcprof2 = np.mean(rec_skysub, axis=0)
        objl, objr, bckl, bckr = find_obj_minima(trcprof2, triml=triml, trimr=trimr, nsmooth=nsmooth)
    elif tracefind == 'standard':
        objl, objr, bckl, bckr = find_objects(trcprof, bgreg, mad)
    else:
        msgs.error("Bad object identification algorithm!!")
    # Remove objects within x percent of the slit edge
    xp = (objr+objl)/float(trcprof.size)/2.
    gdobj = (xp < (1-xedge)) * (xp > xedge)
    objl = objl[gdobj]
    objr = objr[gdobj]
    bckl = bckl[:, gdobj]
    bckr = bckr[:, gdobj]
    if np.sum(~gdobj) > 0:
        msgs.warn("Removed objects near the slit edges")
    #
    nobj = objl.size

    # TODO: 'Manual' should be run in a different function
    if manual is not None and not standard:
        msgs.info('Manual extraction desired. Rejecting all automatically detected '
                  'objects for now.')
        # Work on: Instead of rejecting all objects, prepend the manual extraction object?

        _manual = manual if isinstance(manual, list) else [manual]
        nmanual = len(_manual)
        for i in range(nmanual):
            if not isinstance(_manual[i], (ParSet,dict)):
                raise TypeError('Manual extraction must be defined using a dict or ParSet.')

        if _manual[0]['params'][0] == det:
            nobj = 1
            cent_spatial_manual = _manual[0]['params'][1]

            # Entered into .pypeit file in this format:
            #   [det, x_pixel_location, y_pixel_location,[x_range, y_range]]
            # - 1 or x_pixel_location is spatial pixel;
            # - 2 or y_pixel_location is dispersion/spectral pixel
            width_spatial_manual = _manual[0]['params'][3][0]
            objl = np.array([int(cent_spatial_manual - lordloc[cent_spatial_manual])
                                - width_spatial_manual])
            objr = np.array([int(cent_spatial_manual - lordloc[cent_spatial_manual])
                                + width_spatial_manual])
            bckl = np.zeros((trcprof.shape[0], objl.shape[0]))
            bckr = np.zeros((trcprof.shape[0], objl.shape[0]))

            for o in range(nobj):
                for x in range(1, bgreg + 1):
                    if objl[o] - x >= 0:
                        bckl[objl[o] - x, o] = 1
                    if objr[o] + x <= trcprof.shape[0] - 1:
                        bckr[objr[o] + x, o] = 1
        else:
            nobj = 0
    if nobj == 1:
        msgs.info("Found {0:d} object".format(objl.size))
        msgs.info("Tracing {0:d} object".format(objl.size))
    else:
        msgs.info("Found {0:d} objects".format(objl.size))
        msgs.info("Tracing {0:d} objects".format(objl.size))
    # Max obj
    if maxnumber is not None and nobj > maxnumber:
        nobj = maxnumber
        msgs.warn("Restricting to the brightest {:d} objects found".format(nobj))
        objl = objl[:nobj]
        objr = objr[:nobj]
        bckl = bckl[:, :nobj]
        bckr = bckr[:, :nobj]
    # Trace objects
    cval = np.zeros(nobj)
    allsfit = np.array([])
    allxfit = np.array([])
    clip_image2 = sigma_clip(mask_sigframe, axis=0, sigma=4.)
    for o in range(nobj):
        xfit = np.arange(objl[o], objr[o]).reshape((1, -1))/(npix-1.0)
        #cent = np.ma.sum(mask_sigframe[:,objl[o]:objr[o]]*xfit, axis=1)
        #wght = np.ma.sum(mask_sigframe[:,objl[o]:objr[o]], axis=1)
        cent = np.ma.sum(clip_image2[:, objl[o]:objr[o]]*xfit, axis=1)
        wght = np.ma.sum(clip_image2[:, objl[o]:objr[o]], axis=1)
        cent /= wght
        centfit = cent.filled(maskval)
        specfit = np.linspace(-1.0, 1.0, skysub.shape[0])
        w = np.where(centfit != maskval)
        specfit = specfit[w]
        centfit = centfit[w]
        mskbad, coeffs = utils.robust_polyfit(specfit, centfit, traceorder, function=tracefunc, minv=-1.0, maxv=1.0)
        cval[o] = utils.func_val(coeffs, np.array([0.0]), tracefunc, minv=-1.0, maxv=1.0)[0]
        w = np.where(mskbad == 0.0)
        if w[0].size != 0:
            allxfit = np.append(allxfit, specfit[w])
            allsfit = np.append(allsfit, centfit[w]-cval[o])
    if nobj == 0:
        msgs.warn("No objects detected in slit")
        return trace_object_dict(0, None)
    # Tracing
    msgs.info("Performing global trace to all objects")
    mskbad, coeffs = utils.robust_polyfit(allxfit, allsfit, traceorder, function=tracefunc, minv=-1.0, maxv=1.0)
    trcfunc = utils.func_val(coeffs, np.linspace(-1.0, 1.0, skysub.shape[0]), tracefunc, minv=-1.0, maxv=1.0)
    msgs.info("Constructing a trace for all objects")
    trcfunc = trcfunc.reshape((-1, 1)).repeat(nobj, axis=1)
    trccopy = trcfunc.copy()
    for o in range(nobj):
        trcfunc[:, o] += cval[o]
    if nobj == 1:
        msgs.info("Converting object trace to detector pixels")
    else:
        msgs.info("Converting object traces to detector pixels")
    ofst = lordloc[:, slitn].reshape((-1, 1)).repeat(nobj, axis=1) + triml
    diff = (rordloc[:, slitn].reshape((-1, 1)).repeat(nobj, axis=1)
            - lordloc[:, slitn].reshape((-1, 1)).repeat(nobj, axis=1))
    # Convert central trace
    traces = ofst + (diff-triml-trimr)*trcfunc
    # Convert left object trace
    for o in range(nobj):
        trccopy[:, o] = trcfunc[:, o] - cval[o] + objl[o]/(npix-1.0)
    # TODO -- Consider bringing back the next line, as needed
    #trobjl = ofst + (diff-triml-trimr)*trccopy
    # Convert right object trace
    for o in range(nobj):
        trccopy[:, o] = trcfunc[:, o] - cval[o] + objr[o]/(npix-1.0)
    #trobjr = ofst + (diff-triml-trimr)*trccopy
    # Generate an image of pixel weights for each object
    rec_obj_img, rec_bg_img = trace_objbg_image(lordloc, skysub, slitn,
                                                [objl, objr], [bckl, bckr],
                                                triml=triml, trimr=trimr)
    # Check object traces in ginga
    if debug:
        viewer, ch = ginga.show_image(skysub)
        for ii in range(nobj):
            ginga.show_trace(viewer, ch, traces[:, ii], '{:d}'.format(ii), clear=(ii == 0))
        debugger.set_trace()

    # Trace dict
    tracelist = trace_object_dict(nobj, traces, object=rec_obj_img, background=rec_bg_img,
                                  params=tracepar)

    # Return
    return tracelist


def obj_trace_qa(slf, frame, det, tracelist, root='trace', normalize=True, desc=""):
    """ Generate a QA plot for the object traces in a single slit

    Parameters
    ----------
    frame : ndarray
      image
    ltrace : ndarray
      Left edge traces
    rtrace : ndarray
      Right edge traces
    objids : list
    det : int
    slit : int
    desc : str, optional
      Title
    root : str, optional
      Root name for generating output file, e.g. msflat_01blue_000.fits
    normalize : bool, optional
      Normalize the flat?  If not, use zscale for output
    """

    plt.rcdefaults()

    # Grab the named of the method
    method = inspect.stack()[0][3]
    # Outfile name
    outfile = qa.set_qa_filename(slf._basename, method, det=det)
    #
    ycen = np.arange(frame.shape[0])
    # Normalize flux in the traces
    if normalize:
        ntrc = ltrace.shape[1]
        debugger.set_trace() # NO LONGER SUPPORTED
        nrm_frame = np.zeros_like(frame)
        for ii in range(ntrc):
            xtrc = (ltrace[:,ii] + rtrace[:,ii])/2.
            ixtrc = np.round(xtrc).astype(int)
            # Simple 'extraction'
            dumi = np.zeros( (frame.shape[0],3) )
            for jj in range(3):
                dumi[:,jj] = frame[ycen,ixtrc-1+jj]
            trc = np.median(dumi, axis=1)
            # Find portion of the image and normalize
            for yy in ycen:
                xi = max(0, int(ltrace[yy,ii])-3)
                xe = min(frame.shape[1],int(rtrace[yy,ii])+3)
                # Fill + normalize
                nrm_frame[yy, xi:xe] = frame[yy,xi:xe] / trc[yy]
        sclmin, sclmax = 0.4, 1.1
    else:
        nrm_frame = frame.copy()
        sclmin, sclmax = plot.zscale(nrm_frame)

    # Plot
    plt.clf()
    fig = plt.figure(dpi=1200)

    plt.rcParams['font.family'] = 'times new roman'
    ticks_font = font_manager.FontProperties(family='times new roman',
       style='normal', size=16, weight='normal', stretch='normal')
    ax = plt.gca()
    for label in ax.get_yticklabels() :
        label.set_fontproperties(ticks_font)
    for label in ax.get_xticklabels() :
        label.set_fontproperties(ticks_font)
    cmm = cm.Greys_r
    mplt = plt.imshow(nrm_frame, origin='lower', cmap=cmm, extent=(0., frame.shape[1], 0., frame.shape[0]))
    mplt.set_clim(vmin=sclmin, vmax=sclmax)

    # Axes
    plt.xlim(0., frame.shape[1])
    plt.ylim(0., frame.shape[0])
    plt.tick_params(axis='both', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labelleft='off')

    # Traces
    #iy_mid = int(frame.shape[0] / 2.)
    #iy = np.linspace(iy_mid,frame.shape[0]*0.9,ntrc).astype(int)
    for slit, tlist in enumerate(tracelist):
        if 'nobj' not in tlist.keys():
            continue
        for obj_idx in range(tlist['nobj']):
            obj_img = tlist['object'][:,:,obj_idx].astype(int)
            objl = np.argmax(obj_img, axis=1)
            objr = obj_img.shape[1]-np.argmax(np.rot90(obj_img,2), axis=1)-1  # The -1 is for Python indexing
            #
            # Left
            plt.plot(objl, ycen, 'r:', alpha=0.7, lw=0.1)
            # Right
            plt.plot(objr, ycen, 'c:', alpha=0.7, lw=0.1)
        #if objids is not None:
        #    # Label
        #    # plt.text(ltrace[iy,ii], ycen[iy], '{:d}'.format(ii+1), color='red', ha='center')
        #    lbl = 'O{:03d}'.format(objids[ii])
        #    plt.text((ltrace[iy[ii], ii]+rtrace[iy[ii], ii])/2., ycen[iy[ii]],
        #        lbl, color='green', ha='center')
    # Title
    tstamp = qa.gen_timestamp()
    if desc == "":
        plt.suptitle(tstamp)
    else:
        plt.suptitle(desc+'\n'+tstamp)

    plt.savefig(outfile, dpi=800)
    plt.close()

    plt.rcdefaults()


def find_objects(profile, bgreg, stddev):
    """
    Find significantly detected objects in the profile array
    For all objects found, the background regions will be defined.
    """
    # Input profile must be a vector
    if len(profile.shape) != 1:
        raise ValueError('Input profile array must be 1D.')
    # Get the length of the array
    sz_x = profile.size
    # Copy it to a masked array for processing
    # TODO: Assumes input is NOT a masked array
    _profile = np.ma.MaskedArray(profile.copy())

    # Peaks are found as having flux at > 5 sigma and background is
    # where the flux is not greater than 3 sigma
    gt_5_sigma = profile > 5*stddev
    not_gt_3_sigma = np.invert(profile > 3*stddev)

    # Define the object centroids array
    objl = np.zeros(sz_x, dtype=int)
    objr = np.zeros(sz_x, dtype=int)
    has_obj = np.zeros(sz_x, dtype=bool)

    obj = 0
    while np.any(gt_5_sigma & np.invert(_profile.mask)):
        # Find next maximum flux point
        imax = np.ma.argmax(_profile)
        #  IF one is recovering an object mask with *all* pixels triggered,
        #    then consider uncommenting the next 3 lines -- JXP on 20 Apr 2018
        #if not gt_5_sigma[imax]:
        #    break
        #_profile[imax] = np.ma.masked  # Mask the peak
        # Find the valid source pixels around the peak
        f = np.arange(sz_x)[np.roll(not_gt_3_sigma, -imax)]
        # TODO: the ifs below feel like kludges to match old
        # find_objects function.  In particular, should objr be treated
        # as exclusive or inclusive?
        objl[obj] = imax-sz_x+f[-1] if imax-sz_x+f[-1] > 0 else 0
        objr[obj] = f[0]+imax if f[0]+imax < sz_x else sz_x-1
        #        print('object found: ', imax, f[-1], objl[obj], f[0], objr[obj], sz_x)
        # Mask source pixels and increment for next iteration
        has_obj[objl[obj]:objr[obj]+1] = True
        _profile[objl[obj]:objr[obj]+1] = np.ma.masked
        obj += 1

    # The background is the region away from sources up to the provided
    # region size.  Starting pixel for the left limit...
    s = objl[:obj]-bgreg
    s[s < 0] = 0
    # ... and ending pixel for the right limit.
    e = objr[:obj]+1+bgreg
    e[e > sz_x] = sz_x
    # Flag the possible background regions
    bgl = np.zeros((sz_x,obj), dtype=bool)
    bgr = np.zeros((sz_x,obj), dtype=bool)
    for i in range(obj):
        bgl[s[i]:objl[i],i] = True
        bgr[objl[i]+1:e[i],i] = True

    # Return source region limits and background regions that do not
    # have sources
    return objl[:obj], objr[:obj], (bgl & np.invert(has_obj)[:,None]).astype(int), \
           (bgr & np.invert(has_obj)[:,None]).astype(int)


# USE BoxcarFilter instead!
#def mean_weight(array, weight, rejhilo, maskval):
#    _a = array if rejhilo == 0 else np.sort(array)
#    sumw = np.sum(weight[rejhilo:-rejhilo])
#    sumwa = np.sum(weight[rejhilo:-rejhilo]*_a[rejhilo:-rejhilo])
#    return maskval if sumw == 0.0 else sumwa/sumw
#
#
## Weighted boxcar smooothing with rejection
#def new_smooth_x(array, weight, fact, rejhilo, maskval):
#    hf = fact // 2
#
#    sz_x, sz_y = array.shape
#
#    smtarr = np.zeros((sz_x,sz_y), dtype=float)
#    medarr = np.zeros((fact+1), dtype=float)
#    wgtarr = np.zeros((fact+1), dtype=float)
#
#    for y in range(sz_y):
#        for x in range(sz_x):
#            for b in range(fact+1):
#                if (x+b-hf < 0) or (x+b-hf >= sz_x):
#                    wgtarr[b] = 0.0
#                else:
#                    medarr[b] = array[x+b-hf,y]
#                    wgtarr[b] = weight[x+b-hf,y]
#            smtarr[x,y] = mean_weight(medarr, wgtarr, rejhilo, maskval)
#    return smtarr


def tilts_image(tilts, lordloc, rordloc, pad, sz_y):
    """
    Using the tilt (assumed to be fit with a first order polynomial)
    generate an image of the tilts for each slit.

    Parameters
    ----------
    tilts : ndarray
      An (m x n) 2D array specifying the tilt (i.e. gradient) at each
      pixel along the spectral direction (m) for each slit (n).
    lordloc : ndarray
      Location of the left slit edges
    rordloc : ndarray
      Location of the right slit edges
    pad : int
      Set the tilts within each slit, and extend a number of pixels
      outside the slit edges (this number is set by pad).
    sz_y : int
      Number of detector pixels in the spatial direction.

    Returns
    -------
    tiltsimg : ndarray
      An image the same size as the science frame, containing values from 0-1.
      0/1 corresponds to the bottom/top of the detector (in the spectral direction),
      and constant wavelength is represented by a single value from 0-1.

    """
    sz_x, sz_o = tilts.shape
    dszx = (sz_x-1.0)

    tiltsimg = np.zeros((sz_x,sz_y), dtype=float)
    for o in range(sz_o):
        for x in range(sz_x):
            ow = (rordloc[x,o]-lordloc[x,o])/2.0
            oc = (rordloc[x,o]+lordloc[x,o])/2.0
            ymin = int(oc-ow) - pad
            ymax = int(oc+ow) + 1 + pad
            # Check we are in bounds
            if ymin < 0:
                ymin = 0
            elif ymax < 0:
                continue
            if ymax > sz_y-1:
                ymax = sz_y-1
            elif ymin > sz_y-1:
                continue
            # Set the tilt value at each pixel in this row
            for y in range(ymin, ymax):
                yv = (y-lordloc[x, o])/ow - 1.0
                tiltsimg[x,y] = (tilts[x,o]*yv + x)/dszx
    return tiltsimg


def trace_tilt(slf, det, msarc, slitnum, censpec=None, maskval=-999999.9,
               trthrsh=1000.0, nsmth=0, method = "fweight", wv_calib=None):
    """
    This function performs a PCA analysis on the arc tilts for a single spectrum (or order)
               trthrsh=1000.0, nsmth=0):

    Parameters
    ----------
    slf
    det
    msarc
    slitnum : int
      Slit number, here indexed from 0
    censpec
    maskval
    trthrsh
    nsmth
    method : str (fweight or cc)

    Returns
    -------
    trcdict : dict

    """
    def pad_dict(indict):
        """ If an arc line is considered bad, fill the
        dictionary arrays with null values
        """
        indict["xtfit"].append(None)
        indict["ytfit"].append(None)
        indict["wmask"].append(None)
        return indict

    dnum = settings.get_dnum(det)

    msgs.work("Detecting lines for slit {0:d}".format(slitnum+1))
    ordcen = slf._pixcen[det-1].copy()
    tampl, tcent, twid, _, w, _ = arc.detect_lines(censpec)
    satval = settings.spect[dnum]['saturation']*settings.spect[dnum]['nonlinear']
    # Order of the polynomials to be used when fitting the tilts.
    arcdet = (tcent[w]+0.5).astype(np.int)
    ampl = tampl[w]

    # Determine the best lines to use to trace the tilts
    ncont = 15
    aduse = np.zeros(arcdet.size, dtype=np.bool)  # Which lines should be used to trace the tilts
    w = np.where(ampl >= trthrsh)
    aduse[w] = 1
    # Remove lines that are within ncont pixels
    nuse = np.sum(aduse)
    detuse = arcdet[aduse]
    idxuse = np.arange(arcdet.size)[aduse]
    olduse = aduse.copy()
    for s in range(nuse):
        w = np.where((np.abs(arcdet-detuse[s]) <= ncont) & (np.abs(arcdet-detuse[s]) >= 1.0))[0]
        for u in range(w.size):
            if ampl[w[u]] > ampl[olduse][s]:
                aduse[idxuse[s]] = False
                break
    # Restricted to ID lines? [introduced to avoid LRIS ghosts]
    if settings.argflag['trace']['slits']['tilts']['idsonly']:
        ids_pix = np.round(np.array(wv_calib[str(slitnum)]['xfit'])*(msarc.shape[0]-1))
        idxuse = np.arange(arcdet.size)[aduse]
        for s in idxuse:
            if np.min(np.abs(arcdet[s]-ids_pix)) > 2:
                msgs.info("Ignoring line at row={:d}".format(arcdet[s]))
                aduse[s] = False

    # Divide the detector into Nseg segments,
    # and find the brightest lines in each segment.
    # The total number of lines used to trace the tilts will be  = Nseg*Nuse + Nadd
    # Nseg = 4
    # Nuse = 2
    # Nadd = 8
    # segsz = msarc.shape[0]/float(Nseg)
    # aduse = np.zeros(arcdet.size, dtype=np.bool)  # Which lines should be used to trace the tilts
    # for s in range(Nseg):
    #     w = np.where((arcdet > s*segsz) & (arcdet <= (s+1)*segsz))[0]
    #     segampl = tampl[w]
    #     asrt = np.argsort(segampl)[::-1]
    #     for u in range(Nuse):
    #         aduse[w[asrt[u]]] = True
    # # Now include some additional bright lines
    # asrt = np.argsort(tampl)[::-1]
    # s, u = 0, 0
    # while u < Nadd:
    #     if not aduse[asrt[s]]:
    #         aduse[asrt[s]] = True
    #         u += 1
    #     s += 1

    # Setup the trace dictionary
    trcdict = {"xtfit":[], "ytfit":[], "wmask":[], "arcdet":arcdet, "aduse":aduse, "badlines":0}

    msgs.info("Modelling arc line tilts with {0:d} arc lines".format(np.sum(aduse)))
    if np.sum(aduse) == 0:
        msgs.warn("No arc lines were deemed usable in slit {0:d} for spectral tilt".format(slitnum))
        return None
    # Go along each order and trace the tilts
    # Start by masking every row, then later unmask the rows with usable arc lines
    msgs.work("This next step could be multiprocessed to speed up the reduction")
    nspecfit = 3
    badlines = 0
    for j in range(arcdet.size):
        # For each detection in this order
        #msgs.info("Tracing tilt of arc line {0:d}/{1:d}".format(j+1, arcdet.size))
        # Check if this is a saturated line
        ysat = msarc[arcdet[j]-nspecfit:arcdet[j]+nspecfit+1, ordcen[arcdet[j], slitnum]-nsmth:ordcen[arcdet[j], slitnum]+nsmth+1]
        if np.where(ysat > satval)[0].size != 0:
            aduse[j] = False
            badlines += 1
            trcdict = pad_dict(trcdict)
            continue
        # Get the size of the slit
        sz = int(np.floor(np.abs(slf._rordloc[det-1][arcdet[j], slitnum]-slf._lordloc[det-1][arcdet[j], slitnum])/2.0)) - 2
        xtfit = np.zeros(2*sz+1)
        ytfit = np.ones(2*sz+1)*maskval  # Fitted centroid
        etfit = np.zeros(2*sz+1)  # Fitted centroid error
        mtfit = np.ones(2*sz+1, dtype=np.int)   # Mask of bad fits
        #apfit = np.zeros(2*sz+1)  # Fitted Amplitude
        xfit = np.arange(-nspecfit, nspecfit+1, 1.0)
        wfit = np.ones(xfit.size, dtype=np.float)
        tstcc = True  # A boolean to tell the loop once a good set of pixels has been found to cross-correlate with
        # Fit up
        pcen = arcdet[j]
        if (pcen < nspecfit) or (pcen > msarc.shape[0]-(nspecfit+1)):
            # Too close to the end of the spectrum
            aduse[j] = False
            badlines += 1
            trcdict = pad_dict(trcdict)
            continue
        offchip = False
        centv = None
        for k in range(0, sz+1-nsmth):
            if (pcen < nspecfit) or (pcen > msarc.shape[0]-(nspecfit+1)):
                offchip = True
                break
            if ordcen[pcen, slitnum]+k >= msarc.shape[1]:
                offchip = True
                break
            # yfit = msarc[pcen-nspecfit:pcen+nspecfit+1,ordcen[arcdet[j],0]+k]
            yfit = msarc[pcen-nspecfit:pcen+nspecfit+1, ordcen[arcdet[j], slitnum]+k-nsmth:ordcen[arcdet[j], slitnum]+k+nsmth+1]
            if np.size(yfit) == 0:
                offchip = True
                break
            if len(yfit.shape) == 2:
                yfit = np.median(yfit, axis=1)
            # wgd = np.where((yfit<satval)&(yfit!=maskval))
            wgd = np.where(yfit == maskval)
            if wgd[0].size != 0:
                continue
            if method == "fweight":
                if centv is None:
                    centv = np.sum(yfit * (pcen+xfit))/np.sum(yfit)
                wfit[0] = 0.5 + (pcen-centv)
                wfit[-1] = 0.5 - (pcen-centv)
                sumxw = yfit * (pcen+xfit) * wfit
                sumw = yfit * wfit
                centv = np.sum(sumxw)/np.sum(sumw)
                fail = False
            elif method == "cc":
                # Get a copy of the array that will be used to cross-correlate
                if tstcc:
                    # ccyfit = msarc[arcdet[j]-nspecfit:arcdet[j]+nspecfit+1,ordcen[arcdet[j],0]+k]
                    ccyfit = msarc[arcdet[j]-nspecfit:arcdet[j]+nspecfit+1,ordcen[arcdet[j],slitnum]+k-nsmth:ordcen[arcdet[j],slitnum]+k+nsmth+1]
                    if len(ccyfit.shape) == 2:
                        ccyfit = np.median(ccyfit, axis=1)
                    wgd = np.where(ccyfit == maskval)
                    if wgd[0].size != 0:
                        continue
                    ccval = arcdet[j] + np.sum(xfit*ccyfit)/np.sum(ccyfit)
                    tstcc = False  # Once we have an array, there's no need to keep looking
                cc = np.correlate(ccyfit, yfit, mode='same')
                # TODO: This must have been obsolete even before I (KBW)
                # removed the cython functionality!
                params, fail = utils.gauss_lsqfit(xfit, cc, 0.0)
                centv = ccval + pcen - arcdet[j] - params[1]
            xtfit[k+sz] = ordcen[arcdet[j], slitnum] + k
            ytfit[k+sz] = centv
            etfit[k+sz] = 0.02
            #apfit[k+sz] = params[0]
            if fail:
                mtfit[k+sz] = 1
            else:
                pcen = int(0.5 + centv)
                mtfit[k+sz] = 0
        if offchip:
            # Don't use lines that go off the chip (could lead to a bad trace)
            aduse[j] = False
            badlines += 1
            trcdict = pad_dict(trcdict)
            continue
        for k in range(sz+1-nsmth, sz+1):
            xtfit[k+sz] = ordcen[arcdet[j], slitnum]+k
        # Fit down
        pcen = arcdet[j]
        centv = None
        for k in range(1,sz+1-nsmth):
            if (pcen < nspecfit) or (pcen > msarc.shape[0]-(nspecfit+1)):
                offchip = True
                break
            if ordcen[pcen, slitnum]-k < 0:
                offchip = True
                break
            # yfit = msarc[pcen-nspecfit:pcen+nspecfit+1,ordcen[arcdet[j],0]-k]
            yfit = msarc[pcen - nspecfit:pcen + nspecfit + 1,
                   ordcen[arcdet[j], slitnum] - k - nsmth:ordcen[arcdet[j], slitnum] - k + nsmth + 1]
            if len(yfit.shape) == 2:
                yfit = np.median(yfit, axis=1)
            if np.size(yfit) == 0:
                offchip = True
                break
            wgd = np.where(yfit == maskval)
            if wgd[0].size != 0:
                continue
            if method == "fweight":
                if centv is None:
                    centv = np.sum(yfit * (pcen+xfit))/np.sum(yfit)
                wfit[0] = 0.5 + (pcen-centv)
                wfit[-1] = 0.5 - (pcen-centv)
                sumxw = yfit * (pcen+xfit) * wfit
                sumw = yfit * wfit
                centv = np.sum(sumxw)/np.sum(sumw)
                fail = False
            elif method == "cc":
                # Get a copy of the array that will be used to cross-correlate
                # (testcc is probably already False from the Fit Up part of the loop, but it's best to be sure)
                if tstcc:
                    # ccyfit = msarc[arcdet[j]-nspecfit:arcdet[j]+nspecfit+1,ordcen[arcdet[j],0]-k]
                    ccyfit = msarc[arcdet[j]-nspecfit:arcdet[j]+nspecfit+1,ordcen[arcdet[j],slitnum]-k-nsmth:ordcen[arcdet[j],slitnum]-k+nsmth+1]
                    if len(ccyfit.shape) == 2:
                        ccyfit = np.median(ccyfit, axis=1)
                    wgd = np.where(ccyfit == maskval)
                    if wgd[0].size != 0:
                        continue
                    ccval = arcdet[j] + np.sum(xfit*ccyfit)/np.sum(ccyfit)
                    tstcc = False  # Once we have an array, there's no need to keep looking
                cc = np.correlate(ccyfit, yfit, mode='same')
                # TODO: This must have been obsolete even before I (KBW)
                # removed the cython functionality!
                params, fail = utils.gauss_lsqfit(xfit, cc, 0.0)
                centv = ccval + pcen - arcdet[j] - params[1]
            xtfit[sz-k] = ordcen[arcdet[j], slitnum] - k
            ytfit[sz-k] = centv
            etfit[sz-k] = 0.02
            #apfit[sz-k] = params[0]
            if fail:
                mtfit[sz-k] = 1
            else:
                #from IPython import embed
                if np.isfinite(centv) is False: debugger.set_trace() #embed()
                pcen = int(0.5 + centv)
                mtfit[sz-k] = 0
        '''
        jxp_fix = False
        if jxp_fix:
            from desispec.bootcalib import trace_crude_init
            from desispec.bootcalib import trace_fweight as dbtf
            from pypeit import ginga
            pcen = arcdet[j]
            img = msarc[pcen-nspecfit:pcen+nspecfit+1, ordcen[arcdet[j], slitnum]-sz:ordcen[arcdet[j], slitnum]+sz+1]
            rot_img = np.rot90(img,3) - 980.  # Bias!
            #
            xcen = np.array([6.6]) # 8.5, 173
            ypass = 173
            # Tune-up first
            for ii in range(4):
                xcen, xsig = dbtf(rot_img, xcen, ycen=np.array([ypass]).astype(int), invvar=None, radius=2.)
            xset, xerr = trace_crude_init(rot_img, xcen, ypass, invvar=None, radius=3., maxshift0=0.5, maxshift=0.15, maxerr=0.2)
            #xcen, xsig = dbtf(rot_img, np.array([6.5,6.5]), ycen=np.array([173,173]).astype(int), invvar=None, radius=2.)
            # Convert back
            y0 = pcen+nspecfit
            ycrude = y0 - xset[:,0]
            #debugger.set_trace()
            #cytfit = ytfit.copy()
            #cytfit[np.where(ytfit < 0)] = np.median(cytfit)
            #debugger.xplot(np.arange(img.shape[1]), ycrude, cytfit)
            #
            #trcdict['save_yt'] = ytfit.copy()
            assert ycrude.size == ytfit.size
            ytfit = ycrude
            #debugger.set_trace()
        '''

        if offchip:
            # Don't use lines that go off the chip (could lead to a bad trace)
            aduse[j] = False
            badlines += 1
            trcdict = pad_dict(trcdict)
            continue
        for k in range(sz+1-nsmth, sz+1):
            xtfit[sz-k] = ordcen[arcdet[j], slitnum]-k

        wmask = np.where(mtfit == 0)
        ytfit[np.where(mtfit == 1)] = maskval

        # Append the trace information into the dictionary
        trcdict["xtfit"].append(xtfit.copy())
        trcdict["ytfit"].append(ytfit.copy())
        trcdict["wmask"].append(wmask[0].copy())
    trcdict["aduse"] = aduse
    trcdict["badlines"] = badlines
    msgs.info("Completed spectral tilt tracing".format(np.sum(aduse)))
    return trcdict


def trace_weighted(frame, ltrace, rtrace, mask=None, wght="flux"):
    """ Estimate the trace of an object in a single slit,
    weighted by the specified method.

    Parameters:
    -----------
    frame : 2D ndarray
      Image for tracing
    ltrace : ndarray
      Left slit edge trace
    rtrace : ndarray
      Right slit edge trace
    mask : ndarray, optional
      Mask of pixels to ignore while tracing
    wght : str
      Method that should be used to weight the tracing:
        wght="flux"  will weight the trace by the flux of each pixel
        wght="uniform"  will use uniform weights

    Returns:
    --------
    trace : ndarray
      array containing the trace of the object along the slit
    error : ndarray
      the associated error of the object trace
    """

    nspec, nspat = frame.shape
    lidx = int(np.ceil(np.min(ltrace)))-1
    ridx = int(np.floor(np.max(rtrace)))+1
    if lidx < 0:
        lidx = 0
    elif lidx >= nspat:
        msgs.info("Slit is off the detector - not tracing")
        return None, None
    if ridx >= nspat:
        ridx = nspat-1
    elif ridx <= 0:
        msgs.info("Slit is off the detector - not tracing")
        return None, None
    extfram = frame[:, lidx:ridx+1]
    if isinstance(wght, (str, ustr)):
        if wght == "flux":
            extwght = extfram.copy()
        elif wght == "uniform":
            extwght = np.ones(extfram.shape, dtype=np.float)
    else:
        # Assume extwght is a numpy array
        extwght = wght
    # Calculate the weights
    idxarr = np.outer(np.ones(nspec), np.arange(lidx, ridx+1))
    ldiff = idxarr - ltrace.reshape(nspec, 1)
    rdiff = idxarr - rtrace.reshape(nspec, 1) - 1.0
    msgs.work("Not sure why -1 is needed here, might be a bug about how ltrace and rtrace are defined")
    lwhr = np.where(ldiff <= -1.0)
    rwhr = np.where(rdiff >= +0.0)
    extwght[lwhr] = 0.0
    extwght[rwhr] = 0.0
    lwhr = np.where((-1.0 < ldiff) & (ldiff <= 0.0))
    rwhr = np.where((-1.0 <= rdiff) & (rdiff < 0.0))
    extwght[lwhr] *= 1.0+ldiff[lwhr]
    extwght[rwhr] *= -rdiff[rwhr]
    # Apply mask
    if mask is not None:
        extwght[np.where(mask[:, lidx:ridx+1] != 0)] *= 0.0
    # Determine the weighted trace
    wghtsum = np.sum(extwght, axis=1)
    wfact = 1.0/(wghtsum + (wghtsum == 0.0))
    trace = np.sum(idxarr*extwght, axis=1) * wfact * (wghtsum != 0)
    error = np.sqrt(np.sum((idxarr-trace.reshape(nspec, 1))**2 * extwght, axis=1) * wfact)
    return trace, error

# This should be deprecated, since there is another version of trace_fweight in core.artraceslits

#def trace_fweight(fimage, xinit, ltrace=None, rtraceinvvar=None, radius=3.):
def trace_fweight_deprecated(fimage, xinit, ltrace=None, rtraceinvvar=None, radius=3.):
    """ Python port of trace_fweight.pro from IDLUTILS

    Parameters:
    -----------
    fimage: 2D ndarray
      Image for tracing
    xinit: ndarray
      Initial guesses for x-trace
    invvar: ndarray, optional
      Inverse variance array for the image
    radius: float, optional
      Radius for centroiding; default to 3.0
    """

    # Init
    nx = fimage.shape[1]
    ny = fimage.shape[0]
    ncen = len(xinit)
    xnew = copy.deepcopy(xinit)
    xerr = np.zeros(ncen) + 999.

    ycen = np.arange(ny, dtype=int)
    invvar = 0. * fimage + 1.
    x1 = xinit - radius + 0.5
    x2 = xinit + radius + 0.5
    ix1 = np.floor(x1).astype(int)
    ix2 = np.floor(x2).astype(int)

    fullpix = int(np.maximum(np.min(ix2-ix1)-1, 0))
    sumw = np.zeros(ny)
    sumxw = np.zeros(ny)
    sumwt = np.zeros(ny)
    sumsx1 = np.zeros(ny)
    sumsx2 = np.zeros(ny)
    qbad = np.array([False]*ny) 

    if invvar is None: 
        invvar = np.zeros_like(fimage) + 1. 

    # Compute
    for ii in range(0,fullpix+3):
        spot = ix1 - 1 + ii
        ih = np.clip(spot,0,nx-1)
        xdiff = spot - xinit
        #
        wt = np.clip(radius - np.abs(xdiff) + 0.5,0,1) * ((spot >= 0) & (spot < nx))
        sumw = sumw + fimage[ycen,ih] * wt
        sumwt = sumwt + wt
        sumxw = sumxw + fimage[ycen,ih] * xdiff * wt
        var_term = wt**2 / (invvar[ycen,ih] + (invvar[ycen,ih] == 0))
        sumsx2 = sumsx2 + var_term
        sumsx1 = sumsx1 + xdiff**2 * var_term
        #qbad = qbad or (invvar[ycen,ih] <= 0)
        qbad = np.any([qbad, invvar[ycen,ih] <= 0], axis=0)

    # Fill up
    good = (sumw > 0) & (~qbad)
    if np.sum(good) > 0:
        delta_x = sumxw[good]/sumw[good]
        xnew[good] = delta_x + xinit[good]
        xerr[good] = np.sqrt(sumsx1[good] + sumsx2[good]*delta_x**2)/sumw[good]

    bad = np.any([np.abs(xnew-xinit) > radius + 0.5, xinit < radius - 0.5, xinit > nx - 0.5 - radius], axis=0)
    if np.sum(bad) > 0:
        xnew[bad] = xinit[bad]
        xerr[bad] = 999.0

    # Return
    return xnew, xerr

'''
def echelle_tilt(slf, msarc, det, pcadesc="PCA trace of the spectral tilts", maskval=-999999.9):
    """ Determine the spectral tilts in each echelle order

    Parameters
    ----------
    slf : Class instance
      An instance of the Science Exposure class
    msarc : numpy ndarray
      Wavelength calibration frame that will be used to trace constant wavelength
    det : int
      Index of the detector
    pcadesc : str (optional)
      A description of the tilts
    maskval : float (optional)
      Mask value used in numpy arrays

    Returns
    -------
    tilts : ndarray
      Locations of the left slit edges (in physical pixel coordinates)
    satmask : ndarray
      Locations of the right slit edges (in physical pixel coordinates)
    extrapord : ndarray
      A boolean mask indicating if an order was extrapolated (True = extrapolated)
    """
    arccen, maskslit, satmask = arc.get_censpec(slf._lordloc[det-1], slf._rordloc[det-1],
                                                      slf._pixlocn[det-1], msarc, det, settings.spect, gen_satmask=True)
    # If the user sets no tilts, return here
    if settings.argflag['trace']['slits']['tilts']['method'].lower() == "zero":
        # Assuming there is no spectral tilt
        tilts = np.outer(np.linspace(0.0, 1.0, msarc.shape[0]), np.ones(msarc.shape[1]))
        return tilts, satmask, None
    norders = maskslit.size

    # Now model the tilt for each slit
    tiltang, centval = None, None
    slitcnt = 0
    for o in range(norders):
        if maskslit[o] == 1:
            continue
        # Determine the tilts for this slit
        trcdict = trace_tilt(slf, det, msarc, o, censpec=arccen[:, o])
        slitcnt += 1
        if trcdict is None:
            # No arc lines were available to determine the spectral tilt
            continue
        aduse = trcdict["aduse"]
        arcdet = trcdict["arcdet"]
        if tiltang is None:
            tiltang = maskval * np.ones((aduse.size, norders))
            centval = maskval * np.ones((aduse.size, norders))
            totnum = aduse.size
        else:
            if aduse.size > totnum:
                tiltang = np.append(tiltang, maskval * np.ones((aduse.size-totnum, norders)), axis=0)
                centval = np.append(centval, maskval * np.ones((aduse.size-totnum, norders)), axis=0)
                totnum = aduse.size
        for j in range(aduse.size):
            if not aduse[j]:
                continue
            lrdiff = slf._rordloc[det-1][arcdet[j], o] - slf._lordloc[det-1][arcdet[j], o]
            xtfit = 2.0 * (trcdict["xtfit"][j]-slf._lordloc[det-1][arcdet[j], o])/lrdiff - 1.0
            ytfit = trcdict["ytfit"][j]
            wmask = trcdict["wmask"][j]
            if wmask.size < settings.argflag['trace']['slits']['tilts']['order']+2:
                maskslit[o] = 1
                continue
            null, tcoeff = utils.robust_polyfit(xtfit[wmask], ytfit[wmask],
                                                  settings.argflag['trace']['slits']['tilts']['order'], sigma=2.0)
            # Save the tilt angle
            tiltang[j, o] = tcoeff[1]  # tan(tilt angle)
            centval[j, o] = tcoeff[0]  # centroid of arc line

    msgs.info("Fitting tilt angles")
    tcoeff = np.ones((settings.argflag['trace']['slits']['tilts']['disporder'] + 1, tiltang.shape[1]))
    maskord = np.where(maskslit == 1)[0]
    extrap_ord = np.zeros(norders)
    for o in range(norders):
        if o in maskord:
            extrap_ord[o] = 1
            continue
        w = np.where(tiltang[:, o] != maskval)
        if np.size(w[0]) <= settings.argflag['trace']['slits']['tilts']['disporder'] + 2:
            extrap_ord[o] = 1.0
            maskord = np.append(maskord, o)
        else:
            null, tempc = utils.robust_polyfit(centval[:, o][w], tiltang[:, o][w],
                                                 settings.argflag['trace']['slits']['tilts']['disporder'],
                                                 function=settings.argflag['trace']['slits']['function'], sigma=2.0,
                                                 minv=0.0, maxv=msarc.shape[0] - 1)
            tcoeff[:, o] = tempc
    # Sort which orders are masked
    maskord.sort()
    xv = np.arange(msarc.shape[0])
    tiltval = utils.func_val(tcoeff, xv, settings.argflag['trace']['slits']['function'],
                               minv=0.0, maxv=msarc.shape[0] - 1).T
    ofit = settings.argflag['trace']['slits']['tilts']['params']
    lnpc = len(ofit) - 1
    if np.sum(1.0 - extrap_ord) > ofit[0] + 1:  # Only do a PCA if there are enough good orders
        # Perform a PCA on the tilts
        msgs.info("Performing a PCA on the spectral tilts")
        ordsnd = np.arange(norders) + 1.0
        xcen = xv[:, np.newaxis].repeat(norders, axis=1)
        fitted, outpar = arpca.basis(xcen, tiltval, tcoeff, lnpc, ofit, x0in=ordsnd, mask=maskord, skipx0=False,
                                     function=settings.argflag['trace']['slits']['function'])
        if not msgs._debug['no_qa']:
            #pcadesc = "Spectral Tilt PCA"
#            arqa.pca_plot(slf, outpar, ofit, 'Arc', pcadesc=pcadesc, addOne=False)
            pca.pca_plot(slf, outpar, ofit, 'Arc', pcadesc=pcadesc, addOne=False)
        # Extrapolate the remaining orders requested
        orders = 1.0 + np.arange(norders)
        extrap_tilt, outpar = pca.extrapolate(outpar, orders, function=settings.argflag['trace']['slits']['function'])
        tilts = extrap_tilt
#        arqa.pca_arctilt(slf, tiltang, centval, tilts)
        pca.pca_arctilt(slf, tiltang, centval, tilts)
    else:
        outpar = None
        msgs.warn("Could not perform a PCA when tracing the order tilts" + msgs.newline() +
                  "Not enough well-traced orders")
        msgs.info("Attempting to fit tilts by assuming the tilt is order-independent")
        xtiltfit = np.array([])
        ytiltfit = np.array([])
        for o in range(tiltang.shape[1]):
            w = np.where(tiltang[:, o] != -999999.9)
            if np.size(w[0]) != 0:
                xtiltfit = np.append(xtiltfit, centval[:, o][w])
                ytiltfit = np.append(ytiltfit, tiltang[:, o][w])
        if np.size(xtiltfit) > settings.argflag['trace']['slits']['tilts']['disporder'] + 2:
            tcoeff = utils.func_fit(xtiltfit, ytiltfit, settings.argflag['trace']['slits']['function'],
                                      settings.argflag['trace']['slits']['tilts']['disporder'],
                                      minv=0.0, maxv=msarc.shape[0] - 1)
            tiltval = utils.func_val(tcoeff, xv, settings.argflag['trace']['slits']['function'], minv=0.0,
                                       maxv=msarc.shape[0] - 1)
            tilts = tiltval[:, np.newaxis].repeat(tiltang.shape[1], axis=1)
        else:
            msgs.warn("There are still not enough detections to obtain a reliable tilt trace")
            msgs.info("Assuming there is no tilt")
            tilts = np.zeros_like(slf._lordloc)

    # Generate tilts image
    tiltsimg = tilts_image(tilts, slf._lordloc[det-1], slf._rordloc[det-1],
                           settings.argflag['trace']['slits']['pad'], msarc.shape[1])

    return tiltsimg, satmask, outpar

def multislit_tilt(slf, msarc, det, maskval=-999999.9, doqa=False):
    """ Determine the spectral tilt of each slit in a multislit image

    Parameters
    ----------
    slf : Class instance
      An instance of the Science Exposure class
    msarc : numpy ndarray
      Wavelength calibration frame that will be used to trace constant wavelength
    det : int
      Index of the detector
    maskval : float (optional)
      Mask value used in numpy arrays
    doqa : bool, optional
      Output QA files.  These can be many files and slow for
      lots of slits

    Returns
    -------
    tilts : ndarray
      Locations of the left slit edges (in physical pixel coordinates)
    satmask : ndarray
      Locations of the right slit edges (in physical pixel coordinates)
    extrapord : ndarray
      A boolean mask indicating if an order was extrapolated (True = extrapolated)
    """
    arccen, maskslit, _ = arc.get_censpec(slf._lordloc[det-1], slf._rordloc[det-1],
                                                      slf._pixlocn[det-1], msarc, det, settings.spect,
                                                      gen_satmask=False)
    satmask = np.zeros_like(slf._pixcen)
    # If the user sets no tilts, return here
    if settings.argflag['trace']['slits']['tilts']['method'].lower() == "zero":
        # Assuming there is no spectral tilt
        tilts = np.outer(np.linspace(0.0, 1.0, msarc.shape[0]), np.ones(msarc.shape[1]))
        return tilts, satmask, None

    ordcen = slf._pixcen[det - 1].copy()
    fitxy = [settings.argflag['trace']['slits']['tilts']['order'], 1]

    # maskslit
    if slf._maskslits[det-1] is not None:
        mask = slf._maskslits[det-1] & (maskslit==1)
    else:
        mask = maskslit
    slf._maskslits[det-1] = mask
    gdslits = np.where(mask == 0)[0]

    # Final tilts image
    final_tilts = np.zeros_like(msarc)

    # Now trace the tilt for each slit
    #for  o in range(arccen.shape[1]):
    for slit in gdslits:
        # Determine the tilts for this slit
        trcdict = trace_tilt(slf, det, msarc, slit, censpec=arccen[:, slit], nsmth=3,
                             wv_calib=wv_calib)
        if trcdict is None:
            # No arc lines were available to determine the spectral tilt
            continue
        if msgs._debug['tilts']:
            debugger.chk_arc_tilts(msarc, trcdict, sedges=(slf._lordloc[det-1][:,slit], slf._rordloc[det-1][:,slit]))
            debugger.set_trace()
        # Extract information from the trace dictionary
        aduse = trcdict["aduse"]
        arcdet = trcdict["arcdet"]
        xtfits = trcdict["xtfit"]
        ytfits = trcdict["ytfit"]
        wmasks = trcdict["wmask"]
        badlines = trcdict["badlines"]
        # Initialize some arrays
        maskrows = np.ones(msarc.shape[0], dtype=np.int)
        tcoeff = np.ones((settings.argflag['trace']['slits']['tilts']['order'] + 1, msarc.shape[0]))
        xtilt = np.ones((msarc.shape[1], arcdet.size)) * maskval
        ytilt = np.ones((msarc.shape[1], arcdet.size)) * maskval
        ztilt = np.ones((msarc.shape[1], arcdet.size)) * maskval
        mtilt = np.ones((msarc.shape[1], arcdet.size)) * maskval
        wtilt = np.ones((msarc.shape[1], arcdet.size)) * maskval
        # Analyze each spectral line
        for j in range(arcdet.size):
            if not aduse[j]:
                continue
            xtfit = xtfits[j]
            ytfit = ytfits[j]
            wmask = wmasks[j]
            xint = int(xtfit[0])
            sz = (xtfit.size-1)//2

            # Trim if we are off the detector
            lastx = min(xint + 2 * sz + 1, msarc.shape[1])
            if (lastx-xint) < xtfit.size: # Cut down
                dx = (lastx-xint)-xtfit.size
                xtfit = xtfit[:dx]
                ytfit = ytfit[:dx]
                wmask = wmask[np.where(wmask < (xtfit.size+dx))]

            # Perform a scanning polynomial fit to the tilts
            wmfit = np.where(ytfit != maskval)
            if wmfit[0].size > settings.argflag['trace']['slits']['tilts']['order'] + 1:
                cmfit = utils.func_fit(xtfit[wmfit], ytfit[wmfit],
                                         settings.argflag['trace']['slits']['function'],
                                         settings.argflag['trace']['slits']['tilts']['order'],
                                         minv=0.0, maxv=msarc.shape[1] - 1.0)
                model = utils.func_val(cmfit, xtfit, settings.argflag['trace']['slits']['function'],
                                         minv=0.0, maxv=msarc.shape[1] - 1.0)
            else:
                aduse[j] = False
                badlines += 1
                continue

            if maskval in model:
                # Model contains masked values
                aduse[j] = False
                badlines += 1
                continue

            # Perform a robust polynomial fit to the traces
            if settings.argflag['trace']['slits']['tilts']['method'].lower() == "spca":
                yfit = ytfit[wmask] / (msarc.shape[0] - 1.0)
            else:
                yfit = (2.0 * model[sz] - ytfit[wmask]) / (msarc.shape[0] - 1.0)
            wmsk, mcoeff = utils.robust_polyfit(xtfit[wmask], yfit,
                                                  settings.argflag['trace']['slits']['tilts']['order'],
                                                  function=settings.argflag['trace']['slits']['function'],
                                                  sigma=2.0, minv=0.0, maxv=msarc.shape[1] - 1.0)
            # Update the mask
            wmask = wmask[np.where(wmsk == 0)]

            # Save the tilt angle, and unmask the row
            factr = (msarc.shape[0] - 1.0) * utils.func_val(mcoeff, ordcen[arcdet[j], slit],
                                                              settings.argflag['trace']['slits']['function'],
                                                              minv=0.0, maxv=msarc.shape[1] - 1.0)
            idx = int(factr + 0.5)
            if (idx > 0) and (idx < msarc.shape[0]):
                maskrows[idx] = 0
                tcoeff[:, idx] = mcoeff.copy()
            # Restrict to good IDs?
            if settings.argflag['trace']['slits']['tilts']['idsonly']:
                if not aduse[j]:
                    maskrows[idx] = 1

            xtilt[xint:lastx, j] = xtfit / (msarc.shape[1] - 1.0)
            ytilt[xint:lastx, j] = arcdet[j] / (msarc.shape[0] - 1.0)
            ztilt[xint:lastx, j] = ytfit / (msarc.shape[0] - 1.0)
            if settings.argflag['trace']['slits']['tilts']['method'].lower() == "spline":
                mtilt[xint:lastx, j] = model / (msarc.shape[0] - 1.0)
            elif settings.argflag['trace']['slits']['tilts']['method'].lower() == "interp":
                mtilt[xint:lastx, j] = (2.0 * model[sz] - model) / (msarc.shape[0] - 1.0)
            else:
                mtilt[xint:lastx, j] = (2.0 * model[sz] - model) / (msarc.shape[0] - 1.0)
            wbad = np.where(ytfit == maskval)[0]
            ztilt[xint + wbad, j] = maskval
            if wmask.size != 0:
                sigg = max(1.4826 * np.median(np.abs(ytfit - model)[wmask]) / np.sqrt(2.0), 1.0)
                wtilt[xint:lastx, j] = 1.0 / sigg
            # Extrapolate off the slit to the edges of the chip
            nfit = 6  # Number of pixels to fit a linear function to at the end of each trace
            xlof, xhif = np.arange(xint, xint + nfit), np.arange(lastx - nfit, lastx)
            xlo, xhi = np.arange(xint), np.arange(lastx, msarc.shape[1])
            glon = np.mean(xlof * mtilt[xint:xint + nfit, j]) - np.mean(xlof) * np.mean(mtilt[xint:xint + nfit, j])
            glod = np.mean(xlof ** 2) - np.mean(xlof) ** 2
            clo = np.mean(mtilt[xint:xint + nfit, j]) - (glon / glod) * np.mean(xlof)
            yhi = mtilt[lastx - nfit:lastx, j]
            try:
                ghin = np.mean(xhif * yhi) - np.mean(xhif) * np.mean(yhi)
            except:
                debugger.set_trace()
            ghid = np.mean(xhif ** 2) - np.mean(xhif) ** 2
            chi = np.mean(yhi) - (ghin / ghid) * np.mean(xhif)
            mtilt[0:xint, j] = (glon / glod) * xlo + clo
            mtilt[lastx:, j] = (ghin / ghid) * xhi + chi
        if badlines != 0:
            msgs.warn("There were {0:d} additional arc lines that should have been traced".format(badlines) +
                      msgs.newline() + "(perhaps lines were saturated?). Check the spectral tilt solution")

        # Masking
        maskrw = np.where(maskrows == 1)[0]
        maskrw.sort()
        extrap_row = maskrows.copy()
        xv = np.arange(msarc.shape[1])
        # Tilt values
        tiltval = utils.func_val(tcoeff, xv, settings.argflag['trace']['slits']['function'],
                                   minv=0.0, maxv=msarc.shape[1] - 1.0).T
        msgs.work("May need to do a check here to make sure ofit is reasonable")
        ofit = settings.argflag['trace']['slits']['tilts']['params']
        lnpc = len(ofit) - 1
        # Only do a PCA if there are enough good orders
        if np.sum(1.0 - extrap_row) > ofit[0] + 1:
            # Perform a PCA on the tilts
            msgs.info("Performing a PCA on the tilts")
            ordsnd = np.linspace(0.0, 1.0, msarc.shape[0])
            xcen = xv[:, np.newaxis].repeat(msarc.shape[0], axis=1)
            fitted, outpar = arpca.basis(xcen, tiltval, tcoeff, lnpc, ofit, weights=None,
                                         x0in=ordsnd, mask=maskrw, skipx0=False,
                                         function=settings.argflag['trace']['slits']['function'])
#            qa.pca_plot(slf, outpar, ofit, 'Arc', pcadesc="Spectral Tilt PCA", addOne=False)
            arpca.pca_plot(slf, outpar, ofit, 'Arc', pcadesc="Spectral Tilt PCA", addOne=False)
            # Extrapolate the remaining orders requested
            orders = np.linspace(0.0, 1.0, msarc.shape[0])
            extrap_tilt, outpar = arpca.extrapolate(outpar, orders, function=settings.argflag['trace']['slits']['function'])
            polytilts = extrap_tilt.T
        else:
            # Fit the model with a 2D polynomial
            msgs.warn("Could not perform a PCA when tracing the spectral tilt" + msgs.newline() +
                      "Not enough well-traced arc lines")
            msgs.info("Fitting tilts with a low order, 2D polynomial")
            wgd = np.where(xtilt != maskval)
            coeff = arutils.polyfit2d_general(xtilt[wgd], ytilt[wgd], mtilt[wgd], fitxy)
            polytilts = arutils.polyval2d_general(coeff, np.linspace(0.0, 1.0, msarc.shape[1]),
                                                  np.linspace(0.0, 1.0, msarc.shape[0]))

        if settings.argflag['trace']['slits']['tilts']['method'].lower() == "interp":
            msgs.info("Interpolating and Extrapolating the tilts")
            xspl = np.linspace(0.0, 1.0, msarc.shape[1])
            # yspl = np.append(0.0, np.append(arcdet[np.where(aduse)]/(msarc.shape[0]-1.0), 1.0))
            # yspl = np.append(0.0, np.append(polytilts[arcdet[np.where(aduse)], msarc.shape[1]/2], 1.0))
            ycen = np.diag(polytilts[arcdet[np.where(aduse)], ordcen[arcdet[np.where(aduse), slit]]])
            yspl = np.append(0.0, np.append(ycen, 1.0))
            zspl = np.zeros((msarc.shape[1], np.sum(aduse) + 2))
            zspl[:, 1:-1] = mtilt[:, np.where(aduse)[0]]
            # zspl[:, 1:-1] = polytilts[arcdet[np.where(aduse)[0]], :].T
            zspl[:, 0] = zspl[:, 1] + polytilts[0, :] - polytilts[arcdet[np.where(aduse)[0][0]], :]
            zspl[:, -1] = zspl[:, -2] + polytilts[-1, :] - polytilts[arcdet[np.where(aduse)[0][-1]], :]
            # Make sure the endpoints are set to 0.0 and 1.0
            zspl[:, 0] -= zspl[ordcen[0, 0], slit]
            zspl[:, -1] = zspl[:, -1] - zspl[ordcen[-1, slit], -1] + 1.0
            # Prepare the spline variables
            # if False:
            #     pmin = 0
            #     pmax = -1
            # else:
            #     pmin = int(max(0, np.min(slf._lordloc[det-1])))
            #     pmax = int(min(msarc.shape[1], np.max(slf._rordloc[det-1])))
            # xsbs = np.outer(xspl, np.ones(yspl.size))
            # ysbs = np.outer(np.ones(xspl.size), yspl)
            # zsbs = zspl[wgd]
            # Restrict to good portion of the image
            tiltspl = interpolate.RectBivariateSpline(xspl, yspl, zspl, kx=3, ky=3)
            yval = np.linspace(0.0, 1.0, msarc.shape[0])
            tilts = tiltspl(xspl, yval, grid=True).T
        elif settings.argflag['trace']['slits']['tilts']['method'].lower() == "spline":
            msgs.info("Performing a spline fit to the tilts")
            wgd = np.where((ytilt != maskval) & (ztilt != maskval))
            txsbs = xtilt[wgd]
            tysbs = ytilt[wgd]
            tzsbs = ztilt[wgd]
            twsbs = wtilt[wgd]
            # Append the end points
            wlo = np.where((ytilt == np.min(tysbs)) & (ytilt != maskval) & (ztilt != maskval))
            whi = np.where((ytilt == np.max(tysbs)) & (ytilt != maskval) & (ztilt != maskval))
            xlo = (xtilt[wlo] * (msarc.shape[1] - 1.0)).astype(np.int)
            xhi = (xtilt[whi] * (msarc.shape[1] - 1.0)).astype(np.int)
            xsbs = np.append(xtilt[wlo], np.append(txsbs, xtilt[whi]))
            ysbs = np.append(np.zeros(wlo[0].size), np.append(tysbs, np.ones(whi[0].size)))
            zlo = ztilt[wlo] + polytilts[0, xlo] - polytilts[arcdet[np.where(aduse)[0][0]], xlo]
            zhi = ztilt[whi] + polytilts[-1, xhi] - polytilts[arcdet[np.where(aduse)[0][-1]], xhi]
            zsbs = np.append(zlo, np.append(tzsbs, zhi))
            wsbs = np.append(wtilt[wlo], np.append(twsbs, wtilt[whi]))
            # Generate the spline curve
            tiltspl = interpolate.SmoothBivariateSpline(xsbs, zsbs, ysbs, w=wsbs, kx=3, ky=3,
                                                        s=xsbs.size, bbox=[0.0, 1.0,
                                                        min(zsbs.min(),0.0), max(zsbs.max(),1.0)])
            xspl = np.linspace(0.0, 1.0, msarc.shape[1])
            yspl = np.linspace(0.0, 1.0, msarc.shape[0])
            tilts = tiltspl(xspl, yspl, grid=True).T
            # QA
            if msgs._debug['tilts']:
                tiltqa = tiltspl(xsbs, zsbs, grid=False)
                plt.clf()
                # plt.imshow((zsbs-tiltqa)/zsbs, origin='lower')
                # plt.imshow((ysbs-tiltqa)/ysbs, origin='lower')
                plt.plot(xsbs, (ysbs - tiltqa) / ysbs, 'bx')
                plt.plot(xsbs, 1.0 / (wsbs * ysbs), 'r-')
                plt.ylim(-5e-2, 5e-2)
                # plt.colorbar()
                plt.show()
                debugger.set_trace()
        elif settings.argflag['trace']['slits']['tilts']['method'].lower() == "spca":
            # Slit position
            xspl = np.linspace(0.0, 1.0, msarc.shape[1])
            # Trace positions down center of the order
            ycen = np.diag(polytilts[arcdet[np.where(aduse)], ordcen[arcdet[np.where(aduse)], slit:slit+1]])
            yspl = np.append(0.0, np.append(ycen, 1.0))
            # Trace positions as measured+modeled
            zspl = np.zeros((msarc.shape[1], np.sum(aduse) + 2))
            zspl[:, 1:-1] = polytilts[arcdet[np.where(aduse)[0]], :].T
            zspl[:, 0] = zspl[:, 1] + polytilts[0, :] - polytilts[arcdet[np.where(aduse)[0][0]], :]
            zspl[:, -1] = zspl[:, -2] + polytilts[-1, :] - polytilts[arcdet[np.where(aduse)[0][-1]], :]
            # Make sure the endpoints are set to 0.0 and 1.0
            zspl[:, 0] -= zspl[ordcen[0, slit], 0]
            zspl[:, -1] = zspl[:, -1] - zspl[ordcen[-1, slit], -1] + 1.0
            # Prepare the spline variables
            if False:
                pmin = 0
                pmax = -1
            else:
                pmin = int(max(0, np.min(slf._lordloc[det - 1])))
                pmax = int(min(msarc.shape[1], np.max(slf._rordloc[det - 1])))
            xsbs = np.outer(xspl, np.ones(yspl.size))[pmin:pmax, :]
            ysbs = np.outer(np.ones(xspl.size), yspl)[pmin:pmax, :]
            zsbs = zspl[pmin:pmax, :]
            # Spline
            msgs.work('Consider adding weights to SmoothBivariate in spca')
            tiltspl = interpolate.SmoothBivariateSpline(xsbs.flatten(), zsbs.flatten(),
                                                        ysbs.flatten(), kx=3, ky=3, s=xsbs.size)
            # Finish
            yval = np.linspace(0.0, 1.0, msarc.shape[0])
            tilts = tiltspl(xspl, yval, grid=True).T
            if False:
                tiltqa = tiltspl(xsbs.flatten(), zsbs.flatten(), grid=False).reshape(xsbs.shape)
                plt.clf()
                # plt.imshow((zsbs-tiltqa)/zsbs, origin='lower')
                plt.imshow((ysbs - tiltqa) / ysbs, origin='lower')
                plt.colorbar()
                plt.show()
                debugger.set_trace()
        elif settings.argflag['trace']['slits']['tilts']['method'].lower() == "pca":
            tilts = polytilts.copy()
        # Save into final_tilts
        word = np.where(slf._slitpix[det - 1] == slit+1)
        final_tilts[word] = tilts[word]

        # Now do the QA
        if doqa:
            msgs.info("Preparing arc tilt QA data")
            tiltsplot = tilts[arcdet, :].T
            tiltsplot *= (msarc.shape[0] - 1.0)
            # Shift the plotted tilts about the centre of the slit
            ztilto = ztilt.copy()
            adj = np.diag(tilts[arcdet, ordcen[arcdet, slit:slit+1]])
            zmsk = np.where(ztilto == maskval)
            ztilto = 2.0 * np.outer(np.ones(ztilto.shape[0]), adj) - ztilto
            ztilto[zmsk] = maskval
            ztilto[np.where(ztilto != maskval)] *= (msarc.shape[0] - 1.0)
            for i in range(arcdet.size):
                w = np.where(ztilto[:, i] != maskval)
                if w[0].size != 0:
                    twa = (xtilt[w[0], i] * (msarc.shape[1] - 1.0) + 0.5).astype(np.int)
                    # fitcns = np.polyfit(w[0], ztilt[w[0], i] - tiltsplot[twa, i], 0)[0]
                    fitcns = np.median(ztilto[w[0], i] - tiltsplot[twa, i])
                    # if abs(fitcns) > 1.0:
                    #     msgs.warn("The tilt of Arc Line {0:d} might be poorly traced".format(i+1))
                    tiltsplot[:, i] += fitcns

            xdat = xtilt.copy()
            xdat[np.where(xdat != maskval)] *= (msarc.shape[1] - 1.0)

            msgs.info("Plotting arc tilt QA")
        #    arqa.plot_orderfits(slf, tiltsplot, ztilto, xdata=xdat, xmodl=np.arange(msarc.shape[1]),
        #                        textplt="Arc line", maxp=9, desc="Arc line spectral tilts", maskval=maskval)
            artracewave.plot_orderfits(slf.setup, tiltsplot, ztilto, xdata=xdat, xmodl=np.arange(msarc.shape[1]),
                           textplt="Arc line", maxp=9, desc="Arc line spectral tilts", maskval=maskval, slit=slit)
    # Finish
    return final_tilts, satmask, outpar
'''


def slit_image(scitrace, obj, tilts):
    """ Generate slit image for a given object
    Ignores changing plate scale (for now)
    The slit is approximated as a straight line in this calculation
    which should be reasonably accurate.  Better, the object profile
    generated from this approximation is applied in the same fashion
    so that the 'error' is compensated for.
    Parameters
    ----------
    scitrace
    obj
    Returns
    -------
    slit_img : ndarray
    """
    # Setup
    ximg = np.outer(np.ones(tilts.shape[0]), np.arange(tilts.shape[1]))
    dypix = 1./tilts.shape[0]
    #  Trace
    xtrc = np.round(scitrace['traces'][:, obj]).astype(int)
    wch = np.where((xtrc >= 0) & (xtrc <= tilts.shape[1]-1))
    msgs.work("Use 2D spline to evaluate tilts")
    trc_tilt = np.zeros(tilts.shape[0], dtype=np.float)
    trc_tilt[wch] = tilts[np.arange(tilts.shape[0])[wch], xtrc[wch]]
    trc_tilt_img = np.outer(trc_tilt, np.ones(tilts.shape[1]))
    # Slit image
    msgs.work("Should worry about changing plate scale")
    dy = (tilts - trc_tilt_img)/dypix  # Pixels
    dx = ximg - np.outer(scitrace['traces'][:, obj], np.ones(tilts.shape[1]))
    slit_img = np.sqrt(dx**2 - dy**2)
    neg = dx < 0.
    slit_img[neg] *= -1
    # Return
    return slit_img


def find_obj_minima(trcprof, fwhm=3., nsmooth=3, nfind=8, xedge=0.03,
        sig_thresh=5., peakthresh=None, triml=2, trimr=2, debug=False):
    ''' Find objects using a ported version of nminima from IDL (idlutils)
    
    Parameters
    ----------
    trcprof : ndarray
      1D array of the slit profile including objects
    fwhm : float, optional
      FWHM estimate of the seeing
    nsmooth : int, optional
      Kernel in pixels for Gaussian smoothing of trcprof
    nfind : int, optional
      Number of sources to identify
    xedge : float, optional
      Fraction of the slit on each edge to ignore peaks within
    sig_thresh : float, optional
      Number of sigma to reject when fitting
    peakthresh : float, optional
      Include objects whose peak exceeds this fraction of the highest peak
    triml : int, optional
    trimr : int, optional
    debug : bool, optional

    Returns
    -------
    objl : ndarray  (npeak)
      Left edges for each peak found
    objr : ndarray  (npeak)
      Right edges for each peak found
    bckl : ndarray (npeak, npix)
      Background regions in the slit, left of the peak
    bckr : ndarray (npeak, npix)
      Background regions in the slit, right of the peak
    '''
    npix = trcprof.size
    # Smooth
    if nsmooth > 0:
        yflux = convolve(-1*trcprof, Gaussian1DKernel(nsmooth))
    else:
        yflux = -1*trcprof
    #
    xvec = np.arange(len(yflux))
    # Find peaks
    peaks, sigmas, ledges, redges = utils.find_nminima(yflux, xvec, minsep=fwhm, nfind=nfind, width=int(fwhm))
    fint = interpolate.interp1d(xvec, yflux, bounds_error=False, fill_value=0.)
    ypeaks = -1.*fint(peaks)
    # Sky background (for significance)
    imask = xvec == xvec
    for ipeak in peaks:  # Mask out for sky determination
        ibad = np.abs(xvec-ipeak) < fwhm
        imask[ibad] = False
    clip_yflux = sigma_clip(yflux[imask], axis=0, sigma=2.5)
    sky_std = np.std(clip_yflux)
    # Toss out peaks near edges
    xp = peaks/float(npix)/2.
    gdobj = (xp < (1-xedge)) * (xp > xedge)
    # Keep those above the threshold
    if np.sum(gdobj) > 1:
        threshold = sig_thresh*sky_std
        if peakthresh is not None:
            threshold = max(threshold, peakthresh*max(ypeaks[gdobj]))
        # Parse again
        gdthresh = (ypeaks > threshold)
        if np.any(gdthresh):  # Parse if at least one is bright enough
            gdobj = gdobj & gdthresh
    # Setup for finish
    nobj = np.sum(gdobj)
    objl = ledges[gdobj]  # Might need to set these by threshold
    objr = redges[gdobj]
    bck_mask = np.ones_like(trcprof).astype(int)
    # Mask out objects
    for ledge, redge in zip(objl,objr):
        bck_mask[ledge:redge+1] = 0
    # Mask out edges
    bck_mask[0:triml] = 0
    bck_mask[-trimr:] = 0
    bckl = np.outer(bck_mask, np.ones(nobj))
    bckr = bckl.copy()
    idx = np.arange(npix).astype(int)
    for kk,iobjr in enumerate(objr):
        pos = idx > iobjr
        bckl[pos,kk] = 0
        neg = idx < objl[kk]
        bckr[neg,kk] = 0
    # Return
    return objl, objr, bckl, bckr
