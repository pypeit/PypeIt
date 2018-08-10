from __future__ import (print_function, absolute_import, division, unicode_literals)

import time

import numpy as np
import os

from scipy import signal, ndimage, interpolate

from astropy import units
from astropy.io import fits

from pypit import msgs

from pypit.core import arextract
from pypit.core import arprocimg
from pypit.core import arflat
from pypit.core import arwave
from pypit import artrace
from pypit import arutils
from pypit import arparse as settings
from pypit import arspecobj
from pypit import arpca

from pypit import ardebug as debugger

try:
    basestring
except NameError:
    basestring = str

def background_subtraction(slf, sciframe, varframe, slitn, det, refine=0.0, doqa=True):
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
        if doqa:
            plt_bspline_sky(tilts, sciframe, bgf_flat, gdp)
            debugger.set_trace()
    else:
        msgs.error('Not ready for this method for skysub {:s}'.format(
                settings.argflag['reduce']['skysub']['method'].lower()))
    if np.any(np.isnan(bgframe)):
        msgs.warn("NAN in bgframe.  Replacing with 0")
        bad = np.isnan(bgframe)
        bgframe[bad] = 0.
    return bgframe, nl, nr


def bg_subtraction(slf, tilts, det, sciframe, varframe, bpix, crpix, **kwargs):
    """ Wrapper to run the background subtraction on a series of slits
    Parameters
    ----------
    slf
    det : int
    sciframe : ndarray
    varframe : ndarray
    crpix : ndarray
    kwargs
       Passed to bg_subtraction_slit

    Returns
    -------
    bgframe : ndarray

    """
    # Setup
    bgframe = np.zeros_like(sciframe)
    gdslits = np.where(~slf._maskslits[det-1])[0]

    for slit in gdslits:
        msgs.info("Working on slit: {:d}".format(slit))
        # TODO -- Replace this try/except when a more stable b-spline is used..
        try:
            slit_bgframe = bg_subtraction_slit(slf, det, slit, tilts, sciframe, varframe, bpix, crpix, **kwargs)
        except ValueError:  # Should have been bspline..
            msgs.warn("B-spline sky subtraction failed.  Slit {:d} will no longer be processed..".format(slit))
            #msgs.warn("Continue if you wish..")
            slf._maskslits[det-1][slit] = True
        else:
            bgframe += slit_bgframe
    # Return
    return bgframe



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
        medpix = flxfr[groups == mm]
        if medpix.size == 0:
            profile[mm - 1] = 0.0
        else:
            profile[mm - 1] = np.median(medpix)
    return xedges, profile


def reduce_prepare(slf, sciframe, bpix, datasec_img, scidx, fitsdict, det,
                   mspixelflatnrm=None, standard=False, slitprof=None):
    """ Prepare the Run standard extraction steps on a frame

    Parameters
    ----------
    sciframe : image
      Bias subtracted image (using arload.load_frame)
    bpix : image
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
    dnum = settings.get_dnum(det)
    namp = settings.spect[dnum]['numamplifiers']
    gain_list = settings.spect[dnum]['gain']
    sciframe *= arprocimg.gain_frame(datasec_img, namp, gain_list)
    # Mask
    slf._scimask[det-1] = np.zeros_like(sciframe).astype(int)
    #msgs.info("Masking bad pixels")
    #slf.update_sci_pixmask(det, bpix, 'BadPix')
    # Variance
    msgs.info("Generate raw variance frame (from detected counts [flat fielded])")
    rawvarframe = arprocimg.variance_frame(datasec_img, det, sciframe, scidx,
                                           settings.spect[dnum], fitsdict=fitsdict)

    ###############
    # Subtract off the scattered light from the image
    msgs.work("Scattered light subtraction is not yet implemented...")
    ###############
    # Flat field the science frame (and variance)
    if settings.argflag['reduce']['flatfield']['perform']:
        msgs.info("Flat fielding the science frame")
        # JXP -- I think it is a bad idea to modify the rawvarframe
        #sciframe, rawvarframe = flatfield(slf, sciframe, slf._mspixelflatnrm[det-1], det, varframe=rawvarframe, slitprofile=slf._slitprof[det-1])
        sciframe = arflat.flatfield(sciframe, mspixelflatnrm, bpix, slitprofile=slitprof)
    else:
        msgs.info("Not performing a flat field calibration")
    if not standard:
        slf._sciframe[det-1] = sciframe
        slf._rawvarframe[det-1] = rawvarframe
    ###############
    # Identify cosmic rays
    msgs.work("Include L.A.Cosmic arguments in the settings files")
    if True: crmask = arprocimg.lacosmic(datasec_img, fitsdict, det, sciframe, scidx, settings.spect[dnum], grow=1.5)
    else: crmask = np.zeros(sciframe.shape)
    # Mask
    slf.update_sci_pixmask(det, crmask, 'CR')
    return sciframe, rawvarframe, crmask


def reduce_echelle(slf, sciframe, scidx, fitsdict, det,
                   standard=False, triml=1, trimr=1,
                   mspixelflatnrm=None, doqa=True):
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
    sciframe, rawvarframe, crmask = reduce_prepare(slf, sciframe, scidx, fitsdict, det,
                                                   mspixelflatnrm=mspixelflatnrm,
                                                   standard=standard, slitprof=slitprof)
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
            dnum = settings.get_dnum(det)
            modelvarframe = arprocimg.variance_frame(datasec_img, det, sciframe, scidx,
                settings.spect[dnum], fitsdict=fitsdict, skyframe=bgframe)
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
            if doqa:
#                arqa.pca_plot(slf, outpar, ofit, "Object_Trace", pcadesc="PCA of object trace")
                arpca.pca_plot(slf.setup, outpar, ofit, "Object_Trace", pcadesc="PCA of object trace")
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
    if doqa:
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
        modelvarframe = arprocimg.variance_frame(datasec_img, det, sciframe, scidx,
            settings.spect[dnum], fitsdict=fitsdict, skyframe=bgframe)

    # Perform an optimal extraction
    return reduce_frame(slf, sciframe, rawvarframe, modelvarframe, bgframe, scidx, fitsdict, det, crmask,
                        scitrace=scitrace, standard=standard)


def reduce_multislit(slf, tilts, sciframe, bpix, datasec_img, scidx, fitsdict, det, mswave,
                     mspixelflatnrm=None, standard=False, slitprof=None, debug=False):
    """ Run standard extraction steps on an echelle frame

    Parameters
    ----------
    sciframe : image
      Bias subtracted image (using arload.load_frame)
    bpix : ndarray
      Bad pixel mask
    scidx : int
      Index of the frame
    fitsdict : dict
      Contains relevant information from fits header files
    det : int
      Detector index
    standard : bool, optional
      Standard star frame?
    """
    #
    dnum = settings.get_dnum(det)
    sciframe, rawvarframe, crmask = reduce_prepare(slf, sciframe, bpix, datasec_img,
                                                   scidx, fitsdict, det,
                                                   mspixelflatnrm=mspixelflatnrm,
                                                   slitprof=slitprof)

    # Save sciframe
    slf._sciframe[det-1] = sciframe.copy()

    ###############
    # Estimate Sky Background
    if settings.argflag['reduce']['skysub']['perform']:
        # Perform an iterative background/science extraction
        if debug:
            debugger.set_trace()  # JXP says THIS MAY NOT WORK AS EXPECTED
            msgs.warn("Reading background from 2D image on disk")
            datfil = settings.argflag['run']['directory']['science']+'/spec2d_{:s}.fits'.format(slf._basename.replace(":","_"))
            hdu = fits.open(datfil)
            bgframe = hdu[1].data - hdu[2].data
        else:
            msgs.info("First estimate of the sky background")
            bgframe = bg_subtraction(slf, tilts, det, sciframe, rawvarframe, bpix, crmask)
        modelvarframe = arprocimg.variance_frame(datasec_img, det, sciframe, scidx,
                                       settings.spect[dnum], fitsdict=fitsdict, skyframe=bgframe)
    else:
        modelvarframe = rawvarframe.copy()
        bgframe = np.zeros_like(sciframe)
    if not standard:  # Need to save
        slf._modelvarframe[det - 1] = modelvarframe
        slf._bgframe[det - 1] = bgframe

    ###############
    # Find objects and estimate their traces
    scitrace = artrace.trace_objects_in_slits(slf, det, sciframe-bgframe, modelvarframe, crmask,
                                    bgreg=20, doqa=False, standard=standard)
    if scitrace is None:
        msgs.info("Not performing extraction for science frame"+msgs.newline()+fitsdict['filename'][scidx[0]])
        debugger.set_trace()
        #continue

    # Make sure that there are objects
    noobj = True
    for sl in range(len(scitrace)):
        if 'nobj' in scitrace[sl].keys():  # There can be empty dict's  (skipped slits)
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
            if 'nobj' in scitrace[sl].keys():
                if scitrace[sl]['nobj'] > 0:
                    trcmask += scitrace[sl]['object'].sum(axis=2)
        trcmask[np.where(trcmask > 0.0)] = 1.0
        # Do it
        bgframe = bg_subtraction(slf, tilts, det, sciframe, modelvarframe, bpix, crmask, tracemask=trcmask)
        # Redetermine the variance frame based on the new sky model
        modelvarframe = arprocimg.variance_frame(datasec_img, det, sciframe, scidx,
                                       settings.spect[dnum], fitsdict=fitsdict, skyframe=bgframe)
        # Save
        if not standard:
            slf._modelvarframe[det-1] = modelvarframe
            slf._bgframe[det-1] = bgframe

    ###############
    # Flexure down the slit? -- Not currently recommended
    if settings.argflag['reduce']['flexure']['method'] == 'slitcen':
        flex_dict = arwave.flexure_slit(slf, det)
        arwave.flexure_qa(slf, det, flex_dict, slit_cen=True)

    # Perform an optimal extraction
    msgs.work("For now, perform extraction -- really should do this after the flexure+heliocentric correction")
    return reduce_frame(slf, sciframe, rawvarframe, modelvarframe, bpix, datasec_img, bgframe,
                        scidx, fitsdict, det, crmask, tilts, mswave, standard=standard)


def reduce_frame(slf, sciframe, rawvarframe, modelvarframe, bpix, datasec_img,
                 bgframe, scidx, fitsdict, det, crmask, tilts, mswave,
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
    dnum = settings.get_dnum(det)
    ###############
    # Determine the final trace of the science objects
    if scitrace is None:
        msgs.info("Performing final object trace")
        scitrace = artrace.trace_objects_in_slits(slf, det, sciframe-bgframe, modelvarframe, crmask,
                                        bgreg=20, doqa=(not standard), standard=standard)
    if standard:
    #    slf._msstd[det-1]['trace'] = scitrace
    #    specobjs = arspecobj.init_exp(slf, scidx, det, fitsdict, scitrace, objtype='standard')
    #    slf._msstd[det-1]['spobjs'] = specobjs
        specobjs = arspecobj.init_exp(slf, scidx, det, fitsdict, scitrace, objtype='standard')
    else:
        # Generate SpecObjExp list
        specobjs = arspecobj.init_exp(slf, scidx, det, fitsdict, scitrace, objtype='science')

    slf._scitrace[det-1] = scitrace
    slf._specobjs[det-1] = specobjs

    ###############
    # Extract
    noobj = True
    for sl in range(len(scitrace)):
        if 'nobj' in scitrace[sl].keys():
            if scitrace[sl]['nobj'] != 0:
                noobj = False
    if noobj is True:
        msgs.warn("No objects to extract for science frame"+msgs.newline()+fitsdict['filename'][scidx])
        return True

    # Boxcar
    msgs.info("Performing boxcar extraction")
    bgcorr_box = arextract.boxcar(slf, det, specobjs, sciframe-bgframe, rawvarframe, bpix,
                                  bgframe, crmask, scitrace, mswave)

    # Optimal
    if not standard:

        # KBW: Using variance_frame() in arextract leads to a circular
        # import.  I've changed the arextract.optimal_extract() function
        # to return the object model, then the last step of generating
        # the new variance image is done here.

        msgs.info("Attempting optimal extraction with model profile")
        arextract.obj_profiles(slf, det, specobjs, sciframe-bgframe-bgcorr_box,
                               modelvarframe, bgframe+bgcorr_box, crmask, scitrace, tilts, doqa=False)
#        newvar = arextract.optimal_extract(slf, det, specobjs, sciframe-bgframe-bgcorr_box,
#                                           modelvarframe, bgframe+bgcorr_box, crmask, scitrace)
        obj_model = arextract.optimal_extract(slf, det, slf._specobjs[det-1], sciframe-bgframe-bgcorr_box,
                                              modelvarframe, bgframe+bgcorr_box, crmask, scitrace, tilts,
                                              mswave)
        newvar = arprocimg.variance_frame(datasec_img, det, sciframe-bgframe-bgcorr_box, -1,
                                settings.spect[dnum], skyframe=bgframe+bgcorr_box, objframe=obj_model)
        msgs.work("Should update variance image (and trace?) and repeat")
        #
        arextract.obj_profiles(slf, det, slf._specobjs[det-1], sciframe-bgframe-bgcorr_box,
                               newvar, bgframe+bgcorr_box, crmask, scitrace, tilts)
#        finalvar = arextract.optimal_extract(slf, det, specobjs, sciframe-bgframe-bgcorr_box,
#                                             newvar, bgframe+bgcorr_box, crmask, scitrace)
        obj_model = arextract.optimal_extract(slf, det, specobjs, sciframe-bgframe-bgcorr_box,
                                              newvar, bgframe+bgcorr_box, crmask, scitrace, tilts,
                                              mswave)
        finalvar = arprocimg.variance_frame(datasec_img, det, sciframe-bgframe-bgcorr_box, -1,
                                  settings.spect[dnum], skyframe=bgframe+bgcorr_box, objframe=obj_model)
        slf._modelvarframe[det-1] = finalvar.copy()

    # Flexure correction?
    if settings.argflag['reduce']['flexure']['perform'] and (not standard):
        if settings.argflag['reduce']['flexure']['method'] is not None:
            flex_list = arwave.flexure_obj(slf, det)
            arwave.flexure_qa(slf, det, flex_list)

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

