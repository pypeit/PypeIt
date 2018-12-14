""" Module for sky subtraction
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
import sys, os
from matplotlib import pyplot as plt

#from pydl.pydlutils.bspline import bspline
from pypeit import ginga
from pypeit.core import pydl
from pypeit import msgs
from pypeit import utils
from pypeit import debugger
from pypeit.core import pixels
from pypeit.core import extract
from matplotlib import pyplot as plt

from scipy.special import ndtr



def global_skysub(image, ivar, tilts, thismask, slit_left, slit_righ, inmask = None, bsp=0.6, sigrej=3., trim_edg = (3,3),
                  pos_mask=True, show_fit=False, no_poly=False, npoly = None):
    """
    Perform global sky subtraction on an input slit

    Parameters
    ----------
    image: float ndarray, shape (nspec, nspat)
          Frame to be sky subtracted

    ivar: float ndarray, shape (nspec, nspat)
          Inverse variance image

    tilts: float ndarray, shape (nspec, nspat)
          Tilgs indicating how wavelengths move across the slit

    thismask : numpy boolean array, shape (nspec, nspat)
      Specifies pixels in the slit in question

    slit_left: ndarray of shape (nspec, 1) or (nspec)
      Left slit boundary in floating point pixels.

    slit_righ: ndarray of shape (nspec, 1) or (nspec)
      Right slit boundary in floating point pixels.


    Optional Parameters
    --------------------

    inmask: boolean ndarray, shape (nspec, nspat), default inmask = None
      Input mask for pixels not to be included in sky subtraction fits. True = Good (not masked), False = Bad (masked)

    bsp: float, default bsp = 0.6
      break point spacing in pixel units

    sigrej : float, default sigrej = 3.0
      sigma rejection threshold

    trim_edg: tuple of floats  (left_edge, right_edge), default (3,3)
      indicates how many pixels to trim from left and right slit edges for creating the edgemask. These pixels are
      excluded from sky subtraction fits.

    pos_mask: boolean, defualt pos_mask = True
      First do a prelimnary fit to the log of the sky (i.e. positive pixels only). Then use this fit to create an input
      mask from the residuals lmask = (res < 5.0) & (res > -4.0) for the full fit.
      NOTE: pos_mask should be False for near-IR sky residual subtraction, since fitting the log(sky) requires that the
      counts are positive which will not be the case for i.e. an A-B image. Thus the routine will fail if pos_mask is not
      set to False.

    show_fit: boolean, default show_fit = False
       Plot a fit of the sky pixels and model fit to the screen. This feature will block further execution until the screen is closed.

    Returns
    -------
    bgframe : ndarray
      Returns the model sky background at the pixels where thismask is True.

     >>>  skyframe = np.zeros_like(image)
     >>>  thismask = slitpix == thisslit
     >>>  skyframe[thismask] = global_skysub(image,ivar, tilts, thismask, slit_left, slit_righ)

    """

    # Synthesize ximg, and edgmask  from slit boundaries. Doing this outside this
    # routine would save time. But this is pretty fast, so we just do it here to make the interface simpler.

    # TESTING!!!!
    #no_poly=True
    #show_fit=True

    ximg, edgmask = pixels.ximg_and_edgemask(slit_left, slit_righ, thismask, trim_edg=trim_edg)


    # Init
    (nspec, nspat) = image.shape
    piximg = tilts * (nspec-1)
    if inmask is None:
        inmask = np.copy(thismask)


    # Sky pixels for fitting
    inmask_in = (thismask == True) & (ivar > 0.0) & (inmask == True) & (edgmask == False)
    isrt = np.argsort(piximg[thismask])
    pix = piximg[thismask][isrt]
    sky = image[thismask][isrt]
    sky_ivar = ivar[thismask][isrt]
    ximg_fit = ximg[thismask][isrt]
    inmask_fit = inmask_in[thismask][isrt]
    #spatial = spatial_img[fit_sky][isrt]

    # Restrict fit to positive pixels only and mask out large outliers via a pre-fit to the log.
    if (pos_mask is True):
        pos_sky = (sky > 1.0) & (sky_ivar > 0.0)
        if np.sum(pos_sky) > nspec:
            lsky = np.log(sky[pos_sky])
            lsky_ivar = inmask_fit[pos_sky].astype(float)/3.0** 2  # set errors to just be 3.0 in the log
            #lsky_ivar = np.full(lsky.shape, 0.1)
            # Init bspline to get the sky breakpoints (kludgy)
            #tmp = pydl.bspline(wsky[pos_sky], nord=4, bkspace=bsp)
            lskyset, outmask, lsky_fit, red_chi = utils.bspline_profile(pix[pos_sky], lsky, lsky_ivar, np.ones_like(lsky),
                                                                        inmask = inmask_fit[pos_sky],
                                                                        upper=sigrej, lower=sigrej,
                                                                        kwargs_bspline={'bkspace':bsp},
                                                                        kwargs_reject={'groupbadpix': True, 'maxrej': 10})
            res = (sky[pos_sky] - np.exp(lsky_fit)) * np.sqrt(sky_ivar[pos_sky])
            lmask = (res < 5.0) & (res > -4.0)
            sky_ivar[pos_sky] = sky_ivar[pos_sky] * lmask
            inmask_fit[pos_sky]=(sky_ivar[pos_sky] > 0.0) & lmask

    if no_poly:
        poly_basis= np.ones_like(sky)
    else:
        npercol = np.fmax(np.floor(np.sum(thismask) / nspec), 1.0)
        # Demand at least 10 pixels per row (on average) per degree of the polynomial
        if npoly is None:
            #npoly_in = 7
            #npoly = np.fmax(np.fmin(npoly_in, (np.ceil(npercol / 10.)).astype(int)), 1)
            if npercol > 100:
                npoly = 3
            elif npercol > 40:
                npoly = 2
            else:
                npoly = 1
        poly_basis = pydl.flegendre(2.0*ximg_fit - 1.0, npoly).T

    # Full fit now
    #full_bspline = pydl.bspline(wsky, nord=4, bkspace=bsp, npoly = npoly)
    #skyset, outmask, yfit, _ = utils.bspline_profile(wsky, sky, sky_ivar, poly_basis,
    #                                                   fullbkpt=full_bspline.breakpoints,upper=sigrej, lower=sigrej,
    #                                                   kwargs_reject={'groupbadpix':True, 'maxrej': 10})


    # Perform the full fit now
    skyset, outmask, yfit, _ = utils.bspline_profile(pix, sky, sky_ivar,poly_basis,inmask = inmask_fit, nord = 4,
                                                               upper=sigrej, lower=sigrej,
                                                               kwargs_bspline = {'bkspace':bsp},
                                                               kwargs_reject={'groupbadpix':True, 'maxrej': 10})

    sky_frame = np.zeros_like(image)
    ythis = np.zeros_like(yfit)
    ythis[isrt] = yfit
    sky_frame[thismask] = ythis

    #skyset.funcname ='legendre'
    #skyset.xmin = spat_min
    #skyset.xmax = spat_max

    # Evaluate and save
    #bgframe, _ = skyset.value(piximg[thismask],x2=spatial_img[thismask])

    # Debugging/checking

    # ToDo This QA ceases to make sense I think for 2-d fits. I need to think about what the best QA would be here, but I think
    # probably looking at residuals as a function of spectral and spatial position like in the flat fielding code.
    if show_fit:
        goodbk = skyset.mask
        # This is approximate
        yfit_bkpt = np.interp(skyset.breakpoints[goodbk], pix,yfit)
        plt.clf()
        ax = plt.gca()
        was_fit_and_masked = inmask_fit & ~outmask
        ax.plot(pix, sky, color='k', marker='o', markersize=0.4, mfc='k', fillstyle='full', linestyle='None')
        ax.plot(pix[was_fit_and_masked], sky[was_fit_and_masked], color='red', marker='+', markersize=1.5, mfc='red', fillstyle='full', linestyle='None')
        ax.plot(pix, yfit, color='cornflowerblue')
        ax.plot(skyset.breakpoints[goodbk], yfit_bkpt, color='lawngreen', marker='o', markersize=4.0, mfc='lawngreen', fillstyle='full', linestyle='None')
        ax.set_ylim((0.99*yfit.min(),1.01*yfit.max()))
        plt.show()

    # Return
    # ToDO worth thinking about whether we want to return a mask here. It makese no sense to return outmask
    # in its present form though since that does not refer to the whole image.
    # return bgframe, outmask
    return ythis



# Utility routine used by local_bg_subtraction_slit
def skyoptimal(wave,data,ivar, oprof, sortpix, sigrej = 3.0, npoly = 1, spatial = None, fullbkpt = None):


    nx = data.size
    nc = oprof.shape[0]
    nobj = int(oprof.size / nc)
    if nc != nx:
        raise ValueError('Object profile should have oprof.shape[0] equal to nx')

    msgs.info('Iter     Chi^2     Rejected Pts')
    xmin = 0.0
    xmax = 1.0

    if ((npoly == 1) | (spatial is None)):
        profile_basis = np.column_stack((oprof, np.ones(nx)))
    else:
        xmin = spatial.min()
        xmax = spatial.max()
        x2 = 2.0 * (spatial - xmin) / (xmax - xmin) - 1
        poly_basis = pydl.flegendre(x2, npoly).T
        profile_basis = np.column_stack((oprof, poly_basis))

    relative_mask = (np.sum(oprof, axis=1) > 1e-10)

    indx, = np.where(ivar[sortpix] > 0.0)
    ngood = indx.size
    good = sortpix[indx]
    good = good[wave[good].argsort()]
    relative, = np.where(relative_mask[good])



    sset1, outmask_good1, yfit1, red_chi1 = utils.bspline_profile(wave[good], data[good], ivar[good], profile_basis[good, :],
                                                              fullbkpt=fullbkpt, upper=sigrej, lower=sigrej,
                                                              relative=relative,
                                                              kwargs_reject={'groupbadpix': True, 'maxrej': 5})

    chi2 = (data[good] - yfit1) ** 2 * ivar[good]
    chi2_srt = np.sort(chi2)
    gauss_prob = 1.0 - 2.0 * ndtr(-1.2 * sigrej)
    sigind = int(np.fmin(np.rint(gauss_prob * float(ngood)), ngood - 1))
    chi2_sigrej = chi2_srt[sigind]
    mask1 = (chi2 < chi2_sigrej)
    msgs.info('2nd round....')
    msgs.info('Iter     Chi^2     Rejected Pts')

    sset, outmask_good, yfit, red_chi = utils.bspline_profile(wave[good], data[good], ivar[good] * mask1,
                                                          profile_basis[good, :],
                                                          fullbkpt=fullbkpt, upper=sigrej, lower=sigrej,
                                                          relative=relative,
                                                          kwargs_reject={'groupbadpix': True, 'maxrej': 1})

    ncoeff = npoly + nobj
    skyset = pydl.bspline(None, fullbkpt=sset.breakpoints, nord=sset.nord, npoly=npoly)
    # Set coefficients for the sky.
    # The rehshape below deals with the different sizes of the coeff for npoly = 1 vs npoly > 1
    # and mirrors similar logic in the bspline.py
    skyset.coeff = sset.coeff[nobj:, :].reshape(skyset.coeff.shape)

    skyset.mask = sset.mask
    skyset.xmin = xmin
    skyset.xmax = xmax

    sky_bmodel, _ = skyset.value(wave, x2=spatial)

    obj_bmodel = np.zeros(sky_bmodel.shape)
    objset = pydl.bspline(None, fullbkpt=sset.breakpoints, nord=sset.nord)
    objset.mask = sset.mask
    for i in range(nobj):
        objset.coeff = sset.coeff[i, :]
        obj_bmodel1, _ = objset.value(wave)
        obj_bmodel = obj_bmodel + obj_bmodel1 * profile_basis[:, i]

    outmask = np.zeros(wave.shape, dtype=bool)
    outmask[good] = outmask_good

    return (sky_bmodel, obj_bmodel, outmask)

def local_skysub_extract(sciimg, sciivar, tilts, waveimg, global_sky, rn2_img, thismask, slit_left, slit_righ, sobjs,
                         bsp = 0.6, inmask = None, extract_maskwidth = 3.0, trim_edg = (3,3), std = False, prof_nsigma = None, niter=4,
                         box_rad = 7, sigrej = 3.5,skysample = False, sn_gauss = 3.0, coadd_2d = False, show_profile=False,
                         show_resids=False):

    """Perform local sky subtraction and  extraction

     Parameters
     ----------
     sciimg : numpy float 2-d array (nspec, nspat)
         sky-subtracted image
     sciivar : numpy float 2-d array (nspec, nspat)
         inverse variance of sky-subtracted image
     tilts: ndarray, (nspec, nspat)
         spectral tilts
     waveimg numpy float 2-d array (nspec, nspat)
         2-d wavelength map
     global_sky : ndarray (nspec, nspat)
         Global sky model
     rn2_img:
         Image with the read noise squared per pixel
         object trace


    Optional Parameters
    ----------
    extract_maskwidth: float, default = 3.0
        This parameter determines the initial size of the region in units of fwhm that will be used for local sky subtraction. This
        maskwidth is defined in the obfjind code, but is then updated here as the profile fitting improves the fwhm estimates

     Returns
     -------
     :func:`tuple`

     """


    ximg, edgmask = pixels.ximg_and_edgemask(slit_left, slit_righ, thismask, trim_edg = trim_edg)

    nspat = sciimg.shape[1]
    nspec = sciimg.shape[0]
    piximg = tilts * (nspec-1)

    # Copy the specobjs that will be the output
    nobj = len(sobjs)
    # specobjs = copy.deepcopy(specobjs_in)

    if (prof_nsigma is None):
        prof_nsigma1 = np.full(len(sobjs), None)
    elif len(prof_nsigma) == 1:
        prof_nsigma1 = np.full(nobj, prof_nsigma)
    elif len(prof_nsigma) == nobj:
        prof_nsigma1 = prof_nsigma
    else:
        raise ValueError('Invalid size for prof_nsigma.')

    # Set some rejection parameters based on whether this is a standard or not. Only reject extreme outliers for standards
    # since super high S/N and low order profile models imply we will always have large outliers
    if std is True:
        chi2_sigrej = 100.0
        sigrej_ceil = 1e10
    else:
        chi2_sigrej = 6.0
        sigrej_ceil = 10.0
    # We will use this number later
    gauss_prob = 1.0 - 2.0 * ndtr(-sigrej)

    for iobj in range(nobj):
        sobjs[iobj].prof_nsigma = prof_nsigma1[iobj]

    nspat = sciimg.shape[1]
    nspec = sciimg.shape[0]

    if inmask is None:
        # These values are hard wired for the case where no inmask is provided
        FULLWELL = 5e5
        MINWELL = -1000.0,
        inmask = (sciivar > 0.0) & thismask & np.isfinite(sciimg) & np.isfinite(sciivar) & (sciimg < FULLWELL) & (sciimg > MINWELL)

    # Create the images that will be returned
    outmask = np.copy(inmask)
    modelivar = np.copy(sciivar)
    objimage = np.zeros_like(sciimg)
    skyimage = np.copy(global_sky)
    varnoobj = np.abs(skyimage - np.sqrt(2.0) * np.sqrt(rn2_img)) + rn2_img

    xarr = np.outer(np.ones(nspec), np.arange(nspat))
    yarr = np.outer(np.arange(nspec), np.ones(nspat))

    xa_min = xarr[thismask].min()
    xa_max = xarr[thismask].max()
    ya_min = yarr[thismask].min()
    ya_max = yarr[thismask].max()

    xsize = slit_righ - slit_left
    spatial_img = thismask * ximg * (np.outer(xsize, np.ones(nspat)))

    # Loop over objects and group them
    i1 = 0
    while i1 < nobj:
        group = []
        group.append(i1)
        # The default value of maskwidth = 3.0 * FWHM = 7.05 * sigma in objfind with a log(S/N) correction for bright objects
        mincols = np.maximum(sobjs[i1].trace_spat - sobjs[i1].maskwidth - 1, slit_left)
        maxcols = np.minimum(sobjs[i1].trace_spat + sobjs[i1].maskwidth + 1, slit_righ)
        for i2 in range(i1 + 1, nobj):
            left_edge = sobjs[i2].trace_spat - sobjs[i2].maskwidth - 1
            righ_edge = sobjs[i2].trace_spat + sobjs[i2].maskwidth + 1
            touch = (left_edge < maxcols) & (sobjs[i2].trace_spat > slit_left) & (righ_edge > mincols)
            if touch.any():
                maxcols = np.minimum(np.maximum(righ_edge, maxcols), slit_righ)
                mincols = np.maximum(np.minimum(left_edge, mincols), slit_left)
                group.append(i2)
        # Keep for next iteration
        i1 = max(group) + 1
        # Some bookeeping to define the sub-image and make sure it does not land off the mask
        objwork = len(group)
        scope = np.sum(thismask, axis=0)
        iscp, = np.where(scope)
        imin = min(iscp)
        imax = max(iscp)
        mincol = np.fmax(np.floor(min(mincols)), imin)
        maxcol = np.fmin(np.ceil(max(maxcols)), imax)
        nc = int(maxcol - mincol + 1)
        rows = np.arange(nspec, dtype=np.intp)
        columns = np.arange(mincol, mincol + nc, dtype=np.intp)
        ipix = np.ix_(rows, columns)
        skymask = outmask & ~edgmask
        if nc > 100:
            npoly = 3
        elif nc > 40:
            npoly = 2
        else:
            npoly = 1
        obj_profiles = np.zeros((nspec, nspat, objwork), dtype=float)
        sigrej_eff = sigrej
#        proc_list = [] # List of processes for interactive plotting, these are terminated at the end of each iteration sequence
        for iiter in range(1, niter + 1):
            msgs.info('--------------------------REDUCING: Iteration # ' + '{:2d}'.format(iiter) + ' of ' +
                      '{:2d}'.format(niter) + '---------------------------------------------------')
            img_minsky = sciimg - skyimage
            for ii in range(objwork):
                iobj = group[ii]
                if iiter == 1:
                    # If this is the first iteration, print status message. Initiate profile fitting with a simple
                    # boxcar extraction.
                    msgs.info("----------------------------------- PROFILE FITTING --------------------------------------------------------")
                    msgs.info("Fitting profile for obj # " + "{:}".format(sobjs[iobj].objid) + " of {:}".format(nobj))
                    msgs.info("At x = {:5.2f}".format(sobjs[iobj].spat_pixpos) + " on slit # {:}".format(sobjs[iobj].slitid))
                    msgs.info("------------------------------------------------------------------------------------------------------------")
                    flux = extract.extract_boxcar(img_minsky * outmask, sobjs[iobj].trace_spat, box_rad,
                                          ycen=sobjs[iobj].trace_spec)
                    mvarimg = 1.0 / (modelivar + (modelivar == 0))
                    mvar_box = extract.extract_boxcar(mvarimg * outmask, sobjs[iobj].trace_spat, box_rad,
                                              ycen=sobjs[iobj].trace_spec)
                    pixtot = extract.extract_boxcar(0 * mvarimg + 1.0, sobjs[iobj].trace_spat, box_rad,
                                            ycen=sobjs[iobj].trace_spec)
                    mask_box = (extract.extract_boxcar(~outmask, sobjs[iobj].trace_spat, box_rad,
                                               ycen=sobjs[iobj].trace_spec) != pixtot)
                    box_denom = extract.extract_boxcar(waveimg > 0.0, sobjs[iobj].trace_spat, box_rad,
                                               ycen=sobjs[iobj].trace_spec)
                    wave = extract.extract_boxcar(waveimg, sobjs[iobj].trace_spat, box_rad, ycen=sobjs[iobj].trace_spec) / (
                                box_denom + (box_denom == 0.0))
                    fluxivar = mask_box / (mvar_box + (mvar_box == 0.0))
                else:
                    # For later iterations, profile fitting is based on an optimal extraction
                    last_profile = obj_profiles[:, :, ii]
                    trace = np.outer(sobjs[iobj].trace_spat, np.ones(nspat))
                    objmask = ((xarr >= (trace - 2.0 * box_rad)) & (xarr <= (trace + 2.0 * box_rad)))
                    extract.extract_optimal(sciimg, modelivar, (outmask & objmask), waveimg, skyimage, rn2_img, last_profile,
                                    box_rad, sobjs[iobj])
                    # If the extraction is bad do not update
                    if sobjs[iobj].optimal['MASK'].any():
                        flux = sobjs[iobj].optimal['COUNTS']
                        fluxivar = sobjs[iobj].optimal['COUNTS_IVAR']
                        wave = sobjs[iobj].optimal['WAVE']

                obj_string = 'obj # {:}'.format(sobjs[iobj].objid) + ' on slit # {:}'.format(sobjs[iobj].slitid) + ', iter # {:}'.format(iiter) + ':'
                if wave.any():
                    (profile_model, xnew, fwhmfit, med_sn2) = extract.fit_profile(img_minsky[ipix], (modelivar * outmask)[ipix],
                                                                                  waveimg[ipix],
                                                                                  sobjs[iobj].trace_spat - mincol,
                                                                                  wave, flux, fluxivar,
                                                                                  thisfwhm=sobjs[iobj].fwhm,
                                                                                  maskwidth=sobjs[iobj].maskwidth,
                                                                                  prof_nsigma=sobjs[iobj].prof_nsigma,
                                                                                  sn_gauss=sn_gauss, obj_string = obj_string,
                                                                                  show_profile=show_profile)
                    #proc_list.append(show_proc)

                    # Update the object profile and the fwhm and mask parameters
                    obj_profiles[ipix[0], ipix[1], ii] = profile_model
                    sobjs[iobj].trace_spat = xnew + mincol
                    sobjs[iobj].fwhmfit = fwhmfit
                    sobjs[iobj].fwhm = np.median(fwhmfit)
                    mask_fact = 1.0 + 0.5 * np.log10(np.fmax(np.sqrt(np.fmax(med_sn2, 0.0)), 1.0))
                    maskwidth = extract_maskwidth*np.median(fwhmfit) * mask_fact
                    if sobjs[iobj].prof_nsigma is None:
                        sobjs[iobj].maskwidth = maskwidth
                    else:
                        sobjs[iobj].maskwidth = sobjs[iobj].prof_nsigma * (sobjs[iobj].fwhm / 2.3548)

                else:
                    msgs.warn("Bad extracted wavelengths in local_skysub_extract")
                    msgs.warn("Skipping this profile fit and continuing.....")

            sky_bmodel = np.array(0.0)
            iterbsp = 0
            while (sky_bmodel.any() == False) & (iterbsp <= 5):
                bsp_now = (1.2 ** iterbsp) * bsp
                # if skysample is set, determine optimal break-point spacing
                # directly measuring how well we are sampling of the sky. The
                # bsp in this case correspons to the minimum distance between
                # breakpoints which we allow.
                if skysample:
                    sampmask = (waveimg > 0.0) & (thismask == True)
                    # fullbkpt = skybkpts()
                    # TODO Port long_skybkpts.pro code and put it here.
                else:
                    pixvec = piximg[skymask]
                    srt = pixvec.flatten().argsort()
                    bset0 = pydl.bspline(pixvec.flat[srt], nord=4, bkspace=bsp_now)
                    fullbkpt = bset0.breakpoints
                # check to see if only a subset of the image is used.
                # if so truncate input pixels since this can result in singular matrices
                ibool = (yarr >= ya_min) & (yarr <= ya_max) & (xarr >= xa_min) & (xarr <= xa_max) & (xarr >= mincol) & (
                            xarr <= maxcol) & thismask
                isub, = np.where(ibool.flatten())
                sortpix = (piximg.flat[isub]).argsort()
                ithis, = np.where(thismask.flat[isub])
                keep = (fullbkpt >= piximg.flat[isub[ithis]].min()) & (fullbkpt <= piximg.flat[isub[ithis]].max())
                fullbkpt = fullbkpt[keep]
                obj_profiles_flat = obj_profiles.reshape(nspec * nspat, objwork)
                (sky_bmodel, obj_bmodel, outmask_opt) = skyoptimal(piximg.flat[isub], sciimg.flat[isub],
                                                                   (modelivar * skymask).flat[isub],
                                                                   obj_profiles_flat[isub, :], sortpix,
                                                                   spatial=spatial_img.flat[isub],
                                                                   fullbkpt=fullbkpt, sigrej=sigrej_eff, npoly=npoly)
                iterbsp = iterbsp + 1
                if (sky_bmodel.any() is False) & (iterbsp <= 4):
                    msgs.warn('***************************************')
                    msgs.warn('WARNING: bspline sky-subtraction failed')
                    msgs.warn('Increasing bkpt spacing by 20%. Retry')
                    msgs.warn(
                        'Old bsp = {:5.2f}'.format(bsp_now) + '; New bsp = {:5.2f}'.format(1.2 ** (iterbsp) * bsp))
                    msgs.warn('***************************************')

            if (sky_bmodel.any() == True):
                skyimage.flat[isub] = sky_bmodel
                objimage.flat[isub] = obj_bmodel
                img_minsky.flat[isub] = sciimg.flat[isub] - sky_bmodel
                var = np.abs(sky_bmodel + obj_bmodel - np.sqrt(2.0) * np.sqrt(rn2_img.flat[isub])) + rn2_img.flat[isub]
                var_no = np.abs(sky_bmodel - np.sqrt(2.0) * np.sqrt(rn2_img.flat[isub])) + rn2_img.flat[isub]
                igood1 = skymask.flat[isub]
                #  update the outmask for only those pixels that were fit. This prevents masking of slit edges in outmask
                outmask.flat[isub[igood1]] = outmask_opt[igood1]
                #  For weighted co-adds, the variance of the image is no longer equal to the image, and so the modelivar
                #  eqn. below is not valid. However, co-adds already have the model noise propagated correctly in sciivar,
                #  so no need to re-model the variance
                if coadd_2d is False:
                    modelivar.flat[isub] = (var > 0.0) / (var + (var == 0.0))
                    varnoobj.flat[isub] = var_no
                # Now do some masking based on this round of model fits
                chi2 = (img_minsky.flat[isub] - obj_bmodel) ** 2 * modelivar.flat[isub]
                igood = (skymask.flat[isub]) & (chi2 <= chi2_sigrej ** 2)
                ngd = np.sum(igood)
                if ngd > 0:
                    chi2_good = chi2[igood]
                    chi2_srt = np.sort(chi2_good)
                    sigind = np.fmin(int(np.rint(gauss_prob * float(ngd))), ngd - 1)
                    chi2_sigrej = chi2_srt[sigind]
                    sigrej_eff = np.fmax(np.sqrt(chi2_sigrej), sigrej)
                    #  Maximum sigrej is sigrej_ceil (unless this is a standard)
                    sigrej_eff = np.fmin(sigrej_eff, sigrej_ceil)
                    msgs.info('Measured effective rejection from distribution of chi^2')
                    msgs.info('Instead of rejecting sigrej = {:5.2f}'.format(sigrej) +
                              ', use threshold sigrej_eff = {:5.2f}'.format(sigrej_eff))
                    # Explicitly mask > sigrej outliers using the distribution of chi2 but only in the region that was actually fit.
                    # This prevents e.g. excessive masking of slit edges
                    outmask.flat[isub[igood1]] = outmask.flat[isub[igood1]] & (chi2[igood1] < chi2_sigrej) & (
                                sciivar.flat[isub[igood1]] > 0.0)
                    nrej = outmask.flat[isub[igood1]].sum()
                    msgs.info(
                        'Iteration = {:d}'.format(iiter) + ', rejected {:d}'.format(nrej) + ' of ' + '{:d}'.format(
                            igood1.sum()) + 'fit pixels')

            else:
                msgs.warn('ERROR: Bspline sky subtraction failed after 4 iterations of bkpt spacing')
                msgs.warn('       Moving on......')
                obj_profiles = np.zeros_like(obj_profiles)
                # Just replace with the global sky
                skyimage.flat[isub] = global_sky.flat[isub]

        # Now that iterations are complete, clear the windows if show_profile was set
#        if show_profile:
#            for proc in proc_list:
#                proc.terminate()
#                proc.join()

        # Now that the iterations of profile fitting and sky subtraction are completed,
        # loop over the objwork objects in this grouping and perform the final extractions.
        for ii in range(objwork):
            iobj = group[ii]
            msgs.info('Extracting for obj # {:d}'.format(iobj + 1) + ' of {:d}'.format(nobj) +
                      ' on slit # {:d}'.format(sobjs[iobj].slitid) + ' at x = {:5.2f}'.format(sobjs[iobj].spat_pixpos))
            this_profile = obj_profiles[:, :, ii]
            trace = np.outer(sobjs[iobj].trace_spat, np.ones(nspat))
            objmask = ((xarr >= (trace - 2.0 * box_rad)) & (xarr <= (trace + 2.0 * box_rad)))
            extract.extract_optimal(sciimg, modelivar * thismask, (outmask & objmask), waveimg, skyimage, rn2_img, this_profile,
                            box_rad, sobjs[iobj])
            sobjs[iobj].mincol = mincol
            sobjs[iobj].maxcol = maxcol



    # If requested display the model fits for this slit
    if show_resids:
        viewer, ch = ginga.show_image((sciimg - skyimage - objimage) * np.sqrt(modelivar) * thismask)
        # TODO add error checking here to see if ginga exists
        canvas = viewer.canvas(ch._chname)
        out1 = canvas.clear()
        out2 = ch.cut_levels(-5.0, 5.0)
        out3 = ch.set_color_algorithm('linear')
        # Overplot the traces
        for spec in sobjs:
            if spec.hand_extract_flag is False:
                color = 'magenta'
            else:
                color = 'orange'
            ginga.show_trace(viewer, ch, spec.trace_spat, spec.idx, color=color)

        # These are the pixels that were masked by the extraction
        spec_mask, spat_mask = np.where((outmask == False) & (inmask == True))
        nmask = len(spec_mask)
        # note: must cast numpy floats to regular python floats to pass the remote interface
        points_mask = [dict(type='point', args=(float(spat_mask[i]), float(spec_mask[i]), 2),
                            kwargs=dict(style='plus', color='red')) for i in range(nmask)]

        # These are the pixels that were originally masked
        spec_omask, spat_omask = np.where((inmask == False) & (thismask == True))
        nomask = len(spec_omask)
        # note: must cast numpy floats to regular python floats to pass the remote interface
        points_omask = [dict(type='point', args=(float(spat_omask[i]), float(spec_omask[i]), 2),
                             kwargs=dict(style='plus', color='cyan')) for i in range(nomask)]

        # Labels for the points
        text_mask = [dict(type='text', args=(nspat / 2, nspec / 2, 'masked by extraction'),
                          kwargs=dict(color='red', fontsize=20))]
        text_omask = [dict(type='text', args=(nspat / 2, nspec / 2 + 30, 'masked initially'),
                           kwargs=dict(color='cyan', fontsize=20))]

        canvas_list = points_mask + points_omask + text_mask + text_omask
        canvas.add('constructedcanvas', canvas_list)

    # Clean up any profile plots
#    if show_profile:
#        try:
#            show_proc.terminate()
#            show_proc.join()
#        except:
#            pass

    return (skyimage[thismask], objimage[thismask], modelivar[thismask], outmask[thismask])


