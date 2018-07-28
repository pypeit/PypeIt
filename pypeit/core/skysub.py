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

from scipy.special import ndtr



# ToDO Fix masking logic. This code should also take an ivar for consistency with rest of extraction
def bg_subtraction_slit(slit, slitpix, edge_mask, sciframe, varframe, tilts,
                        bpm=None, crmask=None, tracemask=None, bsp=0.6, sigrej=3., POS_MASK=True,
                        PLOT_FIT=False):
    """
    Perform sky subtraction on an input slit

    Parameters
    ----------
    slit : int
      Slit number; indexed 1, 2,
    slitpix : ndarray
      Specifies pixels in the slits
    edgemask : ndarray
      Mask edges of the slit
    sciframe : ndarray
      science frame
    varframe : ndarray
      Variance array
    tilts : ndarray
      Tilts of the wavelengths
    bpm : ndarray, optional
      Bad pixel mask
    crmask : ndarray
      Cosmic ray mask
    tracemask : ndarray
      Object mask
    bsp : float
      Break point spacing
    sigrej : float
      rejection

    Returns
    -------
    bgframe : ndarray
      Sky background image

    """

    # Init
    bgframe = np.zeros_like(sciframe)
    ivar = utils.calc_ivar(varframe)
    ny = sciframe.shape[0]
    piximg = tilts * (ny-1)

    #
    ordpix = slitpix.copy()
    # Masks
    if bpm is not None:
        ordpix *= 1-bpm.astype(np.int)
    if crmask is not None:
        ordpix *= 1-crmask.astype(np.int)
    if tracemask is not None:
        ordpix *= (1-tracemask.astype(np.int))

    # Sky pixels for fitting
    fit_sky = (ordpix == slit) & (ivar > 0.) & (~edge_mask)
    isrt = np.argsort(piximg[fit_sky])
    wsky = piximg[fit_sky][isrt]
    sky = sciframe[fit_sky][isrt]
    sky_ivar = ivar[fit_sky][isrt]

    # All for evaluation
    all_slit = (slitpix == slit) & (~edge_mask)

    # Restrict fit to positive pixels only and mask out large outliers via a pre-fit to the log
    if (POS_MASK==True):
        pos_sky = (sky > 1.0) & (sky_ivar > 0.)
        if np.sum(pos_sky) > ny:
            lsky = np.log(sky[pos_sky])
            lsky_ivar = lsky * 0. + 0.1

            # Init bspline to get the sky breakpoints (kludgy)
            tmp = pydl.bspline(wsky[pos_sky], nord=4, bkspace=bsp)

            #skybkpt = bspline_bkpts(wsky[pos_sky], nord=4, bkspace=bsp $
            #, / silent)
            if False:
                plt.clf()
                ax = plt.gca()
                ax.scatter(wsky[pos_sky], lsky)
                #ax.scatter(wsky[~full_out], sky[~full_out], color='red')
                #ax.plot(wsky, yfit, color='green')
                plt.show()
                #debugger.set_trace()
            lskyset, outmask, lsky_fit, red_chi = utils.bspline_profile(
                wsky[pos_sky], lsky, lsky_ivar, np.ones_like(lsky),
                fullbkpt = tmp.breakpoints, upper=sigrej, lower=sigrej,
                kwargs_reject={'groupbadpix':True})
            res = (sky[pos_sky] - np.exp(lsky_fit)) * np.sqrt(sky_ivar[pos_sky])
            lmask = (res < 5.0) & (res > -4.0)
            sky_ivar[pos_sky] = sky_ivar[pos_sky] * lmask

    # Full fit now
    full_bspline = pydl.bspline(wsky, nord=4, bkspace=bsp)
    skyset, full_out, yfit, _ = utils.bspline_profile(
        wsky, sky, sky_ivar, np.ones_like(sky),
        fullbkpt=full_bspline.breakpoints,
        upper=sigrej, lower=sigrej, kwargs_reject={'groupbadpix':True, 'maxrej': 10})
    # Evaluate and save
    bgframe[all_slit] = skyset.value(piximg[all_slit])[0] #, skyset)

    # Debugging/checking
    if PLOT_FIT:
        goodbk = skyset.mask
        yfit_bkpt = np.interp(skyset.breakpoints[goodbk], wsky,yfit)
        plt.clf()
        ax = plt.gca()
        was_fit = (sky_ivar > 0.0)
        was_fit_and_masked = (was_fit == True) & (full_out == False)
        ax.plot(wsky[was_fit], sky[was_fit], color='k', marker='o', markersize=0.4, mfc='k', fillstyle='full', linestyle='None')
        ax.plot(wsky[was_fit_and_masked], sky[was_fit_and_masked], color='red', marker='+', markersize=1.5, mfc='red', fillstyle='full', linestyle='None')
        ax.plot(wsky, yfit, color='cornflowerblue')
        ax.plot(skyset.breakpoints[goodbk], yfit_bkpt, color='lawngreen', marker='o', markersize=2.0, mfc='lawngreen', fillstyle='full', linestyle='None')
        ax.set_ylim((0.99*yfit.min(),1.01*yfit.max()))
        plt.show()

    # Return
    return bgframe


# ToDO Fix masking logic. This code should also take an ivar for consistency with rest of extraction
def global_skysub(image, ivar, tilts, thismask, slit_left, slit_righ, inmask = None, bsp=0.6, sigrej=3., TRIM_EDG = (3,3), POS_MASK=True, PLOT_FIT=False):
    """
    Perform global sky subtraction on an input slit

    Parameters
    ----------
    image : ndarray
          Frame to be sky subtracted

    thismask : numpy boolean array
      Specifies pixels in the slit in question
    edgemask : ndarray
      Mask edges of the slit
    sciframe : ndarray
      science frame
    varframe : ndarray
      Variance array
    tilts : ndarray
      Tilts of the wavelengths
    bpm : ndarray, optional
      Bad pixel mask
    crmask : ndarray
      Cosmic ray mask
    tracemask : ndarray
      Object mask
    bsp : float
      Break point spacing
    sigrej : float
      rejection

    Returns
    -------
    bgframe : ndarray
      Sky background image

    """

    # Synthesize ximg, and edgmask  from slit boundaries. Doing this outside this
    # routine would save time. But this is pretty fast, so we just do it here to make the interface simpler.

    ximg, edgmask = pixels.ximg_and_edgemask(slit_left, slit_righ, thismask, trim_edg=TRIM_EDG)

    # Init
    nspec = image.shape[0]
    piximg = tilts * (nspec-1)
    if inmask is None:
        inmask = np.copy(thismask)

    # Sky pixels for fitting
    fit_sky = (thismask == True) & (ivar > 0.0) & (inmask == True) & (edgmask == False)
    isrt = np.argsort(piximg[fit_sky])
    wsky = piximg[fit_sky][isrt]
    sky = image[fit_sky][isrt]
    sky_ivar = ivar[fit_sky][isrt]

    # Restrict fit to positive pixels only and mask out large outliers via a pre-fit to the log
    if (POS_MASK==True):
        pos_sky = (sky > 1.0) & (sky_ivar > 0.0)
        if np.sum(pos_sky) > nspec:
            lsky = np.log(sky[pos_sky])
            lsky_ivar = lsky * 0. + 0.1

            # Init bspline to get the sky breakpoints (kludgy)
            tmp = pydl.bspline(wsky[pos_sky], nord=4, bkspace=bsp)

            #skybkpt = bspline_bkpts(wsky[pos_sky], nord=4, bkspace=bsp $
            #, / silent)
            if False:
                plt.clf()
                ax = plt.gca()
                ax.scatter(wsky[pos_sky], lsky)
                #ax.scatter(wsky[~full_out], sky[~full_out], color='red')
                #ax.plot(wsky, yfit, color='green')
                plt.show()
                #debugger.set_trace()
            lskyset, outmask, lsky_fit, red_chi = utils.bspline_profile(
                wsky[pos_sky], lsky, lsky_ivar, np.ones_like(lsky),
                fullbkpt = tmp.breakpoints, upper=sigrej, lower=sigrej,
                kwargs_reject={'groupbadpix': True, 'maxrej': 10})
            res = (sky[pos_sky] - np.exp(lsky_fit)) * np.sqrt(sky_ivar[pos_sky])
            lmask = (res < 5.0) & (res > -4.0)
            sky_ivar[pos_sky] = sky_ivar[pos_sky] * lmask


    # Full fit now
    full_bspline = pydl.bspline(wsky, nord=4, bkspace=bsp)
    skyset, outmask, yfit, _ = utils.bspline_profile(wsky, sky, sky_ivar, np.ones_like(sky),
                                                       fullbkpt=full_bspline.breakpoints,upper=sigrej, lower=sigrej,
                                                       kwargs_reject={'groupbadpix':True, 'maxrej': 10})
    # Evaluate and save
    bgframe, _ = skyset.value(piximg[thismask])

    # Debugging/checking
    if PLOT_FIT:
        goodbk = skyset.mask
        yfit_bkpt = np.interp(skyset.breakpoints[goodbk], wsky,yfit)
        plt.clf()
        ax = plt.gca()
        was_fit_and_masked = (outmask == False)
        ax.plot(wsky, sky, color='k', marker='o', markersize=0.4, mfc='k', fillstyle='full', linestyle='None')
        ax.plot(wsky[was_fit_and_masked], sky[was_fit_and_masked], color='red', marker='+', markersize=1.5, mfc='red', fillstyle='full', linestyle='None')
        ax.plot(wsky, yfit, color='cornflowerblue')
        ax.plot(skyset.breakpoints[goodbk], yfit_bkpt, color='lawngreen', marker='o', markersize=4.0, mfc='lawngreen', fillstyle='full', linestyle='None')
        ax.set_ylim((0.99*yfit.min(),1.01*yfit.max()))
        plt.show()

    # Return
    # ToDO worth thinking about whether we want to return a mask here. It makese no sense to return outmask
    # in its present form though since that does not refer to the whole image.
#    return bgframe, outmask
    return bgframe



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

    (sset1, outmask_good1, yfit1, red_chi1) = utils.bspline_profile(wave[good], data[good], ivar[good], profile_basis[good, :],
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

    (sset, outmask_good, yfit, red_chi) = utils.bspline_profile(wave[good], data[good], ivar[good] * mask1,
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

def local_skysub_extract(sciimg, sciivar, tilts, waveimg, global_sky, rn2_img, thismask, slit_left, slit_righ, sobjs, bsp,
    TRIM_EDG = (3,3), STD = False, PROF_NSIGMA = None, niter=4, box_rad = 7, sigrej = 3.5, skysample = False,
    FULLWELL = 5e5,MINWELL = -1000.0, SN_GAUSS = 3.0, COADD_2D = False, SHOW_RESIDS=False):


    ximg, edgmask = pixels.ximg_and_edgemask(slit_left, slit_righ, thismask, trim_edg = TRIM_EDG)

    nspat = sciimg.shape[1]
    nspec = sciimg.shape[0]
    piximg = tilts * (nspec-1)

    # Copy the specobjs that will be the output
    nobj = len(sobjs)
    # specobjs = copy.deepcopy(specobjs_in)

    if (PROF_NSIGMA is None):
        prof_nsigma1 = np.full(len(sobjs), None)
    elif len(PROF_NSIGMA) == 1:
        prof_nsigma1 = np.full(nobj, PROF_NSIGMA)
    elif len(PROF_NSIGMA) == nobj:
        prof_nsigma1 = PROF_NSIGMA
    else:
        raise ValueError('Invalid size for PROF_NSIGMA.')

    # Set some rejection parameters based on whether this is a STD or not. Only reject extreme outliers for standards
    # since super high S/N and low order profile models imply we will always have large outliers
    if STD is True:
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

    # Create the imagews that will be returned
    outmask = (sciivar > 0.0) & thismask & np.isfinite(sciimg) & (sciimg < FULLWELL) & (sciimg > MINWELL)
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
        # The default value of maskwidth = 3.0 * FWHM = 7.05 * sigma in long_objfind with a log(S/N) correction for bright objects
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
        for iiter in range(1, niter):
            msgs.info("Iteration # " + "{:2d}".format(iiter) + " of " + "{:2d}".format(niter))
            img_minsky = sciimg - skyimage
            for ii in range(objwork):
                iobj = group[ii]
                if iiter == 1:
                    # If this is the first iteration, print status message. Initiate profile fitting with a simple
                    # boxcar extraction.
                    msgs.info("-------------------REDUCING-------------------")
                    msgs.info("Fitting profile for obj #: " + "{:}".format(sobjs[iobj].objid) + " of {:}".format(nobj))
                    msgs.info(
                        "At x = {:5.2f}".format(sobjs[iobj].spat_pixpos) + " on slit # {:}".format(sobjs[iobj].slitid))
                    msgs.info("----------------------------------------------")
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
                    if sobjs[iobj].optimal['MASK_OPT'].any():
                        flux = sobjs[iobj].optimal['FLUX_OPT']
                        fluxivar = sobjs[iobj].optimal['IVAR_OPT']
                        wave = sobjs[iobj].optimal['WAVE_OPT']

                if wave.any():
                    (profile_model, xnew, fwhmfit, med_sn2) = extract.fit_profile(img_minsky[ipix], (modelivar * outmask)[ipix],
                                                                          waveimg[ipix],
                                                                          sobjs[iobj].trace_spat - mincol,
                                                                          wave, flux, fluxivar,
                                                                          thisfwhm=sobjs[iobj].fwhm,
                                                                          hwidth=sobjs[iobj].maskwidth,
                                                                          PROF_NSIGMA=sobjs[iobj].prof_nsigma,
                                                                          SN_GAUSS=SN_GAUSS)
                    # Update the object profile and the fwhm and mask parameters
                    obj_profiles[ipix[0], ipix[1], ii] = profile_model
                    sobjs[iobj].trace_spat = xnew + mincol
                    sobjs[iobj].fwhmfit = fwhmfit
                    sobjs[iobj].fwhm = np.median(fwhmfit)
                    mask_fact = 1.0 + 0.5 * np.log10(np.fmax(np.sqrt(np.fmax(med_sn2, 0.0)), 1.0))
                    maskwidth = 3.0 * np.median(fwhmfit) * mask_fact
                    if sobjs[iobj].prof_nsigma is None:
                        sobjs[iobj].maskwidth = maskwidth
                    else:
                        sobjs[iobj].maskwidth = sobjs[iobj].prof_nsigma * (sobjs[iobj].fwhm / 2.3548)

                else:
                    msgs.warn("Bad extracted wavelengths in local_skysub")
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
                if COADD_2D is False:
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


        '''   # If requested display the model fits for this grouping
        if SHOW_2D == True:
            viewer, ch = ginga.show_image((sciimg - skyimage) * np.sqrt(modelivar))
            # TODO figure out a way to overplot the pixels that were masked in red like as a scatter plot
            for ii in range(objwork):
                iobj = group[ii]
                if sobjs[iobj].HAND_FLAG == False:
                    color = 'green'
                else:
                    color = 'orange'
                ginga.show_trace(viewer, ch, sobjs[iobj].trace_spat, sobjs[iobj].idx, color=color)
        '''
    # If requested display the model fits for this slit
    if SHOW_RESIDS == True:
        # TODO add error checking here to see if ginga exists
        viewer, ch = ginga.show_image((sciimg - skyimage - objimage) * np.sqrt(modelivar)*thismask)
        # TODO figure out a way to overplot the pixels that were masked in red like as a scatter plot
        for spec in sobjs
            if spec.HAND_FLAG == False:
                color = 'green'
            else:
                color = 'orange'
            ginga.show_trace(viewer, ch, spec.trace_spat, spec.idx, color=color)
        # TODO figure out a way to set the cuts of the ginga viewer to go from -5 to 5


    return (skyimage[thismask], objimage[thismask], modelivar[thismask], outmask[thismask])

def order_pixels(pixlocn, lord, rord):
    """
    Based on physical pixel locations, determine which pixels are within the orders
    """

    sz_x, sz_y, _ = pixlocn.shape
    sz_o = lord.shape[1]

    outfr = np.zeros((sz_x, sz_y), dtype=int)

    for y in range(sz_y):
        for o in range(sz_o):
            indx = (lord[:,o] < rord[:,o]) & (pixlocn[:,y,1] > lord[:,o]) \
                   & (pixlocn[:,y,1] < rord[:,o])
            indx |= ( (lord[:,o] > rord[:,o]) & (pixlocn[:,y,1] < lord[:,o])
                      & (pixlocn[:,y,1] > rord[:,o]) )
            if np.any(indx):
                # Only assign a single order to a given pixel
                outfr[indx,y] = o+1
                break

    return outfr


# This code is deprecated and replaced by bg_subtraction_slit
def orig_bg_subtraction_slit(tslits_dict, pixlocn,
                        slit, tilts, sciframe, varframe, bpix, crpix,
                        settings,
                        tracemask=None,
                        rejsigma=3.0, maskval=-999999.9,
                        method='bspline'):
    """ Extract a science target and background flux
    :param slf:
    :param sciframe:
    :param varframe:
    :return:
    """
    # Unpack tslits
    lordloc = tslits_dict['lcen']
    rordloc = tslits_dict['rcen']
    slitpix = tslits_dict['slitpix']
    # Init
    bgframe = np.zeros_like(sciframe)

    # Begin the algorithm
    # Find which pixels are within the order edges
    msgs.info("Identifying pixels within each order")
    ordpix = order_pixels(pixlocn, lordloc*0.95+rordloc*0.05, lordloc*0.05+rordloc*0.95)

    msgs.info("Applying bad pixel mask")
    ordpix *= (1-bpix.astype(np.int)) * (1-crpix.astype(np.int))
    if tracemask is not None: ordpix *= (1-tracemask.astype(np.int))

    # Construct an array of pixels to be fit with a spline
    whord = np.where(ordpix == slit+1)
    xvpix  = tilts[whord]
    scipix = sciframe[whord]
    xargsrt = np.argsort(xvpix, kind='mergesort')
    sxvpix  = xvpix[xargsrt]
    sscipix = scipix[xargsrt]

    # Reject deviant pixels -- step through every 1.0/sciframe.shape[0] in sxvpix and reject significantly deviant pixels
    edges = np.linspace(min(0.0,np.min(sxvpix)),max(1.0,np.max(sxvpix)),sciframe.shape[0])
    fitcls = np.zeros(sciframe.shape[0])

    # Identify science target
    maskpix = np.zeros(sxvpix.size)
    msgs.info("Identifying pixels containing the science target")
    msgs.work("Speed up this step with multi-processing")
    for i in range(sciframe.shape[0]-1):
        wpix = np.where((sxvpix>=edges[i]) & (sxvpix<=edges[i+1]))
        if (wpix[0].size>5):
            txpix = sxvpix[wpix]
            typix = sscipix[wpix]
            msk, cf = utils.robust_polyfit(txpix, typix, 0, sigma=rejsigma)
            maskpix[wpix] = msk
            #fitcls[i] = cf[0]
            wgd=np.where(msk == 0)
            szt = np.size(wgd[0])
            if szt > 8:
                fitcls[i] = np.mean(typix[wgd][szt//2-3:szt//2+4]) # Average the 7 middle pixels
                #fitcls[i] = np.mean(np.random.shuffle(typix[wgd])[:5]) # Average the 5 random pixels
            else:
                fitcls[i] = cf[0]
    '''
    else:
        debugger.set_trace()
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
        scicopy[np.where(ordpix == slit)] = maskval
        scitilts, _ = artrace.model_tilt(det, scicopy, guesstilts=tilts.copy(),
                                         censpec=fitcls, maskval=maskval, plotQA=True)
        xvpix  = scitilts[whord]
        scipix = sciframe[whord]
        varpix = varframe[whord]
        mskpix = tracemask[whord]
        xargsrt = np.argsort(xvpix, kind='mergesort')
        sxvpix  = xvpix[xargsrt]
        sscipix = scipix[xargsrt]
        svarpix = varpix[xargsrt]
        maskpix = mskpix[xargsrt]
    '''

    # Check the mask is reasonable
    scimask = sciframe.copy()
    rxargsrt = np.argsort(xargsrt, kind='mergesort')
    scimask[whord] *= (1.0-maskpix)[rxargsrt]

    # Now trace the sky lines to get a better estimate of the spectral tilt during the observations
    scifrcp = scimask.copy()
    scifrcp[whord] += (maskval*maskpix)[rxargsrt]
    scifrcp[np.where(ordpix != slit+1)] = maskval

    # Check tilts? -- Can also be error in flat fielding or slit illumination
    if False:
        idx = 1893
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
    if method == 'bspline':
        msgs.info("Using bspline sky subtraction")
        gdp = (scifrcp != maskval) & (ordpix == slit+1) & (varframe > 0.)
        srt = np.argsort(tilts[gdp])
        #bspl = utils.func_fit(tilts[gdp][srt], scifrcp[gdp][srt], 'bspline', 3,
        #                        **settings.argflag['reduce']['skysub']['bspline'])
        ivar = utils.calc_ivar(varframe)
        mask, bspl = utils.robust_polyfit(tilts[gdp][srt], scifrcp[gdp][srt], 3,
                                            function='bspline',
                                            weights=np.sqrt(ivar)[gdp][srt],
                                            sigma=5., maxone=False,
                                            bspline_par=settings['skysub']['bspline'])
        # Just those in the slit
        in_slit = np.where(slitpix == slit+1)
        bgf_flat = utils.func_val(bspl, tilts[in_slit].flatten(), 'bspline')
        #bgframe = bgf_flat.reshape(tilts.shape)
        bgframe[in_slit] = bgf_flat
        if msgs._debug['sky_sub']:
            plt_bspline_sky(tilts, scifrcp, bgf_flat, in_slit, gdp)
            debugger.set_trace()
    else:
        msgs.error('Not ready for this method for skysub {:s}'.format(method))

    if np.sum(np.isnan(bgframe)) > 0:
        msgs.warn("NAN in bgframe.  Replacing with 0")
        bad = np.isnan(bgframe)
        bgframe[bad] = 0.

    # Plot to make sure that the result is good
    return bgframe



def plt_bspline_sky(tilts, scifrcp, bgf_flat, inslit, gdp):
    # Setup
    srt = np.argsort(tilts[inslit].flatten())
    # Plot
    plt.close()
    plt.clf()
    ax = plt.gca()
    ax.scatter(tilts[gdp]*tilts.shape[0], scifrcp[gdp], marker='o')
    ax.plot(tilts[inslit].flatten()[srt]*tilts.shape[0], bgf_flat[srt], 'r-')
    plt.show()


