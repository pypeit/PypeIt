""" Module for sky subtraction

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""
import sys
import os

import numpy as np

from scipy import ndimage
from scipy.special import ndtr

from matplotlib import pyplot as plt

from IPython import embed

from pypeit import msgs, utils, ginga, bspline
from pypeit.images import maskimage
from pypeit.core import pixels, extract, basis
from pypeit.core.moment import moment1d

def skysub_npoly(thismask):
    """
    Utility routine used by global_skysub and local_skysub_extract. Determine the order for the spatial
    polynomial for global sky subtraction and local sky subtraction.

    Args:
        thismask : ndarray, bool, shape (nspec, nspat)
            Specifies pixels in the slit in question

    Returns:
        int: Order of polynomial
    """
    slit_width = np.sum(thismask,axis=1)
    med_slit_width = np.median(slit_width[slit_width > 0])
    nspec_eff = np.sum(slit_width > 0.5*med_slit_width)
    npercol = np.fmax(np.floor(np.sum(thismask)/nspec_eff), 1.0)
    # Demand at least 10 pixels per row (on average) per degree of the polynomial
    if npercol > 100:
        npoly = 3
    elif npercol > 40:
        npoly = 2
    else:
        npoly = 1

    return npoly


def global_skysub(image, ivar, tilts, thismask, slit_left, slit_righ, inmask = None, bsp=0.6, sigrej=3.0, maxiter=35,
                  trim_edg=(3,3), pos_mask=True, show_fit=False, no_poly=False, npoly=None):
    """
    Perform global sky subtraction on an input slit

    Args:
        image: float ndarray, shape (nspec, nspat)
            Frame to be sky subtracted
        ivar: float ndarray, shape (nspec, nspat)
            Inverse variance image
        tilts: float ndarray, shape (nspec, nspat)
            Tilts indicating how wavelengths move across the slit
        thismask : numpy boolean array, shape (nspec, nspat)
            Specifies pixels in the slit in question
        slit_left: ndarray of shape (nspec, 1) or (nspec)
            Left slit boundary in floating point pixels.
        slit_righ: ndarray of shape (nspec, 1) or (nspec)
            Right slit boundary in floating point pixels.
        inmask: boolean ndarray, shape (nspec, nspat), default inmask = None
            Input mask for pixels not to be included in sky subtraction
            fits. True = Good (not masked), False = Bad (masked)
        bsp: float, default bsp = 0.6
            break point spacing in pixel units
        sigrej : float, default sigrej = 3.0
            sigma rejection threshold
        no_poly: bool, optional
            Do not incldue polynomial basis
        trim_edg: tuple of floats  (left_edge, right_edge), default (3,3)
            indicates how many pixels to trim from left and right slit
            edges for creating the edgemask. These pixels are excluded
            from sky subtraction fits.
        pos_mask: boolean, defualt pos_mask = True
            First do a prelimnary fit to the log of the sky (i.e.
            positive pixels only). Then use this fit to create an input
            mask from the residuals lmask = (res < 5.0) & (res > -4.0)
            for the full fit.  NOTE: pos_mask should be False for
            near-IR sky residual subtraction, since fitting the log(sky)
            requires that the counts are positive which will not be the
            case for i.e. an A-B image. Thus the routine will fail if
            pos_mask is not set to False.

        show_fit: boolean, default show_fit = False
            Plot a fit of the sky pixels and model fit to the screen.
            This feature will block further execution until the screen
            is closed.

    Returns:
        `numpy.ndarray`_: Returns the model sky background at the pixels
        where thismask is True::

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
        inmask = (ivar > 0.0) & thismask & np.isfinite(image) & np.isfinite(ivar)

    # Sky pixels for fitting
    inmask_in = thismask & (ivar > 0.0) & inmask & np.invert(edgmask)
    isrt = np.argsort(piximg[thismask])
    pix = piximg[thismask][isrt]
    sky = image[thismask][isrt]
    sky_ivar = ivar[thismask][isrt]
    ximg_fit = ximg[thismask][isrt]
    inmask_fit = inmask_in[thismask][isrt]
    #spatial = spatial_img[fit_sky][isrt]

    # Restrict fit to positive pixels only and mask out large outliers via a pre-fit to the log.
    if pos_mask:
        pos_sky = (sky > 1.0) & (sky_ivar > 0.0)
        if np.sum(pos_sky) > nspec:
            lsky = np.log(sky[pos_sky])
            lsky_ivar = inmask_fit[pos_sky].astype(float)/3.0**2  # set errors to just be 3.0 in the log
            #lsky_ivar = np.full(lsky.shape, 0.1)
            # Init bspline to get the sky breakpoints (kludgy)
            #tmp = pydl.bspline(wsky[pos_sky], nord=4, bkspace=bsp)
            lskyset, outmask, lsky_fit, red_chi, exit_status \
                    = utils.bspline_profile(pix[pos_sky], lsky, lsky_ivar, np.ones_like(lsky),
                                            ingpm=inmask_fit[pos_sky], upper=sigrej, lower=sigrej,
                                            kwargs_bspline={'bkspace':bsp},
                                            kwargs_reject={'groupbadpix': True, 'maxrej': 10})
            res = (sky[pos_sky] - np.exp(lsky_fit)) * np.sqrt(sky_ivar[pos_sky])
            lmask = (res < 5.0) & (res > -4.0)
            sky_ivar[pos_sky] = sky_ivar[pos_sky] * lmask
            inmask_fit[pos_sky]=(sky_ivar[pos_sky] > 0.0) & lmask

    # Include a polynomial basis?
    if no_poly:
        poly_basis = np.ones_like(sky)
        npoly_fit = 1
    else:
        npoly_fit = skysub_npoly(thismask) if npoly is None else npoly
        poly_basis = basis.flegendre(2.0*ximg_fit - 1.0, npoly_fit)

    # Full fit now
    #full_bspline = pydl.bspline(wsky, nord=4, bkspace=bsp, npoly = npoly)
    #skyset, outmask, yfit, _ = utils.bspline_profile(wsky, sky, sky_ivar, poly_basis,
    #                                                   fullbkpt=full_bspline.breakpoints,upper=sigrej, lower=sigrej,
    #                                                   kwargs_reject={'groupbadpix':True, 'maxrej': 10})


    # Perform the full fit now
    msgs.info("Full fit in global sky sub.")
    skyset, outmask, yfit, _, exit_status \
            = utils.bspline_profile(pix, sky, sky_ivar, poly_basis, ingpm=inmask_fit, nord=4,
                                    upper=sigrej, lower=sigrej, maxiter=maxiter,
                                    kwargs_bspline={'bkspace':bsp},
                                    kwargs_reject={'groupbadpix':True, 'maxrej': 10})
    # TODO JFH This is a hack for now to deal with bad fits for which iterations do not converge. This is related
    # to the groupbadpix behavior requested for the djs_reject rejection. It would be good to
    # better understand what this functionality is doing, but it makes the rejection much more quickly approach a small
    # chi^2
    if exit_status == 1:
        msgs.warn('Maximum iterations reached in bspline_profile global sky-subtraction for npoly={:d}.'.format(npoly_fit) +
                  msgs.newline() +
                  'Redoing sky-subtraction without polynomial degrees of freedom')
        poly_basis = np.ones_like(sky)
        # Perform the full fit now
        skyset, outmask, yfit, _, exit_status \
                = utils.bspline_profile(pix, sky, sky_ivar, poly_basis, ingpm=inmask_fit, nord=4,
                                        upper=sigrej, lower=sigrej, maxiter=maxiter,
                                        kwargs_bspline={'bkspace': bsp},
                                        kwargs_reject={'groupbadpix': False, 'maxrej': 10})

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
        ax.plot(pix[inmask_fit], sky[inmask_fit], color='k', marker='o', markersize=0.4, mfc='k', fillstyle='full', linestyle='None')
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



# TODO -- This needs JFH docs, desperately
def skyoptimal(wave, data, ivar, oprof, sortpix, sigrej=3.0, npoly=1, spatial=None, fullbkpt=None):
    """
    Utility routine used by local_bg_subtraction_slit

    Args:
        wave:
        data:
        ivar:
        oprof (ndarray): Flattened object profile in this slit
        sortpix:
        sigrej:
        npoly:
        spatial:
        fullbkpt:

    Returns:
        ndarray, ndarray, ndarray:

    """


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
        poly_basis = basis.flegendre(x2, npoly)
        profile_basis = np.column_stack((oprof, poly_basis))

    relative_mask = (np.sum(oprof, axis=1) > 1e-10)

    indx, = np.where(ivar[sortpix] > 0.0)
    ngood = indx.size
    good = sortpix[indx]
    good = good[wave[good].argsort()]
    relative, = np.where(relative_mask[good])

    outmask = np.zeros(wave.shape, dtype=bool)

    if ngood > 0:
        sset1, outmask_good1, yfit1, red_chi1, exit_status \
                = utils.bspline_profile(wave[good], data[good], ivar[good], profile_basis[good, :],
                                        fullbkpt=fullbkpt, upper=sigrej, lower=sigrej,
                                        relative=relative,
                                        kwargs_reject={'groupbadpix': True, 'maxrej': 5})
    else:
        msgs.warn('All pixels are masked in skyoptimal. Not performing local sky subtraction.')
        return np.zeros_like(wave), np.zeros_like(wave), outmask

    chi2 = (data[good] - yfit1) ** 2 * ivar[good]
    chi2_srt = np.sort(chi2)
    gauss_prob = 1.0 - 2.0 * ndtr(-1.2 * sigrej)
    sigind = int(np.fmin(np.rint(gauss_prob * float(ngood)), ngood - 1))
    chi2_sigrej = chi2_srt[sigind]
    mask1 = (chi2 < chi2_sigrej)

    msgs.info('2nd round....')
    msgs.info('Iter     Chi^2     Rejected Pts')
    if np.any(mask1):
        sset, outmask_good, yfit, red_chi, exit_status \
                = utils.bspline_profile(wave[good], data[good], ivar[good], profile_basis[good,:],
                                        ingpm=mask1, fullbkpt=fullbkpt, upper=sigrej, lower=sigrej,
                                        relative=relative,
                                        kwargs_reject={'groupbadpix': True, 'maxrej': 1})
    else:
        msgs.warn('All pixels are masked in skyoptimal after first round of rejection. Not performing local sky subtraction.')
        return np.zeros_like(wave), np.zeros_like(wave), outmask

    ncoeff = npoly + nobj
#    skyset = pydl.bspline(None, fullbkpt=sset.breakpoints, nord=sset.nord, npoly=npoly)
    skyset = bspline.bspline(None, fullbkpt=sset.breakpoints, nord=sset.nord, npoly=npoly)
    # Set coefficients for the sky.
    # The rehshape below deals with the different sizes of the coeff for npoly = 1 vs npoly > 1
    # and mirrors similar logic in the bspline.py
    skyset.coeff = sset.coeff[nobj:, :].reshape(skyset.coeff.shape)

    skyset.mask = sset.mask
    skyset.xmin = xmin
    skyset.xmax = xmax

    sky_bmodel, _ = skyset.value(wave, x2=spatial)

    obj_bmodel = np.zeros(sky_bmodel.shape)
#    objset = pydl.bspline(None, fullbkpt=sset.breakpoints, nord=sset.nord)
    objset = bspline.bspline(None, fullbkpt=sset.breakpoints, nord=sset.nord)
    objset.mask = sset.mask
    for i in range(nobj):
        objset.coeff = sset.coeff[i, :]
        obj_bmodel1, _ = objset.value(wave)
        obj_bmodel = obj_bmodel + obj_bmodel1 * profile_basis[:, i]

    outmask[good] = outmask_good

    return sky_bmodel, obj_bmodel, outmask


def optimal_bkpts(bkpts_optimal, bsp_min, piximg, sampmask, samp_frac=0.80,
                  skyimage = None, min_spat=None, max_spat=None, debug=False):
    """

    Args:
        bsp_min: float
           Desired B-spline breakpoint spacing in pixels
        piximg: ndarray float, shape = (nspec, nspat)
           Image containing the pixel sampling, i.e. (nspec-1)*tilts
        sampmask: ndarray, bool
           Boolean array indicating the pixels for which the B-spline fit will actually be evaluated. True = Good, False=Bad
    Optional Args:
        samp_frac: float, default = 0.8
           The fraction of spectral direction pixels required to have a sampling difference < bsp_min in order to instead
           adopt a uniform break point spacing, rather adopting the optimally spaced breakpoints.
        skyimage: ndarray, shape = (nspec, nspat), default = None
           Sky model image used only for QA.
        min_spat: float, default = None
           Minimum spatial pixel used for local sky subtraction fitting. Only used for title of QA plot.
        max_spat: float, defualt = None
           Maximum spatial pixel used for local sky subtraction fitting. Only used for title of QA plot.
        debug: bool, default = False
           Show QA plot to debug breakpoint placing.

    Returns:
        fullbkpt: ndarray, float
           Locations of the optimally sampled breakpoints

    """

    pix = piximg[sampmask]
    isrt = pix.argsort()
    pix = pix[isrt]
    piximg_min = pix.min()
    piximg_max = pix.max()
#    bset0 = pydl.bspline(pix, nord=4, bkspace=bsp_min)
    bset0 = bspline.bspline(pix, nord=4, bkspace=bsp_min)
    fullbkpt_grid = bset0.breakpoints
    keep = (fullbkpt_grid >= piximg_min) & (fullbkpt_grid <= piximg_max)
    fullbkpt_grid = fullbkpt_grid[keep]
    used_grid = False
    if not bkpts_optimal:
        msgs.info('bkpts_optimal = False --> using uniform bkpt spacing spacing: bsp={:5.3f}'.format(bsp_min))
        fullbkpt = fullbkpt_grid
        used_grid = True
    else:
        piximg_temp = np.ma.array(np.copy(piximg))
        piximg_temp.mask = np.invert(sampmask)
        samplmin = np.ma.min(piximg_temp,fill_value=np.inf,axis=1)
        samplmin = samplmin[np.invert(samplmin.mask)].data
        samplmax = np.ma.max(piximg_temp,fill_value=-np.inf,axis=1)
        samplmax = samplmax[np.invert(samplmax.mask)].data
        if samplmax.size != samplmin.size:
            msgs.error('This should not happen')
        nbkpt = samplmax.size
        # Determine the sampling. dsamp represents the gap in spectral pixel (wavelength) coverage between
        # subsequent spectral direction pixels in the piximg, i.e. it is the difference between the minimum
        # value of the piximg at spectral direction pixel i+1, and the maximum value of the piximg at spectral
        # direction pixel i. A negative value dsamp < 0 implies continuous sampling with no gaps, i.e. the
        # the arc lines are sufficiently tilted that there is no sampling gap.
        dsamp_init = np.roll(samplmin, -1) - samplmax
        dsamp_init[nbkpt - 1] = dsamp_init[nbkpt - 2]
        kernel_size = int(np.fmax(np.ceil(dsamp_init.size*0.01)//2*2 + 1,15))  # This ensures kernel_size is odd
        dsamp_med = ndimage.filters.median_filter(dsamp_init, size=kernel_size, mode='reflect')
        boxcar_size = int(np.fmax(np.ceil(dsamp_med.size*0.005)//2*2 + 1,5))
        # Boxcar smooth median dsamp
        kernel = np.ones(boxcar_size)/ float(boxcar_size)
        dsamp = ndimage.convolve(dsamp_med, kernel, mode='reflect')
        # if more than samp_frac of the pixels have dsamp < bsp_min than just use a uniform breakpoint spacing
        if np.sum(dsamp <= bsp_min) > samp_frac*nbkpt:
            msgs.info('Sampling of wavelengths is nearly continuous.')
            msgs.info('Using uniform bkpt spacing: bsp={:5.3f}'.format(bsp_min))
            fullbkpt = fullbkpt_grid
            used_grid = True
        else:
            fullbkpt_orig = samplmax + dsamp/2.0
            fullbkpt_orig.sort()
            # Compute the distance between breakpoints
            dsamp_bkpt = fullbkpt_orig-np.roll(fullbkpt_orig, 1)
            dsamp_bkpt[0] = dsamp_bkpt[1]
            # Good breakpoints are those that are at least separated by our original desired bkpt spacing
            igood = dsamp_bkpt >= bsp_min
            if np.any(igood):
                fullbkpt_orig = fullbkpt_orig[igood]
            fullbkpt = fullbkpt_orig.copy()
            # Recompute the distance between breakpoints
            dsamp_bkpt = fullbkpt_orig-np.roll(fullbkpt_orig, 1)
            dsamp_bkpt[0] = dsamp_bkpt[1]
            nbkpt = fullbkpt_orig.size
            for ibkpt in range(nbkpt):
                dsamp_eff = np.fmax(dsamp_bkpt[ibkpt], bsp_min)
                # can we fit in another bkpt?
                if dsamp_bkpt[ibkpt] > 2*dsamp_eff:
                    nsmp = int(np.fmax(np.floor(dsamp_bkpt[ibkpt]/dsamp_eff),2))
                    bkpt_new = fullbkpt_orig[ibkpt - 1] + (np.arange(nsmp - 1) + 1)*dsamp_bkpt[ibkpt]/float(nsmp)
                    indx_arr = np.where(fullbkpt == fullbkpt_orig[ibkpt-1])[0]
                    if len(indx_arr) > 0:
                        indx_bkpt = indx_arr[0]
                        if indx_bkpt == 0:
                            fullbkpt = np.hstack((fullbkpt[0], bkpt_new, fullbkpt[indx_bkpt + 1:]))
                        elif indx_bkpt == (fullbkpt.size-2):
                            fullbkpt = np.hstack((fullbkpt[0:indx_bkpt], bkpt_new, fullbkpt[indx_bkpt + 1]))
                        else:
                            fullbkpt = np.hstack((fullbkpt[0:indx_bkpt], bkpt_new, fullbkpt[indx_bkpt + 1:]))

            fullbkpt.sort()
            keep = (fullbkpt >= piximg_min) & (fullbkpt <= piximg_max)
            fullbkpt = fullbkpt[keep]


    if debug:
        plt.figure(figsize=(14, 6))
        sky = skyimage[sampmask]
        sky = sky[isrt]
        # This is approximate and only for the sake of visualization:
        spat_samp_vec = np.sum(sampmask, axis=1)  # spatial sampling per spectral direction pixel
        spat_samp_med = np.median(spat_samp_vec[spat_samp_vec > 0])
        window_size = int(np.ceil(5 * spat_samp_med))
        sky_med_filt = utils.fast_running_median(sky, window_size)
        sky_bkpt_grid = np.interp(fullbkpt_grid, pix, sky_med_filt)
        sky_bkpt = np.interp(fullbkpt, pix, sky_med_filt)
        plt.clf()
        ax = plt.gca()
        ax.plot(pix, sky, color='k', marker='o', markersize=0.4, mfc='k', fillstyle='full', linestyle='None')
        # ax.plot(pix, sky_med_filt, color='cornflowerblue', label='median sky', linewidth=1.2)
        if used_grid == False:
            ax.plot(fullbkpt_grid, sky_bkpt_grid, color='lawngreen', marker='o', markersize=2.0, mfc='lawngreen',
                    fillstyle='full', linestyle='None', label='uniform bkpt grid')
            color = 'red'
            title_str = ''
        else:
            color = 'lawngreen'
            title_str = 'Used Grid: '
        ax.plot(fullbkpt, sky_bkpt, color=color, marker='o', markersize=4.0, mfc=color,
                fillstyle='full', linestyle='None', label='optimal bkpts')

        ax.set_ylim((0.99 * sky_med_filt.min(), 1.01 * sky_med_filt.max()))
        if min_spat is not None:
            plt.title(title_str + 'bkpt sampling spat pixels {:7.1f}-{:7.1f}'.format(min_spat, max_spat))
        plt.legend()
        plt.show()

    return fullbkpt


def local_skysub_extract(sciimg, sciivar, tilts, waveimg, global_sky, rn2_img, thismask, slit_left, slit_righ, sobjs,
                         spat_pix=None, adderr=0.01, bsp=0.6, inmask=None, extract_maskwidth=4.0, trim_edg=(3,3),
                         std=False, prof_nsigma=None, niter=4, box_rad=7, sigrej=3.5, bkpts_optimal=True,
                         debug_bkpts=False,sn_gauss=4.0, model_full_slit=False, model_noise=True, show_profile=False,
                         show_resids=False):
    """Perform local sky subtraction and  extraction

     Args:
        sciimg : numpy float 2-d array (nspec, nspat)
            sky-subtracted image
        sciivar : numpy float 2-d array (nspec, nspat)
            inverse variance of sky-subtracted image
        tilts: ndarray, (nspec, nspat)
            spectral tilts
        waveimg numpy float 2-d array (nspec, nspat)
            2-d wavelength map
        global_sky : ndarray (nspec, nspat)
            Global sky model produced by global_skysub
        rn2_img:
            Image with the read noise squared per pixel
        thismask : numpy boolean array, shape (nspec, nspat)
            Specifies pixels in the slit in question
        slit_left: ndarray of shape (nspec, 1) or (nspec)
            Left slit boundary in floating point pixels.
        slit_righ: ndarray of shape (nspec, 1) or (nspec)
            Right slit boundary in floating point pixels.
        sobjs:   SpecoObjs object
            Object containing the information about the objects found on
            the slit/order from objfind or ech_objfind
        spat_pix: float ndarray, shape (nspec, nspat), default = None
            Image containing the spatial location of pixels. If not
            input, it will be computed from ``spat_img =
            np.outer(np.ones(nspec), np.arange(nspat))``. This option
            should generally not be used unless one is extracting 2d
            coadds for which a rectified image contains sub-pixel
            spatial information.
        adderr: float, default = 0.01
            Error floor. The quantity adderr**2*sciframe**2 is added in
            qudarature to the variance to ensure that the S/N is never >
            1/adderr, effectively setting a floor on the noise or a
            ceiling on the S/N.
        bsp: float, default = 0.6
            Break point spacing in pixels for the b-spline sky subtraction.
        inmask: ndarray, bool, (nspec, nspat)
            Input mask with True=Good, False=Bad
        extract_maskwidth: float, default = 4.0
            Determines the initial size of the region in units of fwhm
            that will be used for local sky subtraction. This maskwidth
            is defined in the obfjind code, but is then updated here as
            the profile fitting improves the fwhm estimates
        trim_edg: tuple of ints of floats, default = (3,3)
            Number of pixels to be ignored on the (left,right) edges of
            the slit in object/sky model fits.
        std: bool, default = False
            This should be set to True if the object being extracted is
            a standards star so that the reduction parameters can be
            adjusted accordingly.
        prof_nsigma: int or float, default = None
            Number of sigmas that the object profile will be fit, i.e.
            the region extending from -prof_nsigma to +prof_nsigma will
            be fit where sigma = FWHM/2.35. This option should only be
            used for bright large extended source with tails in their
            light profile like elliptical galaxies. If prof_nsigma is
            set then the profiles will no longer be apodized by an
            exponential at large distances from the trace.
        niter: int, default = 4
            Number of iterations for successive profile fitting and local sky-subtraction
        box_rad: int or float, default = 7
            Boxcar radius in *pixels* used for boxcar extraction.
        sigrej:
            Outlier rejection threshold for sky and object fitting
        bkpts_optimal = bool, default = True
            Parameter governing whether spectral direction breakpoints
            for b-spline sky/object modeling are determined optimally.
            If ``bkpts_optima=True``, the optimal break-point spacing
            will be determined directly using the optimal_bkpts function
            by measuring how well we are sampling the sky using ``piximg
            = (nspec-1)*yilyd``. The bsp parameter in this case
            corresponds to the minimum distance between breakpoints
            which we allow.  If ``bkpts_optimal = False``, the
            break-points will be chosen to have a uniform spacing in
            pixel units sets by the bsp parameter, i.e.  using the
            bkspace functionality of the bspline class::

              bset = bspline.bspline(piximg_values, nord=4, bkspace=bsp)
              fullbkpt = bset.breakpoints

        debug_bkpts: bool, default=False
            Make an interactive plot to the screen to indicate how the
            breakpoints are being chosen.
        sn_gauss: int or float, default = 4.0
            The signal to noise threshold above which optimal extraction
            with non-parametric b-spline fits to the objects spatial
            profile will be performed. For objects with median S/N <
            sn_gauss, a Gaussian profile will simply be assumed because
            there is not enough S/N to justify performing a more
            complicated fit.
        model_full_slit: bool, default = False
            Set the maskwidth of the objects to be equal to the slit
            width/2 such that the entire slit will be modeled by the
            local skysubtraction. This mode is recommended for echelle
            spectra with reasonably narrow slits.
        model_noise: bool, default = True
            If True, the model of the object, sky will be combined with
            the rn2img to create (and iteratively update) a model
            inverse variance image. If False, a variance model will not
            be created and instead the input sciivar will always be
            taken to be the inverse variance. Note that in order for the
            noise model to make any sense one needs to be subtracting
            the sky and *not* the sky residuals. In other words, for
            near-IR reductions where difference imaging has been
            performed and this algorithm is used to fit out the sky
            residuals (but not the sky itself) one should definitely set
            model_noise=False since otherwise the code will attempt to
            create a noise model using sky residuals instead of the sky,
            which is incorrect (does not have the right count levels).
            In principle this could be improved if the user could pass
            in a model of what the sky is for near-IR difference imaging
            + residual subtraction
        show_profile: bool, default=False
            Show QA for the object profile fitting to the screen. Note
            that this will show interactive matplotlib plots which will
            block the execution of the code until the window is closed.
        show_resids:
            Show the

    Returns:
        :obj:`tuple`:  Returns (skyimage[thismask], objimage[thismask],
        modelivar[thismask], outmask[thismask])
    """
    # TODO Force traces near edges to always be extracted with a Gaussian profile.


    if inmask is None:
        inmask = (sciivar > 0.0) & thismask & np.isfinite(sciimg) & np.isfinite(sciivar)

    # Adjust maskwidths of the objects such that we will apply the local_skysub_extract to the entire slit
    if model_full_slit:
        max_slit_width = np.max(slit_righ - slit_left)
        for spec in sobjs:
            spec.maskwidth = max_slit_width/2.0

    ximg, edgmask = pixels.ximg_and_edgemask(slit_left, slit_righ, thismask, trim_edg=trim_edg)

    nspat = sciimg.shape[1]
    nspec = sciimg.shape[0]
    piximg = tilts * (nspec-1)

    # Copy the specobjs that will be the output
    nobj = len(sobjs)

    # Set up the prof_nsigma
    if (prof_nsigma is None):
        prof_nsigma1 = np.full(len(sobjs), None)
    elif len(prof_nsigma) == 1:
        prof_nsigma1 = np.full(nobj, prof_nsigma)
    elif len(prof_nsigma) == nobj:
        prof_nsigma1 = prof_nsigma
    else:
        raise ValueError('Invalid size for prof_nsigma.')

    for iobj in range(nobj):
        sobjs[iobj].prof_nsigma = prof_nsigma1[iobj]

    # Set some rejection parameters based on whether this is a standard or not. Only reject extreme outliers for standards
    # since super high S/N and low order profile models imply we will always have large outliers
    if std is True:
        chi2_sigrej = 100.0
        #sigrej_ceil = 1e10
        sigrej = 50.0  # 25 wasn't enough for MagE 2x2 binning (probably undersampled)
    else:
        # TODO Why is this not an input parameter
        chi2_sigrej = 6.0
        #sigrej_ceil = 10.0
    # We will use this number later
    gauss_prob = 1.0 - 2.0 * ndtr(-sigrej)

    # Create the images that will be returned
    outmask = np.copy(inmask)
    modelivar = np.copy(sciivar)
    objimage = np.zeros_like(sciimg)
    skyimage = np.copy(global_sky)
    #varnoobj = np.abs(skyimage - np.sqrt(2.0) * np.sqrt(rn2_img)) + rn2_img

    # TODO Add a line of code here that updates the modelivar using the global sky if nobj = 0 and simply returns
    spat_img = np.outer(np.ones(nspec), np.arange(nspat))
    if spat_pix is None:
        spat_pix = spat_img

    xsize = slit_righ - slit_left
    # TODO Can this be simply replaced with spat_img above (but not spat_pix since that could have holes)
    spatial_img = thismask * ximg * (np.outer(xsize, np.ones(nspat)))

    # Loop over objects and group them
    i1 = 0
    while i1 < nobj:
        group = np.array([], dtype=np.int)
        group = np.append(group, i1)
        # The default value of maskwidth = 3.0 * FWHM = 7.05 * sigma in objfind with a log(S/N) correction for bright objects
        min_spat1 = np.maximum(sobjs[i1].TRACE_SPAT - sobjs[i1].maskwidth - 1, slit_left)
        max_spat1 = np.minimum(sobjs[i1].TRACE_SPAT + sobjs[i1].maskwidth + 1, slit_righ)
        for i2 in range(i1 + 1, nobj):
            left_edge = sobjs[i2].TRACE_SPAT - sobjs[i2].maskwidth - 1
            righ_edge = sobjs[i2].TRACE_SPAT + sobjs[i2].maskwidth + 1
            touch = (left_edge < max_spat1) & (sobjs[i2].TRACE_SPAT > slit_left) & (righ_edge > min_spat1)
            if touch.any():
                max_spat1 = np.minimum(np.maximum(righ_edge, max_spat1), slit_righ)
                min_spat1 = np.maximum(np.minimum(left_edge, min_spat1), slit_left)
                group = np.append(group, i2)
        # Create the local mask which defines the pixels that will be updated by local sky subtraction
        min_spat_img = np.outer(min_spat1, np.ones(nspat))
        max_spat_img = np.outer(max_spat1, np.ones(nspat))
        localmask = (spat_img > min_spat_img) & (spat_img < max_spat_img) & thismask
        npoly = skysub_npoly(localmask)
        # Keep for next iteration
        i1 = group.max() + 1
        # Some bookeeping to define the sub-image and make sure it does not land off the mask
        objwork = len(group)
        scope = np.sum(thismask, axis=0)
        iscp, = np.where(scope)
        imin = iscp.min()
        imax = iscp.max()
        min_spat = np.fmax(np.floor(min(min_spat1)), imin)
        max_spat = np.fmin(np.ceil(max(max_spat1)), imax)
        nc = int(max_spat - min_spat + 1)
        spec_vec = np.arange(nspec, dtype=np.intp)
        spat_vec = np.arange(min_spat, min_spat + nc, dtype=np.intp)
        ipix = np.ix_(spec_vec, spat_vec)
        obj_profiles = np.zeros((nspec, nspat, objwork), dtype=float)
        sigrej_eff = sigrej
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
                    msgs.info("Fitting profile for obj # " + "{:}".format(sobjs[iobj].OBJID) + " of {:}".format(nobj))
                    msgs.info("At x = {:5.2f}".format(sobjs[iobj].SPAT_PIXPOS) + " on slit # {:}".format(sobjs[iobj].slit_order))
                    msgs.info("------------------------------------------------------------------------------------------------------------")

                    # TODO -- Use extract_specobj_boxcar to avoid code duplication
                    extract.extract_boxcar(sciimg, modelivar, outmask, waveimg,
                                           skyimage, rn2_img, box_rad, sobjs[iobj])
                    flux = sobjs[iobj].BOX_COUNTS
                    fluxivar = sobjs[iobj].BOX_COUNTS_IVAR * sobjs[iobj].BOX_MASK
                    wave = sobjs[iobj].BOX_WAVE
                    #flux = moment1d(img_minsky * outmask, sobjs[iobj].TRACE_SPAT, 2*box_rad,
                    #                row=sobjs[iobj].trace_spec)[0]
                    #mvarimg = 1.0 / (modelivar + (modelivar == 0))
                    #mvar_box = moment1d(mvarimg * outmask, sobjs[iobj].TRACE_SPAT, 2*box_rad,
                    #                    row=sobjs[iobj].trace_spec)[0]
                    #pixtot = moment1d(0 * mvarimg + 1.0, sobjs[iobj].TRACE_SPAT, 2*box_rad,
                    #                  row=sobjs[iobj].trace_spec)[0]
                    #mask_box = moment1d(np.invert(outmask), sobjs[iobj].TRACE_SPAT, 2*box_rad,
                    #                    row=sobjs[iobj].trace_spec)[0] != pixtot
                    #box_denom = moment1d(waveimg > 0.0, sobjs[iobj].TRACE_SPAT, 2*box_rad,
                    #                     row=sobjs[iobj].trace_spec)[0]
                    #wave = moment1d(waveimg, sobjs[iobj].TRACE_SPAT, 2*box_rad,
                    #                row=sobjs[iobj].trace_spec)[0] \
                    #            / (box_denom + (box_denom == 0.0))
                    #fluxivar = mask_box / (mvar_box + (mvar_box == 0.0))
                else:
                    # For later iterations, profile fitting is based on an optimal extraction
                    last_profile = obj_profiles[:, :, ii]
                    trace = np.outer(sobjs[iobj].TRACE_SPAT, np.ones(nspat))
                    objmask = ((spat_img >= (trace - 2.0 * box_rad)) & (spat_img <= (trace + 2.0 * box_rad)))
                    # Boxcar
                    extract.extract_boxcar(sciimg, modelivar, (outmask & objmask),
                                                   waveimg, skyimage, rn2_img, box_rad,
                                                   sobjs[iobj])
                    # Optimal
                    extract.extract_optimal(sciimg, modelivar, (outmask & objmask), waveimg, skyimage, rn2_img, thismask,
                                            last_profile, box_rad, sobjs[iobj])
                    # If the extraction is bad do not update
                    if 'OPT_MASK' in sobjs[iobj].keys():
                        if sobjs[iobj].OPT_MASK.any():
                            flux = sobjs[iobj].OPT_COUNTS
                            fluxivar = sobjs[iobj].OPT_COUNTS_IVAR*sobjs[iobj].OPT_MASK
                            wave = sobjs[iobj].OPT_WAVE

                obj_string = 'obj # {:}'.format(sobjs[iobj].OBJID) + ' on slit # {:}'.format(sobjs[iobj].slit_order) + ', iter # {:}'.format(iiter) + ':'
                if wave.any():
                    sign = sobjs[iobj].sign
                    # TODO This is "sticky" masking. Do we want it to be?
                    profile_model, trace_new, fwhmfit, med_sn2 = extract.fit_profile(
                        sign*img_minsky[ipix], (modelivar * outmask)[ipix],waveimg[ipix], thismask[ipix], spat_pix[ipix], sobjs[iobj].TRACE_SPAT,
                        wave, sign*flux, fluxivar, inmask = outmask[ipix],
                        thisfwhm=sobjs[iobj].FWHM, maskwidth=sobjs[iobj].maskwidth,
                        prof_nsigma=sobjs[iobj].prof_nsigma, sn_gauss=sn_gauss, obj_string=obj_string,
                        show_profile=show_profile)
                    #proc_list.append(show_proc)

                    # Update the object profile and the fwhm and mask parameters
                    obj_profiles[ipix[0], ipix[1], ii] = profile_model
                    sobjs[iobj].TRACE_SPAT = trace_new
                    sobjs[iobj].FWHMFIT = fwhmfit
                    sobjs[iobj].FWHM = np.median(fwhmfit)
                    mask_fact = 1.0 + 0.5 * np.log10(np.fmax(np.sqrt(np.fmax(med_sn2, 0.0)), 1.0))
                    maskwidth = extract_maskwidth*np.median(fwhmfit) * mask_fact
                    if sobjs[iobj].prof_nsigma is None:
                        sobjs[iobj].maskwidth = maskwidth
                    else:
                        sobjs[iobj].maskwidth = sobjs[iobj].prof_nsigma * (sobjs[iobj].FWHM / 2.3548)

                else:
                    msgs.warn("Bad extracted wavelengths in local_skysub_extract")
                    msgs.warn("Skipping this profile fit and continuing.....")

            sky_bmodel = np.array(0.0)
            iterbsp = 0
            while (not sky_bmodel.any()) & (iterbsp <= 4):
                bsp_now = (1.2 ** iterbsp) * bsp
                fullbkpt = optimal_bkpts(bkpts_optimal, bsp_now, piximg, localmask, debug=(debug_bkpts & (iiter == niter)),
                                         skyimage=skyimage, min_spat=min_spat, max_spat=max_spat)
                # check to see if only a subset of the image is used.
                # if so truncate input pixels since this can result in singular matrices
                isub, = np.where(localmask.flatten())
                sortpix = (piximg.flat[isub]).argsort()
                obj_profiles_flat = obj_profiles.reshape(nspec * nspat, objwork)

                skymask = outmask & np.invert(edgmask)
                sky_bmodel, obj_bmodel, outmask_opt = skyoptimal(piximg.flat[isub], sciimg.flat[isub],
                                                                 (modelivar * skymask).flat[isub],
                                                                 obj_profiles_flat[isub, :], sortpix,
                                                                 spatial=spatial_img.flat[isub],
                                                                 fullbkpt=fullbkpt, sigrej=sigrej_eff, npoly=npoly)
                iterbsp = iterbsp + 1
                if (not sky_bmodel.any()) & (iterbsp <= 3):
                    msgs.warn('***************************************')
                    msgs.warn('WARNING: bspline sky-subtraction failed')
                    msgs.warn('Increasing bkpt spacing by 20%. Retry')
                    msgs.warn(
                        'Old bsp = {:5.2f}'.format(bsp_now) + '; New bsp = {:5.2f}'.format(1.2 ** (iterbsp) * bsp))
                    msgs.warn('***************************************')

            if sky_bmodel.any():
                skyimage.flat[isub] = sky_bmodel
                objimage.flat[isub] = obj_bmodel
                img_minsky.flat[isub] = sciimg.flat[isub] - sky_bmodel
                #var_no = np.abs(sky_bmodel - np.sqrt(2.0) * np.sqrt(rn2_img.flat[isub])) + rn2_img.flat[isub]
                igood1 = skymask.flat[isub]
                #  update the outmask for only those pixels that were fit. This prevents masking of slit edges in outmask
                outmask.flat[isub[igood1]] = outmask_opt[igood1]
                #from pypeit.ginga import show_image
                #show_image(sciimg*outmask)
                #debugger.set_trace()
                #  For weighted co-adds, the variance of the image is no longer equal to the image, and so the modelivar
                #  eqn. below is not valid. However, co-adds already have the model noise propagated correctly in sciivar,
                #  so no need to re-model the variance.
                if model_noise:
                    var = np.abs(sky_bmodel + obj_bmodel - np.sqrt(2.0) * np.sqrt(rn2_img.flat[isub])) + rn2_img.flat[isub]
                    var = var + adderr**2*(np.abs(sky_bmodel + obj_bmodel))**2
                    modelivar.flat[isub] = (var > 0.0) / (var + (var == 0.0))
                    #varnoobj.flat[isub] = var_no
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
                    #sigrej_eff = np.fmin(sigrej_eff, sigrej_ceil)
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

        # Now that the iterations of profile fitting and sky subtraction are completed,
        # loop over the objwork objects in this grouping and perform the final extractions.
        for ii in range(objwork):
            iobj = group[ii]
            msgs.info('Extracting obj # {:d}'.format(iobj + 1) + ' of {:d}'.format(nobj) +
                      ' with objid = {:d}'.format(sobjs[iobj].OBJID) + ' on slit # {:d}'.format(sobjs[iobj].slit_order) +
                      ' at x = {:5.2f}'.format(sobjs[iobj].SPAT_PIXPOS))
            this_profile = obj_profiles[:, :, ii]
            trace = np.outer(sobjs[iobj].TRACE_SPAT, np.ones(nspat))
            # Optimal
            objmask = ((spat_img >= (trace - 2.0 * box_rad)) & (spat_img <= (trace + 2.0 * box_rad)))
            extract.extract_optimal(sciimg, modelivar * thismask, (outmask & objmask), waveimg, skyimage, rn2_img, thismask, this_profile,
                            box_rad, sobjs[iobj])
            # Boxcar
            extract.extract_boxcar(sciimg, modelivar*thismask, (outmask & objmask),
                                           waveimg, skyimage, rn2_img, box_rad, sobjs[iobj])
            sobjs[iobj].min_spat = min_spat
            sobjs[iobj].max_spat = max_spat


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
            ginga.show_trace(viewer, ch, spec.TRACE_SPAT, spec.name, color=color)

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

    return (skyimage[thismask], objimage[thismask], modelivar[thismask], outmask[thismask])


def ech_local_skysub_extract(sciimg, sciivar, mask, tilts, waveimg, global_sky, rn2img, norders,
                             left, right, slitmask, sobjs, order_vec, spat_pix=None,
                             fit_fwhm=False, min_snr=2.0,bsp=0.6, extract_maskwidth=4.0,
                             trim_edg=(3,3), std=False, prof_nsigma=None, niter=4, box_rad_order=7,
                             sigrej=3.5, bkpts_optimal=True, sn_gauss=4.0, model_full_slit=False,
                             model_noise=True, debug_bkpts=False, show_profile=False,
                             show_resids=False, show_fwhm=False):
    """
    Perform local sky subtraction, profile fitting, and optimal extraction slit by slit

    Args:
        sciimg:
        sciivar:
        mask:
        tilts:
        waveimg:
        global_sky:
        rn2img:
        norders (:obj:`int`):
            Number of echelle orders
        left (`numpy.ndarray`_):
            Spatial-pixel coordinates for the left edges of each
            order.
        right (`numpy.ndarray`_):
            Spatial-pixel coordinates for the right edges of each
            order.
        slitmask (`numpy.ndarray`_):
            Image identifying the 0-indexed order associated with
            each pixel. Pixels with -1 are not associatead with any
            order.
        sobjs:
        order_vec:
        spat_pix:
        fit_fwhm:
        min_snr:
        bsp:
        extract_maskwidth:
        trim_edg:
        std:
        prof_nsigma:
        niter:
        box_rad_order (int???):
            Code assumes an np.ndarray even though the default value is int!!
        sigrej:
        bkpts_optimal:
        sn_gauss:
        model_full_slit:
        model_noise:
        debug_bkpts:
        show_profile:
        show_resids:
        show_fwhm:

    Returns:
        skymodel, objmodel, ivarmodel, outmask, sobjs

    """

    bitmask = maskimage.ImageBitMask()

    # Allocate the images that are needed
    # Initialize to mask in case no objects were found
    outmask = np.copy(mask)
    extractmask = (mask == 0)
    # TODO case of no objects found should be properly dealt with by local_skysub_extract
    # Initialize to zero in case no objects were found
    objmodel = np.zeros_like(sciimg)
    # Set initially to global sky in case no objects were found
    skymodel  = np.copy(global_sky)
    # Set initially to sciivar in case no obects were found.
    ivarmodel = np.copy(sciivar)
    sobjs = sobjs.copy()

    slit_vec = np.arange(norders)

    if (np.sum(sobjs.sign > 0) % norders) == 0:
        nobjs = int((np.sum(sobjs.sign > 0)/norders))
    else:
        msgs.error('Number of specobjs in sobjs is not an integer multiple of the number or ordres!')

    order_snr = np.zeros((norders, nobjs))
    uni_objid = np.unique(sobjs[sobjs.sign > 0].ECH_OBJID)
    for iord in range(norders):
        for iobj in range(nobjs):
            ind = (sobjs.ECH_ORDERINDX == iord) & (sobjs.ECH_OBJID == uni_objid[iobj])
            order_snr[iord,iobj] = sobjs[ind].ech_snr

    # Compute the average SNR and find the brightest object
    snr_bar = np.mean(order_snr,axis=0)
    srt_obj = snr_bar.argsort()[::-1]
    ibright = srt_obj[0] # index of the brightest object
    # Now extract the orders in descending order of S/N for the brightest object
    srt_order_snr = order_snr[:,ibright].argsort()[::-1]
    fwhm_here = np.zeros(norders)
    fwhm_was_fit = np.zeros(norders,dtype=bool)
    # Print out a status message
    str_out = ''
    for iord in srt_order_snr:
        str_out += '{:<8d}{:<8d}{:>10.2f}'.format(slit_vec[iord], order_vec[iord], order_snr[iord,ibright]) + msgs.newline()
    dash = '-'*27
    dash_big = '-'*40
    msgs.info(msgs.newline() + 'Reducing orders in order of S/N of brightest object:' + msgs.newline() + dash +
              msgs.newline() + '{:<8s}{:<8s}{:>10s}'.format('slit','order','S/N') + msgs.newline() + dash +
              msgs.newline() + str_out)
    # Loop over orders in order of S/N ratio (from highest to lowest) for the brightest object
    for iord in srt_order_snr:
        order = order_vec[iord]
        msgs.info("Local sky subtraction and extraction for slit/order: {:d}/{:d}".format(iord,order))
        other_orders = (fwhm_here > 0) & np.invert(fwhm_was_fit)
        other_fit    = (fwhm_here > 0) & fwhm_was_fit
        # Loop over objects in order of S/N ratio (from highest to lowest)
        for iobj in srt_obj:
            if (order_snr[iord, iobj] <= min_snr) & (np.sum(other_orders) >= 3):
                if iobj == ibright:
                    # If this is the brightest object then we extrapolate the FWHM from a fit
                    #fwhm_coeffs = np.polyfit(order_vec[other_orders], fwhm_here[other_orders], 1)
                    #fwhm_fit_eval = np.poly1d(fwhm_coeffs)
                    #fwhm_fit = fwhm_fit_eval(order_vec[iord])
                    fwhm_was_fit[iord] = True
                    # Either perform a linear fit to the FWHM or simply take the median
                    if fit_fwhm:
                        minx = 0.0
                        maxx = fwhm_here[other_orders].max()
                        # ToDO robust_poly_fit needs to return minv and maxv as outputs for the fits to be usable downstream
                        fit_mask, fwhm_coeffs = utils.robust_polyfit_djs(order_vec[other_orders], fwhm_here[other_orders],1,
                                                                        function='polynomial',maxiter=25,lower=2.0, upper=2.0,
                                                                        maxrej=1,sticky=False, minx=minx, maxx=maxx)
                        fwhm_this_ord = utils.func_val(fwhm_coeffs, order_vec[iord], 'polynomial', minx=minx, maxx=maxx)
                        fwhm_all = utils.func_val(fwhm_coeffs, order_vec, 'polynomial', minx=minx, maxx=maxx)
                        fwhm_str = 'linear fit'
                    else:
                        fit_mask = np.ones_like(order_vec[other_orders],dtype=bool)
                        fwhm_this_ord = np.median(fwhm_here[other_orders])
                        fwhm_all = np.full(norders,fwhm_this_ord)
                        fwhm_str = 'median '
                    indx = (sobjs.ECH_OBJID == uni_objid[iobj]) & (sobjs.ECH_ORDERINDX == iord)
                    for spec in sobjs[indx]:
                        spec.FWHM = fwhm_this_ord

                    str_out = ''
                    for slit_now, order_now, snr_now, fwhm_now in zip(slit_vec[other_orders], order_vec[other_orders],order_snr[other_orders,ibright], fwhm_here[other_orders]):
                        str_out += '{:<8d}{:<8d}{:>10.2f}{:>10.2f}'.format(slit_now, order_now, snr_now, fwhm_now) + msgs.newline()
                    msgs.info(msgs.newline() + 'Using' +  fwhm_str + ' for FWHM of object={:d}'.format(uni_objid[iobj]) +
                              ' on slit/order: {:d}/{:d}'.format(iord,order) + msgs.newline() + dash_big +
                              msgs.newline() + '{:<8s}{:<8s}{:>10s}{:>10s}'.format('slit', 'order','SNR','FWHM') +
                              msgs.newline() + dash_big +
                              msgs.newline() + str_out[:-8] +
                              fwhm_str.upper() +  ':{:<8d}{:<8d}{:>10.2f}{:>10.2f}'.format(iord, order, order_snr[iord,ibright], fwhm_this_ord) +
                              msgs.newline() + dash_big)
                    if show_fwhm:
                        plt.plot(order_vec[other_orders][fit_mask], fwhm_here[other_orders][fit_mask], marker='o', linestyle=' ',
                        color='k', mfc='k', markersize=4.0, label='orders informing fit')
                        if np.any(np.invert(fit_mask)):
                            plt.plot(order_vec[other_orders][np.invert(fit_mask)],
                                     fwhm_here[other_orders][np.invert(fit_mask)], marker='o', linestyle=' ',
                                     color='magenta', mfc='magenta', markersize=4.0, label='orders rejected by fit')
                        if np.any(other_fit):
                            plt.plot(order_vec[other_fit], fwhm_here[other_fit], marker='o', linestyle=' ',
                            color='lawngreen', mfc='lawngreen',markersize=4.0, label='fits to other low SNR orders')
                        plt.plot([order_vec[iord]], [fwhm_this_ord], marker='o', linestyle=' ',color='red', mfc='red', markersize=6.0,label='this order')
                        plt.plot(order_vec, fwhm_all, color='cornflowerblue', zorder=10, linewidth=2.0, label=fwhm_str)
                        plt.legend()
                        plt.show()
                else:
                    # If this is not the brightest object then assign it the FWHM of the brightest object
                    indx     = np.where((sobjs.ECH_OBJID == uni_objid[iobj]) & (sobjs.ECH_ORDERINDX == iord))[0][0]
                    indx_bri = np.where((sobjs.ECH_OBJID == uni_objid[ibright]) & (sobjs.ECH_ORDERINDX == iord))[0][0]
                    spec = sobjs[indx]
                    spec.FWHM = sobjs[indx_bri].FWHM

        thisobj = (sobjs.ECH_ORDERINDX == iord) # indices of objects for this slit
        thismask = (slitmask == iord) # pixels for this slit
        # True  = Good, False = Bad for inmask
        inmask = (mask == 0) & thismask
        # Local sky subtraction and extraction
        skymodel[thismask], objmodel[thismask], ivarmodel[thismask], extractmask[thismask] = local_skysub_extract(
            sciimg, sciivar, tilts, waveimg, global_sky,rn2img, thismask,
            left[:,iord], right[:,iord], sobjs[thisobj], spat_pix=spat_pix,
            inmask=inmask,std = std, bsp=bsp, extract_maskwidth=extract_maskwidth, trim_edg=trim_edg,
            prof_nsigma=prof_nsigma, niter=niter, box_rad=box_rad_order[iord], sigrej=sigrej, bkpts_optimal=bkpts_optimal,
            sn_gauss=sn_gauss, model_full_slit=model_full_slit, model_noise=model_noise, debug_bkpts=debug_bkpts,
            show_resids=show_resids, show_profile=show_profile)

        # update the FWHM fitting vector for the brighest object
        indx = (sobjs.ECH_OBJID == uni_objid[ibright]) & (sobjs.ECH_ORDERINDX == iord)
        fwhm_here[iord] = np.median(sobjs[indx].FWHMFIT)
        # Did the FWHM get updated by the profile fitting routine in local_skysub_extract? If so, include this value
        # for future fits
        if np.abs(fwhm_here[iord] - sobjs[indx].FWHM) >= 0.01:
            fwhm_was_fit[iord] = False

    # Set the bit for pixels which were masked by the extraction.
    # For extractmask, True = Good, False = Bad
    iextract = (mask == 0) & (extractmask == False)
    # Undefined inverse variances
    outmask[iextract] = bitmask.turn_on(outmask[iextract], 'EXTRACT')

    # Return
    return skymodel, objmodel, ivarmodel, outmask, sobjs


