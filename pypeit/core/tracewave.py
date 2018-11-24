""" Module for methods related to tracing arc/sky lines across a slit/order
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import inspect
import copy

import numpy as np

from scipy import interpolate

import matplotlib.pyplot as plt

from pypeit import msgs
from pypeit.core import arc
from pypeit import utils
from pypeit.core import parse
from pypeit.core import pca
from pypeit.core import qa
from pypeit.core import trace_slits
from pypeit.core import extract
from pypeit import debugger
from astropy.stats import sigma_clipped_stats

try:
    from pypeit import ginga
except ImportError:
    pass


def trace_tilts_guess(arcimg, lines_spec, lines_spat, slit_width, thismask, inmask=None, tilts_guess=None, fwhm=4.0,
                      ncoeff=3, maxdev_fit=0.1,percentile_reject=0.10, max_badpix_frac=0.20,
                      maxerr=1.0, maxshift=3.0, maxshift0=3.0,nave=5, show_fits=False):

    slit_widp2 = slit_width + 2
    slit_width_even = slit_widp2 if slit_widp2 % 2 == 0 else slit_widp2 + 1
    trace_int = slit_width_even//2

    nspec, nspat = arcimg.shape
    do_crude = True if tilts_guess is None else False
    nlines = len(lines_spec)

    nsub = 2 * trace_int + 1

    lines_spat_int = np.round(lines_spat).astype(int)

    spat_min = np.zeros(nlines, dtype=int)
    spat_max = np.zeros(nlines, dtype=int)

    if inmask is None:
        inmask = thismask

    tilts_sub = np.zeros((nsub, nlines))
    tilts_sub_fit = np.zeros((nsub, nlines))
    tilts_sub_err = np.zeros((nsub, nlines))
    tilts_sub_mask = np.zeros((nsub, nlines), dtype=bool)
    tilts_sub_spat = np.outer(np.arange(nsub), np.ones(nlines))
    tilts_sub_spec = np.outer(np.ones(nsub), lines_spec)

    tilts = np.zeros((nspat, nlines))
    tilts_fit = np.zeros((nspat, nlines))
    tilts_err = np.zeros((nspat, nlines))
    tilts_mask = np.zeros((nspat, nlines), dtype=bool)  # This is true if the pixel was in a region traced
    tilts_spat = np.outer(np.arange(nspat), np.ones(nlines))
    tilts_spec = np.outer(np.ones(nspat), lines_spec)

    # Transposed image and masks for traceing
    arcimg_trans = (arcimg * thismask).T
    inmask_trans = (inmask * thismask).T.astype(float)

    # 1) Trace the tilts from a guess. If no guess is provided from a previous iteration use trace_crude
    for iline in range(nlines):
        spat_min[iline] = lines_spat_int[iline] - trace_int
        spat_max[iline] = lines_spat_int[iline] + trace_int + 1
        sub_img = arcimg_trans[np.fmax(spat_min[iline], 0):np.fmin(spat_max[iline], nspat - 1), :]
        sub_inmask = inmask_trans[np.fmax(spat_min[iline], 0):np.fmin(spat_max[iline], nspat - 1), :]
        if do_crude:
            tilts_guess_now, err_now = trace_slits.trace_crude_init(
                sub_img, np.array([lines_spec[iline]]), (sub_img.shape[0] - 1) // 2, invvar=sub_inmask, radius=fwhm,
                nave=nave, maxshift0=maxshift0, maxshift=maxshift, maxerr=maxerr)
        else:
            tilts_guess_now = tilts_guess[:, iline]
        # Do iterative flux weighted tracing and polynomial fitting to refine these traces. This must also be done in a loop
        # since the sub image is different for every aperture, i.e. each aperature has its own image
        tilts_sub_fit_fw, tilts_sub_fw, tilts_sub_fw_err, tset_fw = extract.iter_tracefit(
            sub_img, tilts_guess_now, ncoeff, inmask=sub_inmask, fwhm=fwhm, maxdev=maxdev_fit, niter=6, idx=str(iline),
            show_fits=show_fits)
        tilts_sub_fit_gw, tilts_sub_gw, tilts_sub_gw_err, tset_gw = extract.iter_tracefit(
            sub_img, tilts_sub_fit_fw, ncoeff, inmask=sub_inmask, fwhm=fwhm, maxdev=maxdev_fit, niter=3, idx=str(iline),
            show_fits=show_fits)

        # Pack the results into arrays, accounting for possibly falling off the image
        # Deal with possibly falling off the chip
        if spat_min[iline] < 0:
            tilts_sub[-spat_min[iline]:, iline], tilts_sub_fit[-spat_min[iline]:,
                                                 iline] = tilts_sub_gw.flatten(), tilts_sub_fit_gw.flatten()
            tilts_sub_err[-spat_min[iline]:, iline] = tilts_sub_gw_err.flatten()
            tilts_sub_mask[-spat_min[iline]:, iline] = True
            tilts[np.fmax(spat_min[iline], 0):np.fmin(spat_max[iline], nspat - 1), iline] = tilts_sub[-spat_min[iline]:,
                                                                                            iline]
            tilts_fit[np.fmax(spat_min[iline], 0):np.fmin(spat_max[iline], nspat - 1), iline] = tilts_sub_fit[
                                                                                                -spat_min[iline]:,
                                                                                                iline]
            tilts_err[np.fmax(spat_min[iline], 0):np.fmin(spat_max[iline], nspat - 1), iline] = tilts_sub_err[
                                                                                                -spat_min[iline]:,
                                                                                                iline]
            tilts_mask[np.fmax(spat_min[iline], 0):np.fmin(spat_max[iline], nspat - 1), iline] = tilts_sub_mask[
                                                                                                 -spat_min[iline]:,
                                                                                                 iline]
        elif spat_max[iline] > (nspat - 1):
            tilts_sub[:-(spat_max[iline] - nspat + 1), iline], tilts_sub_fit[:-(spat_max[iline] - nspat + 1),
                                                               iline] = tilts_sub_gw.flatten(), tilts_sub_fit_gw.flatten()
            tilts_sub_err[:-(spat_max[iline] - nspat + 1), iline] = tilts_sub_gw_err.flatten()
            tilts_sub_mask[:-(spat_max[iline] - nspat + 1), iline] = True
            tilts[np.fmax(spat_min[iline], 0):np.fmin(spat_max[iline], nspat - 1), iline] = tilts_sub[:-(
                        spat_max[iline] - nspat + 1), iline]
            tilts_fit[np.fmax(spat_min[iline], 0):np.fmin(spat_max[iline], nspat - 1), iline] = tilts_sub_fit[:-(
                        spat_max[iline] - nspat + 1), iline]
            tilts_err[np.fmax(spat_min[iline], 0):np.fmin(spat_max[iline], nspat - 1), iline] = tilts_sub_err[:-(
                        spat_max[iline] - nspat + 1), iline]
            tilts_mask[np.fmax(spat_min[iline], 0):np.fmin(spat_max[iline], nspat - 1), iline] = tilts_sub_mask[:-(
                        spat_max[iline] - nspat + 1), iline]
        else:
            tilts_sub[:, iline], tilts_sub_fit[:, iline] = tilts_sub_gw.flatten(), tilts_sub_fit_gw.flatten()
            tilts_sub_err[:, iline] = tilts_sub_gw_err.flatten()
            tilts_sub_mask[:, iline] = True
            tilts[np.fmax(spat_min[iline], 0):np.fmin(spat_max[iline], nspat - 1), iline] = tilts_sub[:, iline]
            tilts_fit[np.fmax(spat_min[iline], 0):np.fmin(spat_max[iline], nspat - 1), iline] = tilts_sub_fit[:, iline]
            tilts_err[np.fmax(spat_min[iline], 0):np.fmin(spat_max[iline], nspat - 1), iline] = tilts_sub_err[:, iline]
            tilts_mask[np.fmax(spat_min[iline], 0):np.fmin(spat_max[iline], nspat - 1), iline] = tilts_sub_mask[:,
                                                                                                 iline]

    # Create the mask for the bad lines. Define the error on the bad tilt as being the
    bad_mask = (tilts_sub_err > 900) | (tilts_sub_mask == False)
    bad_pixel_count = np.sum(bad_mask, 0)
    dev_mean, dev_median, dev_sig = sigma_clipped_stats(np.abs(tilts_sub - tilts_sub_fit), mask=bad_mask, sigma=4.0,
                                                        axis=0)
    dev_mad = 1.4826*dev_median.data
    use_tilt = (dev_mad < np.quantile(dev_mad, 1.0 - percentile_reject)) & (bad_pixel_count < max_badpix_frac * nsub)
    msgs.info('Number of usable arc lines for tilts: {:d}/{:d}'.format(np.sum(use_tilt),nlines))

    # Tighten it up with Gaussian weighted centroiding
    trc_tilt_dict = dict(nspec = nspec, nspat = nspat, nsub = nsub, nlines = nlines, spat_min=spat_min, spat_max=spat_max, do_crude=do_crude, use_tilt=use_tilt,
                         tilts_sub_spec=tilts_sub_spec, tilts_sub_spat=tilts_sub_spat,
                         tilts_sub=tilts_sub, tilts_sub_fit=tilts_sub_fit, tilts_sub_err=tilts_sub_err,
                         tilts_sub_mask=tilts_sub_mask,
                         tilts_spec=tilts_spec, tilts_spat=tilts_spat,
                         tilts=tilts, tilts_fit=tilts_fit, tilts_err=tilts_err, tilts_mask=tilts_mask)

    return trc_tilt_dict


def coeff2tilts(coeff2, shape, func2D, max_tilt=1.2, min_tilt=-0.2):
    # Compute the tilts image
    polytilts = utils.polyval2d_general(coeff2, np.linspace(0.0, 1.0, shape[1]),
                                        np.linspace(0.0, 1.0, shape[0]),
                                        minx=0., maxx=1., miny=0., maxy=1., function=func2D)
    # y normalization and subtract
    ynorm = np.outer(np.linspace(0., 1., shape[0]), np.ones(shape[1]))
    polytilts = ynorm - polytilts / (shape[0] - 1)
    # JFH Added this to ensure that tilts are never crazy values due to extrapolation of fits which can break
    # wavelength solution fitting
    polytilts = np.fmax(np.fmin(polytilts, max_tilt), min_tilt)

    return polytilts


# TODO: Change yorder to "dispaxis_order"?
def fit_tilts(msarc, slit, all_tilts, order=2, yorder=4, func2D='legendre2d', maskval=-999999.9,
              setup=None, doqa=True, show_QA=False, out_dir=None):
    # ToDO please add some docs.
    """

    Args:
        msarc:
        slit:
        all_tilts:
        order:
        yorder:
        func2D:
        maskval:
        setup:
        doqa:
        show_QA:
        out_dir:

    Returns:

    """
    # Unpack
    xtilt, ytilt, mtilt, wtilt = all_tilts
    #
    # ToDO change variable labels spec_order, spat_order, get rid of this +1
    fitxy = [order + 1, yorder]

    # Fit the inverted model with a 2D polynomial
    msgs.info("Fitting tilts with a low order, 2D {:s}".format(func2D))

    inmask = xtilt != maskval
    fitmask, coeff2 = utils.robust_polyfit_djs(xtilt[inmask], mtilt[inmask] - ytilt[inmask], fitxy,
                                               x2=mtilt[inmask] / (msarc.shape[0] - 1.0),
                                               function='legendre2d', maxiter=30, maxrej=10, lower=3.0, upper=3.0,
                                               minx=0.0, maxx=1.0, minx2=0.0, maxx2=1.0,
                                               use_mad=True, sticky=True)
    plt.figure(figsize=(8, 20))
    deltay_fit = utils.func_val(coeff2, xtilt[inmask], 'legendre2d', x2=mtilt[inmask] / (msarc.shape[0] - 1.0),
                                minx=0.0, maxx=1.0, minx2=0.0, maxx2=1.0)
    plt.plot(xtilt[inmask], ytilt[inmask] + deltay_fit, 'go', mfc='g', markersize=5.0, markeredgewidth=1.0,
             label='2D Fit')
    plt.plot(xtilt[inmask][fitmask], mtilt[inmask][fitmask], 'ko', mfc='None', markersize=7.0, markeredgewidth=1.0,
             label='Good Points')
    plt.plot(xtilt[inmask][~fitmask], mtilt[inmask][~fitmask], 'ro', mfc='None', markersize=7.0, markeredgewidth=1.0,
             label='Rejected Points')
    plt.legend()
    plt.show()
    res2 = (mtilt[inmask][fitmask] - ytilt[inmask][fitmask]) - deltay_fit[fitmask]
    msgs.info("RMS (pixels): {}".format(np.std(res2)))

    #    plt.figure(figsize=(8,20))
    #    deltay_fit = utils.func_val(coeff2, xtilt[inmask], func2D, x2=mtilt[inmask]/(msarc.shape[0]-1.0), minx=0.0, maxx=1.0, minx2=0.0, maxx2=1.0)
    #    plt.plot(xtilt[inmask], ytilt[inmask] + deltay_fit, 'go', mfc='g',markersize=7.0, markeredgewidth=1.0, label = '2D Fit')
    #    plt.plot(xtilt[inmask][fitmask], mtilt[inmask][fitmask], 'ko', mfc='None',markersize=5.0, markeredgewidth=1.0, label = 'Good Points')
    #    plt.plot(xtilt[inmask][~fitmask], mtilt[inmask][~fitmask], 'ro', mfc='None',markersize=5.0, markeredgewidth=1.0, label = 'Rejected Points')
    #    plt.legend()
    #    plt.show()

    # wgd = np.where(xtilt != maskval)
    # Invert

    """
        inmask = xtilt != maskval
        coeff2 = utils.polyfit2d_general(xtilt[inmask], mtilt[inmask]/(msarc.shape[0]-1),
                                           mtilt[inmask]-ytilt[inmask], fitxy,
                                           minx=0., maxx=1., miny=0., maxy=1., function=func2D)

        # TODO -- Add a rejection iteration (or two)

        # Residuals
        xv2 = utils.scale_minmax(xtilt[inmask], minx=0., maxx=1)
        yv2 = utils.scale_minmax(mtilt[inmask]/(msarc.shape[0]-1), minx=0., maxx=1)
        deltay_fit = np.polynomial.legendre.legval2d(xv2, yv2, coeff2)
        res2 = (mtilt[inmask]-ytilt[inmask]) - deltay_fit
        msgs.info("RMS (pixels): {}".format(np.std(res2)))

        plt.figure(figsize=(8, 20))
        plt.plot(xtilt[inmask], mtilt[inmask] - deltay_fit, 'go', mfc='g', markersize=7.0, markeredgewidth=1.0,
                 label='2D Fit')
        plt.plot(xtilt[inmask], ytilt[inmask], 'ko', mfc='None', markersize=5.0, markeredgewidth=1.0,label='Good Points')
        plt.legend()
        plt.show()
    """

    # QA
    if doqa:
        plot_tiltres(setup, mtilt[inmask], ytilt[inmask], deltay_fit, slit=slit, show_QA=show_QA, out_dir=out_dir)

    # Compute the tilts image
    polytilts = coeff2tilts(coeff2, msarc.shape, func2D)

    # Return
    outpar = None
    return polytilts, coeff2, outpar
