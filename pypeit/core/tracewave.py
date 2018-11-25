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


def tilts_find_lines(arc_spec, slit_cen, tracethresh=10.0, sigdetect=5.0, npix_neigh=5.0,
                    only_these_lines=None, fwhm=4.0, nonlinear_counts=1e10, fit_frac_fwhm=1.25, cont_frac_fwhm=1.0,
                    max_frac_fwhm=2.0, cont_samp=30, niter_cont=3, debug = False):


    nspec = arc_spec.size
    spec_vec = np.arange(nspec)
    # Find peaks with a liberal threshold of sigdetect = 5.0
    tampl_tot, tampl_cont_tot, tcent_tot, twid_tot, _, wgood, _, nsig_tot = arc.detect_lines(
        arc_spec, sigdetect=sigdetect, fwhm=fwhm, fit_frac_fwhm=fit_frac_fwhm, cont_frac_fwhm=cont_frac_fwhm,
        max_frac_fwhm=max_frac_fwhm, cont_samp=cont_samp, niter_cont=niter_cont, nonlinear_counts=nonlinear_counts,
        debug=debug)
    # Good lines
    arcdet = tcent_tot[wgood]
    nsig = nsig_tot[wgood]

    # Determine the best lines to use to trace the tilts
    aduse = np.zeros(arcdet.size, dtype=np.bool)  # Which lines should be used to trace the tilts
    w = np.where(nsig >= tracethresh)
    aduse[w] = 1
    # Remove lines that are within npix_neigh pixels
    nuse = np.sum(aduse)
    detuse = arcdet[aduse]
    idxuse = np.arange(arcdet.size)[aduse]
    olduse = aduse.copy()
    for s in range(nuse):
        w = np.where((np.abs(arcdet - detuse[s]) <= npix_neigh) & (np.abs(arcdet - detuse[s]) >= 1.0))[0]
        for u in range(w.size):
            if nsig[w[u]] > nsig[olduse][s]:
                aduse[idxuse[s]] = False
                break

    # Restricted to ID lines? [introduced to avoid LRIS ghosts]
    if only_these_lines is not None:
        ids_pix = np.array(only_these_lines)
        idxuse = np.arange(arcdet.size)[aduse]
        for s in idxuse:
            if np.min(np.abs(arcdet[s] - ids_pix)) > 2.0:
                msgs.info("Ignoring line at spectral position={:6.1f} which was not identified".format(arcdet[s]))
                aduse[s] = False

    # Final spectral positions of arc lines we will trace
    lines_spec = arcdet[aduse]
    nlines = len(lines_spec)
    if nlines == 0:
        msgs.warn('No arc lines were deemed usable on this slit. Cannot compute tilts. Try lowering tracethresh.')
        return None
    else:
        msgs.info('Modelling arc line tilts with {:d} arc lines'.format(nlines))

    # Spatial position of line, i.e. the central trace interpolated onto the spectral pixel of the line
    lines_spat = np.interp(lines_spec, spec_vec, slit_cen)

    return lines_spec, lines_spat


def trace_tilts(arcimg, lines_spec, lines_spat, slit_width, thismask, inmask=None, tilts_guess=None, fwhm=4.0,
                ncoeff=3, maxdev_fit=0.1,percentile_reject=0.10, max_badpix_frac=0.20,
                maxerr=1.0, maxshift=3.0, maxshift0=3.0,nave=5, show_fits=False, debug = False):

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
    thismask_trans = thismask.T

    # 1) Trace the tilts from a guess. If no guess is provided from a previous iteration use trace_crude
    for iline in range(nlines):
        spat_min[iline] = lines_spat_int[iline] - trace_int
        spat_max[iline] = lines_spat_int[iline] + trace_int + 1
        min_spat = np.fmax(spat_min[iline], 0)
        max_spat = np.fmin(spat_max[iline], nspat - 1)
        sub_img = arcimg_trans[min_spat:max_spat, :]
        sub_inmask = inmask_trans[min_spat:max_spat,:]
        sub_thismask = thismask_trans[min_spat:max_spat,:]
        if do_crude: # First time tracing, do a trace crude
            tilts_guess_now, err_now = trace_slits.trace_crude_init(
                sub_img, np.array([lines_spec[iline]]), (sub_img.shape[0] - 1) // 2, invvar=sub_inmask, radius=fwhm,
                nave=nave, maxshift0=maxshift0, maxshift=maxshift, maxerr=maxerr)
        else:  # A guess was provided, use that as the crutch, but determine if it is a full trace or a sub-trace
            if tilts_guess.shape[0] == nspat:
                # This is full image size tilt trace, sub-window it
                tilts_guess_now = tilts_guess[min_spat:max_spat, iline]
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
        # Boxcar extract the thismask to have a mask indicating whether a tilt is defined along the spatial direction
        tilts_sub_mask_box = (extract.extract_boxcar(sub_thismask, tilts_sub_fit_gw, fwhm/2.0) > 0.99*fwhm)
        # Pack the results into arrays, accounting for possibly falling off the image
        # Deal with possibly falling off the chip
        if spat_min[iline] < 0:
            tilts_sub[     -spat_min[iline]:,iline] = tilts_sub_gw.flatten()
            tilts_sub_fit[ -spat_min[iline]:,iline] = tilts_sub_fit_gw.flatten()
            tilts_sub_err[ -spat_min[iline]:,iline] = tilts_sub_gw_err.flatten()
            tilts_sub_mask[-spat_min[iline]:,iline] = tilts_sub_mask_box
            tilts[     min_spat:max_spat,iline] = tilts_sub[     -spat_min[iline]:,iline]
            tilts_fit[ min_spat:max_spat,iline] = tilts_sub_fit[ -spat_min[iline]:,iline]
            tilts_err[ min_spat:max_spat,iline] = tilts_sub_err[ -spat_min[iline]:,iline]
            tilts_mask[min_spat:max_spat,iline] = tilts_sub_mask[-spat_min[iline]:,iline]
        elif spat_max[iline] > (nspat - 1):
            tilts_sub[     :-(spat_max[iline] - nspat + 1),iline] = tilts_sub_gw.flatten()
            tilts_sub_fit[ :-(spat_max[iline] - nspat + 1),iline] = tilts_sub_fit_gw.flatten()
            tilts_sub_err[ :-(spat_max[iline] - nspat + 1),iline] = tilts_sub_gw_err.flatten()
            tilts_sub_mask[:-(spat_max[iline] - nspat + 1),iline] = tilts_sub_mask_box
            tilts[     min_spat:max_spat,iline] = tilts_sub[     :-(spat_max[iline] - nspat + 1),iline]
            tilts_fit[ min_spat:max_spat,iline] = tilts_sub_fit[ :-(spat_max[iline] - nspat + 1),iline]
            tilts_err[ min_spat:max_spat,iline] = tilts_sub_err[ :-(spat_max[iline] - nspat + 1),iline]
            tilts_mask[min_spat:max_spat,iline] = tilts_sub_mask[:-(spat_max[iline] - nspat + 1),iline]
        else:
            tilts_sub[     :,iline] = tilts_sub_gw.flatten()
            tilts_sub_fit[ :,iline] = tilts_sub_fit_gw.flatten()
            tilts_sub_err[ :,iline] = tilts_sub_gw_err.flatten()
            tilts_sub_mask[:,iline] = tilts_sub_mask_box
            tilts[     min_spat:max_spat,iline] = tilts_sub[     :,iline]
            tilts_fit[ min_spat:max_spat,iline] = tilts_sub_fit[ :,iline]
            tilts_err[ min_spat:max_spat,iline] = tilts_sub_err[ :,iline]
            tilts_mask[min_spat:max_spat,iline] = tilts_sub_mask[:,iline]

    # Create the mask for the bad lines. Define the error on the bad tilt as being the
    bad_mask = (tilts_sub_err > 900) | (tilts_sub_mask == False)
    bad_pixel_count = np.sum(bad_mask, 0)
    dev_mean, dev_median, dev_sig = sigma_clipped_stats(np.abs(tilts_sub - tilts_sub_fit), mask=bad_mask, sigma=4.0,
                                                        axis=0)
    dev_mad = 1.4826*dev_median.data
    use_tilt = (dev_mad < np.quantile(dev_mad, 1.0 - percentile_reject)) & (bad_pixel_count < max_badpix_frac * nsub)
    nuse = np.sum(use_tilt)
    msgs.info('Number of usable arc lines for tilts: {:d}/{:d}'.format(nuse,nlines))


    # Tighten it up with Gaussian weighted centroiding
    trc_tilt_dict = dict(nspec = nspec, nspat = nspat, nsub = nsub, nlines = nlines, spat_min=spat_min, spat_max=spat_max, do_crude=do_crude, use_tilt=use_tilt,
                         tilts_sub_spec=tilts_sub_spec, tilts_sub_spat=tilts_sub_spat,
                         tilts_sub=tilts_sub, tilts_sub_fit=tilts_sub_fit, tilts_sub_err=tilts_sub_err,
                         tilts_sub_mask=tilts_sub_mask,
                         tilts_spec=tilts_spec, tilts_spat=tilts_spat,
                         tilts=tilts, tilts_fit=tilts_fit, tilts_err=tilts_err, tilts_mask=tilts_mask)

    return trc_tilt_dict



# TODO: Change yorder to "dispaxis_order"?
def fit_tilts(trc_tilt_dict, spat_order=3, spec_order=4, func2D='legendre2d', maskval=-999999.9,
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

    nspec = trc_tilt_dict['nspec']
    nspat = trc_tilt_dict['nspat']
    tilts = trc_filt_dict['tilts']
    tilts_spat = trc_tilt_dict['tilts_spat']
    tilts_spec = trc_tilt_dict['tilts_spec']
    tilts_mask = trc_tilt_dict['tilts_spec']
    # Unpack
    xtilt, ytilt, mtilt, wtilt = all_tilts
    #
    # ToDO change variable labels spec_order, spat_order, get rid of this +1
    fitxy = [spat_order, spec_order]

    # Fit the inverted model with a 2D polynomial
    msgs.info("Fitting tilts with a low order, 2D {:s}".format(func2D))

    #                                               spatial position
    fitmask, coeff2 = utils.robust_polyfit_djs(tilts_spat[tilts_mask], tilts[tilts_mask] - tilts_spec[inmask], fitxy,x2=mtilt[inmask],
                                               function='legendre2d', maxiter=30, maxrej=10, lower=3.0, upper=3.0,
                                               minx=0.0, maxx=float(nspat-1), minx2=0.0, maxx2=float(nspec-1),
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


    # QA
    if doqa:
        plot_tiltres(setup, mtilt[inmask], ytilt[inmask], deltay_fit, slit=slit, show_QA=show_QA, out_dir=out_dir)

    # Compute the tilts image
    polytilts = coeff2tilts(coeff2, msarc.shape, func2D)

    # Return
    outpar = None
    return polytilts, coeff2, outpar



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
def fit_tilts_old(msarc, slit, all_tilts, order=2, yorder=4, func2D='legendre2d', maskval=-999999.9,
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
