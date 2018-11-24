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

try:
    from pypeit import ginga
except ImportError:
    pass


def trace_tilts_guess(arcimg, lines_spec, lines_spat, trace_int, thismask, inmask=None, tilts_guess=None, fwhm=4.0,
                      ncoeff=3, maxdev_fit=0.1,percentile_reject=0.10, max_badpix_frac=0.20,
                      maxerr=1.0, maxshift=3.0, maxshift0=3.0,nave=5, show_fits=False):

    nspec, nspat = arcimg.shape
    do_crude = True if tilts_guess is None else False
    nlines = len(lines_spec)
    if trace_int % 2 != 0:
        msgs.error('The trace_int parameter must be an even integer')
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
    trc_tilt_dict = dict(spat_min=spat_min, spat_max=spat_max, do_crude=do_crude, use_tilt=use_tilt,
                         tilts_sub_spec=tilts_sub_spec, tilts_sub_spat=tilts_sub_spat,
                         tilts_sub=tilts_sub, tilts_sub_fit=tilts_sub_fit, tilts_sub_err=tilts_sub_err,
                         tilts_sub_mask=tilts_sub_mask,
                         tilts_spec=tilts_spec, tilts_spat=tilts_spat,
                         tilts=tilts, tilts_fit=tilts_fit, tilts_err=tilts_err, tilts_mask=tilts_mask)

    return trc_tilt_dict
