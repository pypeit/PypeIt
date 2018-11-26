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


def tilts_find_lines(arc_spec, slit_cen, tracethresh=10.0, sigdetect=5.0, nfwhm_neigh=5.0,
                    only_these_lines=None, fwhm=4.0, nonlinear_counts=1e10, fit_frac_fwhm=1.25, cont_frac_fwhm=1.0,
                    max_frac_fwhm=2.0, cont_samp=30, niter_cont=3, debug_lines=False, debug_peaks=False):


    nspec = arc_spec.size
    spec_vec = np.arange(nspec)
    # Find peaks with a liberal threshold of sigdetect = 5.0
    tampl_tot, tampl_cont_tot, tcent_tot, twid_tot, _, wgood, arc_cont_sub, nsig_tot = arc.detect_lines(
        arc_spec, sigdetect=sigdetect, fwhm=fwhm, fit_frac_fwhm=fit_frac_fwhm, cont_frac_fwhm=cont_frac_fwhm,
        max_frac_fwhm=max_frac_fwhm, cont_samp=cont_samp, niter_cont=niter_cont, nonlinear_counts=nonlinear_counts,
        debug=debug_peaks)
    # Good lines
    arcdet = tcent_tot[wgood]
    arc_ampl = tampl_cont_tot[wgood]
    nsig = nsig_tot[wgood]

    npix_neigh = nfwhm_neigh*fwhm
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


    if debug_lines:
        xrng = np.arange(nspec)
        plt.figure(figsize=(14, 6))
        plt.plot(xrng, arc_cont_sub, color='black', drawstyle='steps-mid', lw=3, label='arc', linewidth=1.0)
        plt.plot(arcdet[~aduse], arc_ampl[~aduse], 'r+', markersize=6.0, label='bad for tilts')
        plt.plot(arcdet[aduse], arc_ampl[aduse], 'g+', markersize=6.0, label='good for tilts')
        if nonlinear_counts < 1e9:
            plt.hlines(nonlinear_counts, xrng.min(), xrng.max(), color='orange', linestyle='--', linewidth=2.0,
            label='nonlinear', zorder=10)
        plt.title('Good Lines = {:d}'.format(np.sum(aduse)) + ',  Bad Lines = {:d}'.format(np.sum(~aduse)))
        plt.ylim(arc_cont_sub.min(), 1.5 * arc_cont_sub.max())
        plt.legend()
        plt.show()

    # Spatial position of line, i.e. the central trace interpolated onto the spectral pixel of the line
    lines_spat = np.interp(lines_spec, spec_vec, slit_cen)

    return lines_spec, lines_spat


def trace_tilts_work(arcimg, lines_spec, lines_spat, thismask, inmask=None, tilts_guess=None, fwhm=4.0,
                     spat_order=3, maxdev_tracefit=1.0,sigrej_trace=3.0, max_badpix_frac=0.20,
                     tcrude_maxerr=1.0, tcrude_maxshift=3.0, tcrude_maxshift0=3.0,tcrude_nave=5, show_fits=False, debug = False):

    # Smash the
    slit_width = (np.sum(thismask, axis=1)).max()
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
                nave=tcrude_nave, maxshift0=tcrude_maxshift0, maxshift=tcrude_maxshift, maxerr=tcrude_maxerr)
        else:  # A guess was provided, use that as the crutch, but determine if it is a full trace or a sub-trace
            if tilts_guess.shape[0] == nspat:
                # This is full image size tilt trace, sub-window it
                tilts_guess_now = tilts_guess[min_spat:max_spat, iline]
            else:
                # If it is a sub-trace, deal with falling off the image
                if spat_min[iline] < 0:
                    tilts_guess_now = tilts_guess[-spat_min[iline]:,iline]
                elif spat_max[iline] > (nspat-1):
                    tilts_guess_now = tilts_guess[:-(spat_max[iline] - nspat + 1),iline]
                else:
                    tilts_guess_now = tilts_guess[:, iline]

        # Do iterative flux weighted tracing and polynomial fitting to refine these traces. This must also be done in a loop
        # since the sub image is different for every aperture, i.e. each aperature has its own image
        tilts_sub_fit_fw, tilts_sub_fw, tilts_sub_fw_err, tset_fw = extract.iter_tracefit(
            sub_img, tilts_guess_now, spat_order, inmask=sub_inmask, fwhm=fwhm, maxdev=maxdev_tracefit, niter=6, idx=str(iline),
            show_fits=show_fits)
        tilts_sub_fit_gw, tilts_sub_gw, tilts_sub_gw_err, tset_gw = extract.iter_tracefit(
            sub_img, tilts_sub_fit_fw, spat_order, inmask=sub_inmask, fwhm=fwhm, maxdev=maxdev_tracefit, niter=3, idx=str(iline),
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
    good_line = np.any(bad_mask == False,axis=0) # Is it masked everywhere?
    # Median absolute deviation for each line quantifies the goodnes of tracing
    dev_mad = 1.4826*dev_median.data
    # Now reject outliers from this distribution
    dev_mad_dist_median = np.median(dev_mad[good_line])
    dev_mad_dist_mad = 1.4826*np.median(np.abs(dev_mad[good_line] - dev_mad_dist_median)) # i.e. this is like the sigma
    # Reject lines that are sigrej trace outliers
    mad_keep = np.abs((dev_mad - dev_mad_dist_median)/dev_mad_dist_mad) < sigrej_trace

    use_tilt = (mad_keep) & (bad_pixel_count < max_badpix_frac * nsub) & good_line
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


def trace_tilts(arcimg, lines_spec, lines_spat, thismask, inmask=None,fwhm=4.0,spat_order=5, maxdev_tracefit=1.0,
                sigrej_trace=3.0, max_badpix_frac=0.20, tcrude_nave = 5,
                npca = 1, coeff_npoly_pca = 1, sigrej_pca = 2.0,debug_pca = True, show_tracefits=False):

    """
    Use a PCA model to determine the best object (or slit edge) traces for echelle spectrographs.

    Parameters
    ----------
    arcimg:  ndarray, float (nspec, nspat)
       Image of arc or sky that will be used for tracing tilts.
    lines_spec: ndarray, float (nlines,)
       Array containing arc line centroids along the center of the slit for each arc line that will be traced. This is
       in pixels in image coordinates.
    lines_spat: ndarray, float (nlines,)
       Array contianing the spatial position of the center of the slit along which the arc was extracted. This is is in
       pixels in image coordinates.
    thismask: ndarray, boolean (nspec, nsapt)
        Boolean mask image specifying the pixels which lie on the slit/order to search for objects on.
        The convention is: True = on the slit/order, False  = off the slit/order. This must be the same size as the arcimg.
    Optional Parameters
    -------------------
    inmask: float ndarray, default = None
        Input mask image.
    fwhm: float
       Expected FWHM of the arc lines.
    spat_order: int, default = None
       Order of the legendre polynomial that will be fit to the tilts.
    maxdev_tracefit: float, default = 1.0
       Maximum absolute deviation for the arc tilt fits during iterative trace fitting (flux weighted, then gaussian weighted)
    sigrej_trace: float, default =  3.0
       From each line we compute a median absolute deviation of the trace from the polynomial fit. We then
       analyze the distribution of maximxum absolute deviations (MADs) for all the lines, and reject sigrej_trace outliers
       from that distribution.
    max_badpix_frac: float, default = 0.20
       Maximum fraction of total pixels that can be masked by the trace_gweight algorithm
       (because the residuals are too large) to still be usable for tilt fitting.
    tcrude_nave: int, default = 5
       Trace crude is used to determine the initial arc line tilts, which are then iteratively fit. Trace crude
       can optionally boxcar smooth the image (along the spatial direction of the image, i.e. roughly along the arc line tilts)
       to improve the tracing.
    npca: int, default = 1
       Tilts are initially traced and then a PCA is performed. The PCA is used to determine better crutches for a second
       round of improved tilt tracing. This parameter is the order of that PCA and determined how much the tilts behavior
       is being compressed. npca = 0 would be just using the mean tilt. This PCA is only an intermediate step to
       improve the crutches and is an attempt to make the tilt tracing that goes into the final fit more robust.
    coeff_npoly_pca: int, default = 1
       Order of polynomial fits used for PCA coefficients fitting for the PCA described above.
    sigrej_pca: float, default = 2.0
       Significance threhsold for rejection of outliers from fits to PCA coefficients for the PCA described above.
    show_tracefits: bool, default = False
       If true the fits will be shown to each arc line trace by iter_fitting.py


    Returns:
    --------
    """

    tcrude_maxerr = fwhm/4.0
    tcrude_maxshift = 3.0*fwhm/4.0
    tcrude_maxshift0 = fwhm

    trace_dict0 = trace_tilts_work(
        arcimg, lines_spec, lines_spat, thismask, inmask=inmask,tilts_guess=None, fwhm=fwhm, spat_order=spat_order,
        maxdev_tracefit=maxdev_tracefit,sigrej_trace=sigrej_trace, max_badpix_frac=max_badpix_frac,
        tcrude_maxerr=tcrude_maxerr, tcrude_maxshift=tcrude_maxshift,tcrude_maxshift0=tcrude_maxshift0,
        tcrude_nave=tcrude_nave, show_fits=show_tracefits)

    # TODO THE PCA may not be necessary. It appears to improve the results though for some instruments where the
    # tracing is problematic. We could consider making this optional to speed things up.
    debug_pca_fit = False
    if debug_pca_fit:
        # !!!! FOR TESTING ONLY!!!!  Evaluate the model fit to the tilts for all of our lines
        tilt_fit_dict0 = fit_tilts(trc_tilt_dict0, spat_order=spat_order, spec_order=4, debug=True,
                                   maxdev=1.0, sigrej=3.0,doqa=True, setup='test', slit=0, show_QA=False, out_dir='./')

    # Do a PCA fit, which rejects some outliers
    iuse = trc_tilt_dict0['use_tilt']
    nuse = np.sum(iuse)
    msgs.info('PCA modeling {:d} good tilts'.format(nuse))
    pca_fit, poly_fit_dict, pca_mean, pca_vectors = extract.pca_trace(
        trc_tilt_dict0['tilts_sub_fit'], predict=np.invert(iuse), npca=npca, coeff_npoly=coeff_npoly_pca,
        lower=sigrej_pca, upper=sigrej_pca, order_vec=lines_spec, xinit_mean=lines_spec,
        minv=0.0, maxv=float(trc_tilt_dict0['nsub'] - 1), debug=debug_pca)

    # Now trace again with the PCA predictions as the starting crutches
    trace_dict1 = trace_tilts_work(arcimg, lines_spec, lines_spat, thismask, inmask=inmask,tilts_guess=pca_fit,
                                      fwhm=fwhm, spat_order=spat_order, maxdev_tracefit=maxdev_tracefit,sigrej_trace=sigrej_trace,
                                      max_badpix_frac=max_badpix_frac,show_fits=show_tracefits)

    return trace_dict1


def fit_tilts(trc_tilt_dict, spat_order=3, spec_order=4, maxdev = 1.0, sigrej = 3.0, func2d='legendre2d', doqa=True, setup = 'test',
              slit = '0', show_QA=False, out_dir=None, debug=True):
    """

    Parameters
    ----------
    trc_tilt_dict: dict
        Diciontary containing tilt info

    Optional Parameters
    -------------------
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
    use_tilt = trc_tilt_dict['use_tilt']                 # mask for good/bad tilts, based on aggregate fit, frac good pixels
    tilts = trc_tilt_dict['tilts'][:,use_tilt]           # gaussian weighted centroid
    tilts_fit = trc_tilt_dict['tilts_fit'][:,use_tilt]   # legendre polynomial fit
    tilts_err = trc_tilt_dict['tilts_err'][:,use_tilt]   # gaussian weighted centroidding error
    tilts_spat = trc_tilt_dict['tilts_spat'][:,use_tilt] # spatial pixel position
    tilts_spec = trc_tilt_dict['tilts_spec'][:,use_tilt] # line_spec spectral pixel position
    tilts_mask = trc_tilt_dict['tilts_mask'][:,use_tilt] # Reflects if trace is on the slit

    # Let's just let the code do the rejection
    #delta_tilt  = np.abs(tilts - tilts_fit)

    # Do one last round of rejection here at the pixel level, i.e. we already rejected lines before
    #tot_mask = tilts_mask & (delta_tilt < maxdev) & (tilts_err < 900)
    tot_mask = tilts_mask & (tilts_err < 900)
    fitxy = [spat_order, spec_order]

    # Fit the inverted model with a 2D polynomial
    msgs.info("Fitting tilts with a low order, 2D {:s}".format(func2d))

    # This is a 2-d polynomila, i.e. z(x,y), or following how robust_polyfit works y(x,x2)
    # x = spatial position along image, i.e. np.arange(nspat) for each line
    # y = spectral pixel where arc line was detected, i.e. line_spec replicated everywhere for each line
    # z = (spectral_line_trace - spectral_pos on extracted arc), i.e. tilt - y

    # We fit the function in this way i.e. with independent variable y being the actual tilt, becuase we are trying to
    # determine the mapping which brings the spectral line back to the same spectral_pos on extracted arc for which
    # we determined the wavelength solution

    fitmask, coeff2 = utils.robust_polyfit_djs(tilts_spat[tot_mask], (tilts_fit[tot_mask] - tilts_spec[tot_mask]),
                                               fitxy,x2=tilts_fit[tot_mask],
                                               function='legendre2d', maxiter=100, lower=sigrej, upper=sigrej,
                                               maxdev=maxdev,
                                               minx=0.0, maxx=float(nspat-1), minx2=0.0, maxx2=float(nspec-1),
                                               use_mad=True, sticky=False)
    delta_spec_fit = utils.func_val(coeff2, tilts_spat[tot_mask], func2d, x2=tilts_fit[tot_mask],
                                    minx=0.0, maxx=float(nspat-1), minx2=0.0, maxx2=float(nspec-1))
    res2 = (tilts_fit[tot_mask][fitmask] - tilts_spec[tot_mask][fitmask]) - delta_spec_fit[fitmask]
    msgs.info("RMS (pixels): {}".format(np.std(res2)))


    if debug:
        plt.figure(figsize=(8, 20))
        plt.plot(tilts_spat[tot_mask],tilts_spec[tot_mask] + delta_spec_fit, 'go', mfc='g', markersize=5.0,
                 markeredgewidth=1.0,label='2D Fit')
        plt.plot(tilts_spat[tot_mask][fitmask], tilts_fit[tot_mask][fitmask], 'ko', mfc='None', markersize=7.0,
                 markeredgewidth=1.0,label='Good Points')
        plt.plot(tilts_spat[tot_mask][~fitmask], tilts_fit[tot_mask][~fitmask], 'ro', mfc='None', markersize=7.0,
                 markeredgewidth=1.0,
        label='Rejected Points')
        plt.legend()
        plt.show()

    # QA
    if doqa:
        plot_tiltres(setup, tilts_fit[tot_mask], tilts_spec[tot_mask], delta_spec_fit, slit=slit, show_QA=show_QA, out_dir=out_dir)

    tilt_fit_dict = dict(nspec = nspec, nspat = nspat, ngood_lines=np.sum(use_tilt), npix_fit = np.sum(tot_mask),
                         npix_rej = np.sum(fitmask == False), coeff2=coeff2, spec_order = spec_order, spat_order = spat_order,
                         minx = 0.0, maxx = float(nspat-1), minx2 = 0.0, maxx2 = float(nspec-1), func =func2d)

    return tilt_fit_dict



def fit2piximg(tilt_fit_dict):
    """

    Parameters
    ----------
    tilt_fit_dict: dict
        Tilt fit dictioary produced by fit_tilts

    Returns
    -------
    piximg: ndarray, float
       Image indicating how spectral pixel locations move across the image. This output is used in the pipeline.
    """

    # Compute the tilts image
    nspec = tilt_fit_dict['nspec']
    nspat = tilt_fit_dict['nspat']

    spat_vec = np.arange(tilt_fit_dict['nspat'])
    tilt_vec = np.arange(tilt_fit_dict['nspec'])

    tilt_min_spec_fit = utils.polyval2d_general(
        tilt_fit_dict['coeff2'], spat_vec, tilt_vec, minx=tilt_fit_dict['minx'], maxx=tilt_fit_dict['maxx'],
        miny=tilt_fit_dict['minx2'], maxy=tilt_fit_dict['maxx2'],
        function=tilt_fit_dict['func'])
    # y normalization and subtract
    spec_img = np.outer(np.arange(nspec), np.ones(nspat))
    piximg = spec_img - tilt_min_spec_fit
    # Added this to ensure that tilts are never crazy values due to extrapolation of fits which can break
    # wavelength solution fitting
    piximg = np.fmax(np.fmin(piximg, nspec),-1.0)

    return piximg


def fit2tilts(tilt_fit_dict, spat_vec, tilt_vec):
    """

    Parameters
    ----------
    tilt_fit_dict: dict
        Tilt fit dictioary produced by fit_tilts

    spat_vec: ndarray, float, default = None
        Spatial positions where tilt fit is desired
    tilt_vec: ndarray, float, default = None
        Spectral positions at which the tilt fit is desired. Note that the fit was actually done with the tilt
        as the y-independent variable, see above.

    Returns
    -------
    tilt_min_spec_fit: ndarray, float
       The actualy thing that was fit in computing the tilts, which is the trace of the spectra lines minuse the
       spectral pixel location of the central trace. See above. This output is only really used for debugging and
       internal functionality.

    """

    # Compute the tilts image
    nspec = tilt_fit_dict['nspec']
    nspat = tilt_fit_dict['nspat']

    tilt_min_spec_fit = utils.polyval2d_general(
        tilt_fit_dict['coeff2'], spat_vec, tilt_vec,
        minx=tilt_fit_dict['minx'], maxx=tilt_fit_dict['maxx'],
        miny=tilt_fit_dict['minx2'], maxy=tilt_fit_dict['maxx2'],
        function=tilt_fit_dict['func'])
    return tilt_min_spec_fit





def plot_tiltres(setup, mtilt, ytilt, yfit, slit=None, outfile=None, show_QA=False, out_dir=None):
    """ Generate a QA plot of the residuals for the fit to the tilts
    One slit at a time

    Parameters
    ----------
    """

    plt.rcdefaults()
    plt.rcParams['font.family']= 'times new roman'

    # Outfil
    method = inspect.stack()[0][3]
    if (outfile is None) and (not show_QA):
        outfile = qa.set_qa_filename(setup, method, slit=slit, out_dir=out_dir)

    # Setup
    plt.figure(figsize=(8, 4.0))
    plt.clf()
    ax = plt.gca()

    # Scatter plot
    res = (mtilt-ytilt) - yfit
    ax.scatter(mtilt, res)

    rms = np.std(res)
    ax.text(0.90, 0.90, 'Slit {:d}:  RMS (pix) = {:0.5f}'.format(slit, rms),
            transform=ax.transAxes, size='large', ha='right', color='black')
    # Label
    ax.set_xlabel('Row')
    ax.set_ylabel('Residual (pix)')

    # Finish
    plt.tight_layout(pad=0.2, h_pad=0.0, w_pad=0.0)
    if show_QA:
        plt.show()
    else:
        plt.savefig(outfile, dpi=400)
    plt.close()

    plt.rcdefaults()

    return





def coeff2tilts_old(coeff2, shape, func2D, max_tilt=1.2, min_tilt=-0.2):
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
