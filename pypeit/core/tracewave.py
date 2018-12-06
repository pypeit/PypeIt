""" Module for methods related to tracing arc/sky lines across a slit/order
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import inspect
import copy

import numpy as np
import matplotlib.pyplot as plt

from pypeit import msgs
from pypeit.core import arc
from pypeit import utils
from pypeit.core import qa
from pypeit.core import trace_slits
from pypeit.core import extract
from astropy.stats import sigma_clipped_stats



try:
    from pypeit import ginga
except ImportError:
    pass


def tilts_find_lines(arc_spec, slit_cen, tracethresh=10.0, sig_neigh=5.0, nfwhm_neigh=5.0,
                    only_these_lines=None, fwhm=4.0, nonlinear_counts=1e10, fit_frac_fwhm=1.25, cont_frac_fwhm=1.0,
                    max_frac_fwhm=2.0, cont_samp=30, niter_cont=3, debug_lines=False, debug_peaks=False):


    nspec = arc_spec.size
    spec_vec = np.arange(nspec)
    # Find peaks with a liberal threshold of sigdetect = 5.0
    tampl_tot, tampl_cont_tot, tcent_tot, twid_tot, _, wgood, arc_cont_sub, nsig_tot = arc.detect_lines(
        arc_spec, sigdetect=np.min([sig_neigh,tracethresh]), fwhm=fwhm, fit_frac_fwhm=fit_frac_fwhm, cont_frac_fwhm=cont_frac_fwhm,
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
    # Remove lines that are within npix_neigh pixels.
    # #ToDO replce this with a near-neighbor based approach, where
    # we identify groups and take the brightest line in a given group?
    nuse = np.sum(aduse)
    detuse = arcdet[aduse]
    idxuse = np.arange(arcdet.size)[aduse]
    olduse = aduse.copy()
    for s in range(nuse):
        w = np.where((np.abs(arcdet - detuse[s]) <= npix_neigh) & (np.abs(arcdet - detuse[s]) >= 1.0) & (nsig > sig_neigh))[0]
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


def trace_tilts_work(arcimg, lines_spec, lines_spat, thismask, inmask=None, gauss=False, tilts_guess=None, fwhm=4.0,
                     spat_order=3, maxdev_tracefit=0.2,sigrej_trace=3.0, max_badpix_frac=0.20,
                     tcrude_maxerr=1.0, tcrude_maxshift=3.0, tcrude_maxshift0=3.0,tcrude_nave=5, show_tracefits=False):


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
    gauss: bool, default = False
        If true the code will trace the arc lines usign Gaussian weighted centroiding (trace_gweight) instead of the default,
        which is flux weighted centroiding (trace_fweight)
    tilts_guess: float ndarray, default = None
        A guess for the tilts used for running this tilt tracing in an iterative manner. If the tilts_guess is not None, it
        should be an array containing the tilts from a previous iteration which will be used as a crutch for the tracing of
         the tilts. The default is None, which is how this code is run on a first iteration. In that case the crutces are
         determined via trace_crude, and then the flux (or Gaussian) weighted tracing is performed.
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
    tcrude_maxerr: float, default = 1.0
         maxerr parameter for trace crude
    tcrude_maxshift: float, default = 3.0
         maxshift parameter for trace crude
    tcrude_maxshift0: float, default = 3.0
         maxshift0 parameter for trace crude
    tcrude_nave: int, default = 5
       Trace crude is used to determine the initial arc line tilts, which are then iteratively fit. Trace crude
       can optionally boxcar smooth the image (along the spatial direction of the image, i.e. roughly along the arc line tilts)
       to improve the tracing.
    show_tracefits: bool, default = False
       If true the fits will be shown to each arc line trace by iterative_fitting
    """


    nspec, nspat = arcimg.shape
    spec_vec = np.arange(nspec)
    spat_vec = np.arange(nspat)
    slit_widp2 = int(np.ceil((np.sum(thismask,axis=1)).max()) + 2)
    slit_width_even = slit_widp2 if slit_widp2 % 2 == 0 else slit_widp2 + 1
    trace_int = slit_width_even//2

    maxdev = maxdev_tracefit*fwhm     # maxdev is fraction of fwhm
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
    tilts_sub_spec = np.outer(np.ones(nsub), lines_spec)
    tilts_sub_spat = np.outer(np.arange(nsub), np.ones(nlines))
    tilts_sub_dspat = np.zeros_like(tilts_sub_spat)


    tilts = np.zeros((nspat, nlines))
    tilts_fit = np.zeros((nspat, nlines))
    tilts_err = np.zeros((nspat, nlines))
    tilts_mask = np.zeros((nspat, nlines), dtype=bool)  # This is true if the pixel was in a region traced
    tilts_spec = np.outer(np.ones(nspat), lines_spec)
    tilts_spat = np.outer(np.arange(nspat), np.ones(nlines))
    tilts_dspat= np.zeros_like(tilts_spat)
    # interpolate traces onto this location and compute the pixel offset from the center of the slit
    for iline in range(nlines):
        tilts_dspat[:,iline] = (spat_vec - lines_spat[iline])

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
        tilts_sub_fit_out, tilts_sub_out, tilts_sub_err_out, tset_out = extract.iter_tracefit(
            sub_img, tilts_guess_now, spat_order, inmask=sub_inmask, fwhm=fwhm, maxdev=maxdev, niter=6, idx=str(iline),
            show_fits=show_tracefits)
        if gauss: # If gauss is set, do a Gaussian refinement to the flux weighted tracing
            tilts_sub_fit_gw, tilts_sub_gw, tilts_sub_err_gw, tset_gw = extract.iter_tracefit(
                sub_img, tilts_sub_fit_out, spat_order, inmask=sub_inmask, fwhm=fwhm, maxdev=maxdev, niter=3, idx=str(iline),
                show_fits=show_tracefits)
            tilts_sub_fit_out = tilts_sub_fit_gw
            tilts_sub_out = tilts_sub_gw
            tilts_sub_err_out = tilts_sub_err_gw
        # Boxcar extract the thismask to have a mask indicating whether a tilt is defined along the spatial direction
        tilts_sub_mask_box = (extract.extract_boxcar(sub_thismask, tilts_sub_fit_out, fwhm/2.0) > 0.99*fwhm)
        # Pack the results into arrays, accounting for possibly falling off the image
        # Deal with possibly falling off the chip
        if spat_min[iline] < 0:
            tilts_sub[      -spat_min[iline]:,iline] = tilts_sub_out.flatten()
            tilts_sub_fit[  -spat_min[iline]:,iline] = tilts_sub_fit_out.flatten()
            tilts_sub_err[  -spat_min[iline]:,iline] = tilts_sub_err_out.flatten()
            tilts_sub_mask[ -spat_min[iline]:,iline] = tilts_sub_mask_box
            tilts_sub_dspat[-spat_min[iline]:,iline] = tilts_dspat[min_spat:max_spat,iline]
            tilts[     min_spat:max_spat,iline] = tilts_sub[     -spat_min[iline]:,iline]
            tilts_fit[ min_spat:max_spat,iline] = tilts_sub_fit[ -spat_min[iline]:,iline]
            tilts_err[ min_spat:max_spat,iline] = tilts_sub_err[ -spat_min[iline]:,iline]
            tilts_mask[min_spat:max_spat,iline] = tilts_sub_mask[-spat_min[iline]:,iline]
        elif spat_max[iline] > (nspat - 1):
            tilts_sub[      :-(spat_max[iline] - nspat + 1),iline] = tilts_sub_out.flatten()
            tilts_sub_fit[  :-(spat_max[iline] - nspat + 1),iline] = tilts_sub_fit_out.flatten()
            tilts_sub_err[  :-(spat_max[iline] - nspat + 1),iline] = tilts_sub_err_out.flatten()
            tilts_sub_mask[ :-(spat_max[iline] - nspat + 1),iline] = tilts_sub_mask_box
            tilts_sub_dspat[:-(spat_max[iline] - nspat + 1),iline] = tilts_dspat[min_spat:max_spat,iline]

            tilts[     min_spat:max_spat,iline] = tilts_sub[     :-(spat_max[iline] - nspat + 1),iline]
            tilts_fit[ min_spat:max_spat,iline] = tilts_sub_fit[ :-(spat_max[iline] - nspat + 1),iline]
            tilts_err[ min_spat:max_spat,iline] = tilts_sub_err[ :-(spat_max[iline] - nspat + 1),iline]
            tilts_mask[min_spat:max_spat,iline] = tilts_sub_mask[:-(spat_max[iline] - nspat + 1),iline]
        else:
            tilts_sub[      :,iline] = tilts_sub_out.flatten()
            tilts_sub_fit[  :,iline] = tilts_sub_fit_out.flatten()
            tilts_sub_err[  :,iline] = tilts_sub_err_out.flatten()
            tilts_sub_mask[ :,iline] = tilts_sub_mask_box
            tilts_sub_dspat[:,iline] = tilts_dspat[min_spat:max_spat,iline]
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

    if (nuse < 0.05*nlines):
        msgs.warn('This slit/order would have too many rejected lines.' + msgs.newline() +
                  'We should be rejecting nuse ={:d} out of nlines = {:d} total lines'.format(nuse,nlines) +  msgs.newline() +
                  'We are will proceed without rejecting anything but something is probably very wrong with this slit/order.')
        use_tilt = np.ones(nlines,dtype=bool)
        nuse = nlines

    # Tighten it up with Gaussian weighted centroiding
    trc_tilt_dict = dict(nspec = nspec, nspat = nspat, nsub = nsub, nlines = nlines, nuse = nuse,
                         spat_min=spat_min, spat_max=spat_max, do_crude=do_crude, fwhm = fwhm, use_tilt=use_tilt,
                         tilts_sub_spec=tilts_sub_spec, tilts_sub_spat=tilts_sub_spat, tilts_sub_dspat=tilts_sub_dspat,
                         tilts_sub=tilts_sub, tilts_sub_fit=tilts_sub_fit, tilts_sub_err=tilts_sub_err,
                         tilts_sub_mask=tilts_sub_mask,
                         tilts_spec=tilts_spec, tilts_spat=tilts_spat, tilts_dspat=tilts_dspat,
                         tilts=tilts, tilts_fit=tilts_fit, tilts_err=tilts_err, tilts_mask=tilts_mask)

    return trc_tilt_dict


def trace_tilts(arcimg, lines_spec, lines_spat, thismask, inmask=None, gauss=False, fwhm=4.0,spat_order=5, maxdev_tracefit=0.2,
                sigrej_trace=3.0, max_badpix_frac=0.20, tcrude_nave = 5,
                npca = 1, coeff_npoly_pca = 2, sigrej_pca = 2.0,debug_pca = False, show_tracefits=False):

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
    gauss: bool, default = False
        If true the code will trace the arc lines usign Gaussian weighted centroiding (trace_gweight) instead of the default,
        which is flux weighted centroiding (trace_fweight)
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
        arcimg, lines_spec, lines_spat, thismask, inmask=inmask, gauss=gauss, tilts_guess=None, fwhm=fwhm, spat_order=spat_order,
        maxdev_tracefit=maxdev_tracefit,sigrej_trace=sigrej_trace, max_badpix_frac=max_badpix_frac,
        tcrude_maxerr=tcrude_maxerr, tcrude_maxshift=tcrude_maxshift,tcrude_maxshift0=tcrude_maxshift0,
        tcrude_nave=tcrude_nave, show_tracefits=show_tracefits)

    # TODO THE PCA may not be necessary. It appears to improve the results though for some instruments where the
    # tracing is problematic. We could consider making this optional to speed things up.
    debug_pca_fit = False
    if debug_pca_fit:
        # !!!! FOR TESTING ONLY!!!!  Evaluate the model fit to the tilts for all of our lines
        msgs.info('TESTING: Performing an initial fit before PCA.')
        # JFH Note spec_order is hard wired here as we don't pass it in
        tilt_fit_dict0 = fit_tilts(trace_dict0, spat_order=spat_order, spec_order=6, debug=True,
                                   maxdev=0.2, sigrej=3.0,doqa=True, setup='test', slit=0, show_QA=True)

    # Do a PCA fit, which rejects some outliers
    iuse = trace_dict0['use_tilt']
    nuse = np.sum(iuse)
    msgs.info('PCA modeling {:d} good tilts'.format(nuse))
    pca_fit, poly_fit_dict, pca_mean, pca_vectors = extract.pca_trace(
        trace_dict0['tilts_sub_fit'], predict=np.invert(iuse), npca=npca, coeff_npoly=coeff_npoly_pca,
        lower=sigrej_pca, upper=sigrej_pca, order_vec=lines_spec, xinit_mean=lines_spec,
        minv=0.0, maxv=float(trace_dict0['nsub'] - 1), debug=debug_pca)

    # Now trace again with the PCA predictions as the starting crutches
    trace_dict1 = trace_tilts_work(arcimg, lines_spec, lines_spat, thismask, inmask=inmask, gauss=gauss, tilts_guess=pca_fit,
                                      fwhm=fwhm, spat_order=spat_order, maxdev_tracefit=maxdev_tracefit,sigrej_trace=sigrej_trace,
                                      max_badpix_frac=max_badpix_frac,show_tracefits=show_tracefits)

    return trace_dict1


def fit_tilts(trc_tilt_dict, spat_order=3, spec_order=4, maxdev = 0.2, sigrej = 3.0, func2d='legendre2d', doqa=True, setup = 'test',
              slit = 0, show_QA=False, out_dir=None, debug=False):
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

    import matplotlib as mpl
    from matplotlib.lines import Line2D

    nspec = trc_tilt_dict['nspec']
    nspat = trc_tilt_dict['nspat']
    fwhm = trc_tilt_dict['fwhm']
    maxdev_pix = maxdev*fwhm
    xnspecmin1 = float(nspec-1)
    xnspatmin1 = float(nspat-1)
    nspat = trc_tilt_dict['nspat']
    use_tilt = trc_tilt_dict['use_tilt']                 # mask for good/bad tilts, based on aggregate fit, frac good pixels
    nuse = np.sum(use_tilt)
    tilts = trc_tilt_dict['tilts'][:,use_tilt]   # legendre polynomial fit
#   JFH Before we were fitting the fits. Now we fit the actual flux weighted centroided tilts.
#   tilts_fit = trc_tilt_dict['tilts_fit'][:,use_tilt]   # legendre polynomial fit
    tilts_err = trc_tilt_dict['tilts_err'][:,use_tilt]   # gaussian weighted centroidding error
    tilts_dspat = trc_tilt_dict['tilts_dspat'][:,use_tilt] # spatial offset from the central trace
    tilts_spec = trc_tilt_dict['tilts_spec'][:,use_tilt] # line_spec spectral pixel position
    tilts_mask = trc_tilt_dict['tilts_mask'][:,use_tilt] # Reflects if trace is on the slit

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

    # Fits are done in dimensionless coordinates to allow for different binnings between i.e. the arc and the science frame
    fitmask, coeff2 = utils.robust_polyfit_djs(tilts_dspat.flatten()/xnspatmin1,
                                               (tilts.flatten() - tilts_spec.flatten())/xnspecmin1,
                                               fitxy,x2=tilts.flatten()/xnspecmin1, inmask = tot_mask.flatten(),
                                               function=func2d, maxiter=100, lower=sigrej, upper=sigrej,
                                               maxdev=maxdev_pix/xnspecmin1,minx=-1.0, maxx=1.0, minx2=0.0, maxx2=1.0,
                                               use_mad=True, sticky=False)

    fitmask = fitmask.reshape(tilts_dspat.shape)
    delta_spec_fit1 = xnspecmin1*utils.func_val(coeff2, tilts_dspat[tot_mask]/xnspatmin1, func2d, x2=tilts[tot_mask]/xnspecmin1,
                                               minx=-1.0, maxx=1.0, minx2=0.0, maxx2=1.0)
    delta_spec_fit = np.zeros_like(tilts_dspat)
    delta_spec_fit[tot_mask] = delta_spec_fit1
    # Residuals in pixels
    res = (tilts[fitmask] - tilts_spec[fitmask]) - delta_spec_fit[fitmask]
    rms = np.std(res)
    msgs.info("RMS (pixels): {}".format(rms))
    msgs.info("RMS/FWHM: {}".format(rms/fwhm))
    # These are locations that were fit but were rejected
    rej_mask = tot_mask & np.invert(fitmask)

    # TODO: Make these plots a part of the standard QA and write to a file
    if debug:

        # QA1
        xmin = 1.1*tilts_dspat[tot_mask].min()
        xmax = 1.1*tilts_dspat[tot_mask].max()

        # Show the fit
        fig, ax = plt.subplots(figsize=(12, 18))
        # dummy mappable shows the spectral pixel
        #dummie_cax = ax.scatter(lines_spec, lines_spec, c=lines_spec, cmap=cmap)
        ax.cla()
        ax.plot(tilts_dspat[tot_mask], tilts[tot_mask], color='black', linestyle=' ', mfc ='None', marker='o',
                markersize = 9.0, markeredgewidth=1.0,zorder=4, label='Good Tilt')
        ax.plot(tilts_dspat[rej_mask], tilts[rej_mask], color ='red',linestyle=' ', mfc = 'None', marker='o',
                markersize=9.0, markeredgewidth=2.0, zorder=5, label='Rejected')
        ax.plot(tilts_dspat[tot_mask], tilts_spec[tot_mask] + delta_spec_fit[tot_mask], color='black', linestyle=' ', marker='o',
                markersize = 2.0,markeredgewidth=1.0,zorder=1)

        ax.set_xlim((xmin,xmax))
        ax.set_xlabel('Spatial Offset from Central Trace (pixels)', fontsize=15)
        ax.set_ylabel('Spectral Pixel',fontsize=15)
        ax.legend()
        ax.set_title('Tilts vs Fit (spat_order, spec_order)=({:d},{:d}) for slit={:d}: RMS = {:5.3f}, '
                     'RMS/FWHM={:5.3f}'.format(spat_order,spec_order,slit,rms, rms/fwhm),fontsize=15)
        plt.show()

        # QA2
        # Show the fit residuals as a function of spatial position
        line_indx = np.outer(np.ones(nspat), np.arange(nuse))
        lines_spec = tilts_spec[0,:]
        cmap = mpl.cm.get_cmap('coolwarm', nuse)

        fig, ax = plt.subplots(figsize=(14, 12))
        # dummy mappable shows the spectral pixel
        dummie_cax = ax.scatter(lines_spec, lines_spec, c=lines_spec, cmap=cmap)
        ax.cla()

        for iline in range(nuse):
            iall = (line_indx == iline) & tot_mask
            irej = (line_indx == iline) & tot_mask & rej_mask
            this_color = cmap(iline)
            # plot the residuals
            ax.plot(tilts_dspat[iall], (tilts[iall] - tilts_spec[iall]) - delta_spec_fit[iall], color=this_color,
                    linestyle='-', linewidth=3.0, marker='None', alpha=0.5)
            ax.plot(tilts_dspat[irej],(tilts[irej] - tilts_spec[irej]) - delta_spec_fit[irej],linestyle=' ',
                    marker='o', color = 'limegreen', mfc='limegreen', markersize=5.0)

        ax.hlines(0.0, xmin, xmax, linestyle='--', linewidth=2.0, color='k', zorder=10)

        legend_elements = [Line2D([0], [0], color='cornflowerblue', linestyle='-', linewidth=3.0, label='residual'),
                           Line2D([0], [0], color='limegreen', linestyle=' ', marker='o', mfc='limegreen', markersize=7.0, label='rejected')]

        ax.set_xlim((xmin,xmax))
        ax.set_xlabel('Spatial Offset from Central Trace (pixels)')
        ax.set_ylabel('Arc Line Tilt Residual (pixels)')
        ax.legend(handles=legend_elements)
        cb = fig.colorbar(dummie_cax, ticks=lines_spec)
        cb.set_label('Spectral Pixel')
        plt.show()

    # TODO Add QA where we overlay the final model of the tilts on the image using ginga!

    # QA
    if doqa:
        plot_tiltres(setup, tilts[tot_mask], tilts_spec[tot_mask], delta_spec_fit[tot_mask], fitmask[tot_mask], fwhm, slit=slit, show_QA=show_QA, out_dir=out_dir)

    tilt_fit_dict = dict(nspec = nspec, nspat = nspat, ngood_lines=np.sum(use_tilt), npix_fit = np.sum(tot_mask),
                         npix_rej = np.sum(fitmask == False), coeff2=coeff2, spec_order = spec_order, spat_order = spat_order,
                         minx = -1.0, maxx = 1.0, minx2 = 0.0, maxx2 = 1.0, func=func2d)

    return tilt_fit_dict



def fit2tilts(shape, slit_cen, coeff2, func):
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
    nspec, nspat = shape
    xnspecmin1 = float(nspec-1)
    xnspatmin1 = float(nspat-1)

    spat_img = np.outer(np.ones(nspec), np.arange(nspat)) # spatial position everywhere along image
    slit_cen_img = np.outer(slit_cen, np.ones(nspat))   # center of the slit replicated spatially
    # Normalized spatial offset image (from central trace)
    dspat_img_nrm = (spat_img - slit_cen_img)/xnspatmin1
    # normalized spec image
    spec_img_nrm = np.outer(np.arange(nspec), np.ones(nspat))/xnspecmin1

    delta_spec_fit_nrm = utils.func_val(coeff2, dspat_img_nrm, func, x2=spec_img_nrm, minx=-1.0, maxx=1.0, minx2=0.0, maxx2=1.0)
    tilts = (spec_img_nrm - delta_spec_fit_nrm)
    # Added this to ensure that tilts are never crazy values due to extrapolation of fits which can break
    # wavelength solution fitting
    tilts = np.fmax(np.fmin(tilts, 1.2),-0.2)

    return tilts


def plot_tiltres(setup, mtilt, ytilt, yfit, fitmask, fwhm, slit=None, outfile=None, show_QA=False, out_dir=None):
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
    plt.figure(figsize=(8, 4))
    plt.clf()
    ax = plt.gca()

    # Scatter plot
    res = (mtilt - ytilt) - yfit
    xmin = 0.95*np.asarray([ytilt,mtilt]).min()
    xmax = 1.05*np.asarray([ytilt,mtilt]).max()
    ax.hlines(0.0, xmin, xmax,linestyle='--', color='green')
    ax.plot(ytilt[fitmask], (res[fitmask]), 'ko', mfc='k', markersize=5.0, label='Good Points')
    ax.plot(ytilt[~fitmask],(res[~fitmask]), 'ro', mfc='r', markersize=5.0, label='Rejected Points')
    #ax.scatter(mtilt, res)
    rms = np.std(res[fitmask])
    ax.text(0.90, 0.90, 'Slit {:d}:  RMS (pixels) = {:0.5f}'.format(slit, rms),
            transform=ax.transAxes, size='large', ha='right', color='black',fontsize=16)
    ax.text(0.90, 0.80, ' Slit {:d}:  RMS/FWHM = {:0.5f}'.format(slit, rms/fwhm),
            transform=ax.transAxes, size='large', ha='right', color='black',fontsize=16)
    # Label
    ax.set_xlabel('Spectral Pixel')
    ax.set_ylabel('RMS (pixels)')
    ax.set_title('RMS of Each Arc Line Traced')
    ax.set_xlim((xmin,xmax))
    # Finish
    plt.tight_layout(pad=0.2, h_pad=0.0, w_pad=0.0)
    if show_QA:
        plt.show()
    else:
        plt.savefig(outfile, dpi=400)
    plt.close()

    plt.rcdefaults()

    return




#TODO this should be deleted, but I'm saving some it for now becasue I'm still testing.
'''


        # JFH Code block starts here
        ########
        # def get_tilts(npca = 2, fwhm=4.0, ncoeff=5, maxdev_tracefit=0.1,percentile_reject=0.10, max_badpix_frac=0.20, tcrude_maxerr=1.0,
        # tcrude_maxshift=3.0, tcrude_maxshift0=3.0, tcrude_nave=5,)

        # DEBUGGIN
        slit = 4

        from pypeit.core import pixels
        from pypeit.core import extract


        nspat = self.msarc.shape[1]
        nspec = self.msarc.shape[0]
        arcimg = self.msarc
        arc_spec = self.arccen[:, slit]
        slit_left = self.tslits_dict['lcen'][:,slit].copy()
        slit_righ = self.tslits_dict['rcen'][:,slit].copy()
        inmask = (self.bpm == False)
        # Center of the slit
        slit_cen = (slit_left + slit_righ)/2.0

        slitmask = self.spectrograph.slitmask(self.tslits_dict)
        thismask = slitmask == slit

        # Tilt specific Optional parameters
        tracethresh = 10.0 # significance threshold for an arc line to be traced
        sigdetect = 5.0 # This is the significance threshold for finding neighboring lines. The number of these neighboring lines
        #  determines how many nsig > tracethresh lines that may be removed because they are too close.
        nfwhm_neigh = 3.0
        only_these_lines = None # These are lines from the wavelength solution that are known to be good. If set code only uses these
        # identified lines to avoid problems with ghosts (is this still necessary with new tracing code?)
        maxdev_tracefit = 1.0 # maximum deviation for iterative trace fitting (flux weighted, then gaussian weighted)
        sigrej_trace = 3.0 # From each line we compute a median absolute deviation of the trace from the polynomial fit. We then
        # analyze the distribution of mads for all the lines, and reject sigrej_trace outliers from that distribution.
        npca = 1 # Order of PCA for iteration to improve tracing
        coeff_npoly_pca = 1 # Polynomial order for fit to PCA coefficients to improve tracing
        sigrej_pca = 2.0 # Outlier rejection significance for PCA fit to arc line coefficients
        ncoeff = 5 # Order of legendre polynomial fits to the tilts

        max_badpix_frac = 0.20 # Maximum fraction of total pixels masked by the trace_gweight algorithm (because the residuals are too large)
        # Trace Crude parameters
        tcrude_maxerr = 1.0 #
        tcrude_maxshift = 3.0
        tcrude_maxshift0 = 3.0
        tcrude_nave = 5
        show_tracefits = False # show the trace fits

        # These are parameters for 2D fitting
        spec_order = 4
        maxdev_2d = 1.0 # Maximum absolute deviation used in rejection for 2d legendre polynomial fit to ensemble of all spectral line tilts
        sigrej_2d = 3.0 # Outlier significance threshold for rejection for 2d legendre polynomial fit to ensemble of all spectral line tilts
        show_tilts = True
        debug = True

        trc_tilt_dict0 = tracewave.trace_tilts_work(
            arcimg, lines_spec, lines_spat, thismask, inmask=inmask,
            tilts_guess=None,fwhm=fwhm,ncoeff=ncoeff, maxdev_tracefit=maxdev_tracefit,
            sigrej_trace=sigrej_trace, max_badpix_frac=max_badpix_frac,
            tcrude_maxerr=tcrude_maxerr, tcrude_maxshift=tcrude_maxshift,
            tcrude_maxshift0=tcrude_maxshift0,
            tcrude_nave=tcrude_nave,show_fits=show_tracefits)
        # Now evaluate the model of the tilts for all of our lines
        # Testing!!!!
        # Now perform a fit to the tilts
        tilt_fit_dict0 = tracewave.fit_tilts(trc_tilt_dict0, spat_order=ncoeff, spec_order=spec_order,debug=True,
                                             maxdev=maxdev_2d, sigrej=sigrej_2d,
                                             doqa=True,setup='test',slit=0, show_QA=False, out_dir='./')

        # Do a PCA fit, which rejects some outliers
        iuse = trc_tilt_dict0['use_tilt']
        nuse =np.sum(iuse)
        msgs.info('PCA modeling {:d} good tilts'.format(nuse))
        #pca_explained_var = None
        # TODO Should we truncate this PCA by hand, or just let it explain variance
        # PCA fit good orders, and predict bad orders
        pca_fit, poly_fit_dict, pca_mean, pca_vectors = extract.pca_trace(
            trc_tilt_dict0['tilts_sub_fit'], predict=np.invert(iuse), npca = npca, coeff_npoly=coeff_npoly_pca,
            lower=sigrej_pca, upper=sigrej_pca, order_vec=lines_spec,xinit_mean=lines_spec,
            minv=0.0, maxv=float(trc_tilt_dict0['nsub'] - 1), debug=True)

        # Now trace again with the PCA predictions as the starting crutches
        trc_tilt_dict1 = tracewave.trace_tilts_work(
            arcimg, lines_spec, lines_spat, thismask, inmask=inmask,
            tilts_guess=pca_fit, fwhm=fwhm,ncoeff=ncoeff, maxdev_tracefit=maxdev_tracefit,
            sigrej_trace=sigrej_trace, max_badpix_frac=max_badpix_frac,
            show_fits=show_tracefits)


        # Now perform a fit to the tilts
        tilt_fit_dict1 = tracewave.fit_tilts(trc_tilt_dict1, spat_order=ncoeff, spec_order=spec_order, debug=True,
                                             maxdev=maxdev_2d, sigrej=sigrej_2d,
                                             doqa=True,setup='test',slit=0, show_QA=False, out_dir='./')
        # Now evaluate the model of the tilts for all of our lines
        piximg1 = tracewave.fit2piximg(tilt_fit_dict1)

        # Now trace again with the piximg model as the starting crutches


        """
        # Since the y independent variable is the tilts in the way we do the 2d fit, and soem tilts are spurios, it does
        # no good to evaluate the global fit at these spurious tilts to get the new tracing crutches. The following is a
        # hack to determine this from the piximg. Perhaps there is a better way.
        spec_img = np.outer(np.arange(nspec), np.ones(nspat))
        spec_vec = np.arange(nspec)
        nlines=len(lines_spec)
        interp_flag = np.ones(nlines,dtype=bool)
        tilts_crutch = np.zeros((nspat, nlines))
        spat_min = trc_tilt_dict1['spat_min']
        spat_max = trc_tilt_dict1['spat_max']
        piximg1[np.invert(thismask)] = 1e10
        # Is there a faster more clever way to do this?
        for iline in range(nlines):
            min_spat = np.fmax(spat_min[iline], 0)
            max_spat = np.fmin(spat_max[iline], nspat - 1)
            min_spec = int(np.fmax(np.round(lines_spec[iline]) - np.round(0.01*nspec),0))
            max_spec = int(np.fmin(np.round(lines_spec[iline]) + np.round(0.01*nspec),nspec-1))
            piximg_sub = piximg1[min_spec:max_spec,:]
            for ispat in range(min_spat,max_spat):
                if np.any(np.diff(piximg_sub[:,ispat] < 0)):
                    # If we ever encounter an unsorted piximg_sub the logic below makes no sense so just use the
                    # previous polynomial fit as the crutch and exit this loop
                    tilts_crutch[ispat,iline] = trc_tilt_dict1['tilts_fit'][:,iline]
                    msgs.warn('piximg is not monotonic. Your tilts are probably bad')
                    break
                else:
                    tilts_crutch[ispat,iline] = np.interp(lines_spec[iline],piximg_sub[:,ispat],spec_vec[min_spec:max_spec])


        trc_tilt_dict1['tilts_crutch'] = tilts_crutch
        # Now trace again with the PCA predictions as the starting crutches
        trc_tilt_dict2 = tracewave.trace_tilts(arcimg, lines_spec, lines_spat, slit_width, thismask, inmask=inmask,
                                               tilts_guess=tilts_crutch, fwhm=fwhm,ncoeff=ncoeff, maxdev_tracefit=maxdev_tracefit,
                                               sigrej_trace=sigrej_trace, max_badpix_frac=max_badpix_frac,
                                               show_fits=show_tracefits)

        # Now perform a second fit to the tilts
        tilt_fit_dict2 = tracewave.fit_tilts(trc_tilt_dict2, spat_order=ncoeff, spec_order=spec_order, debug=True,
                                             maxdev=maxdev_2d, sigrej=sigrej_2d,
                                             doqa=True,setup='test',slit=0, show_QA=False, out_dir='./')
        # Now evaluate the model of the tilts for all of our lines
        piximg2 = tracewave.fit2piximg(tilt_fit_dict2)
        """

#        from matplotlib import pyplot as plt
#        tilt_mask = trc_tilt_dict1['tilts_mask']
#        plt.plot(trc_tilt_dict1['tilts_spat'][tilt_mask], tilts_crutch[tilt_mask], 'ko', markersize=2.0)
#        plt.plot(trc_tilt_dict1['tilts_spat'][tilt_mask], trc_tilt_dict1['tilts_fit'][tilt_mask], 'ro', mfc='none', markersize=2.0)


        # Now do a fit


        """
        for iline in range(nlines):
            min_spat = np.fmax(spat_min[iline], 0)
            max_spat = np.fmin(spat_max[iline], nspat - 1)
            sub_img = (piximg*thismask)[:, min_spat:max_spat]
            spec_img_sub = spec_img[  :, min_spat:max_spat]
            ispec_min = np.argmin(np.abs(sub_img - lines_spec[iline]), axis=0)
            tilts_crutch[min_spat:max_spat,iline] = spec_img_sub[ispec_min,np.arange(sub_img.shape[1])]
        """



'''