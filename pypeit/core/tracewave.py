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
import matplotlib as mpl
from matplotlib.lines import Line2D
import scipy
from pypeit.core import pydl

try:
    from pypeit import ginga
except ImportError:
    pass


def tilts_find_lines(arc_spec, slit_cen, tracethresh=10.0, sig_neigh=5.0, nfwhm_neigh=3.0,
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


def trace_tilts_work(arcimg, lines_spec, lines_spat, thismask, slit_cen, inmask=None, gauss=False, tilts_guess=None, fwhm=4.0,
                     spat_order=3, maxdev_tracefit=0.02,sigrej_trace=3.0, max_badpix_frac=0.30,
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
    maxdev_tracefit: float, default = 0.2
       Maximum absolute deviation for the arc tilt fits during iterative trace fitting expressed in units of the fwhm.
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
    slit_width_even = np.fmax(slit_widp2 if slit_widp2 % 2 == 0 else slit_widp2 + 1, nspec-1)
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

    # The sub arrays hold the sub-imaged tilts
    #tilts_sub = np.zeros((nsub, nlines))       # Thee trace_fweight (or gweighed) tilts

    #tilts_sub_err = np.zeros((nsub, nlines))   # errors on the tilts (only used for masking but not weighted fitting)
    #tilts_sub_mask = np.zeros((nsub, nlines), dtype=bool)  # mask indicating where the tilts are actually covered in the sub image, i.e. where thismask != False
    #tilts_sub_spec = np.outer(np.ones(nsub), lines_spec)   # spectral coordinate of each tilt, which is the arc line spectral pixel location
    #tilts_sub_spec_fit = np.zeros((nsub, nlines))   # spectral position determined by evaluating the tilt fit at the center of the slit

    #tilts_sub_dspat = np.zeros_like(tilts_sub_spat) # delta position of the tilt in pixels, i.e. difference between slitcen and the spatial coordinate above

    # PCA fitting uses the sub-imaged fits, so we need them
    tilts_sub_fit = np.zeros((nsub, nlines))   # legendre polynomial fits to the tilt traces
    tilts_sub_spat = np.outer(np.arange(nsub), np.ones(nlines)) # spatial coordinate along each tilt

    tilts = np.zeros((nspat, nlines))      # The trace_fweight (or gweighed) tilts
    tilts_fit = np.zeros((nspat, nlines))  # legendre polynomial fits to the tilt traces
    tilts_err = np.zeros((nspat, nlines))   # errors on the tilts (only used for masking but not weighted fitting)
    tilts_mask = np.zeros((nspat, nlines), dtype=bool)  # This is true if the pixel was in a region traced
    tilts_spec = np.zeros((nspat, nlines)) # spectral position determined by evaluating the tilt fit at the center of the slit
    tilts_spat = np.outer(np.arange(nspat), np.ones(nlines))  # spatial coordinate along each tilt
    tilts_dspat= np.zeros_like(tilts_spat)  # delta position of the tilt in pixels, i.e. difference between slitcen and the spatial coordinate above

    # Transposed image and masks for traceing
    arcimg_trans = (arcimg * thismask).T
    inmask_trans = (inmask * thismask).T.astype(float)
    thismask_trans = thismask.T

    # 1) Trace the tilts from a guess. If no guess is provided from a previous iteration use trace_crude
    for iline in range(nlines):
        # We sub-image each tilt using a symmetric window about the (integer) spatial location of each line,
        # which is the slitcen evaluated at the line spectral position.
        spat_min[iline] = lines_spat_int[iline] - trace_int  # spat_min is the minium location of the sub-image
        spat_max[iline] = lines_spat_int[iline] + trace_int + 1  # spat_max is the maximum location of the sub-image
        min_spat = np.fmax(spat_min[iline], 0)  # These min_spat and max_spat are to prevent leaving the image
        max_spat = np.fmin(spat_max[iline], nspat - 1)
        sub_img = arcimg_trans[min_spat:max_spat, :]
        sub_inmask = inmask_trans[min_spat:max_spat,:]
        sub_thismask = thismask_trans[min_spat:max_spat,:]
        if do_crude: # First time tracing, do a trace crude
            tilts_guess_now, err_now = trace_slits.trace_crude_init(
                sub_img, np.array([lines_spec[iline]]), (sub_img.shape[0] - 1) // 2, invvar=sub_inmask, radius=fwhm,
                nave=tcrude_nave, maxshift0=tcrude_maxshift0, maxshift=tcrude_maxshift, maxerr=tcrude_maxerr)
            tilts_guess_now=tilts_guess_now.flatten()
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
        # Boxcar extract the thismask to have a mask indicating whether a tilt is defined along the spatial direction
        tilts_sub_mask_box = (extract.extract_boxcar(sub_thismask, tilts_guess_now, fwhm/2.0) > 0.99*fwhm)
        # If more than 80% of the pixels are masked, then don't mask at all. This happens when the traces leave the good
        # part of the slit. If we proceed with everything masked the iter_tracefit fitting will crash.
        if (np.sum(tilts_sub_mask_box) < 0.8*nsub):
            tilts_sub_mask_box = np.ones_like(tilts_sub_mask_box)
        # Do iterative flux weighted tracing and polynomial fitting to refine these traces. This must also be done in a loop
        # since the sub image is different for every aperture, i.e. each aperature has its own image
        tilts_sub_fit_out, tilts_sub_out, tilts_sub_err_out, tset_out = extract.iter_tracefit(
            sub_img, tilts_guess_now, spat_order, inmask=sub_inmask, trc_inmask = tilts_sub_mask_box, fwhm=fwhm,
            maxdev=maxdev, niter=6, idx=str(iline),show_fits=show_tracefits, xmin=0.0,xmax=float(nsub-1))
        tilts_sub_mask_box = (extract.extract_boxcar(sub_thismask, tilts_sub_fit_out, fwhm/2.0) > 0.99*fwhm)
        if gauss: # If gauss is set, do a Gaussian refinement to the flux weighted tracing
            if (np.sum(tilts_sub_mask_box) < 0.8 * nsub):
                tilts_sub_mask_box = np.ones_like(tilts_sub_mask_box)
            tilts_sub_fit_gw, tilts_sub_gw, tilts_sub_err_gw, tset_gw = extract.iter_tracefit(
                sub_img, tilts_sub_fit_out, spat_order, inmask=sub_inmask, trc_inmask = tilts_sub_mask_box, fwhm=fwhm,
                maxdev=maxdev, niter=3, idx=str(iline),show_fits=show_tracefits, xmin=0.0, xmax=float(nsub-1))
            tilts_sub_fit_out = tilts_sub_fit_gw
            tilts_sub_out = tilts_sub_gw
            tilts_sub_err_out = tilts_sub_err_gw
        tilts_sub_mask_box = (extract.extract_boxcar(sub_thismask, tilts_sub_fit_out, fwhm/2.0) > 0.99*fwhm)

        # Pack the results into arrays, accounting for possibly falling off the image
        # Deal with possibly falling off the chip

        # This is the same for all cases since it is the evaluation of a fit
        tilts_sub_fit[:, iline] = tset_out.xy(tilts_sub_spat[:, iline].reshape(1, nsub))[1]

        # We use the tset_out.xy to evaluate the trace across the whole sub-image even for pixels off the slit. This
        # guarantees that the fits are always evaluated across the whole sub-image which is required for the PCA step.
        if spat_min[iline] < 0:
            #tilts_sub[      -spat_min[iline]:,iline] = tilts_sub_out.flatten()
            #tilts_sub_err[  -spat_min[iline]:,iline] = tilts_sub_err_out.flatten()
            #tilts_sub_mask[ -spat_min[iline]:,iline] = tilts_sub_mask_box.flatten()
            #tilts_sub_dspat[-spat_min[iline]:,iline] = tilts_dspat[min_spat:max_spat,iline]
            tilts[     min_spat:max_spat,iline] = tilts_sub_out.flatten() # tilts_sub[     -spat_min[iline]:,iline]
            tilts_fit[ min_spat:max_spat,iline] = tilts_sub_fit[ -spat_min[iline]:,iline]
            tilts_err[ min_spat:max_spat,iline] = tilts_sub_err_out.flatten() # tilts_sub_err[ -spat_min[iline]:,iline]
            tilts_mask[min_spat:max_spat,iline] = tilts_sub_mask_box.flatten() #tilts_sub_mask[-spat_min[iline]:,iline]
        elif spat_max[iline] > (nspat - 1):
            #tilts_sub[      :-(spat_max[iline] - nspat + 1),iline] = tilts_sub_out.flatten()
            #tilts_sub_err[  :-(spat_max[iline] - nspat + 1),iline] = tilts_sub_err_out.flatten()
            #tilts_sub_mask[ :-(spat_max[iline] - nspat + 1),iline] = tilts_sub_mask_box.flatten()
            #tilts_sub_dspat[:-(spat_max[iline] - nspat + 1),iline] = tilts_dspat[min_spat:max_spat,iline]
            tilts[     min_spat:max_spat,iline] = tilts_sub_out.flatten() # tilts_sub[     :-(spat_max[iline] - nspat + 1),iline]
            tilts_fit[ min_spat:max_spat,iline] = tilts_sub_fit[ :-(spat_max[iline] - nspat + 1),iline]
            tilts_err[ min_spat:max_spat,iline] = tilts_sub_err_out.flatten() # tilts_sub_err[ :-(spat_max[iline] - nspat + 1),iline]
            tilts_mask[min_spat:max_spat,iline] = tilts_sub_mask_box.flatten()  #tilts_sub_mask[:-(spat_max[iline] - nspat + 1),iline]
        else:
            #tilts_sub[      :,iline] = tilts_sub_out.flatten()
            #tilts_sub_err[  :,iline] = tilts_sub_err_out.flatten()
            #tilts_sub_mask[ :,iline] = tilts_sub_mask_box.flatten()
            #tilts_sub_dspat[:,iline] = tilts_dspat[min_spat:max_spat,iline]
            tilts[     min_spat:max_spat,iline] = tilts_sub_out.flatten() #tilts_sub[     :,iline]
            tilts_fit[ min_spat:max_spat,iline] = tilts_sub_fit[ :,iline]
            tilts_err[ min_spat:max_spat,iline] = tilts_sub_err_out.flatten() # tilts_sub_err[ :,iline]
            tilts_mask[min_spat:max_spat,iline] = tilts_sub_mask_box.flatten() #tilts_sub_mask[:,iline]

        # Now use these fits to the traces to get a more robust value of the tilt spectral position and spatial
        # offset from the trace than what was initially determined from the 1d arc line spectrum. This is technically
        # where the slit_cen cross the tilts_fit, but it is a tricky since they are parameterized by different
        # independent variables (slit_cen uses spec_vec, whereas tilts_fit uses spat_vec). This code uses
        # a trick of interpolating the slit_cen onto the arc pixels. If we find it fails, then replace with something
        # simpler that simply iterates to zero on where the two cross.

        # ToDO Fix this later with an iterative thing that also updates spatial reference position of the tilt
        #imask = tilts_mask[:, iline]
        #slit_cen_spat_ontilt = np.interp(tilts_fit[imask,iline],spec_vec, slit_cen)
        #delta_spat = slit_cen_spat_ontilt - tilts_spat[imask,iline]
        # Grab the monotonic indices
        #ediff = np.ediff1d(delta_spat,to_begin=0.0)
        #mono_ind = np.sign(ediff) == np.sign(np.median(ediff))
        #zero_cross_spat = (scipy.interpolate.interp1d(delta_spat[mono_ind],(tilts_spat[imask,iline])[mono_ind],assume_sorted=False))(0.0)
        #spec_fit_now = np.interp(zero_cross_spat,tilts_spat[imask,iline], tilts_fit[imask,iline])
        #spat_fit_now = np.interp(spec_fit_now,spec_vec, slit_cen)
        #tilts_spec[:, iline] = np.full(nspat, spec_fit_now)

        tilts_dspat[:, iline] = (spat_vec - lines_spat[iline])
        imask = tilts_mask[:,iline]
        try:
            spec_fit_now = np.interp(0.0, tilts_dspat[imask, iline], tilts_fit[imask, iline])
        except ValueError:
            spec_fit_now = lines_spec[iline]
        tilts_spec[:,iline] = np.full(nspat, spec_fit_now)


    # Create the mask for the bad lines. Define the error on the bad tilt as being the
    bad_mask = (tilts_err > 900) | (tilts_mask == False)
    on_slit = np.sum(tilts_mask,0)
    on_slit_bad = np.sum((tilts_mask & (tilts_err > 900)),0)
    bad_frac = on_slit_bad/on_slit

    dev_mean, dev_median, dev_sig = sigma_clipped_stats(np.abs(tilts - tilts_fit), mask=bad_mask, sigma=4.0,axis=0)
    good_line = np.any(bad_mask == False,axis=0) # Is it masked everywhere?
    # Median absolute deviation for each line quantifies the goodnes of tracing
    dev_mad = 1.4826*dev_median
    # Now reject outliers from this distribution
    dev_mad_dist_median = np.median(dev_mad[good_line])
    dev_mad_dist_mad = 1.4826*np.median(np.abs(dev_mad[good_line] - dev_mad_dist_median)) # i.e. this is like the sigma
    # Reject lines that are sigrej trace outliers
    mad_rej = np.abs((dev_mad - dev_mad_dist_median)/dev_mad_dist_mad) < sigrej_trace

    # Do we need this dev_mad < maxdev step?
    use_tilt = (mad_rej) & (bad_frac < max_badpix_frac) & good_line & (dev_mad < maxdev)
    nuse = np.sum(use_tilt)
    msgs.info('Number of usable arc lines for tilts: {:d}/{:d}'.format(nuse,nlines))

    tilts_mad = np.outer(np.ones(nspat),dev_mad)

    if (nuse < 0.05*nlines):
        msgs.warn('This slit/order would have too many rejected lines.' + msgs.newline() +
                  'We should be rejecting nuse ={:d} out of nlines = {:d} total lines'.format(nuse,nlines) +  msgs.newline() +
                  'We are will proceed without rejecting anything but something is probably very wrong with this slit/order.')
        use_tilt = np.ones(nlines,dtype=bool)
        nuse = nlines

    # Tighten it up with Gaussian weighted centroiding
    trc_tilt_dict = dict(nspec = nspec, nspat = nspat, nsub = nsub, nlines = nlines, nuse = nuse,
                         spat_min=spat_min, spat_max=spat_max, do_crude=do_crude, fwhm = fwhm, use_tilt=use_tilt,
                         tilts_sub_spat=tilts_sub_spat, tilts_sub_fit=tilts_sub_fit,
                         tilts_mad = tilts_mad,tilts_spec=tilts_spec, tilts_spat=tilts_spat, tilts_dspat=tilts_dspat,
                         tilts=tilts, tilts_fit=tilts_fit, tilts_err=tilts_err, tilts_mask=tilts_mask)

    return trc_tilt_dict


def trace_tilts(arcimg, lines_spec, lines_spat, thismask, slit_cen, inmask=None, gauss=False, fwhm=4.0,spat_order=5, maxdev_tracefit=0.2,
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
       Maximum absolute deviation for the arc tilt fits during iterative trace fitting expressed in units of the fwhm.
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
        arcimg, lines_spec, lines_spat, thismask, slit_cen, inmask=inmask, gauss=gauss, tilts_guess=None, fwhm=fwhm, spat_order=spat_order,
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
    trace_dict1 = trace_tilts_work(arcimg, lines_spec, lines_spat, thismask, slit_cen, inmask=inmask, gauss=gauss, tilts_guess=pca_fit,
                                      fwhm=fwhm, spat_order=spat_order, maxdev_tracefit=maxdev_tracefit,sigrej_trace=sigrej_trace,
                                      max_badpix_frac=max_badpix_frac,show_tracefits=show_tracefits)

    return trace_dict1


def fit_tilts(trc_tilt_dict, thismask, slit_cen, spat_order=3, spec_order=4, maxdev = 0.2,
              maxrej = None,
              maxiter = 100, sigrej = 3.0, pad_spec = 30, pad_spat =5,
              func2d='legendre2d', doqa=True, master_key='test',
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

    nspec = trc_tilt_dict['nspec']
    nspat = trc_tilt_dict['nspat']
    fwhm = trc_tilt_dict['fwhm']
    maxdev_pix = maxdev*fwhm
    xnspecmin1 = float(nspec-1)
    xnspatmin1 = float(nspat-1)
    nspat = trc_tilt_dict['nspat']
    use_tilt = trc_tilt_dict['use_tilt']                 # mask for good/bad tilts, based on aggregate fit, frac good pixels
    nuse = np.sum(use_tilt)
    tilts = trc_tilt_dict['tilts']   # legendre polynomial fit
    #JFH Before we were fitting the fits. Now we fit the actual flux weighted centroided tilts.
    tilts_err = trc_tilt_dict['tilts_err']   # flux weighted centroidding error
    tilts_dspat = trc_tilt_dict['tilts_dspat'] # spatial offset from the central trace
    #tilts_spat = trc_tilt_dict['tilts_dspat'][:,use_tilt] # spatial offset from the central trace
    tilts_spec = trc_tilt_dict['tilts_spec'] # line spectral pixel position from legendre fit evaluated at slit center
    tilts_mask = trc_tilt_dict['tilts_mask'] # Reflects if trace is on the slit
    tilts_mad = trc_tilt_dict['tilts_mad']   # quantitfies aggregate error of this tilt

    use_mask = np.outer(np.ones(nspat,dtype=bool),use_tilt)
    tot_mask = tilts_mask & (tilts_err < 900) & use_mask
    fitxy = [spec_order, spat_order]

    # Fit the inverted model with a 2D polynomial
    msgs.info("Fitting tilts with a low order, 2D {:s}".format(func2d))

    adderr = 0.03
    tilts_sigma = ((tilts_mad < 100.0) & (tilts_mad > 0.0))*np.sqrt(np.abs(tilts_mad)**2 + adderr**2)

    fitmask, coeff2 = utils.robust_polyfit_djs(tilts_spec.flatten()/xnspecmin1, (tilts.flatten() - tilts_spec.flatten())/xnspecmin1,
                                               fitxy, x2=tilts_dspat.flatten()/xnspatmin1, inmask = tot_mask.flatten(),
                                               sigma=tilts_sigma.flatten()/xnspecmin1,
                                               function=func2d, maxiter=maxiter, lower=sigrej, upper=sigrej,
                                               maxdev=maxdev_pix/xnspecmin1,minx=-0.0, maxx=1.0, minx2=-1.0, maxx2=1.0,
                                               use_mad=False, sticky=False)
    fitmask = fitmask.reshape(tilts_dspat.shape)
    # Compute a rejection  mask that we will use later. These are locations that were fit but were rejected
    rej_mask = tot_mask & np.invert(fitmask)
    # Compute and store the 2d tilts fit
    delta_tilt_1 = xnspecmin1*utils.func_val(coeff2, tilts_spec[tilts_mask]/xnspecmin1, func2d, x2=tilts_dspat[tilts_mask]/xnspatmin1,
                                            minx=0.0, maxx=1.0, minx2=-1.0, maxx2=1.0)
    delta_tilt = np.zeros_like(tilts_dspat)
    tilts_2dfit = np.zeros_like(tilts_dspat)
    delta_tilt[tilts_mask] = delta_tilt_1
    tilts_2dfit[tilts_mask] = tilts_spec[tilts_mask] + delta_tilt[tilts_mask]
    # Add the 2d fit to the tracetilt dictionary
    trc_tilt_dict_out = copy.deepcopy(trc_tilt_dict)
    trc_tilt_dict_out['tilt_2dfit'] = tilts_2dfit

    # Report the residuals in pixels
    res_fit = tilts[fitmask] - tilts_2dfit[fitmask]
    rms_fit = np.std(res_fit)
    msgs.info("Residuals: 2D Legendre Fit")
    msgs.info("RMS (pixels): {}".format(rms_fit))
    msgs.info("RMS/FWHM: {}".format(rms_fit/fwhm))

    msgs.info('Inverting the fit to generate the tilts image')
    spec_vec = np.arange(nspec)
    spat_vec = np.arange(nspat)
    spat_img, spec_img = np.meshgrid(spat_vec, spec_vec)
    # We do some padding here to guarantee that the tilts arc lines falling off the image get tilted onto the image
    spec_vec_pad = np.arange(-pad_spec, nspec + pad_spec)
    spat_vec_pad = np.arange(-pad_spat, nspat + pad_spat)
    spat_img_pad, spec_img_pad = np.meshgrid(np.arange(-pad_spat, nspat + pad_spat),np.arange(-pad_spec, nspec + pad_spec))
    slit_cen_pad = (scipy.interpolate.interp1d(spec_vec, slit_cen, bounds_error=False, fill_value='extrapolate'))(spec_vec_pad)
    thismask_pad = np.zeros_like(spec_img_pad, dtype=bool)
    ind_spec, ind_spat = np.where(thismask)
    slit_cen_img_pad = np.outer(slit_cen_pad, np.ones(nspat + 2 * pad_spat))  # center of the slit replicated spatially
    # Normalized spatial offset image (from central trace)
    dspat_img_nrm = (spat_img_pad - slit_cen_img_pad) / xnspatmin1
    # normalized spec image
    spec_img_nrm = spec_img_pad / xnspecmin1
    # Embed the old thismask in the new larger padded thismask
    thismask_pad[ind_spec + pad_spec, ind_spat + pad_spat] = thismask[ind_spec, ind_spat]
    # Now grow the thismask_pad
    kernel = np.ones((2*pad_spec, 2*pad_spat)) / float(4 * pad_spec * pad_spat)
    thismask_grow = scipy.ndimage.convolve(thismask_pad.astype(float), kernel, mode='nearest') > 0.0
    # Evaluate the tilts on the padded image grid
    tiltpix = spec_img_pad[thismask_grow] + xnspecmin1 * utils.func_val(coeff2, spec_img_nrm[thismask_grow], func2d,
                                                                         x2=dspat_img_nrm[thismask_grow],
                                                                         minx=0.0, maxx=1.0, minx2=-1.0, maxx2=1.0)
    # Now do one last fit to invert the function above to obtain the final tilts model in normalized image coordinates
    inmask = np.isfinite(tiltpix)
    sigma = np.full_like(spec_img_pad, 10.0)
    # JFH What I find confusing is that this last fit was actually what Burles was doing on the raw tilts, so why was that failing?
    fitmask_tilts, coeff2_tilts = utils.robust_polyfit_djs(tiltpix/xnspecmin1, spec_img_pad[thismask_grow]/xnspecmin1,
                                                           fitxy, x2=spat_img_pad[thismask_grow]/xnspatmin1,
                                                           sigma=sigma[thismask_grow]/xnspecmin1,
                                                           upper=5.0, lower=5.0, maxdev=10.0/xnspecmin1,
                                                           inmask=inmask, function=func2d, maxiter=20,
                                                           minx=0.0, maxx=1.0, minx2=0.0, maxx2=1.0, use_mad=False)
    irej = np.invert(fitmask_tilts) & inmask
    msgs.info('Rejected {:d}/{:d} pixels in final inversion tilts image fit'.format(np.sum(irej),np.sum(inmask)))
    # normalized tilts image
    tilts_img = utils.func_val(coeff2_tilts, spec_img/xnspecmin1, func2d, x2=spat_img/xnspatmin1,minx=0.0, maxx=1.0, minx2=0.0, maxx2=1.0)
    tilts_img = np.fmax(np.fmin(tilts_img, 1.2),-0.2)

    tilt_fit_dict = dict(nspec = nspec, nspat = nspat, ngood_lines=np.sum(use_tilt), npix_fit = np.sum(tot_mask),
                         npix_rej = np.sum(fitmask == False), coeff2=coeff2_tilts, spec_order = spec_order, spat_order = spat_order,
                         minx = 0.0, maxx = 1.0, minx2 = 0.0, maxx2 = 1.0, func=func2d)

    # Now do some QA
    if doqa:
        plot_tilt_2d(tilts_dspat, tilts, tilts_2dfit, tot_mask, rej_mask, spat_order, spec_order, rms_fit, fwhm,
                     slit=slit, setup=master_key, show_QA=show_QA, out_dir=out_dir)
        plot_tilt_spat(tilts_dspat, tilts, tilts_2dfit, tilts_spec, tot_mask, rej_mask, spat_order, spec_order, rms_fit, fwhm,
                       slit=slit, setup=master_key, show_QA=show_QA, out_dir=out_dir)
        plot_tilt_spec(tilts_spec, tilts, tilts_2dfit, tot_mask, rej_mask, rms_fit, fwhm, slit=slit,
                       setup = master_key, show_QA=show_QA, out_dir=out_dir)

    return tilts_img, tilt_fit_dict, trc_tilt_dict_out


    #fitmask, coeff2 = fit_tilts_rej(
    #    tilts_dspat, tilts_spec_fit, tilts, tilts_invvar, tot_mask, slit_cen, spat_order, spec_order,
    #    maxdev = 1.0, maxrej=maxrej, sigrej = sigrej, maxiter = maxiter)


    #result = optimize.minimize(fit_tilts_func, coeff2.flatten(), tol=0.01, args=(tilts_dspat, tilts_spec_fit, tilts, tilts_invvar,
    #                                                                             tot_mask, fitmask, slit_cen, spat_order, spec_order))
    #bounds = [(i,j) for i,j in zip(0.8*coeff2.flatten(),1.2*coeff2.flatten())]
    #result_df = optimize.differential_evolution(
    #    fit_tilts_func, bounds, tol=0.01, disp=True, polish=True,
    #    args=(tilts_dspat, tilts_spec_fit, tilts, tilts_invvar,tot_mask, fitmask, slit_cen, spat_order, spec_order))



    # This is a second way of computing the 2d fit at the tilts. Do this just for now as a consistency check that our tilts_img is okay

    # Testing. For the moment do this on exactly the same set of lines
    #tilt_2dfit_piximg_all = eval_2d_at_tilts(trc_tilt_dict['tilts_spec'], trc_tilt_dict['tilts_mask'], trc_tilt_dict['tilts'] (nspec, nspat), thismask, slit_cen, coeff2, func2d)
    #tilts_2dfit_piximg = tilt_2dfit_piximg_all[:, use_tilt]
    #tilts_2dfit_piximg = eval_2d_at_tilts(tilts_spec, tilts_mask, (nspec, nspat), thismask, slit_cen, coeff2, func2d)

    # Actual 2D Model Tilt Residuals
    #res_real = tilts[fitmask] - tilts_2dfit_piximg[fitmask]
    #rms_real = np.std(res_real)
    #msgs.info("Residuals: Actual 2D Tilt Residuals from piximg")
    #msgs.info("RMS (pixels): {}".format(rms_real))
    #msgs.info("RMS/FWHM: {}".format(rms_real/fwhm))



# JFH This is now defunct
def gridtilts(shape, thismask, slit_cen, coeff2, func2d, spec_order, spat_order, pad_spec=30, pad_spat = 5, method='interp'):
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
    spec_vec = np.arange(nspec)
    spat_vec = np.arange(nspat)

    # JFH This histogram method is not preferred, since it basically does NGP. It is however super fast, so for big images
    # it is useful to have it
    if 'hist2d' in method:
        oversamp_spec=5
        oversamp_spat=3
        spec_ind, spat_ind = np.where(thismask)
        min_spec = spec_ind.min() - pad_spec
        max_spec = spec_ind.max() + pad_spec
        num_spec = max_spec - min_spec + 1

        min_spat = spat_ind.min() - pad_spat
        max_spat = spat_ind.max() + pad_spat
        num_spat = max_spat - min_spat + 1

        spec_lin = np.linspace(min_spec,max_spec,num = int(np.round(num_spec*oversamp_spec)))
        spat_lin = np.linspace(min_spat,max_spat,num = int(np.round(num_spat*oversamp_spat)))
        spat_img, spec_img = np.meshgrid(spat_lin, spec_lin)
        # Normalized spatial offset image (from central trace)
        slit_cen_lin = (scipy.interpolate.interp1d(np.arange(nspec),slit_cen,bounds_error=False,fill_value='extrapolate'))(spec_lin)
        slit_cen_img = np.outer(slit_cen_lin, np.ones(spat_img.shape[1]))  # center of the slit replicated spatially
        dspat_img_nrm = (spat_img - slit_cen_img)/xnspatmin1
        spec_img_nrm = spec_img/xnspecmin1

        # normalized spec image
        tracepix = spec_img + xnspecmin1*utils.func_val(coeff2, spec_img_nrm, func2d, x2=dspat_img_nrm,
                                                              minx=0.0, maxx=1.0, minx2=-1.0, maxx2=1.0)
        norm_img, spec_edges, spat_edges = np.histogram2d(tracepix.flatten(), spat_img.flatten(),
                                                          bins=[np.arange(nspec+1), np.arange(nspat+1)], density=False)
        weigh_img, spec_edges, spat_edges = np.histogram2d(tracepix.flatten(), spat_img.flatten(),
                                                           bins=[np.arange(nspec+1), np.arange(nspat+1)],
                                                           weights = spec_img.flatten(),density=False)
        piximg =(norm_img > 0.0)*weigh_img/(norm_img + (norm_img == 0.0))
        inmask = thismask & (norm_img > 0) & (piximg/xnspecmin1 > -0.2) & (piximg/xnspecmin1 < 1.2)
    # This is the defulat method although scipy.interpolate.griddata is a bit slow
    elif 'interp' in method:
        spec_vec_pad = np.arange(-pad_spec,nspec+pad_spec)
        spat_vec_pad = np.arange(-pad_spat,nspat+pad_spat)
        spat_img, spec_img = np.meshgrid(spat_vec, spec_vec)
        spat_img_pad, spec_img_pad = np.meshgrid(np.arange(-pad_spat,nspat+pad_spat),np.arange(-pad_spec,nspec+pad_spec))
        slit_cen_pad = (scipy.interpolate.interp1d(spec_vec,slit_cen,bounds_error=False,fill_value='extrapolate'))(spec_vec_pad)
        thismask_pad = np.zeros_like(spec_img_pad,dtype=bool)
        ind_spec, ind_spat = np.where(thismask)
        slit_cen_img_pad= np.outer(slit_cen_pad, np.ones(nspat + 2*pad_spat))  # center of the slit replicated spatially
        # Normalized spatial offset image (from central trace)
        dspat_img_nrm = (spat_img_pad - slit_cen_img_pad)/xnspatmin1
        # normalized spec image
        spec_img_nrm = spec_img_pad/xnspecmin1
        # Embed the old thismask in the new larger padded thismask
        thismask_pad[ind_spec + pad_spec,ind_spat + pad_spat] = thismask[ind_spec,ind_spat]
        # Now grow the thismask_pad
        kernel = np.ones((2*pad_spec, 2*pad_spat))/float(4*pad_spec*pad_spat)
        thismask_grow = scipy.ndimage.convolve(thismask_pad.astype(float), kernel, mode='nearest') > 0.0
        # Evaluate the tilts on the padded image grid
        tracepix = spec_img_pad[thismask_grow] + xnspecmin1*utils.func_val(coeff2, spec_img_nrm[thismask_grow], func2d, x2=dspat_img_nrm[thismask_grow],
                                                              minx=0.0, maxx=1.0, minx2=-1.0, maxx2=1.0)
        ## TESTING STARTS
        """
        ikeep = np.isfinite(tracepix)
        sigma = np.full_like(spec_img_pad[thismask_grow], 10.0)/xnspecmin1
        fitxy = [spec_order, spat_order]
        fitmask, coeff2_tilts = utils.robust_polyfit_djs(tracepix/xnspecmin1, spec_img_pad[thismask_grow]/xnspecmin1,
                                                         fitxy, x2=spat_img_pad[thismask_grow]/xnspatmin1,
                                                         sigma=sigma,
                                                         upper=5.0, lower=5.0, maxdev=10.0/xnspecmin1,
                                                         inmask=ikeep, function=func2d, maxiter=20,
                                                         minx=0.0, maxx=1.0, minx2=0.0, maxx2=1.0, use_mad=False)
        ## TESTING ENDS
        # values(points) \equiv spec_pos(tilt,spat_pos) which is the piximg that we want to create via griddata interpolation
        """
        ikeep = np.isfinite(tracepix) 
        points = np.stack((tracepix[ikeep], spat_img_pad[thismask_grow][ikeep]), axis=1)
        values =spec_img_pad[thismask_grow][ikeep]
        piximg = scipy.interpolate.griddata(points, values, (spec_img, spat_img), method='cubic')
        inmask = thismask & np.isfinite(piximg) & (piximg/xnspecmin1 > -0.2) & (piximg/xnspecmin1 < 1.2)

    # Now simply do a 2d polynomial fit with just rejection of crazy behavior, i.e. 10 pixels
    fitxy = [spec_order, spat_order]
    sigma = np.full_like(spec_img,10.0)/xnspecmin1
    fitmask, coeff2_tilts = utils.robust_polyfit_djs(spec_img.flatten()/xnspecmin1, piximg.flatten()/xnspecmin1,
                                                     fitxy, x2=spat_img.flatten()/xnspatmin1, sigma = sigma.flatten(),
                                                     upper=5.0, lower=5.0, maxdev = 10.0/xnspecmin1,
                                                     inmask=inmask.flatten(), function=func2d, maxiter=20,
                                                     minx=0.0, maxx=1.0, minx2=0.0,maxx2=1.0,use_mad=False)
    irej = np.invert(fitmask) & inmask.flatten()
    msgs.info('Rejected {:d}/{:d} pixels in final tilts image after gridding'.format(np.sum(irej),np.sum(inmask)))
    # normalized tilts image

    tilts = utils.func_val(coeff2_tilts, spec_img/xnspecmin1, func2d, x2=spat_img/xnspatmin1,minx=0.0, maxx=1.0, minx2=0.0, maxx2=1.0)
    tilts = np.fmax(np.fmin(tilts, 1.2),-0.2)
    # Added this to ensure that tilts are never crazy values due to extrapolation of fits which can break
    # wavelength solution fitting

    return coeff2_tilts, tilts



def fit2tilts(shape, coeff2, func2d):
    """

    Parameters
    ----------
    shape: tuple of ints,
        shape of image
    coeff2: ndarray, float
        result of griddata tilt fit
    func2d: str
        the 2d function used to fit the tilts
    Returns
    -------
    tilts: ndarray, float
       Image indicating how spectral pixel locations move across the image. This output is used in the pipeline.
    """

    # Compute the tilts image
    nspec, nspat = shape
    xnspecmin1 = float(nspec-1)
    xnspatmin1 = float(nspat-1)
    spec_vec = np.arange(nspec)
    spat_vec = np.arange(nspat)
    spat_img, spec_img = np.meshgrid(spat_vec, spec_vec)
    tilts = utils.func_val(coeff2, spec_img/xnspecmin1, func2d, x2=spat_img/xnspatmin1, minx=0.0, maxx=1.0, minx2=0.0, maxx2=1.0)
    # Added this to ensure that tilts are never crazy values due to extrapolation of fits which can break
    # wavelength solution fitting
    tilts = np.fmax(np.fmin(tilts, 1.2),-0.2)
    return tilts



def fit2tilts_backup(shape, slit_cen, coeff2, func):
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

    tilts = utils.func_val(coeff2, dspat_img_nrm, func, x2=spec_img_nrm, minx=-1.0, maxx=1.0, minx2=0.0, maxx2=1.0)
    # Added this to ensure that tilts are never crazy values due to extrapolation of fits which can break
    # wavelength solution fitting
    tilts = np.fmax(np.fmin(tilts, 1.2),-0.2)

    return tilts

def plot_tilt_2d(tilts_dspat, tilts, tilts_model, tot_mask, rej_mask, spat_order, spec_order, rms, fwhm,
                 slit=0, setup='A', outfile=None, show_QA=False, out_dir=None):


    plt.rcdefaults()
    plt.rcParams['font.family']= 'Helvetica'

    # Outfile
    method = inspect.stack()[0][3]
    if (outfile is None):
        outfile = qa.set_qa_filename(setup, method, slit=slit, out_dir=out_dir)


    # Show the fit
    fig, ax = plt.subplots(figsize=(12, 18))
    ax.cla()
    ax.plot(tilts_dspat[tot_mask], tilts[tot_mask], color='black', linestyle=' ', mfc='None', marker='o',
            markersize=9.0, markeredgewidth=1.0, zorder=4, label='Good Tilt')
    ax.plot(tilts_dspat[rej_mask], tilts[rej_mask], color='red', linestyle=' ', mfc='None', marker='o',
            markersize=9.0, markeredgewidth=2.0, zorder=5, label='Rejected')
    ax.plot(tilts_dspat[tot_mask], tilts_model[tot_mask], color='black', linestyle=' ', marker='o',
            markersize=2.0, markeredgewidth=1.0, zorder=1, label='2D Model')

    xmin = 1.1 * tilts_dspat[tot_mask].min()
    xmax = 1.1 * tilts_dspat[tot_mask].max()
    ax.set_xlim((xmin, xmax))
    ax.set_xlabel('Spatial Offset from Central Trace (pixels)', fontsize=15)
    ax.set_ylabel('Spectral Pixel', fontsize=15)
    ax.legend()
    ax.set_title('Tilts vs Fit (spat_order, spec_order)=({:d},{:d}) for slit={:d}: RMS = {:5.3f}, '
                 'RMS/FWHM={:5.3f}'.format(spat_order, spec_order, slit, rms, rms / fwhm), fontsize=15)

    # Finish
    #plt.tight_layout(pad=1.0, h_pad=1.0, w_pad=1.0)

    if outfile is not None:
        plt.savefig(outfile, dpi=400)

    if show_QA:
        plt.show()

    plt.close()
    plt.rcdefaults()



def plot_tilt_spec(tilts_spec_fit, tilts, tilts_model, tot_mask, rej_mask, rms, fwhm,
                   slit=0, setup = 'A', outfile=None, show_QA=False, out_dir=None):
    """ Generate a QA plot of the residuals for the fit to the tilts in the spectral direction one slit at a time

    Parameters
    ----------
    """

    plt.rcdefaults()
    plt.rcParams['font.family']= 'Helvetica'

    # Outfil
    method = inspect.stack()[0][3]
    if (outfile is None):
        outfile = qa.set_qa_filename(setup, method, slit=slit, out_dir=out_dir)

    # Setup
    plt.figure(figsize=(14, 6))
    plt.clf()
    ax = plt.gca()

    # Scatter plot
    res = (tilts - tilts_model)

    nspat, nuse = tilts.shape
    # Show the fit residuals as a function of spatial position
    line_indx = np.outer(np.ones(nspat), np.arange(nuse))

    xmin = 0.90*(tilts_spec_fit.min())
    xmax = 1.10*(tilts_spec_fit.max())

    ax.hlines(0.0, xmin, xmax,linestyle='--', color='green')

    for iline in range(nuse):
        iall = (line_indx == iline) & tot_mask
        igd = (line_indx == iline) & tot_mask & (rej_mask == False)
        irej = (line_indx == iline) & tot_mask & rej_mask

        ax.plot(tilts_spec_fit[igd], (res[igd]), 'ko', mfc='k', markersize=4.0)
        ax.plot(tilts_spec_fit[irej],(res[irej]), 'ro', mfc='r', markersize=4.0)
        # Compute the RMS for this line
        all_rms = np.std(res[iall])
        good_rms = np.std(res[igd])
        # ToDo show the mean here as well
        if np.any(igd):
            ax.plot(tilts_spec_fit[igd][0], all_rms, marker='s',linestyle=' ', color='g', mfc='g', markersize=7.0)
            ax.plot(tilts_spec_fit[igd][0], good_rms, marker='^', linestyle=' ', color='orange', mfc='orange', markersize=7.0)

    ax.text(0.90, 0.90, 'Slit {:d}:  Residual (pixels) = {:0.5f}'.format(slit, rms),
            transform=ax.transAxes, size='large', ha='right', color='black',fontsize=16)
    ax.text(0.90, 0.80, ' Slit {:d}:  RMS/FWHM = {:0.5f}'.format(slit, rms/fwhm),
            transform=ax.transAxes, size='large', ha='right', color='black',fontsize=16)
    # Label
    ax.set_xlabel('Spectral Pixel')
    ax.set_ylabel('RMS (pixels)')
    ax.set_title('RMS of Each Arc Line Traced')
    ax.set_xlim((xmin,xmax))
    ax.set_ylim((-5.0*rms,5.0*rms))
    # Legend
    legend_elements = [Line2D([0], [0],linestyle=' ', color='k', marker='o', mfc='k', markersize=4.0, label='good'),
                       Line2D([0], [0], linestyle=' ', color='r', marker='o', mfc='r', markersize=4.0, label='rejected'),
                       Line2D([0], [0], linestyle=' ', color='g', marker='s', mfc='g', markersize=7.0, label='all RMS'),
                       Line2D([0], [0], linestyle=' ', color='orange', marker='^', mfc='orange', markersize=7.0, label='good RMS')]
    ax.legend(handles=legend_elements)

    # Finish
    plt.tight_layout(pad=0.2, h_pad=0.0, w_pad=0.0)

    if outfile is not None:
        plt.savefig(outfile, dpi=400)

    if show_QA:
        plt.show()

    plt.close()
    plt.rcdefaults()



def plot_tilt_spat(tilts_dspat, tilts, tilts_model, tilts_spec_fit, tot_mask, rej_mask,spat_order, spec_order, rms, fwhm,
                   setup='A', slit=0, outfile=None, show_QA=False, out_dir=None):


    import matplotlib as mpl
    from matplotlib.lines import Line2D

    plt.rcdefaults()
    plt.rcParams['font.family']= 'Helvetica'

    # Outfil
    method = inspect.stack()[0][3]
    if (outfile is None):
        outfile = qa.set_qa_filename(setup, method, slit=slit, out_dir=out_dir)

    nspat, nuse = tilts_dspat.shape
    # Show the fit residuals as a function of spatial position
    line_indx = np.outer(np.ones(nspat), np.arange(nuse))
    lines_spec = tilts_spec_fit[0, :]
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
        ax.plot(tilts_dspat[iall], tilts[iall] - tilts_model[iall], color=this_color,
                linestyle='-', linewidth=3.0, marker='None', alpha=0.5)
        ax.plot(tilts_dspat[irej], tilts[irej] - tilts_model[irej], linestyle=' ',
                marker='o', color='limegreen', mfc='limegreen', markersize=5.0)

    xmin = 1.1 * tilts_dspat[tot_mask].min()
    xmax = 1.1 * tilts_dspat[tot_mask].max()
    ax.hlines(0.0, xmin, xmax, linestyle='--', linewidth=2.0, color='k', zorder=10)

    ax.set_xlim((xmin, xmax))
    ax.set_xlabel('Spatial Offset from Central Trace (pixels)')
    ax.set_ylabel('Arc Line Tilt Residual (pixels)')

    legend_elements = [Line2D([0], [0], color='cornflowerblue', linestyle='-', linewidth=3.0, label='residual'),
                       Line2D([0], [0], color='limegreen', linestyle=' ', marker='o', mfc='limegreen', markersize=7.0,
                              label='rejected')]
    ax.legend(handles=legend_elements)
    ax.set_title('Tilts vs Fit (spat_order, spec_order)=({:d},{:d}) for slit={:d}: RMS = {:5.3f}, '
                 'RMS/FWHM={:5.3f}'.format(spat_order, spec_order, slit, rms, rms / fwhm), fontsize=15)
    cb = fig.colorbar(dummie_cax, ticks=lines_spec)
    cb.set_label('Spectral Pixel')

    # Finish
    plt.tight_layout(pad=0.2, h_pad=0.0, w_pad=0.0)

    if outfile is not None:
        plt.savefig(outfile, dpi=400)

    if show_QA:
        plt.show()

    plt.close()
    plt.rcdefaults()
