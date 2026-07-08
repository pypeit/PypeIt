"""
Refactored spatial profile fitting for optimal extraction.

This module contains :func:`fit_profile_refactor`, a drop-in replacement for
:func:`~pypeit.core.spatialprofile.fit_profile` with the following improvements:

- **Bug 1 fixed**: ``nb = np.round(prof_nsigma > 10)`` replaced with
  ``nb = max(1, round(prof_nsigma / 10))``.
- **Bug 2 fixed**: :func:`_findfwhm` uses linear interpolation for sub-pixel
  accuracy at the half-maximum crossing (was a discrete-sample overestimate).
- **Accuracy**: :func:`_return_gaussian` uses the pixel-integrated
  error-function form instead of a point-sampled Gaussian (important for
  FWHM < 5 px; row sums of the returned profile are now ~1.0, matching the
  B-spline path).
- **Performance**: all ``np.outer(v, np.ones(n))`` replaced with ``v[:, None]``
  broadcasting; the two apodization-limit ``while True`` loops replaced with
  vectorised searches over pre-computed derivative arrays.

The legacy :func:`~pypeit.core.spatialprofile.fit_profile` in
``spatialprofile.py`` is left unchanged so both implementations can be compared
directly.

.. include:: ../include/links.rst
"""

import astropy.stats
from IPython import embed
import numpy as np
import scipy.interpolate
import scipy.ndimage
import scipy.special

from pypeit import log
from pypeit import utils
from pypeit.core.fitting import iterative_bspline_fit
from pypeit.core import pydl
from pypeit.core.spatialprofile import qa_fit_profile


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


def _inverse(arr):
    """
    Calculate the inverse of the provided array.

    This is virtually identical to :func:`~pypeit.utils.inverse`, except that it
    does not force the result to be positive.

    Parameters
    ----------
    arr : :class:`numpy.ndarray`
        Input array

    Returns
    -------
    :class:`numpy.ndarray`
        Calculation of ``1/arr`` that controls for where ``arr`` is identically
        0.
    """
    gpm = arr != 0.0
    return gpm / (arr + np.logical_not(gpm))


def _findfwhm(model, sig_x):
    """
    Calculate the FWHM of a profile with sub-pixel precision.

    Locates the half-maximum crossing on each side of the peak by linear
    interpolation between the two bracketing sample points, eliminating the
    systematic overestimate that arises when the crossing falls between
    adjacent discrete samples (see Bug 2 in the refactoring plan).

    Parameters
    ----------
    model : :class:`numpy.ndarray`
        1-D model profile
    sig_x : :class:`numpy.ndarray`
        Spatial coordinates normalized such that the nominal center is at 0 and
        the coordinates are in units of a previously estimated 1-sigma width.
        Shape must match ``model``.

    Returns
    -------
    peak : :obj:`float`
        Peak value of the profile.
    peak_x : :obj:`float`
        ``sig_x`` value at the peak.
    lwhm : :obj:`float`
        ``sig_x`` value at the left half-maximum (sub-pixel accuracy).
    rwhm : :obj:`float`
        ``sig_x`` value at the right half-maximum (sub-pixel accuracy).
    """
    # Mask pixels beyond 2 sigma
    _model = np.ma.MaskedArray(model, mask=np.abs(sig_x)>2.)
    # Find the peak coordinate and value
    peak_i = np.ma.argmax(_model)
    peak = _model.data[peak_i]
    peak_x = sig_x[peak_i]
    # Mask all the values less than half of the peak
    _model[_model < 0.5 * peak] = np.ma.masked
    # Get the indices of the unmasked pixels that bracket the peak
    # WARNING: This *assumes* there is only one coherent section with values
    # above 0.5 * peak and within |sig_x| < 2.
    lind, rind = np.ma.flatnotmasked_edges(_model)
    # Get the left edge of the FWHM range
    if lind > 0:
        lwhm = utils.linear_interpolate(
            model[lind-1], sig_x[lind-1], model[lind], sig_x[lind], 0.5*peak
        )
    else:
        lwhm = sig_x[0]
    # Get the right edge of the FWHM range
    if rind < model.size-1:
        rwhm = utils.linear_interpolate(
            model[rind], sig_x[rind], model[rind+1], sig_x[rind+1], 0.5*peak
        )
    else:
        rwhm = sig_x[-1]

    return peak, peak_x, lwhm, rwhm


def _gaussian_profile(x, center, sigma):
    r"""
    Construct a pixel-integrated Gaussian spatial profile.

    Computes the fraction of flux in each pixel as the integral of a Gaussian
    across the pixel width using the error function, rather than evaluating the
    Gaussian density at the pixel centre.  The correction is significant for
    narrow profiles: for FWHM = 3 px the peak value is underestimated by ~10 %
    with the point-sampled form; the integrated form reduces this to < 1 %.

    The row sums of the returned profile are approximately 1.0 (the erf
    differences telescope to 1 for wide enough slits), matching the B-spline
    path.

    .. important::

        All inputs *must* be in pixels and the coordinate (`x`) array must
        provide spatial coordinates that increase by 1 per provided array
        element.  That is, the function assumes that the spatial step in `x` is
        always 1.

    Parameters
    ----------
    x : array-like
        Pixel coordinate array.  May be 1-D with shape
        :math:`(N_{\rm spat},)` or 2-D with shape
        :math:`(N_{\rm spec}, N_{\rm spat})`.  For the 2-D case, coordinates
        are expected to vary primarily along the second axis.
    center : :obj:`float`, array-like
        Profile centre in pixels.  Must be a scalar if ``x`` is 1-D.  If
        ``x`` is 2-D, may be a scalar or a 1-D array of length
        :math:`N_{\rm spec}`.
    sigma : :obj:`float`, array-like
        Profile standard deviation in pixels.  Same shape constraints as
        ``center``.

    Returns
    -------
    profile : :class:`numpy.ndarray`
        Pixel-integrated Gaussian profile, same shape as ``x``.
    """
    _x = np.asarray(x)
    if _x.ndim == 2:
        _center = np.atleast_1d(center)
        _sigma = np.atleast_1d(sigma)
        if _center.size > 1:
            _center = _center[:, None]
        if _sigma.size > 1:
            _sigma = _sigma[:, None]
    else:
        _center, _sigma = center, sigma
    s = _inverse(_sigma * np.sqrt(2.0))
    sigma_x = (_x - _center) * s
    delta = 0.5 * s
    profile = 0.5 * (
        scipy.special.erf(sigma_x + delta) - scipy.special.erf(sigma_x - delta)
    )
    # Set pixels beyond 5 sigma to 0.
    profile[sigma_x**2 >= 25/2] = 0.0
    invalid = np.logical_not(np.isfinite(profile))
    if np.any(invalid):
        log.warning('Setting NaN/Inf pixel values in object profile to zero.')
        profile[invalid] = 0.0
    return profile


def _return_gaussian(
    spat_img, norm_obj, center, sigma, fwhm, med_sn2, obj_string, show_profile, ind=None,
    l_limit=None, r_limit=None, xlim=None, xtrunc=1e6
):
    r"""
    Build a Gaussian profile and log fit information and optionally show QA.

    Delegates profile construction to :func:`_gaussian_profile`, then logs the
    FWHM and S/N and optionally displays the QA plot via
    :func:`~pypeit.core.spatialprofile.qa_fit_profile`.

    Parameters
    ----------
    spat_img : :class:`numpy.ndarray`
        Spatial-coordinate image in pixels, shape
        :math:`(N_{\rm spec}, N_{\rm spat})`.
    norm_obj : :class:`numpy.ndarray`
        Normalised object image, used only for the QA plot.
    center : :obj:`float` or :class:`numpy.ndarray`
        Object centre in pixels; scalar or shape :math:`(N_{\rm spec},)`.
    sigma : :obj:`float` or :class:`numpy.ndarray`
        Profile standard deviation in pixels; same shape constraints as
        ``center``.
    fwhm : :obj:`float`
        FWHM for the log message in pixels.
    med_sn2 : :obj:`float`
        Median (S/N)^2, used only for the QA plot label.
    obj_string : :obj:`str`
        Object identifier for the QA plot title.
    show_profile : :obj:`bool`
        If ``True``, display the QA plot.
    ind : :class:`numpy.ndarray`, optional
        Flat indices of the good pixels for the QA plot.
    l_limit, r_limit : :obj:`float`, optional
        Profile limits drawn on the QA plot.
    xlim : :obj:`float`, optional
        x-axis half-range for the QA plot.
    xtrunc : :obj:`float`, optional
        Truncation radius in sigma for the QA plot.

    Returns
    -------
    profile_model : :class:`numpy.ndarray`
        Pixel-integrated Gaussian profile, shape
        :math:`(N_{\rm spec}, N_{\rm spat})`.
    """
    profile_model = _gaussian_profile(spat_img, center, sigma)
    title_string = f'{obj_string}, FWHM={fwhm:.2f}, S/N={np.sqrt(med_sn2):.3f}'
    log.info(title_string)
    if not show_profile:
        return profile_model

    _center = np.atleast_1d(center)
    _sigma = np.atleast_1d(sigma)
    if _center.size > 1:
        _center = _center[:, None]
    if _sigma.size > 1:
        _sigma = _sigma[:, None]
    sigma_x = (spat_img - _center) / _sigma
    qa_fit_profile(
        sigma_x, norm_obj, profile_model, title=title_string, l_limit=l_limit, r_limit=r_limit,
        ind=ind, xlim=xlim, xtrunc=xtrunc,
    )
    return profile_model


def _fit_spectrum_and_normalize(
    wave, flux, fluxivar, waveimg, image, ivar, totmask, spec_img, percentile_sn2, fwhm,
    sn_cap=100.0, good_pix_frac=0.05
):
    r"""
    Fit flux B-splines, estimate S/N, and build the normalised-object image.

    Encapsulates the spectral-fitting and normalisation steps of
    :func:`fit_profile_refactor` (blocks 0 and 2–4).  Constructs the 1-D
    slit-pixel arrays from ``totmask`` (and optionally ``spec_img``), then fits
    and returns them together with ``spec_x`` for downstream use.  Returns a
    success flag; on failure the caller should construct a Gaussian fallback.

    Parameters
    ----------
    wave : :class:`numpy.ndarray`
        Extracted wavelength array, shape :math:`(N_{\rm spec},)`.
    flux : :class:`numpy.ndarray`
        Extracted flux array, shape :math:`(N_{\rm spec},)`.
    fluxivar : :class:`numpy.ndarray`
        Inverse variance of ``flux``, shape :math:`(N_{\rm spec},)`.
    waveimg : :class:`numpy.ndarray`
        Wavelength image, shape :math:`(N_{\rm spec}, N_{\rm spat})`.
    image : :class:`numpy.ndarray`
        Sky-subtracted science image, same shape as ``waveimg``.
    ivar : :class:`numpy.ndarray`
        Inverse variance of ``image``, same shape.
    totmask : :class:`numpy.ndarray`
        Combined boolean slit mask, same shape.  Valid slit pixels are those
        where ``totmask`` is ``True``.
    spec_img : :class:`numpy.ndarray` or None
        Integer image where ``spec_img[i, j] = i``, same shape.  If ``None``,
        spectral row indices are derived from ``totmask`` via
        ``np.where(totmask)[0]``.
    percentile_sn2 : :obj:`float`
        Upper percentile of :math:`(S/N)^2` used to estimate ``med_sn2``.
    fwhm : :obj:`float`
        Current FWHM estimate in pixels, used only in warning messages.
    sn_cap : :obj:`float`, optional
        Maximum S/N per pixel used when clipping ``fluxivar`` via
        :func:`~pypeit.utils.clip_ivar`.  Default is 100.0.
    good_pix_frac : :obj:`float`, optional
        Minimum fraction of in-range spectral pixels that must pass the
        ``indsp`` quality cuts for the fit to proceed.  Default is 0.05.

    Returns
    -------
    success : :obj:`bool`
        ``False`` if spectral fitting failed and a Gaussian should be
        returned; ``True`` otherwise.
    med_sn2 : :obj:`float`
        Median :math:`(S/N)^2` of the object; ``0.0`` on failure.
    norm_obj_x : :class:`numpy.ndarray` or None
        Normalised-object values for each valid slit pixel, shape
        :math:`(N_{\rm pix},)`.
    norm_ivar_x : :class:`numpy.ndarray` or None
        Inverse variance of ``norm_obj_x``, same shape.
    xtemp_x : :class:`numpy.ndarray` or None
        Cumulative S/N-weighted spectral coordinate for each valid slit pixel,
        shape :math:`(N_{\rm pix},)`.
    sn2_x : :class:`numpy.ndarray` or None
        :math:`(S/N)^2` interpolated onto each valid slit pixel, shape
        :math:`(N_{\rm pix},)`.
    spec_x : :class:`numpy.ndarray` or None
        Spectral row index of each valid slit pixel, shape :math:`(N_{\rm pix},)`
        int.  ``None`` on failure.
    """
    # Extract the wavelength image data into a 1D vector and get the valid
    # wavelength range of the 2D image
    wave_x = waveimg[totmask]
    wave_min = wave_x.min()
    wave_max = wave_x.max()
    good_wave = (wave >= wave_min) & (wave <= wave_max)

    # Smooth the 1D spectrum and its errors
    flux_sm = scipy.ndimage.median_filter(flux, size=5, mode='reflect')
    fluxivar_sm = scipy.ndimage.median_filter(fluxivar, size=5, mode='reflect')
    fluxivar_sm[fluxivar <= 0.0] = 0.0

    # Impose a maximum S/N
    fluxivar_sm = utils.clip_ivar(flux_sm, fluxivar_sm, sn_cap)

    # Determine if there are enough pixels to model the 1D spectrum
    indsp = good_wave & np.isfinite(flux_sm) & (flux_sm > -1000.0) & (fluxivar_sm > 0.0)
    eligible_pixels = np.sum(good_wave)
    if eligible_pixels == 0 or np.sum(indsp) < good_pix_frac * eligible_pixels:
        log.warning(
            'There are no pixels eligible to be fit for the object profile.  There is likely an '
            f'issue in local_skysub_extract. Returning a Gaussian with fwhm = {fwhm:.3f}'
        )
        return False, 0.0, None, None, None, None, None

    # Construct a high-order model of the spectrum.  The 2nd fit is used to
    # identify outlier pixels (bmask2) that are ignored through the rest of the
    # calculation.
    _, bmask, spline_flux, *_ = iterative_bspline_fit(
        wave, flux_sm, ivar=fluxivar_sm, gpm=indsp,
        kwargs_knots={'stride': 1.5}, kwargs_reject={'groupbadpix': True, 'maxrej': 1}
    )
    _, bmask2, spline_flux, *_ = iterative_bspline_fit(
        wave, flux_sm, ivar=fluxivar_sm, gpm=bmask,
        kwargs_knots={'stride': 1.5}, kwargs_reject={'groupbadpix': True, 'maxrej': 1}
    )
    try:
        # Construct a relatively low-order model to fall back on when the S/N is
        # low.
        _, cmask, cont_flux, *_ = iterative_bspline_fit(
            wave, flux_sm, ivar=fluxivar_sm, gpm=bmask2,
            kwargs_knots={'stride': 30}, kwargs_reject={'groupbadpix': True, 'maxrej': 1}
        )
    except Exception:
        log.warning(
            'Problem estimating S/N ratio of spectrum.  There is likely an issue in '
            f'local_skysub_extract. Returning a Gaussian with fwhm = {fwhm:.3f}.'
        )
        return False, 0.0, None, None, None, None, None

    # Compute the S/N over the set of good pixels
    fluxivar_sm = np.fmax(fluxivar_sm, 0)
    fluxivar_sm[np.logical_not(bmask2)] = 0.0
    sn2 = np.fmax(spline_flux**2 * fluxivar_sm, 0)
    sn2_indsp = sn2[indsp]
    if np.any(sn2_indsp > 0):
        sn2_percentile = np.percentile(sn2_indsp, percentile_sn2)
        _, med_sn2, _ = astropy.stats.sigma_clipped_stats(
            sn2_indsp[sn2_indsp > sn2_percentile], sigma_lower=3.0, sigma_upper=5.0
        )
    else:
        med_sn2 = 0.0
    log.info(f'sqrt(med(S/N)^2) = {np.sqrt(med_sn2):.2f}')

    # Ensure the wavelengths are sorted for interpolation
    isrt = np.argsort(wave[indsp], kind='stable')

    # Build the cumulative, S/N-weighted spectral coordinates
    nspec = wave.size
    sn2_1 = np.zeros(nspec)
    s2_interp = scipy.interpolate.interp1d(
        wave[indsp][isrt], sn2_indsp[isrt], assume_sorted=False, bounds_error=False, fill_value=0.0
    )
    sn2_1[good_wave] = s2_interp(wave[good_wave])
    cumweight = np.cumsum(4.0 + np.sqrt(np.fmax(sn2_1, 0.0)))

    # Sample the spectral coordinates onto the 2D subsample
    if spec_img is None:
        spec_x = np.where(totmask)[0]
    else:
        spec_x = spec_img[totmask]
    xtemp_x = cumweight[spec_x] / cumweight[-1]

    # Assign an S/N^2 value for each image pixel
    sn2_med_filt = scipy.ndimage.median_filter(sn2_indsp, size=9, mode='reflect')
    sn2_interp = scipy.interpolate.interp1d(
        wave[indsp][isrt], sn2_med_filt[isrt], assume_sorted=False, bounds_error=False,
        fill_value='extrapolate'
    )
    sn2_x = sn2_interp(wave_x)

    # Construct the model spectrum sampled over the 2D subsample
    npix = wave_x.size
    if med_sn2 <= 2.0:
        _, _, sigma1 = astropy.stats.sigma_clipped_stats(
            flux[indsp], sigma_lower=3.0, sigma_upper=5.0
        )
        spline_img_x = np.full(npix, np.fmax(sigma1, 0.0))
    else:
        # Fill masked regions
        spline_flux = np.where(good_wave, spline_flux, 0.0)
        spline_flux = pydl.djs_maskinterp(spline_flux, np.logical_not(bmask2))
        cont_flux = np.where(good_wave, cont_flux, 0.0)
        cont_flux = pydl.djs_maskinterp(cont_flux, np.logical_not(cmask))
        if med_sn2 <= 5.0:
            spline_flux = cont_flux

        # Identify more bad pixels and replace them
        badpix = (spline_flux <= 0.5) | np.logical_not(bmask2)
        goodval = (cont_flux > 0.0) & (cont_flux < 5e5)
        indbad1 = badpix & goodval
        if np.any(indbad1):
            spline_flux[indbad1] = cont_flux[indbad1]
        indbad2 = badpix & np.logical_not(goodval)
        if np.any(indbad2) or np.any(np.logical_not(badpix)):
            spline_flux[indbad2] = np.median(spline_flux[np.logical_not(badpix)])

        # Median filter
        spline_flux = scipy.ndimage.median_filter(spline_flux, size=5, mode='reflect')

        # Interpolate onto the flattened 2D grid
        isrt1 = np.argsort(wave[good_wave], kind='stable')
        spline_img_interp = scipy.interpolate.interp1d(
            wave[good_wave][isrt1], spline_flux[good_wave][isrt1], assume_sorted=False,
            bounds_error=False, fill_value='extrapolate'
        )
        spline_img_x = spline_img_interp(wave_x)

    # Normalize the 2D data to remove the spectral variations
    norm_obj_x = np.where(spline_img_x != 0.0, image[totmask] / spline_img_x, 0.0)
    norm_ivar_x = ivar[totmask] * spline_img_x**2
    ivar_mask_x = (
        (norm_obj_x > -0.2) & (norm_obj_x < 0.7) & np.isfinite(norm_obj_x)
        & np.isfinite(norm_ivar_x)
    )
    norm_ivar_x[np.logical_not(ivar_mask_x)] = 0.0

    return True, med_sn2, norm_obj_x, norm_ivar_x, xtemp_x, sn2_x, spec_x


def _compute_bspline_knots(dspat_x, sigma, spec_x, med_sn2, prof_nsigma, good_x):
    r"""
    Initialise spatial-coordinate arrays and B-spline knot grid.

    Encapsulates Block 6 of :func:`fit_profile_refactor`: computes the
    normalised coordinate ``sigma_x`` for each valid slit pixel, the profile
    half-extent ``limit``, the profile sigma bounds, and the sinh-spaced
    interior B-spline knot positions.

    Parameters
    ----------
    dspat_x : :class:`numpy.ndarray`
        Spatial separation from the trace for each valid slit pixel,
        shape :math:`(N_{\rm pix},)`.
    sigma : :class:`numpy.ndarray`
        Current Gaussian sigma estimate per spectral pixel,
        shape :math:`(N_{\rm spec},)`.
    spec_x : :class:`numpy.ndarray`
        Spectral row index of each valid slit pixel, shape
        :math:`(N_{\rm pix},)` int.
    med_sn2 : :obj:`float`
        Median :math:`(S/N)^2` of the object.
    prof_nsigma : :obj:`float` or None
        If set, fit the profile out to this many sigma (extended objects);
        overrides the automatic sigma bounds computed from ``med_sn2``.
    good_x : :class:`numpy.ndarray`
        Boolean array marking pixels with positive ``norm_ivar_x``,
        shape :math:`(N_{\rm pix},)`.

    Returns
    -------
    sigma_x : :class:`numpy.ndarray`
        Normalised spatial coordinate (units of ``sigma``), shape
        :math:`(N_{\rm pix},)`.
    limit : :obj:`float`
        Profile half-extent in units of ``sigma``.
    min_sigma : :obj:`float`
        Lower bound on ``sigma_x`` for profile fitting.
    max_sigma : :obj:`float`
        Upper bound on ``sigma_x`` for profile fitting.
    bkpt : :class:`numpy.ndarray`
        Interior B-spline knot positions in ``sigma_x`` units.
    """
    sigma_x = dspat_x / sigma[spec_x]
    limit = scipy.special.erfcinv(0.1 / np.sqrt(med_sn2)) * np.sqrt(2.0)
    if prof_nsigma is None:
        sinh_space = 0.25 * np.log10(np.fmax(1000.0 / np.sqrt(med_sn2), 10.0))
        abs_sigma = np.fmin((np.abs(sigma_x[good_x])).max(), 2.0 * limit)
        min_sigma = np.fmax(sigma_x[good_x].min(), -abs_sigma)
        max_sigma = np.fmin(sigma_x[good_x].max(), abs_sigma)
        nb = (np.arcsinh(abs_sigma) / sinh_space).astype(int) + 1
    else:
        log.info(f"Using prof_nsigma={prof_nsigma:.2f} for extended/bright objects")
        # Bug 1 fix: was np.round(prof_nsigma > 10) which evaluates a boolean.
        nb = max(1, round(prof_nsigma / 10))
        max_sigma = prof_nsigma
        min_sigma = -prof_nsigma
        sinh_space = np.arcsinh(prof_nsigma) / nb
    rb = np.sinh((np.arange(nb) + 0.5) * sinh_space)
    bkpt = np.concatenate([(-rb)[::-1], rb])
    keep = (bkpt >= min_sigma) & (bkpt <= max_sigma)
    return sigma_x, limit, min_sigma, max_sigma, bkpt[keep] 


def _deriv_walk(bset, init_limit, sign):
    r"""
    Evaluate the log-derivative B-spline walk for one apodization side.

    Builds a vector of candidate limit positions stepping inward from
    ``init_limit``, evaluates the logarithmic derivative
    :math:`d(\ln P)/d(\ln |\sigma_x|)` of the profile at each position via
    finite differences, then returns the position where that derivative first
    reaches its extreme value (capped at ±1).

    Parameters
    ----------
    bset : :class:`~pypeit.core.bspline.BSpline`
        Final 1-D B-spline profile.
    init_limit : :obj:`float`
        Initial apodization limit (left: negative, right: positive) in
        ``sigma_x`` units, as returned by the threshold walk in
        :func:`_build_profile`.
    sign : :obj:`int`
        ``-1`` for the left side (negative ``sigma_x``); ``+1`` for the
        right side.

    Returns
    -------
    limit : :obj:`float`
        Refined apodization limit in ``sigma_x`` units.
    deriv : :obj:`float`
        Log-derivative at ``limit``.
    fit_val : :obj:`float`
        Profile value at ``limit``, used to match the exponential tail.
    """
    step = -sign * 0.1
    lim_vec = np.arange(init_limit + step, sign * 1.0, step)
    if lim_vec.size == 0:
        lim_vec = np.asarray([sign * 1.1])
    fit1, _ = bset.value(lim_vec)
    fit2, _ = bset.value(lim_vec * 0.9)
    deriv_vec = ((np.ma.log(fit2) - np.ma.log(fit1)) / (0.1 * lim_vec)).filled(0.0)
    if sign < 0:
        extreme = np.fmax(deriv_vec.min(), -1.0)
        cross = np.where(deriv_vec <= extreme)[0]
    else:
        extreme = np.fmin(deriv_vec.max(), 1.0)
        cross = np.where(deriv_vec >= extreme)[0]
    if cross.size > 0:
        idx = cross[0]
        return lim_vec[idx], deriv_vec[idx], fit1[idx]
    return sign * 1.0, deriv_vec[-1], fit1[-1]


def _build_profile(
    bset, sigma_x, min_sigma, max_sigma, apodize, ss, median_fit, min_level, limit
):
    r"""
    Locate apodization limits and apply exponential tails to the profile.

    Encapsulates Blocks 11-12 of :func:`fit_profile_refactor`.  Evaluates the
    final B-spline profile over all good pixels, re-locates the profile peak
    and apodization limits from the logarithmic derivative of the profile,
    then replaces the profile outside those limits with matched exponential
    tails.

    Parameters
    ----------
    bset : BSpline
        Final 1-D B-spline profile.
    sigma_x : :class:`numpy.ndarray`
        Normalised spatial coordinate, shape :math:`(N_{\rm pix},)`.
    min_sigma : :obj:`float`
        Lower bound on ``sigma_x`` for the good-pixel mask.
    max_sigma : :obj:`float`
        Upper bound on ``sigma_x`` for the good-pixel mask.
    apodize : :obj:`bool`
        Apply an exponential apodization function to the edges of the profile
        model.
    ss : :class:`numpy.ndarray`
        Indices that sort ``sigma_x`` in ascending order.  Ignored if
        ``apodize`` is False.
    median_fit : :obj:`float`
        Baseline level of the profile (subtracted before peak-finding).  Ignored
        if ``apodize`` is False.
    min_level : :obj:`float`
        Minimum profile level used to walk the initial limit positions.  Ignored
        if ``apodize`` is False.
    limit : :obj:`float`
        Profile half-extent in units of ``sigma_x``.  Ignored if ``apodize`` is
        False.

    Returns
    -------
    full_bsp : :class:`numpy.ndarray`
        Profile array of length :math:`N_{\rm pix}` with exponential tails
        applied outside ``[l_limit, r_limit]``.
    l_limit : :obj:`float`
        Final left apodization limit in ``sigma_x`` units (0.0 when
        ``prof_nsigma`` is set).
    r_limit : :obj:`float`
        Final right apodization limit in ``sigma_x`` units (0.0 when
        ``prof_nsigma`` is set).
    """
    # Sample the Bspline profile over the valid coordinates
    igood = (sigma_x > min_sigma) & (sigma_x < max_sigma)
    sigma_x_igood = sigma_x[igood]
    full_bsp = np.zeros(sigma_x.size)
    full_bsp_igood, _ = bset.value(sigma_x_igood)
    full_bsp[igood] = full_bsp_igood

    if not apodize:
        return full_bsp, 0.0, 0.0

    # Find the profile peak
    isrt2 = sigma_x_igood.argsort(kind='stable')
    peak, peak_x, lwhm, rwhm = _findfwhm(full_bsp_igood[isrt2] - median_fit, sigma_x_igood[isrt2])

    # Determine the left and right edge of the profile
    sorted_sigma = sigma_x[ss]
    sorted_bsp = full_bsp[ss]
    threshold = min_level + median_fit

    left_cond = (
        (sorted_bsp < threshold) & (sorted_sigma < peak_x)
    ) | (sorted_sigma < peak_x - limit)
    edges_l = np.ma.flatnotmasked_edges(
        np.ma.array(sorted_sigma, mask=np.logical_not(left_cond))
    )
    l_limit = (sorted_sigma[edges_l[1]] if edges_l is not None else sorted_sigma[0]) - 0.1

    right_cond = (
        (sorted_bsp < threshold) & (sorted_sigma > peak_x)
    ) | (sorted_sigma > peak_x + limit)
    edges_r = np.ma.flatnotmasked_edges(
        np.ma.array(sorted_sigma, mask=np.logical_not(right_cond))
    )
    r_limit = (sorted_sigma[edges_r[0]] if edges_r is not None else sorted_sigma[-1]) + 0.1

    l_limit, l_deriv, l_fit_val = _deriv_walk(bset, l_limit, sign=-1)
    r_limit, r_deriv, r_fit_val = _deriv_walk(bset, r_limit, sign=1)

    if l_deriv < 0 and r_deriv > 0:
        # Apply the exponential apodization
        left = sigma_x < l_limit
        full_bsp[left] = l_fit_val * np.exp(-(sigma_x[left] - l_limit) * l_deriv)
        right = sigma_x > r_limit
        full_bsp[right] = r_fit_val * np.exp(-(sigma_x[right] - r_limit) * r_deriv)
    else:
        log.warning('Derivative computation failed.  Object profile not apodized.')

    return full_bsp, l_limit, r_limit


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def fit_profile_refactor(
    image, ivar, waveimg, thismask, spat_img, trace_in, wave, flux, fluxivar, inmask=None,
    spec_img=None, thisfwhm=4.0, max_trace_corr=2.0, sn_gauss=4.0, percentile_sn2=70.0,
    prof_nsigma=None, no_deriv=False, gauss=False, obj_string='', show_profile=False
):
    r"""
    Fit a non-parametric object profile to an object spectrum.

    Refactored drop-in replacement for
    :func:`~pypeit.core.spatialprofile.fit_profile`.  The signature, return
    values, and overall algorithm are identical; see that function's docstring
    for full parameter and return-value documentation.  Improvements relative to
    the legacy implementation are described in the module docstring above.

    Parameters
    ----------
    image : :class:`numpy.ndarray`
        Sky-subtracted science image, shape :math:`(N_{\rm spec}, N_{\rm spat})`.
    ivar : :class:`numpy.ndarray`
        Inverse variance of ``image``, same shape.
    waveimg : :class:`numpy.ndarray`
        Wavelength image, same shape.
    thismask : :class:`numpy.ndarray`
        Boolean slit mask, same shape.
    spat_img : :class:`numpy.ndarray`
        Spatial-coordinate image, same shape.
    trace_in : :class:`numpy.ndarray`
        Object trace, shape :math:`(N_{\rm spec},)`.
    wave : :class:`numpy.ndarray`
        Extracted wavelength array, shape :math:`(N_{\rm spec},)`.
    flux : :class:`numpy.ndarray`
        Extracted flux array, shape :math:`(N_{\rm spec},)`.
    fluxivar : :class:`numpy.ndarray`
        Inverse variance of ``flux``, shape :math:`(N_{\rm spec},)`.
    inmask : :class:`numpy.ndarray`, optional
        Additional boolean mask.
    spec_img : :class:`numpy.ndarray`, optional
        Integer image where ``spec_img[i, j] = i``; same shape as ``image``.
        If ``None`` (default), row indices are derived from ``totmask`` via
        ``np.where(totmask)[0]``.  Construct as
        ``np.broadcast_to(np.arange(nspec, dtype=int)[:, None], image.shape)``.
    thisfwhm : :obj:`float`, optional
        Initial FWHM estimate in pixels.
    max_trace_corr : :obj:`float`, optional
        Maximum trace correction in pixels.
    sn_gauss : :obj:`float`, optional
        S/N threshold below which a Gaussian profile is returned.
    percentile_sn2 : :obj:`float`, optional
        Upper percentile of S/N^2 used to estimate ``med_sn2``.
    prof_nsigma : :obj:`float`, optional
        If set, fit the profile out to this many sigma (extended objects).
    no_deriv : :obj:`bool`, optional
        Disable exponential apodization.  If ``prof_nsigma`` is provided, this
        is always set to True regardless of the input value.
    gauss : :obj:`bool`, optional
        Force the Gaussian fallback regardless of S/N.
    obj_string : :obj:`str`, optional
        Object label for log messages and QA plots.
    show_profile : :obj:`bool`, optional
        Display the QA plot.

    Returns
    -------
    profile_model : :class:`numpy.ndarray`
        2-D normalised spatial profile, shape :math:`(N_{\rm spec}, N_{\rm spat})`.
    xnew : :class:`numpy.ndarray`
        Corrected trace, shape :math:`(N_{\rm spec},)`.
    fwhmfit : :class:`numpy.ndarray`
        FWHM estimate per spectral pixel, shape :math:`(N_{\rm spec},)`.
    med_sn2 : :obj:`float`
        Median S/N^2 of the object.
    """

    # ------------------------------------------------------------------
    # Block 1 — Initialisation
    # ------------------------------------------------------------------
    totmask = (ivar > 0.) & thismask
    if inmask is not None:
        totmask &= inmask

    nspec, nspat = image.shape
    _thisfwhm = np.fmax(thisfwhm, 1.0)
    fwhmfit = np.full(nspec, _thisfwhm)
    sig2fwhm = np.sqrt(8*np.log(2))
    sigma = fwhmfit / sig2fwhm

    # ------------------------------------------------------------------
    # Blocks 0 and 2–5 — Slit-pixel extraction, spectral fitting, normalisation
    # ------------------------------------------------------------------
    (
        success, med_sn2, norm_obj_x, norm_ivar_x, xtemp_x, sn2_x, spec_x
    ) = _fit_spectrum_and_normalize(
        wave=wave, flux=flux, fluxivar=fluxivar,
        waveimg=waveimg, image=image, ivar=ivar, totmask=totmask, spec_img=spec_img,
        percentile_sn2=percentile_sn2, fwhm=_thisfwhm
    )
    if success:
        good_x = norm_ivar_x > 0
        gauss_kwargs = {'ind': good_x, 'xtrunc': 7.0}
        ngood = good_x.sum()
        log.info(f"Gaussian vs b-spline of width {_thisfwhm:.2f} pixels")
    else:
        med_sn2 = 0.0
        show_profile = False
        gauss_kwargs = {}
    if not success or gauss or ngood < 10 or med_sn2 < sn_gauss**2:
        log.info(
            'Returning Gaussian profile for one of the following reasons: Determination of the '
            'S/N failed, a Gaussian profile was specifically requested, there are too few good '
            f'pixels, or the measured S/N is below the provided limit ({sn_gauss:.1f}).'
        )
        profile_model = _return_gaussian(
            spat_img, norm_obj_x, trace_in, sigma, _thisfwhm, med_sn2, obj_string, show_profile,
            **gauss_kwargs
        )
        return profile_model, trace_in, fwhmfit, med_sn2

    # ------------------------------------------------------------------
    # Block 6 — Knot setup
    # ------------------------------------------------------------------
    npix = spec_x.size
    dspat_x = spat_img[totmask] - trace_in[spec_x]
    sigma_x, limit, min_sigma, max_sigma, bkpt = _compute_bspline_knots(
        dspat_x=dspat_x, sigma=sigma, spec_x=spec_x, med_sn2=med_sn2,
        prof_nsigma=prof_nsigma, good_x=good_x
    )

    # ------------------------------------------------------------------
    # Block 7 — Initial profile fit
    # ------------------------------------------------------------------
    valid_err = norm_ivar_x > 0
    good_snr = valid_err & (sn2_x > sn_gauss**2)
    good_pos = valid_err & (sigma_x >= min_sigma) & (sigma_x <= max_sigma)
    if good_snr.sum() >= 0.2 * good_pos.sum():
        inside = np.where(good_snr & good_pos)[0]
    else:
        inside = np.where(good_pos)[0]

    si = inside[np.argsort(sigma_x[inside], kind='stable')]

    bset, *_ = iterative_bspline_fit(
        sigma_x[si], norm_obj_x[si], ivar=norm_ivar_x[si],
        nord=4, kwargs_knots={'interior': bkpt}, maxiter=15, upper=1, lower=1
    )
    mode_fit, _ = bset.value(sigma_x[si])
    median_fit = np.median(norm_obj_x[valid_err])
    if np.abs(median_fit) > 0.01:
        log.info(f"Median flux level in profile is not zero: median = {median_fit:.4f}")
    else:
        median_fit = 0.0

    sorted_sigma = sigma_x[si]
    peak, peak_x, lwhm, rwhm = _findfwhm(mode_fit - median_fit, sorted_sigma)
    trace_corr = np.full(nspec, peak_x)
    min_level = peak * np.exp(-0.5 * limit**2)
    threshold = min_level + median_fit

    bspline_fwhm = (rwhm - lwhm) * _thisfwhm / sig2fwhm
    log.info(
        f'Bspline FWHM: {bspline_fwhm:.4f}, compared to initial object finding FWHM: '
        f'{_thisfwhm:.4f}'
    )
    sigma *= (rwhm - lwhm) / sig2fwhm
    limit *= (rwhm - lwhm) / sig2fwhm

    left_cond = (
        (mode_fit < threshold) & (sorted_sigma < peak_x)
    ) | (sorted_sigma < peak_x - limit)
    m_left = np.ma.array(sorted_sigma, mask=~left_cond)
    edges_l = np.ma.flatnotmasked_edges(m_left)
    l_limit = sorted_sigma[edges_l[1]] if edges_l is not None else min_sigma

    right_cond = (
        (mode_fit < threshold) & (sorted_sigma > peak_x)
    ) | (sorted_sigma > peak_x + limit)
    m_right = np.ma.array(sorted_sigma, mask=~right_cond)
    edges_r = np.ma.flatnotmasked_edges(m_right)
    r_limit = sorted_sigma[edges_r[0]] if edges_r is not None else max_sigma

    log.info(
        f"Trace limits: limit={limit:.4f}, min_level={min_level:.4f}, l_limit={l_limit:.4f}, "
        f"r_limit={r_limit:.4f}"
    )

    mask_x = np.zeros(npix, dtype=bool)
    mask_x[si] = (norm_ivar_x[si] > 0) & (np.abs(norm_obj_x[si] - mode_fit) < 0.1)
    inside, = np.where((sorted_sigma > l_limit) & (sorted_sigma < r_limit) & mask_x[si])
    ninside = inside.size

    if ninside < 10:
        log.info(
            'Returning Gaussian profile because there are too few pixels inside l_limit and '
            'r_limit.'
        )
        profile_model = _return_gaussian(
            spat_img, norm_obj_x, trace_in + trace_corr * sigma, sigma, bspline_fwhm, med_sn2,
            obj_string, show_profile, ind=good_x, l_limit=l_limit, r_limit=r_limit, xlim=7.0,
        )
        return profile_model, trace_in, fwhmfit, med_sn2

    # ------------------------------------------------------------------
    # Block 8 — Iterative trace and width correction (sigma_iter=3)
    # ------------------------------------------------------------------
    sigma_iter = 3
    xx = np.bincount(spec_x, weights=xtemp_x, minlength=nspec) / nspat
    isort = xtemp_x[si[inside]].argsort(kind='stable')
    inside = si[inside[isort]]
    pb = np.ones(inside.size)
    area = np.ones(nspec)

    for iiter in range(sigma_iter):
        mode_zero = pb * bset.value(sigma_x[inside])[0]

        mode_shift = pb * (
            bset.value(sigma_x[inside] - 0.5)[0] - bset.value(sigma_x[inside] + 0.5)[0]
        )
        mode_shift[(sigma_x[inside] <= l_limit + 0.5) | (sigma_x[inside] >= r_limit - 0.5)] = 0.0

        mode_stretch = pb * bset.value(sigma_x[inside] / 1.3)[0] / 1.3 - mode_zero

        nbkpts = int(np.log10(np.fmax(med_sn2, 11.0)))

        profile_basis = np.column_stack((mode_zero, mode_shift))
        mode_shift_bspl, mode_shift_gpm, *_ = iterative_bspline_fit(
            xtemp_x[inside], norm_obj_x[inside], ivar=norm_ivar_x[inside],
            basis=profile_basis, maxiter=1, kwargs_knots={'count': nbkpts}
        )
        if not np.any(mode_shift_gpm):
            log.info(
                'Returning a Gaussian profile because the B-spline fit to trace correction failed '
                f'for ninside = {ninside}.'
            )
            profile_model = _return_gaussian(
                spat_img, norm_obj_x, trace_in + trace_corr * sigma, sigma, bspline_fwhm,
                med_sn2, obj_string, show_profile, ind=good_x, l_limit=l_limit,
                r_limit=r_limit, xlim=7.0,
            )
            return profile_model, trace_in, fwhmfit, med_sn2

        temp_set = mode_shift_bspl.to_1d()
        h0, _ = temp_set.value(xx)
        h1, _ = temp_set.value(xx, coeff=mode_shift_bspl.coeff[:, 1])
        ratio_10 = h1 * _inverse(h0)
        delta_trace_corr = ratio_10 / (1.0 + np.abs(ratio_10) / 0.1)
        trace_corr += delta_trace_corr

        profile_basis = np.column_stack((mode_zero, mode_stretch))
        mode_stretch_bspl, mode_stretch_gpm, *_ = iterative_bspline_fit(
            xtemp_x[inside], norm_obj_x[inside], ivar=norm_ivar_x[inside],
            basis=profile_basis, maxiter=1, kwargs_knots={'full': mode_shift_bspl.breakpoints}
        )
        if not np.any(mode_stretch_gpm):
            log.info(
                'Returning a Gaussian profile because the B-spline fit to the width correction '
                f'failed for ninside = {ninside}.'
            )
            profile_model = _return_gaussian(
                spat_img, norm_obj_x, trace_in + trace_corr * sigma, sigma, bspline_fwhm,
                med_sn2, obj_string, show_profile, ind=good_x, l_limit=l_limit,
                r_limit=r_limit, xlim=7.0,
            )
            return profile_model, trace_in, fwhmfit, med_sn2

        temp_set = mode_stretch_bspl.to_1d()
        h0, _ = temp_set.value(xx)
        h2, _ = temp_set.value(xx, coeff=mode_stretch_bspl.coeff[:, 1])
        h0 = np.fmax(h0 + h2 * mode_stretch.sum() / mode_zero.sum(), 0.1)
        ratio_20 = h2 * _inverse(h0)
        sigma_factor = 0.3 * ratio_20 / (1.0 + np.abs(ratio_20))

        log.info(f"Iteration # {iiter+1}")
        log.info(
            f"Median abs value of trace correction = {np.median(np.abs(delta_trace_corr)):.3f}"
        )
        log.info(f"Median abs value of width correction = {np.median(np.abs(sigma_factor)):.3f}")

        sigma *= (1.0 + sigma_factor)
        area *= h0 / (1.0 + sigma_factor)

        sigma_x = dspat_x / sigma[spec_x] - trace_corr[spec_x]

        # TODO: In the previous version, iiter was iterated from 1 to sigma_iter
        # and this was check was for if iiter < sigma_iter - 1.  When moving to
        # a 0-indexed iterator, this makes this if statement sigma_iter - 2.  Is
        # that the original intention, or was this an index mistake?
        if iiter < sigma_iter - 2:
            ss = sigma_x[inside].argsort(kind='stable')
            pb = area[spec_x[inside]]
            keep = (bkpt >= sigma_x[inside].min()) & (bkpt <= sigma_x[inside].max())
            if keep.sum() == 0:
                keep = np.ones(bkpt.size, dtype=bool)
            bset, bset_gpm, *_ = iterative_bspline_fit(
                sigma_x[inside[ss]], norm_obj_x[inside[ss]],
                ivar=norm_ivar_x[inside[ss]], basis=pb[ss], nord=4,
                kwargs_knots={'interior': bkpt[keep]}, maxiter=2
            )
            if not np.any(bset_gpm):
                log.info(
                    'Returning a Gaussian profile because the B-spline profile fit in '
                    f'trace/width loop failed for ninside = {ninside}'
                )
                profile_model = _return_gaussian(
                    spat_img, norm_obj_x, trace_in + trace_corr * sigma, sigma, bspline_fwhm,
                    med_sn2, obj_string, show_profile, ind=good_x, l_limit=l_limit,
                    r_limit=r_limit, xlim=7.0,
                )
                return profile_model, trace_in, fwhmfit, med_sn2

            bset = bset.to_1d()

    # ------------------------------------------------------------------
    # Block 9 — Final trace
    # ------------------------------------------------------------------
    if np.median(np.abs(trace_corr * sigma)) < max_trace_corr:
        xnew = trace_corr * sigma + trace_in
    else:
        xnew = trace_in
    fwhmfit = sigma * sig2fwhm

    # ------------------------------------------------------------------
    # Block 10 — Final profile fit
    # ------------------------------------------------------------------
    ss = sigma_x.argsort(kind='stable')
    inside, = np.where(
        (sigma_x[ss] >= min_sigma) & (sigma_x[ss] <= max_sigma) & mask_x[ss]
        & np.isfinite(norm_obj_x[ss]) & np.isfinite(norm_ivar_x[ss])
    )
    pb_x = area[spec_x]
    bset, outmask, *_ = iterative_bspline_fit(
        sigma_x[ss[inside]], norm_obj_x[ss[inside]],
        ivar=norm_ivar_x[ss[inside]], basis=pb_x[ss[inside]], nord=4,
        kwargs_knots={'interior': bkpt}, upper=10, lower=10
    )
    bset = bset.to_1d()

    # ------------------------------------------------------------------
    # Blocks 11–12 — Apodization limit search and exponential tails
    # ------------------------------------------------------------------
    full_bsp, l_limit, r_limit = _build_profile(
        bset, sigma_x, min_sigma, max_sigma, (prof_nsigma is None and not no_deriv),
        ss, median_fit, min_level, limit
    )

    # ------------------------------------------------------------------
    # Block 13 — Normalisation, logging, QA, return
    # ------------------------------------------------------------------
    profile_x = full_bsp * pb_x
    profile_model = np.zeros((nspec, nspat))
    profile_model[totmask] = profile_x

    res_mode = (
        (norm_obj_x[ss[inside]] - profile_x[ss[inside]])
        * np.sqrt(norm_ivar_x[ss[inside]])
    )
    chi_med = np.median(res_mode[outmask & (norm_ivar_x[ss[inside]] > 0)] ** 2)

    log.info("--------------------  Results of Profile Fit --------------------")
    log.info(
        f" min(fwhmfit)={fwhmfit.min():.2f}; max(fwhmfit)={fwhmfit.max():.2f};"
        f" median(chi^2)={chi_med:.2f}; nbkpts={bkpt.size}"
    )
    log.info("-----------------------------------------------------------------")

    if not np.all(np.isfinite(xnew)):
        log.warning(
            'NaN/Inf pixel values exist in the trace correction.  Returning the original trace.'
        )
        xnew = trace_in

    invalid = np.logical_not(np.isfinite(profile_model))
    if np.any(invalid):
        log.warning('Setting NaN/Inf pixel values in object profile model to zero.')
        profile_model[invalid] = 0.0

    # Normalise each spectral row to unit sum
    row_sums = profile_model.sum(axis=1)
    indx = row_sums > 0.0
    if np.any(indx):
        profile_model[indx,:] /= row_sums[:,None]
        profile_model[np.logical_not(indx),:] = 0.

    info_string = (
        f"FWHM range:{fwhmfit.min():.2f} - {fwhmfit.max():.2f}, S/N={np.sqrt(med_sn2):.3f}"
        f", median(chi^2)={chi_med:.3f}"
    )
    if show_profile:
        qa_fit_profile(
            sigma_x, norm_obj_x * _inverse(pb_x), full_bsp,
            l_limit=l_limit, r_limit=r_limit,
            ind=ss[inside], xlim=prof_nsigma, title=f'{obj_string} {info_string}'
        )

    return profile_model, xnew, fwhmfit, med_sn2
