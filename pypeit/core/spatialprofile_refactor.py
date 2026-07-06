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
import numpy as np
import scipy.interpolate
import scipy.ndimage
import scipy.special

from pypeit import log
from pypeit import utils
from pypeit.core.bspline import BSpline, Knots
from pypeit.core.fitting import iterative_bspline_fit
from pypeit.core import pydl
from pypeit.core.spatialprofile import qa_fit_profile


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _bspline2d_to_1d(bset2d):
    """
    Convert the first component of a ``BSpline2D`` to a 1-D ``BSpline``.

    :func:`~pypeit.core.fitting.iterative_bspline_fit` returns a ``BSpline2D``
    (``coeff`` shape ``(nknots, npoly)``) whenever a ``basis`` array is
    supplied.  This helper extracts the zeroth coefficient column so that the
    result can be evaluated at arbitrary sigma_x values with ``.value()``.

    Parameters
    ----------
    bset2d : BSpline2D
        Two-dimensional B-spline returned by ``iterative_bspline_fit``.

    Returns
    -------
    BSpline
        One-dimensional B-spline with ``coeff = bset2d.coeff[:, 0]``.
    """
    bset = BSpline(knots=Knots(full=bset2d.breakpoints), nord=bset2d.nord)
    bset.bkpt_gpm = bset2d.bkpt_gpm.copy()
    bset.coeff = bset2d.coeff[:, 0]
    return bset


def _findfwhm(model, sig_x):
    """
    Calculate the FWHM of an profile with sub-pixel precision.

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
    sigma_x = (_x - _center) / _sigma / np.sqrt(2.0)
    delta = 0.5 / _sigma / np.sqrt(2.0)
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
    info_string = f"FWHM={fwhm:6.2f}, S/N={np.sqrt(med_sn2):8.3f}"
    title_string = obj_string + ', ' + info_string
    log.info(title_string)
    if show_profile:
        _center = np.atleast_1d(center)
        _sigma = np.atleast_1d(sigma)
        if _center.size > 1:
            _center = _center[:, None]
        if _sigma.size > 1:
            _sigma = _sigma[:, None]
        sigma_x = (spat_img - _center) / _sigma
        qa_fit_profile(
            sigma_x, norm_obj, profile_model,
            title=title_string, l_limit=l_limit, r_limit=r_limit,
            ind=ind, xlim=xlim, xtrunc=xtrunc,
        )
    return profile_model


def _fit_spectrum_and_normalize(
    wave, flux, fluxivar, waveimg, image, ivar, totmask, percentile_sn2, fwhm
):
    """
    Fit flux B-splines, estimate S/N, and build the normalised-object image.

    Encapsulates the spectral-fitting and normalisation steps of
    :func:`fit_profile_refactor` (blocks 2–4).  Returns a success flag; on
    failure the caller should construct a Gaussian fallback.

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
        Combined boolean mask, same shape as ``image``.
    percentile_sn2 : :obj:`float`
        Upper percentile of :math:`(S/N)^2` used to estimate ``med_sn2``.
    fwhm : :obj:`float`
        Current FWHM estimate in pixels, used only in warning messages.

    Returns
    -------
    success : :obj:`bool`
        ``False`` if spectral fitting failed and a Gaussian should be
        returned; ``True`` otherwise.
    med_sn2 : :obj:`float`
        Median :math:`(S/N)^2` of the object; ``0.0`` on failure.
    norm_obj : :class:`numpy.ndarray` or None
        Normalised-object image, shape :math:`(N_{\rm spec}, N_{\rm spat})`.
    norm_ivar : :class:`numpy.ndarray` or None
        Inverse variance of ``norm_obj``, same shape.
    good : :class:`numpy.ndarray` or None
        Boolean flat array of length :math:`N_{\rm spec} N_{\rm spat}`
        marking pixels with positive ``norm_ivar``.
    xtemp : :class:`numpy.ndarray` or None
        Cumulative S/N-weighted spectral coordinate,
        shape :math:`(N_{\rm spec}, N_{\rm spat})`.
    sn2_img : :class:`numpy.ndarray` or None
        :math:`(S/N)^2` image interpolated onto the 2D grid, same shape as
        ``image``.
    """
    nspec, nspat = image.shape
    sn2_img = np.zeros((nspec, nspat))
    spline_img = np.zeros((nspec, nspat))

    flux_sm = scipy.ndimage.median_filter(flux, size=5, mode='reflect')
    fluxivar_sm0 = scipy.ndimage.median_filter(fluxivar, size=5, mode='reflect')
    fluxivar_sm0[fluxivar <= 0.0] = 0.0
    wave_min = waveimg[totmask].min()
    wave_max = waveimg[totmask].max()

    # ------------------------------------------------------------------
    # Block 2 — Flux B-spline fitting
    # ------------------------------------------------------------------
    sn_cap = 100.0
    fluxivar_sm = utils.clip_ivar(flux_sm, fluxivar_sm0, sn_cap)
    indsp = (
        (wave >= wave_min) & (wave <= wave_max) & np.isfinite(flux_sm) & (flux_sm > -1000.0)
        & (fluxivar_sm > 0.0)
    )
    eligible_pixels = np.sum((wave >= wave_min) & (wave <= wave_max))
    good_pix_frac = 0.05
    if np.sum(indsp) < good_pix_frac * eligible_pixels or eligible_pixels == 0:
        log.warning(
            'There are no pixels eligible to be fit for the object profile.  There is likely an '
            f'issue in local_skysub_extract. Returning a Gaussian with fwhm = {fwhm:.3f}'
        )
        return False, 0.0, None, None, None, None, None

    b_answer, bmask, *_ = iterative_bspline_fit(
        wave[indsp], flux_sm[indsp], ivar=fluxivar_sm[indsp],
        kwargs_knots={'stride': 1.5}, kwargs_reject={'groupbadpix': True, 'maxrej': 1}
    )
    b_answer, bmask2, *_ = iterative_bspline_fit(
        wave[indsp], flux_sm[indsp], ivar=np.where(bmask, fluxivar_sm[indsp], 0.0),
        kwargs_knots={'stride': 1.5}, kwargs_reject={'groupbadpix': True, 'maxrej': 1}
    )
    c_answer, cmask, *_ = iterative_bspline_fit(
        wave[indsp], flux_sm[indsp], ivar=np.where(bmask2, fluxivar_sm[indsp], 0.0),
        kwargs_knots={'stride': 30}, kwargs_reject={'groupbadpix': True, 'maxrej': 1}
    )
    spline_flux, _ = b_answer.value(wave[indsp])
    try:
        cont_flux, _ = c_answer.value(wave[indsp])
    except Exception:
        log.warning(
            'Problem estimating S/N ratio of spectrum\nThere is likely an issue in '
            f'local_skysub_extract. Returning a Gaussian with fwhm = {fwhm:.3f}'
        )
        return False, 0.0, None, None, None, None, None

    # ------------------------------------------------------------------
    # Block 3 — S/N estimation
    # ------------------------------------------------------------------
    sqrt_ivar = np.sqrt(np.fmax(fluxivar_sm[indsp], 0))
    sqrt_ivar[np.logical_not(bmask2)] = 0.0
    sn2 = np.fmax(spline_flux * sqrt_ivar, 0)**2
    ind_nonzero = sn2 > 0
    if ind_nonzero.sum() > 0:
        sn2_percentile = np.percentile(sn2, percentile_sn2)
        _, med_sn2, _ = astropy.stats.sigma_clipped_stats(
            sn2[sn2 > sn2_percentile], sigma_lower=3.0, sigma_upper=5.0
        )
    else:
        med_sn2 = 0.0

    spline_flux1 = np.zeros(nspec)
    cont_flux1 = np.zeros(nspec)
    sn2_1 = np.zeros(nspec)
    ispline = (wave >= wave_min) & (wave <= wave_max)
    spline_tmp, _ = b_answer.value(wave[ispline])
    spline_flux1[ispline] = spline_tmp
    cont_tmp, _ = c_answer.value(wave[ispline])
    cont_flux1[ispline] = cont_tmp
    isrt = np.argsort(wave[indsp], kind='stable')
    s2_interp = scipy.interpolate.interp1d(
        wave[indsp][isrt], sn2[isrt], assume_sorted=False, bounds_error=False, fill_value=0.0
    )
    sn2_1[ispline] = s2_interp(wave[ispline])
    bmask_1d = np.zeros(nspec, dtype=bool)
    bmask_1d[indsp] = bmask2
    spline_flux1 = pydl.djs_maskinterp(spline_flux1, np.logical_not(bmask_1d))
    cmask2_1d = np.zeros(nspec, dtype=bool)
    cmask2_1d[indsp] = cmask
    cont_flux1 = pydl.djs_maskinterp(cont_flux1, np.logical_not(cmask2_1d))
    _, _, sigma1 = astropy.stats.sigma_clipped_stats(
        flux[indsp], sigma_lower=3.0, sigma_upper=5.0
    )

    sn2_med_filt = scipy.ndimage.median_filter(sn2, size=9, mode='reflect')
    if np.any(totmask):
        sn2_interp = scipy.interpolate.interp1d(
            wave[indsp][isrt], sn2_med_filt[isrt], assume_sorted=False, bounds_error=False,
            fill_value='extrapolate'
        )
        sn2_img[totmask] = sn2_interp(waveimg[totmask])
    else:
        log.warning('All pixels are masked')

    log.info(f'sqrt(med(S/N)^2) = {np.sqrt(med_sn2):.2f}')

    # ------------------------------------------------------------------
    # Block 4 — Normalised-object image construction
    # ------------------------------------------------------------------
    if med_sn2 <= 2.0:
        spline_img[totmask] = np.fmax(sigma1, 0)
    else:
        if med_sn2 <= 5.0:
            spline_flux1 = cont_flux1
        badpix = (spline_flux1 <= 0.5) | np.logical_not(bmask_1d)
        goodval = (cont_flux1 > 0.0) & (cont_flux1 < 5e5)
        indbad1 = badpix & goodval
        if indbad1.sum() > 0:
            spline_flux1[indbad1] = cont_flux1[indbad1]
        indbad2 = badpix & np.logical_not(goodval)
        ngood0 = np.logical_not(badpix).sum()
        if indbad2.sum() > 0 or ngood0 > 0:
            spline_flux1[indbad2] = np.median(spline_flux1[np.logical_not(badpix)])
        spline_flux1 = scipy.ndimage.median_filter(spline_flux1, size=5, mode='reflect')

        if np.any(totmask):
            igd = (wave >= wave_min) & (wave <= wave_max)
            isrt1 = np.argsort(wave[igd], kind='stable')
            spline_img_interp = scipy.interpolate.interp1d(
                wave[igd][isrt1], spline_flux1[igd][isrt1], assume_sorted=False,
                bounds_error=False, fill_value='extrapolate'
            )
            spline_img[totmask] = spline_img_interp(waveimg[totmask])
        else:
            spline_img[totmask] = np.fmax(sigma1, 0)

    norm_obj = np.where(spline_img != 0.0, image / spline_img, 0.0)
    norm_ivar = ivar * spline_img ** 2

    ivar_mask = (
        (norm_obj > -0.2) & (norm_obj < 0.7) & totmask & np.isfinite(norm_obj)
        & np.isfinite(norm_ivar)
    )
    norm_ivar[np.logical_not(ivar_mask)] = 0.0
    good = norm_ivar.flatten() > 0.0

    # xtemp: cumulative S/N-weighted spectral coordinate (used for trace correction)
    row_weights = 4.0 + np.sqrt(np.fmax(sn2_1, 0.0))
    xtemp = np.cumsum(np.repeat(row_weights, nspat)).reshape((nspec, nspat))
    xtemp /= xtemp.max()

    return True, med_sn2, norm_obj, norm_ivar, good, xtemp, sn2_img


def _compute_bspline_knots(dspat, sigma, med_sn2, prof_nsigma, good):
    """
    Initialise spatial-coordinate arrays and B-spline knot grid.

    Encapsulates Block 6 of :func:`fit_profile_refactor`: computes the
    spatial-separation image ``dspat``, the normalised coordinate ``sigma_x``,
    the profile half-extent ``limit``, the profile sigma bounds, and the
    sinh-spaced interior B-spline knot positions.

    Parameters
    ----------
    dspat : :class:`numpy.ndarray`
        Spatial separation from the trace, shape
        :math:`(N_{\rm spec}, N_{\rm spat})`.
    sigma : :class:`numpy.ndarray`
        Current Gaussian sigma estimate per spectral pixel,
        shape :math:`(N_{\rm spec},)`.
    med_sn2 : :obj:`float`
        Median :math:`(S/N)^2` of the object.
    prof_nsigma : :obj:`float` or None
        If set, fit the profile out to this many sigma (extended objects);
        overrides the automatic sigma bounds computed from ``med_sn2``.
    good : :class:`numpy.ndarray`
        Boolean flat index array of pixels with positive ``norm_ivar``,
        used to clip the sigma range to the data extent.

    Returns
    -------
    sigma_x : :class:`numpy.ndarray`
        Normalised spatial coordinate (units of ``sigma``), same shape.
    limit : :obj:`float`
        Profile half-extent in units of ``sigma``.
    min_sigma : :obj:`float`
        Lower bound on ``sigma_x`` for profile fitting.
    max_sigma : :obj:`float`
        Upper bound on ``sigma_x`` for profile fitting.
    bkpt : :class:`numpy.ndarray`
        Interior B-spline knot positions in ``sigma_x`` units.
    """
    sigma_x = dspat / sigma[:, None]
    limit = scipy.special.erfcinv(0.1 / np.sqrt(med_sn2)) * np.sqrt(2.0)
    if prof_nsigma is None:
        sinh_space = 0.25 * np.log10(np.fmax(1000.0 / np.sqrt(med_sn2), 10.0))
        abs_sigma = np.fmin((np.abs(sigma_x.flat[good])).max(), 2.0 * limit)
        min_sigma = np.fmax(sigma_x.flat[good].min(), -abs_sigma)
        max_sigma = np.fmin(sigma_x.flat[good].max(), abs_sigma)
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
    bkpt = bkpt[keep]
    return sigma_x, limit, min_sigma, max_sigma, bkpt


def _apodize_profile(
    bset, sigma_x, min_sigma, max_sigma, ss, median_fit, min_level, limit, prof_nsigma, no_deriv
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
        Normalised spatial coordinate image, shape
        :math:`(N_{\rm spec}, N_{\rm spat})`.
    min_sigma : :obj:`float`
        Lower bound on ``sigma_x`` for the good-pixel mask.
    max_sigma : :obj:`float`
        Upper bound on ``sigma_x`` for the good-pixel mask.
    ss : :class:`numpy.ndarray`
        Flat indices that sort ``sigma_x`` in ascending order.
    median_fit : :obj:`float`
        Baseline level of the profile (subtracted before peak-finding).
    min_level : :obj:`float`
        Minimum profile level used to walk the initial limit positions.
    limit : :obj:`float`
        Profile half-extent in units of ``sigma_x``.
    prof_nsigma : :obj:`float` or None
        If set, the fit covers this many sigma (extended objects).  Triggers
        the JXP kludge that zeros the returned limits so QA does not draw
        them, and forces ``no_deriv=True``.
    no_deriv : :obj:`bool`
        If ``True``, skip exponential apodization entirely.

    Returns
    -------
    full_bsp : :class:`numpy.ndarray`
        Flat profile array of length :math:`N_{\rm spec} N_{\rm spat}` with
        exponential tails applied outside ``[l_limit, r_limit]``.
    l_limit : :obj:`float`
        Final left apodization limit in ``sigma_x`` units (0.0 when
        ``prof_nsigma`` is set).
    r_limit : :obj:`float`
        Final right apodization limit in ``sigma_x`` units (0.0 when
        ``prof_nsigma`` is set).
    """
    # ------------------------------------------------------------------
    # Block 11 — Apodization limit search (vectorised)
    # ------------------------------------------------------------------
    igood = (sigma_x.flatten() > min_sigma) & (sigma_x.flatten() < max_sigma)
    full_bsp = np.zeros(sigma_x.size)
    sigma_x_igood = sigma_x.flat[igood]
    yfit_out, _ = bset.value(sigma_x_igood)
    full_bsp[igood] = yfit_out
    isrt2 = sigma_x_igood.argsort(kind='stable')
    peak, peak_x, lwhm, rwhm = _findfwhm(
        yfit_out[isrt2] - median_fit, sigma_x_igood[isrt2]
    )

    left_bool = (
        ((full_bsp[ss] < min_level + median_fit) & (sigma_x.flat[ss] < peak_x))
        | (sigma_x.flat[ss] < peak_x - limit)
    )[::-1]
    ind_left, = np.where(left_bool)
    lp = np.fmax(ind_left.min(), 0) if ind_left.size > 0 else 0

    righ_bool = (
        (full_bsp[ss] < min_level + median_fit) & (sigma_x.flat[ss] > peak_x)
    ) | (sigma_x.flat[ss] > peak_x + limit)
    ind_righ, = np.where(righ_bool)
    rp = np.fmax(ind_righ.min(), 0) if ind_righ.size > 0 else 0

    l_limit = (sigma_x.flat[ss][::-1])[lp] - 0.1
    r_limit = sigma_x.flat[ss[rp]] + 0.1

    # Logarithmic-derivative vectors (pre-computed for both the apodization
    # search and the subsequent limit walk)
    l_lim_vec = np.arange(l_limit + 0.1, -1.0, 0.1)
    l_lim_vec = np.asarray([-1.1]) if l_lim_vec.size == 0 else l_lim_vec
    l_fit1, _ = bset.value(l_lim_vec)
    l_fit2, _ = bset.value(l_lim_vec * 0.9)
    l_deriv_vec = (np.log(l_fit2) - np.log(l_fit1)) / (0.1 * l_lim_vec)
    l_deriv_max = np.fmax(l_deriv_vec.min(), -1.0)

    r_lim_vec = np.arange(r_limit - 0.1, 1.0, -0.1)
    r_lim_vec = np.asarray([1.1]) if r_lim_vec.size == 0 else r_lim_vec
    r_fit1, _ = bset.value(r_lim_vec)
    r_fit2, _ = bset.value(r_lim_vec * 0.9)
    r_deriv_vec = (np.log(r_fit2) - np.log(r_fit1)) / (0.1 * r_lim_vec)
    r_deriv_max = np.fmin(r_deriv_vec.max(), 1.0)

    # Vectorised limit walk (replaces two while-True loops)
    l_cross = np.where(l_deriv_vec <= l_deriv_max)[0]
    if l_cross.size > 0:
        idx = l_cross[0]
        l_limit = l_lim_vec[idx]
        l_deriv = l_deriv_vec[idx]
        l_fit_val = l_fit1[idx]
    else:
        l_limit = -1.0
        l_deriv = l_deriv_vec[-1]
        l_fit_val = l_fit1[-1]

    r_cross = np.where(r_deriv_vec >= r_deriv_max)[0]
    if r_cross.size > 0:
        idx = r_cross[0]
        r_limit = r_lim_vec[idx]
        r_deriv = r_deriv_vec[idx]
        r_fit_val = r_fit1[idx]
    else:
        r_limit = 1.0
        r_deriv = r_deriv_vec[-1]
        r_fit_val = r_fit1[-1]

    # ------------------------------------------------------------------
    # Block 12 — Exponential apodization
    # ------------------------------------------------------------------
    # JXP kludge: zero limits for extended-object fits so QA won't draw them
    if prof_nsigma is not None:
        l_limit = 0.0
        r_limit = 0.0
        no_deriv = True

    if l_deriv < 0 and r_deriv > 0 and not no_deriv:
        left = sigma_x.flatten() < l_limit
        full_bsp[left] = np.exp(-(sigma_x.flat[left] - l_limit) * l_deriv) * l_fit_val
        right = sigma_x.flatten() > r_limit
        full_bsp[right] = np.exp(-(sigma_x.flat[right] - r_limit) * r_deriv) * r_fit_val

    return full_bsp, l_limit, r_limit


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def fit_profile_refactor(
    image, ivar, waveimg, thismask, spat_img, trace_in, wave, flux, fluxivar, inmask=None,
    thisfwhm=4.0, max_trace_corr=2.0, sn_gauss=4.0, percentile_sn2=70.0, prof_nsigma=None,
    no_deriv=False, gauss=False, obj_string='', show_profile=False
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
        Additional boolean mask; defaults to ``(ivar > 0) & thismask``.
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
    # Blocks 2–5 — Spectral fitting and normalisation
    # ------------------------------------------------------------------
    success, med_sn2, norm_obj, norm_ivar, good, xtemp, sn2_img = _fit_spectrum_and_normalize(
        wave, flux, fluxivar, waveimg, image, ivar, totmask, percentile_sn2, _thisfwhm
    )
    if success:
        gauss_kwargs = {'ind': good, 'xtrunc': 7.0}
        ngood = good.sum()
        log.info(f"Gaussian vs b-spline of width {_thisfwhm:.2f} pixels")
    else:
        norm_obj = None
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
            spat_img, norm_obj, trace_in, sigma, _thisfwhm, med_sn2, obj_string, show_profile,
            **gauss_kwargs
        )
        return profile_model, trace_in, fwhmfit, med_sn2

    # ------------------------------------------------------------------
    # Block 6 — Knot setup
    # ------------------------------------------------------------------
    dspat = spat_img - trace_in[:, None]
    sigma_x, limit, min_sigma, max_sigma, bkpt = _compute_bspline_knots(
        dspat, sigma, med_sn2, prof_nsigma, good
    )
    trace_corr = np.zeros(nspec)

    # ------------------------------------------------------------------
    # Block 7 — Initial profile fit
    # ------------------------------------------------------------------
    GOOD_PIX = (sn2_img > sn_gauss**2) & (norm_ivar > 0)
    IN_PIX = (sigma_x >= min_sigma) & (sigma_x <= max_sigma) & (norm_ivar > 0)
    if GOOD_PIX.sum() >= 0.2 * IN_PIX.sum():
        inside, = np.where((GOOD_PIX & IN_PIX).flatten())
    else:
        inside, = np.where(IN_PIX.flatten())

    si = inside[np.argsort(sigma_x.flat[inside], kind='stable')]
    sr = si[::-1]

    bset, bmask, *_ = iterative_bspline_fit(
        sigma_x.flat[si], norm_obj.flat[si], ivar=norm_ivar.flat[si],
        nord=4, kwargs_knots={'interior': bkpt}, maxiter=15, upper=1, lower=1
    )
    mode_fit, _ = bset.value(sigma_x.flat[si])
    median_fit = np.median(norm_obj[norm_ivar > 0.0])
    if np.abs(median_fit) > 0.01:
        log.info(f"Median flux level in profile is not zero: median = {median_fit:.4f}")
    else:
        median_fit = 0.0

    peak, peak_x, lwhm, rwhm = _findfwhm(mode_fit - median_fit, sigma_x.flat[si])
    trace_corr = np.full(nspec, peak_x)
    min_level = peak * np.exp(-0.5 * limit ** 2)

    bspline_fwhm = (rwhm - lwhm) * _thisfwhm / sig2fwhm
    log.info(
        f'Bspline FWHM: {bspline_fwhm:.4f}, compared to initial object finding FWHM: '
        f'{_thisfwhm:.4f}'
    )
    sigma *= (rwhm - lwhm) / sig2fwhm
    limit *= (rwhm - lwhm) / sig2fwhm

    rev_fit = mode_fit[::-1]
    lind, = np.where(
        ((rev_fit < (min_level + median_fit)) & (sigma_x.flat[sr] < peak_x))
        | (sigma_x.flat[sr] < (peak_x - limit))
    )
    l_limit = sigma_x.flat[sr[lind.min()]] if lind.size > 0 else min_sigma

    rind, = np.where(
        ((mode_fit < (min_level + median_fit)) & (sigma_x.flat[si] > peak_x))
        | (sigma_x.flat[si] > (peak_x + limit))
    )
    r_limit = sigma_x.flat[si[rind.min()]] if rind.size > 0 else max_sigma

    log.info(
        f"Trace limits: limit={limit:.4f}, min_level={min_level:.4f}, l_limit={l_limit:.4f}, "
        f"r_limit={r_limit:.4f}"
    )

    mask = np.zeros(nspec * nspat, dtype=bool)
    mask[si] = (norm_ivar.flat[si] > 0) & (np.abs(norm_obj.flat[si] - mode_fit) < 0.1)
    inside, = np.where(
        (sigma_x.flat[si] > l_limit) & (sigma_x.flat[si] < r_limit) & mask[si]
    )
    ninside = inside.size

    if ninside < 10:
        log.info("Too few pixels inside l_limit and r_limit")
        log.info("Returning Gaussian profile")
        profile_model = _return_gaussian(
            spat_img, norm_obj, trace_in + trace_corr * sigma, sigma, bspline_fwhm, med_sn2,
            obj_string, show_profile, ind=good, l_limit=l_limit, r_limit=r_limit, xlim=7.0,
        )
        return profile_model, trace_in, fwhmfit, med_sn2

    # ------------------------------------------------------------------
    # Block 8 — Iterative trace and width correction (sigma_iter=3)
    # ------------------------------------------------------------------
    sigma_iter = 3
    isort = xtemp.flat[si[inside]].argsort(kind='stable')
    inside = si[inside[isort]]
    pb = np.ones(inside.size)
    area = np.ones(nspec)

    for iiter in range(sigma_iter):
        mode_zero, _ = bset.value(sigma_x.flat[inside])
        mode_zero = mode_zero * pb

        mode_min05, _ = bset.value(sigma_x.flat[inside] - 0.5)
        mode_plu05, _ = bset.value(sigma_x.flat[inside] + 0.5)
        mode_shift = (mode_min05 - mode_plu05) * pb
        mode_shift[
            (sigma_x.flat[inside] <= (l_limit + 0.5)) | (sigma_x.flat[inside] >= (r_limit - 0.5))
        ] = 0.0

        mode_by13, _ = bset.value(sigma_x.flat[inside] / 1.3)
        mode_stretch = mode_by13 * pb / 1.3 - mode_zero

        nbkpts = int(np.log10(np.fmax(med_sn2, 11.0)))
        xx = xtemp.sum(axis=1) / nspat
        profile_basis = np.column_stack((mode_zero, mode_shift))

        mode_shift_out = iterative_bspline_fit(
            xtemp.flat[inside], norm_obj.flat[inside], ivar=norm_ivar.flat[inside],
            basis=profile_basis, maxiter=1, kwargs_knots={'count': nbkpts}
        )
        if not np.any(mode_shift_out[1]):
            log.info(
                'Returning a Gaussian profile because the B-spline fit to trace correction failed '
                f'for ninside = {ninside}.'
            )
            profile_model = _return_gaussian(
                spat_img, norm_obj, trace_in + trace_corr * sigma, sigma, bspline_fwhm, med_sn2,
                obj_string, show_profile, ind=good, l_limit=l_limit, r_limit=r_limit, xlim=7.0,
            )
            return profile_model, trace_in, fwhmfit, med_sn2

        mode_shift_set = mode_shift_out[0]
        temp_set = BSpline(knots=Knots(full=mode_shift_set.breakpoints), nord=mode_shift_set.nord)
        temp_set.bkpt_gpm = mode_shift_set.bkpt_gpm
        temp_set.coeff = mode_shift_set.coeff[:, 0]
        h0, _ = temp_set.value(xx)
        temp_set.coeff = mode_shift_set.coeff[:, 1]
        h1, _ = temp_set.value(xx)
        ratio_10 = h1 / (h0 + (h0 == 0.0))
        delta_trace_corr = ratio_10 / (1.0 + np.abs(ratio_10) / 0.1)
        trace_corr = trace_corr + delta_trace_corr

        profile_basis = np.column_stack((mode_zero, mode_stretch))
        mode_stretch_out = iterative_bspline_fit(
            xtemp.flat[inside], norm_obj.flat[inside], ivar=norm_ivar.flat[inside],
            basis=profile_basis, maxiter=1, kwargs_knots={'full': mode_shift_set.breakpoints}
        )
        if not np.any(mode_stretch_out[1]):
            log.info(
                'Returning a Gaussian profile because the B-spline fit to the width correction '
                f'failed for ninside = {ninside}.'
            )
            profile_model = _return_gaussian(
                spat_img, norm_obj, trace_in + trace_corr * sigma, sigma, bspline_fwhm, med_sn2,
                obj_string, show_profile, ind=good, l_limit=l_limit, r_limit=r_limit, xlim=7.0,
            )
            return profile_model, trace_in, fwhmfit, med_sn2

        mode_stretch_set = mode_stretch_out[0]
        temp_set = BSpline(
            knots=Knots(full=mode_stretch_set.breakpoints), nord=mode_stretch_set.nord
        )
        temp_set.bkpt_gpm = mode_stretch_set.bkpt_gpm
        temp_set.coeff = mode_stretch_set.coeff[:, 0]
        h0, _ = temp_set.value(xx)
        temp_set.coeff = mode_stretch_set.coeff[:, 1]
        h2, _ = temp_set.value(xx)
        h0 = np.fmax(h0 + h2 * mode_stretch.sum() / mode_zero.sum(), 0.1)
        ratio_20 = h2 / (h0 + (h0 == 0.0))
        sigma_factor = 0.3 * ratio_20 / (1.0 + np.abs(ratio_20))

        log.info(f"Iteration # {iiter+1}")
        log.info(
            f"Median abs value of trace correction = {np.median(np.abs(delta_trace_corr)):.3f}"
        )
        log.info(f"Median abs value of width correction = {np.median(np.abs(sigma_factor)):.3f}")

        sigma = sigma * (1.0 + sigma_factor)
        area = area * h0 / (1.0 + sigma_factor)

        sigma_x = dspat / sigma[:, None] - trace_corr[:, None]

        # TODO: In the previous version, iiter was iterated from 1 to sigma_iter
        # and this was check was for if iiter < sigma_iter - 1.  When moving to
        # a 0-indexed iterator, this makes this if statement sigma_iter - 2.  Is
        # that the original intention, or was this an index mistake?
        if iiter < sigma_iter - 2:
            ss = sigma_x.flat[inside].argsort(kind='stable')
            pb = np.repeat(area, nspat)[inside]
            keep = (bkpt >= sigma_x.flat[inside].min()) & (bkpt <= sigma_x.flat[inside].max())
            if keep.sum() == 0:
                keep = np.ones(bkpt.size, dtype=bool)
            bset_out = iterative_bspline_fit(
                sigma_x.flat[inside[ss]], norm_obj.flat[inside[ss]],
                ivar=norm_ivar.flat[inside[ss]], basis=pb[ss], nord=4,
                kwargs_knots={'interior': bkpt[keep]}, maxiter=2
            )
            if not np.any(bset_out[1]):
                log.info(
                    'Returning a Gaussian profile because the B-spline profile fit in '
                    f'trace/width loop failed for ninside = {ninside}'
                )
                profile_model = _return_gaussian(
                    spat_img, norm_obj, trace_in + trace_corr * sigma, sigma, bspline_fwhm,
                    med_sn2, obj_string, show_profile, ind=good, l_limit=l_limit, r_limit=r_limit,
                    xlim=7.0,
                )
                return profile_model, trace_in, fwhmfit, med_sn2

            bset = _bspline2d_to_1d(bset_out[0])

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
    ss = sigma_x.flatten().argsort(kind='stable')
    inside, = np.where(
        (sigma_x.flat[ss] >= min_sigma) & (sigma_x.flat[ss] <= max_sigma) & mask[ss]
        & np.isfinite(norm_obj.flat[ss]) & np.isfinite(norm_ivar.flat[ss])
    )
    pb = area[:, None] * np.ones((nspec, nspat))
    bset_out = iterative_bspline_fit(
        sigma_x.flat[ss[inside]], norm_obj.flat[ss[inside]],
        ivar=norm_ivar.flat[ss[inside]], basis=pb.flat[ss[inside]], nord=4,
        kwargs_knots={'interior': bkpt}, upper=10, lower=10
    )
    bset = _bspline2d_to_1d(bset_out[0])
    outmask = bset_out[1]

    # ------------------------------------------------------------------
    # Blocks 11–12 — Apodization limit search and exponential tails
    # ------------------------------------------------------------------
    full_bsp, l_limit, r_limit = _apodize_profile(
        bset, sigma_x, min_sigma, max_sigma, ss, median_fit, min_level,
        limit, prof_nsigma, no_deriv
    )

    # ------------------------------------------------------------------
    # Block 13 — Normalisation, logging, QA, return
    # ------------------------------------------------------------------
    full_bsp = full_bsp.reshape(nspec, nspat)
    profile_model = full_bsp * pb

    res_mode = (
        (norm_obj.flat[ss[inside]] - profile_model.flat[ss[inside]])
        * np.sqrt(norm_ivar.flat[ss[inside]])
    )
    chi_good = outmask & (norm_ivar.flat[ss[inside]] > 0)
    chi_med = np.median(res_mode[chi_good] ** 2)

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
    if row_sums.sum() > 0.0:
        norm_safe = np.where(row_sums[:, None] > 0.0, row_sums[:, None], 1.0)
        profile_model = np.where(row_sums[:, None] > 0.0, profile_model / norm_safe, 0.0)

    info_string = (
        f"FWHM range:{fwhmfit.min():.2f} - {fwhmfit.max():.2f}, S/N={np.sqrt(med_sn2):.3f}"
        f", median(chi^2)={chi_med:.3f}"
    )
    if show_profile:
        qa_fit_profile(
            sigma_x, norm_obj / (pb + (pb == 0.0)), full_bsp, l_limit=l_limit, r_limit=r_limit,
            ind=ss[inside], xlim=prof_nsigma, title=f'{obj_string} {info_string}'
        )

    return profile_model, xnew, fwhmfit, med_sn2
