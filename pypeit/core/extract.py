"""
Module for PypeIt extraction code

.. include:: ../include/links.rst

"""

from IPython import embed
import numpy as np

from pypeit import log
from pypeit import PypeItError
from pypeit import utils
from pypeit.core import procimg
from pypeit.core.moment import moment1d


def extract_optimal(
    imgminsky, ivar, mask, waveimg, skyimg, thismask, oprof, min_frac_use=0.9, fwhmimg=None,
    flatimg=None, base_var=None, count_scale=None, noise_floor=None, box_radius=None,
    trace_spec=None, trace_spat=None
):
    r"""
    Perform optimal extraction `(Horne 1986)
    <https://ui.adsabs.harvard.edu/abs/1986PASP...98..609H/abstract>`__.

    When necessary, the function will fall back to a boxcar extraction algorithm
    for setting the wavelength for a wavelength channel that is fully masked.
    If this happens, the ``box_radius``, ``trace_spat``, and ``trace_spec``
    parameters *must* be provided; an exception is raised otherwise.

    Parameters
    ----------
    imgminsky : `numpy.ndarray`_
        Floating-point science image minus skymodel (i.e., imgminsky = sciimg - skyimg)
        with shape :math:`(N_{\rm spec}, N_{\rm spat})`.
        The first dimension (:math:`N_{\rm spec}`) is spectral, and second dimension
        (:math:`N_{\rm spat}`) is spatial.
    ivar : `numpy.ndarray`_
        Floating-point inverse variance image for the science image.
        It can be a model image, or deduced from ``sciimg``. Shape
        must match ``sciimg``, :math:`(N_{\rm spec}, N_{\rm spat})`.
    mask : `numpy.ndarray`_
        Boolean image representing the good-pixel mask for the science image.
        The pixels that have value of True are good to be used.
        Must have the same shape as ``sciimg``, :math:`(N_{\rm spec}, N_{\rm spat})`.
    waveimg : `numpy.ndarray`_
        Floating-point wavelength image. Must have the same shape as ``sciimg``,
        :math:`(N_{\rm spec}, N_{\rm spat})`.
    skyimg : `numpy.ndarray`_
        Floating-point image containing the modeled sky.
        Must have the same shape as ``sciimg``, :math:`(N_{\rm spec}, N_{\rm spat})`.
    thismask : `numpy.ndarray`_
        Boolean image indicating which pixels are on the slit/order in question.
        Must have the same shape as ``sciimg``, :math:`(N_{\rm spec}, N_{\rm spat})`.
    oprof : `numpy.ndarray`_
         Floating-point image containing the profile of the object that is
         going to be extracted. Must have the same shape as ``sciimg``,
         :math:`(N_{\rm spec}, N_{\rm spat})`.
    min_frac_use : :obj:`float`, optional
        Minimum accepted value for the sum of the normalized object profile across the spatial direction.
        For each spectral pixel, if the majority of the object profile has been masked, i.e.,
        the sum of the normalized object profile across the spatial direction is less than `min_frac_use`,
        the optimal extraction will also be masked. The default value is 0.05.
    fwhmimg : `numpy.ndarray`_, None, optional:
        Floating-point image containing the modeled spectral FWHM (in pixels) at every pixel location.
        Must have the same shape as ``sciimg``, :math:`(N_{\rm spec}, N_{\rm spat})`.
    flatimg : `numpy.ndarray`_, None, optional:
        Floating-point image containing the unnormalized flat-field image. This image
        is used to extract the blaze function. Must have the same shape as ``sciimg``,
        :math:`(N_{\rm spec}, N_{\rm spat})`.
    base_var : `numpy.ndarray`_, optional
        Floating-point "base-level" variance image set by the detector properties and
        the image processing steps. See :func:`~pypeit.core.procimg.base_variance`.
        Must have the same shape as ``sciimg``, :math:`(N_{\rm spec}, N_{\rm spat})`.
    count_scale : :obj:`float` or `numpy.ndarray`_, optional
        A scale factor, :math:`s`, that *has already been applied* to the
        provided science image. It accounts for the number of frames contributing to
        the provided counts, and the relative throughput factors that can be measured
        from flat-field frames. For example, if the image has been flat-field
        corrected, this is the inverse of the flat-field counts.  If None, set
        to 1.  If a single float, assumed to be constant across the full image.
        If an array, the shape must match ``base_var``.  The variance will be 0
        wherever :math:`s \leq 0`, modulo the provided ``adderr``.  This is one
        of the components needed to construct the model variance; see
        ``model_noise``.
    noise_floor : :obj:`float`, optional
        A fraction of the counts to add to the variance, which has the effect of
        ensuring that the S/N is never greater than ``1/noise_floor``; see
        :func:`~pypeit.core.procimg.variance_model`.  If None, no noise floor is
        added.
    box_radius : :obj:`float`, optional
        When necessary, use this boxcar radius in pixels to determine the
        wavelength of a fully masked spectral channel.
    trace_spec : :class:`numpy.ndarray`, optional
        When necessary, use these spectral trace positions to determine the
        wavelength of a fully masked spectral channel.
    trace_spat : :class:`numpy.ndarray`, optional
        When necessary, use these spatial trace positions to determine the
        wavelength of a fully masked spectral channel.

    Returns
    -------
    wave : :class:`numpy.ndarray`
        Wavelength vector
    counts : :class:`numpy.ndarray`
        Extracted counts
    counts_ivar : :class:`numpy.ndarray`
        Inverse variance in the extracted counts
    counts_sig : :class:`numpy.ndarray`
        1-sigma uncertainties in the extracted counts
    counts_nivar : :class:`numpy.ndarray`
        Same as ``counts_nivar``, but excludes shot-noise contributions from the
        extracted object.
    gpm : :class:`numpy.ndarray`
        Good-pixel mask selected valid data in the extracted spectrum.
    fwhm : :class:`numpy.ndarray`
        Estimated FWHM of the wavelength-dependent spectral resolution element.
        This will be None if ``fwhmimg`` is not provided.
    flat : :class:`numpy.ndarray`
        Spectrum extracted from the normalized flat-field image at the same
        location as the object spectrum.  This will be None if ``flatimg`` is
        not provided.
    counts_sky : :class:`numpy.ndarray`
        Sky-only spectrum extraced from the sky image (``skyimg``).
    counts_sig_det : :class:`numpy.ndarray`
        Same as ``counts_sig`` but only includes uncertainties introduced by the
        detector.
    frac_use : :class:`numpy.ndarray`
        Wavelength-dependent fraction of pixels in the object profile subimage
        used for the extraction.
    chi2 : :class:`numpy.ndarray`
        Wavelength-dependent reduced chi-square of the model fit.
    """
    # TODO This makes no sense for difference imaging? Not sure we need NIVAR anyway
    var_no = None if base_var is None \
                else procimg.variance_model(base_var, counts=skyimg, count_scale=count_scale,
                                            noise_floor=noise_floor)

    ispec, ispat = np.where(oprof > 0.0)

    # Exit gracefully if we have no positive object profiles, since that means something was wrong with object fitting
    if not np.any(oprof > 0.0):
        log.warning('Object profile is zero everywhere. This aperture is junk.')
        return

    mincol = np.min(ispat)
    maxcol = np.max(ispat) + 1
    nsub = maxcol - mincol

    mask_sub = mask[:,mincol:maxcol]
    thismask_sub = thismask[:, mincol:maxcol]
    wave_sub = waveimg[:,mincol:maxcol]
    ivar_sub = np.fmax(ivar[:,mincol:maxcol],0.0) # enforce positivity since these are used as weights
    vno_sub = None if var_no is None else np.fmax(var_no[:,mincol:maxcol],0.0)

    base_sub = None if base_var is None else base_var[:,mincol:maxcol]
    img_sub = imgminsky[:,mincol:maxcol]
    sky_sub = skyimg[:,mincol:maxcol]
    oprof_sub = oprof[:,mincol:maxcol]
    if fwhmimg is not None:
        fwhmimg_sub = fwhmimg[:,mincol:maxcol]
    if flatimg is not None:
        flatimg_sub = flatimg[:,mincol:maxcol]
    # enforce normalization and positivity of object profiles
    norm = np.nansum(oprof_sub,axis = 1)
    norm_oprof = np.outer(norm, np.ones(nsub))
    oprof_sub = np.fmax(oprof_sub/norm_oprof, 0.0)

    ivar_denom = np.nansum(mask_sub*oprof_sub, axis=1)
    mivar_num = np.nansum(mask_sub*ivar_sub*oprof_sub**2, axis=1)
    mivar_opt = mivar_num/(ivar_denom + (ivar_denom == 0.0))
    flux_opt = np.nansum(mask_sub*ivar_sub*img_sub*oprof_sub, axis=1)/(mivar_num + (mivar_num == 0.0))
    # Optimally extracted noise variance (sky + read noise) only. Since
    # this variance is not the same as that used for the weights, we
    # don't get the usual cancellation. Additional denom factor is the
    # analog of the numerator in Horne's variance formula. Note that we
    # are only weighting by the profile (ivar_sub=1) because
    # otherwise the result depends on the signal (bad).
    nivar_num = np.nansum(mask_sub*oprof_sub**2, axis=1) # Uses unit weights
    if vno_sub is None:
        nivar_opt = None
    else:
        nvar_opt = ivar_denom * np.nansum(mask_sub * vno_sub * oprof_sub**2, axis=1) \
                            / (nivar_num**2 + (nivar_num**2 == 0.0))
        nivar_opt = 1.0/(nvar_opt + (nvar_opt == 0.0))
    # Optimally extract sky and (read noise)**2 in a similar way
    sky_opt = ivar_denom*(np.nansum(mask_sub*sky_sub*oprof_sub**2, axis=1))/(nivar_num**2 + (nivar_num**2 == 0.0))
    if base_var is None:
        base_opt = None
    else:
        base_opt = ivar_denom * np.nansum(mask_sub * base_sub * oprof_sub**2, axis=1) \
                        / (nivar_num**2 + (nivar_num**2 == 0.0))
        base_opt = np.sqrt(base_opt)
        base_opt[np.isnan(base_opt)]=0.0

    tot_weight = np.nansum(mask_sub*ivar_sub*oprof_sub, axis=1)
    prof_norm = np.nansum(oprof_sub, axis=1)
    # NOTE: Frac_use is also equal to np.nansum(mask_sub * oprof_sub, axis=1)
    frac_use = (prof_norm > 0.0)*np.nansum((mask_sub*ivar_sub > 0.0)*oprof_sub, axis=1)/(prof_norm + (prof_norm == 0.0))

    # Use the same weights = oprof^2*mivar for the wavelenghts as the flux.
    # Note that for the flux, one of the oprof factors cancels which does
    # not for the wavelengths.
    wave_opt = np.nansum(mask_sub*ivar_sub*wave_sub*oprof_sub**2, axis=1)/(mivar_num + (mivar_num == 0.0))
    mask_opt = (tot_weight > 0.0) & (frac_use > min_frac_use) & (mivar_num > 0.0) & (ivar_denom > 0.0) & \
               np.isfinite(wave_opt) & (wave_opt > 0.0)
    fwhm_opt = None
    if fwhmimg is not None:
        fwhm_opt = np.nansum(mask_sub*ivar_sub*fwhmimg_sub*oprof_sub, axis=1) * utils.inverse(tot_weight)
    blaze_opt = None
    if flatimg is not None:
        blaze_opt = np.nansum(mask_sub*ivar_sub*flatimg_sub*oprof_sub, axis=1) * utils.inverse(tot_weight)
    # Interpolate wavelengths over masked pixels
    badwvs = (mivar_num <= 0) | np.logical_not(np.isfinite(wave_opt)) | (wave_opt <= 0.0)
    if badwvs.any():
        oprof_smash = np.nansum(thismask_sub*oprof_sub**2, axis=1)
        # Can we use the profile average wavelengths instead?
        oprof_good = badwvs & (oprof_smash > 0.0)
        if oprof_good.any():
            wave_opt[oprof_good] = np.nansum(
                wave_sub[oprof_good,:]*thismask_sub[oprof_good,:]*oprof_sub[oprof_good,:]**2, axis=1)/\
                                   np.nansum(thismask_sub[oprof_good,:]*oprof_sub[oprof_good,:]**2, axis=1)
        oprof_bad = badwvs & ((oprof_smash <= 0.0) | (np.isfinite(oprof_smash) == False) | (wave_opt <= 0.0) | (np.isfinite(wave_opt) == False))
        if oprof_bad.any():
            # If there are no good profile wavelengths, use boxcar wavelengths for these pixels
            # get boxcar_radius
            if None in [box_radius, trace_spec, trace_spat]:
                raise PypeItError(
                    'Fully masked wavelength channels detected; must provide box_radius, '
                    'trace_spec, and trace_spat for fallback boxcar determination of wavelength.'
                )
            if trace_spec.shape != trace_spat.shape:
                raise PypeItError(
                    'Spectral and spatial locations of the extraction trace must have the same '
                    'length.'
                )
            box_denom_no_mask = moment1d(waveimg > 0.0, trace_spat, 2 * box_radius, row=trace_spec)[0]
            wave_no_mask = moment1d(waveimg, trace_spat, 2 * box_radius, row=trace_spec)[0] / (
                        box_denom_no_mask + (box_denom_no_mask == 0.0))
            wave_opt[oprof_bad] = wave_no_mask[oprof_bad]

    flux_model = np.outer(flux_opt,np.ones(nsub))*oprof_sub
    chi2_num = np.nansum((img_sub - flux_model)**2*ivar_sub*mask_sub,axis=1)
    chi2_denom = np.fmax(np.nansum(ivar_sub*mask_sub > 0.0, axis=1) - 1.0, 1.0)
    chi2 = chi2_num/chi2_denom

    # Calculate the Angstroms/pixel and Spectral FWHM
    if fwhm_opt is not None:
        fwhm_opt *= np.gradient(wave_opt)  # Convert pixel FWHM to Angstroms
    # Normalize the blaze function
    if blaze_opt is not None:
        blaze_opt /= np.nanmax(blaze_opt)

    _ivar = mivar_opt*np.logical_not(badwvs)
    _sig = np.sqrt(utils.inverse(_ivar))
    return (
        wave_opt,
        flux_opt,
        _ivar,
        _sig,
        None if nivar_opt is None else nivar_opt*np.logical_not(badwvs),
        mask_opt*np.logical_not(badwvs),
        fwhm_opt,
        blaze_opt,
        sky_opt,
        base_opt,
        frac_use,
        chi2
    )


def extract_asym_boxcar(sciimg, left_trace, righ_trace, gpm=None, ivar=None):
    r"""
    Perform asymmetric boxcar extraction of the flux between two traces.

    Parameters
    ----------
    sciimg : `numpy.ndarray`_
        Floating-point science image with shape :math:`(N_{\rm spec}, N_{\rm spat})`.
        The first dimension (:math:`N_{\rm spec}`) is spectral, and second dimension
        (:math:`N_{\rm spat}`) is spatial.
    left_trace, right_trace : `numpy.ndarray`_
        Left and right trace boundaries of the extraction region for each aperture.
        They are 2-d floating-point arrays with shape :math:`(N_{\rm spec}, N_{\rm apertures})`.
    gpm : `numpy.ndarray`_, optional
        Boolean image representing the good-pixel mask for the science image.
        The pixels that have value of True are good to be used.
        Must have the same shape as ``sciimg``, :math:`(N_{\rm spec}, N_{\rm spat})`.
    ivar : `numpy.ndarray`_, optional
        Floating-point inverse variance image for the science image.
        It can be a model image, or deduced from ``sciimg``. Shape
        must match ``sciimg``, :math:`(N_{\rm spec}, N_{\rm spat})`.
        If not None, the inverse variance of the boxcar extracted flux
        will be returned.

    Returns
    -------
    flux_out : `numpy.ndarray`_
        2-d floating-point array containing, for each aperture, the boxcar
        extracted flux as a function of spectral position. Shape is
        :math:`(N_{\rm spec}, N_{\rm apertures})`.
    gpm_box : `numpy.ndarray`_
        2-d Boolean-point array representing, for each aperture, the good-pixel
        mask for the boxcar extracted flux. The pixels that have value of True
        are good to be used. Shape is :math:`(N_{\rm spec}, N_{\rm apertures})`.
    box_npix :  `numpy.ndarray`_
        2-d floating-point array containing, for each aperture, the number of pixels in
        each spectral position that contributed to the boxcar sum of the flux.
        Shape is :math:`(N_{\rm spec}, N_{\rm apertures})`.
    ivar_out : `numpy.ndarray`_
        2-d floating-point array containing, for each aperture, the inverse variance
        of the boxcar extracted flux as a function of spectral position. Shape is
        :math:`(N_{\rm spec}, N_{\rm apertures})`. This is  only be returned if
        the input parameter `ivar` is not None.

    """
    ivar1 = np.ones_like(sciimg) if ivar is None else ivar
    gpm1 = ivar1 > 0.0 if gpm is None else gpm

    flux_box = moment1d(sciimg*gpm1, (left_trace+righ_trace)/2.0, (righ_trace-left_trace))[0]
    #box_denom = moment1d(gpm1, (left_trace+righ_trace)/2.0, (righ_trace-left_trace))[0]

    pixtot = moment1d(sciimg*0 + 1.0, (left_trace+righ_trace)/2.0, (righ_trace-left_trace))[0]
    pixmsk = moment1d(ivar1*gpm1 == 0.0, (left_trace+righ_trace)/2.0, (righ_trace-left_trace))[0]

    # If every pixel is masked then mask the boxcar extraction
    gpm_box = (pixmsk != pixtot)

    varimg = 1.0 / (ivar1 + (ivar1 == 0.0))
    var_box = moment1d(varimg * gpm1, (left_trace+righ_trace)/2.0, (righ_trace-left_trace))[0]

    ivar_box = 1.0/(var_box + (var_box == 0.0))

    flux_out = flux_box*gpm_box
    ivar_out = ivar_box*gpm_box
    box_npix = pixtot - pixmsk

    if ivar is None:
        return flux_out, gpm_box, box_npix
    else:
        return flux_out, gpm_box, box_npix, ivar_out


def extract_boxcar(
    box_radius, trace_spat, imgminsky, ivar, mask, waveimg, skyimg, fwhmimg=None,
    flatimg=None, base_var=None, count_scale=None, noise_floor=None, trace_spec=None
):
    r"""
    Perform boxcar extraction.

    Parameters
    ----------
    box_radius : :obj:`float`
        The boxcar radius (half of the full width) to use for the extraction.
    trace_spat : :class:`numpy.ndarray`
        The spatial pixel to use for the center of the extraction aperture as a
        function of the spectral dimension.  That is, these are the pixels in
        the *second* axis of the input image at which to perform the extraction.
        If ``trace_spec`` is None, the shape must be :math:`(N_{\rm spec},)`;
        otherwise, the shape must match ``trace_spec``.
    imgminsky : `numpy.ndarray`_
        Floating-point science image minus skymodel (i.e., imgminsky = sciimg - skyimg)
        with shape :math:`(N_{\rm spec}, N_{\rm spat})`.
        The first dimension (:math:`N_{\rm spec}`) is spectral, and second dimension
        (:math:`N_{\rm spat}`) is spatial.
    ivar : `numpy.ndarray`_
        Floating-point inverse variance image for the science image.
        It can be a model image, or deduced from ``sciimg``. Shape
        must match ``sciimg``, :math:`(N_{\rm spec}, N_{\rm spat})`.
    mask : `numpy.ndarray`_
        Boolean image representing the good-pixel mask for the science image.
        The pixels that have value of True are good to be used.
        Must have the same shape as ``sciimg``, :math:`(N_{\rm spec}, N_{\rm spat})`.
    waveimg : `numpy.ndarray`_
        Floating-point wavelength image. Must have the same shape as ``sciimg``,
        :math:`(N_{\rm spec}, N_{\rm spat})`.
    skyimg : `numpy.ndarray`_
        Floating-point image containing the modeled sky.
        Must have the same shape as ``sciimg``, :math:`(N_{\rm spec}, N_{\rm spat})`.
    fwhmimg : `numpy.ndarray`_, None, optional
        Floating-point image containing the modeled spectral FWHM (in pixels) at every pixel location.
        Must have the same shape as ``sciimg``, :math:`(N_{\rm spec}, N_{\rm spat})`.
    flatimg : `numpy.ndarray`_, None, optional
        Floating-point image containing the normalized flat-field. Must have the same shape as
        ``sciimg``, :math:`(N_{\rm spec}, N_{\rm spat})`.
    base_var : `numpy.ndarray`_, optional
        Floating-point "base-level" variance image set by the detector properties and
        the image processing steps. See :func:`~pypeit.core.procimg.base_variance`.
        Must have the same shape as ``sciimg``, :math:`(N_{\rm spec}, N_{\rm spat})`.
    count_scale : :obj:`float` or `numpy.ndarray`_, optional
        A scale factor, :math:`s`, that *has already been applied* to the
        provided science image. It accounts for the number of frames contributing to
        the provided counts, and the relative throughput factors that can be measured
        from flat-field frames plus a scaling factor applied if the counts of each frame are
        scaled to the mean counts of all frames. For example, if the image has been flat-field
        corrected, this is the inverse of the flat-field counts.  If None, set
        to 1.  If a single float, assumed to be constant across the full image.
        If an array, the shape must match ``base_var``.  The variance will be 0
        wherever :math:`s \leq 0`, modulo the provided ``adderr``.  This is one
        of the components needed to construct the model variance; see
        ``model_noise``.
    noise_floor : :obj:`float`, optional
        A fraction of the counts to add to the variance, which has the effect of
        ensuring that the S/N is never greater than ``1/noise_floor``; see
        :func:`~pypeit.core.procimg.variance_model`.  If None, no noise floor is
        added.
    trace_spec : :class:`numpy.ndarray`, optional
        The pixel locations along the spectral axis at which to place each
        aperture defined by ``trace_spat``.  If provided, the shape must match
        ``trace_spat``.  If None, the code assumes ``trace_spat`` covers the
        full image and runs from 0 to the number of spectral pixels.

    Returns
    -------
    wave : :class:`numpy.ndarray`
        Wavelength vector
    counts : :class:`numpy.ndarray`
        Extracted counts
    counts_ivar : :class:`numpy.ndarray`
        Inverse variance in the extracted counts
    counts_sig : :class:`numpy.ndarray`
        1-sigma uncertainties in the extracted counts
    counts_nivar : :class:`numpy.ndarray`
        Same as ``counts_nivar``, but excludes shot-noise contributions from the
        extracted object.
    gpm : :class:`numpy.ndarray`
        Good-pixel mask selected valid data in the extracted spectrum.
    fwhm : :class:`numpy.ndarray`
        Estimated FWHM of the wavelength-dependent spectral resolution element.
        This will be None if ``fwhmimg`` is not provided.
    flat : :class:`numpy.ndarray`
        Spectrum extracted from the normalized flat-field image at the same
        location as the object spectrum.  This will be None if ``flatimg`` is
        not provided.
    counts_sky : :class:`numpy.ndarray`
        Sky-only spectrum extraced from the sky image (``skyimg``).
    counts_sig_det : :class:`numpy.ndarray`
        Same as ``counts_sig`` but only includes uncertainties introduced by the
        detector.
    npix : :class:`numpy.ndarray`
        Wavelength-dependent (fractional) number of valid pixels included in the
        extraction.
    """
    if trace_spec is None:
        trace_spec = np.arange(imgminsky.shape[0])
    if trace_spec.shape != trace_spat.shape:
        raise PypeItError(
            'Spectral and spatial locations of the extraction trace must have the same length.'
        )

    # TODO This makes no sense for difference imaging? Not sure we need NIVAR anyway
    var_no = None if base_var is None \
                else procimg.variance_model(base_var, counts=skyimg, count_scale=count_scale,
                                            noise_floor=noise_floor)

    # Fill in the boxcar extraction tags
    flux_box = moment1d(imgminsky*mask, trace_spat, 2*box_radius, row=trace_spec)[0]
    # Denom is computed in case the trace goes off the edge of the image
    box_denom = moment1d(waveimg*mask > 0.0, trace_spat, 2*box_radius, row=trace_spec)[0]
    wave_box = moment1d(waveimg*mask, trace_spat, 2*box_radius,
                        row=trace_spec)[0] / (box_denom + (box_denom == 0.0))
    fwhm_box = None
    if fwhmimg is not None:
        fwhm_box = moment1d(fwhmimg*mask, trace_spat, 2*box_radius, row=trace_spec)[0]
    blaze_box = None
    if flatimg is not None:
        blaze_box = moment1d(flatimg*mask, trace_spat, 2*box_radius, row=trace_spec)[0]
    varimg = 1.0/(ivar + (ivar == 0.0))
    var_box = moment1d(varimg*mask, trace_spat, 2*box_radius, row=trace_spec)[0]
    nvar_box = None if var_no is None \
                else moment1d(var_no*mask, trace_spat, 2*box_radius, row=trace_spec)[0]
    sky_box = moment1d(skyimg*mask, trace_spat, 2*box_radius, row=trace_spec)[0]
    if base_var is None:
        base_box = None
    else:
        _base_box = moment1d(base_var*mask, trace_spat, 2*box_radius, row=trace_spec)[0]
        base_posind = (_base_box > 0.0)
        base_box = np.zeros(_base_box.shape, dtype=float)
        base_box[base_posind] = np.sqrt(_base_box[base_posind])
    pixtot = moment1d(ivar*0 + 1.0, trace_spat, 2*box_radius, row=trace_spec)[0]
    pixmsk = moment1d(ivar*mask == 0.0, trace_spat, 2*box_radius, row=trace_spec)[0]
    # If every pixel is masked then mask the boxcar extraction
    mask_box = (pixmsk != pixtot) & np.isfinite(wave_box) & (wave_box > 0.0)
    bad_box = (wave_box <= 0.0) | np.logical_not(np.isfinite(wave_box)) | (box_denom == 0.0)
    # interpolate bad wavelengths over masked pixels
    if bad_box.any():
        box_denom_no_mask = moment1d(waveimg > 0.0, trace_spat, 2 * box_radius, row=trace_spec)[0]
        wave_no_mask = moment1d(waveimg, trace_spat, 2 * box_radius, row=trace_spec)[0] / (box_denom_no_mask + (box_denom_no_mask == 0.0))
        wave_box[bad_box] = wave_no_mask[bad_box]

    ivar_box = 1.0/(var_box + (var_box == 0.0))
    nivar_box = None if nvar_box is None else 1.0/(nvar_box + (nvar_box == 0.0))

    # Calculate the Angstroms/pixel and the final spectral FWHM value
    if fwhm_box is not None:
        ang_per_pix = np.gradient(wave_box)
        fwhm_box *= ang_per_pix * utils.inverse(pixtot - pixmsk)  # Need to divide by total number of unmasked pixels
    # Normalize the blaze function
    if blaze_box is not None:
        blaze_box *= utils.inverse(pixtot - pixmsk)  # Need to divide by total number of unmasked pixels
        blaze_box *= utils.inverse(np.nanmax(blaze_box[mask_box]))  # Now normalize to the peak value

    _ivar = ivar_box*mask_box*np.logical_not(bad_box)
    _sig = np.sqrt(utils.inverse(_ivar))
    return (
        wave_box,
        flux_box*mask_box,
        _ivar,
        _sig,
        None if nivar_box is None else nivar_box*mask_box*np.logical_not(bad_box),
        mask_box*np.logical_not(bad_box),
        fwhm_box,
        blaze_box,
        sky_box,
        base_box,
        pixtot-pixmsk
    )


def extract_hist_spectrum(waveimg, frame, gpm=None, bins=1000):
    """
    Generate a quick spectrum using the nearest grid point (histogram) algorithm.

    Args:
        waveimg (`numpy.ndarray`_):
            A 2D image of the wavelength at each pixel.
        frame (`numpy.ndarray`_):
            The frame to use to extract a spectrum. Shape should be the same as waveimg
        gpm (`numpy.ndarray`_, optional):
            A boolean array indicating the pixels to include in the histogram (True = include)
        bins (`numpy.ndarray`_, int, optional):
            Either a 1D array indicating the bin edges to be used for the histogram,
            or an integer that specifies the number of bin edges to generate

    Returns:
        A tuple containing the wavelength and spectrum at the centre of each histogram bin. Both
        arrays returned in the tuple are `numpy.ndarray`_.
    """
    # Check the inputs
    if waveimg.shape != frame.shape:
        raise PypeItError("Wavelength image is not the same shape as the input frame")
    # Check the GPM
    _gpm = gpm if gpm is not None else waveimg > 0
    if waveimg.shape != _gpm.shape:
        raise PypeItError("Wavelength image is not the same shape as the GPM")
    # Set the bins
    if isinstance(bins, int):
        _bins = np.linspace(np.min(waveimg[_gpm]), np.max(waveimg[_gpm]), bins)
    elif isinstance(bins, np.ndarray):
        _bins = bins
    else:
        raise PypeItError("Argument 'bins' should be an integer or a numpy array")

    # Construct a histogram and the normalisation
    hist, edge = np.histogram(waveimg[gpm], bins=_bins, weights=frame[gpm])
    cntr, edge = np.histogram(waveimg[gpm], bins=_bins)
    # Normalise
    cntr = cntr.astype(float)
    spec = hist * utils.inverse(cntr)
    # Generate the corresponding wavelength array - set it to be the bin centre
    wave = 0.5 * (_bins[1:] + _bins[:-1])
    return wave, spec
