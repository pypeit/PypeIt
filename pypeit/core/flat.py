"""
Core module for methods related to flat fielding.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import inspect
import copy
import os

import numpy as np
import scipy.interpolate
import scipy.ndimage
import matplotlib.pyplot as plt

from IPython import embed

from pypeit import msgs
from pypeit.core import parse
from pypeit.core import pixels
from pypeit.core import tracewave
from pypeit.core import coadd
from pypeit import utils
from pypeit.core import pydl

# TODO: Put this in utils
def linear_interpolate(x1, y1, x2, y2, x):
    r"""
    Interplate or extrapolate between two points.

    Given a line defined two points, :math:`(x_1,y_1)` and
    :math:`(x_2,y_2)`, return the :math:`y` value of a new point on
    the line at coordinate :math:`x`.

    This function is meant for speed. No type checking is performed and
    the only check is that the two provided ordinate coordinates are not
    numerically identical. By definition, the function will extrapolate
    without any warning.

    Args:
        x1 (:obj:`float`):
            First abscissa position
        y1 (:obj:`float`):
            First ordinate position
        x2 (:obj:`float`):
            Second abscissa position
        y3 (:obj:`float`):
            Second ordinate position
        x (:obj:`float`):
            Abcissa for new value

    Returns:
        :obj:`float`: Interpolated/extrapolated value of ordinate at
        :math:`x`.
    """
    return y1 if np.isclose(x1,x2) else y1 + (x-x1)*(y2-y1)/(x2-x1)


# TODO: Make this function more general and put it in utils.
def sorted_flat_data(data, coo, gpm=None):
    """
    Sort a set of data by the provided coordinates.

    Args:
        data (`numpy.ndarray`_):
            Data array with arbirary shape and data type.
        coo (`numpy.ndarray`_):
            Relevant coordinate array. Shape must match ``data``.
        gpm (`numpy.ndarray`_, optional):
            Good-pixel mask for array. Used to select data (where
            ``gpm`` is True) to sort and return. Shape must match
            ``data``.  If None, all data is used.

    Returns:
        tuple: Four `numpy.ndarray`_ objects are returned:

            - A boolean array with the pixels used in the sorting.
              Shape is identical to ``data``. If ``gpm`` is provided,
              this is identicall (i.e., not a copy) of the input
              array; otherwise, it is ``np.ones(data.shape,
              dtype=bool)``.
            - A vector with the length of ``numpy.sum(gpm)`` with the
              indices that sorts the flattened list of good
              coordinate values.
            - A vector with the sorted coordinate data.
            - A vector with the data sorted by the respective
              coordinates.

        To reconstruct the input data array for the good pixels::

            _data = np.zeros(data.shape, dtype=data.dtype)
            _data[gpm] = srt_data[np.argsort(srt)]

        where ``data`` is the input array, ``gpm`` is the first
        returned object, ``srt`` is the second returned object, and
        ``srt_data`` is the last returned object of this method.

    """
    if gpm is None:
        gpm = np.ones(data.shape, dtype=bool)

    # Sort the pixels by their spatial coordinate. NOTE: By default
    # np.argsort sorts the data over the last axis. To avoid coo[gpm]
    # returning an array (which will happen if the gpm is not provided
    # as an argument), all the arrays are explicitly flattened.
    srt = np.argsort(coo[gpm].ravel())
    coo_data = coo[gpm].ravel()[srt]
    flat_data = data[gpm].ravel()[srt]
    return gpm, srt, coo_data, flat_data


def illum_filter(spat_flat_data_raw, med_width):
    """
    Filter the flat data to produce the empirical illumination
    profile.

    This is primarily a convenience method for
    :func:`construct_illum_profile`. The method first median filters
    with a window set by ``med_width`` and then Gaussian-filters the
    result with a kernel sigma set to be the maximum of 0.5 or
    ``med_width``/20.

    Args:
        spat_flat_data_raw (`numpy.ndarray`_);
            Raw flat data collapsed along the spectral direction.
        med_width (:obj:`int`):
            Width of the median filter window.

    Returns:
        `numpy.ndarray`_: Returns the filtered spatial profile of the
        flat data.
    """
    # Median filter the data
    spat_flat_data = utils.fast_running_median(spat_flat_data_raw, med_width)
    # Gaussian filter the data with a kernel that is 1/20th of the
    # median-filter width (or at least 0.5 pixels where here a "pixel"
    # is just the index of the data to fit)
    return scipy.ndimage.gaussian_filter1d(spat_flat_data, np.fmax(med_width/20.0, 0.5),
                                             mode='nearest')


def construct_illum_profile(norm_spec, spat_coo, slitwidth, spat_gpm=None, spat_samp=5,
                            illum_iter=0, illum_rej=None, debug=False):
    """
    Construct the slit illumination profile.

    Provided an image with the spectral response normalized out, this
    iteratively filters and rejects the flat-field data to construct
    the empirical slit illumination profile. The data are collapsed
    spectrally using the provided coordinate array to construct a 1D
    profile. Nominally, the provided spatial coordinates and
    good-pixel mask should be for a single slit.

    The iterations involve constructing the illumination profile
    using :func:`illum_filter` and then rejecting deviant residuals.
    Each rejection iteration recomputes the standard deviation and
    pixels to reject from the full input set (i.e., rejected pixels
    are not kept between iterations). Rejection iterations are only
    performed if ``illum_iter > 0 and illum_rej is not None``.

    Args:
        norm_spec (`numpy.ndarray`_):
            Flat-field image (2D array) with the spectral response
            normalized out.
        spat_coo (`numpy.ndarray`_):
            An image with the slit spatial coordinates, expected to
            be with respect to a single slit and span the full image
            region selected by the good-pixel mask (``spat_gpm``).
            Shape must match ``norm_spec``.
        slitwidth (:obj:`float`):
            Fiducial slit width used to set the median-filter window
            size.
        spat_gpm (`numpy.ndarray`_, optional):
            The good-pixel mask that selects the pixels to include in
            the slit illumination profile calculation. If None, **all
            pixels in the provided images are used**. For virtually
            all practical purposes, this array should be provided.
        spat_samp (:obj:`int`, :obj:`float`, optional):
            Spatial sampling for slit illumination function. This is
            the width of the median filter in detector pixels used to
            determine the slit illumination function, and thus sets
            the minimum scale on which the illumination function will
            have features.
        illum_iter (:obj:`int`, optional):
            Iteratively construct the slit illumination profile and
            reject outliers. To include rejection iterations, this
            must be larger than 0, and you have to provide the sigma
            threshold (``illum_rej``); otherwise, no iterations are
            performed.
        illum_rej (:obj:`float`, optional):
            Sigma rejection threshold for iterations. If None, no
            rejection iterations will be performed, regardless of the
            value of ``illum_iter``.
        debug (:obj:`bool`, optional):
            Construct plots output to the screen that show the result
            of each iteration. Regardless of this flag, no plots are
            shown if there are no iterations.

    Returns:
        tuple: Five `numpy.ndarray`_ objects are returned:

            - A boolean array with the pixels used in the
              construction of the illumination profile. Shape is
              identical to ``norm_spec``.
            - A vector with the length of the number of good pixels
              (sum of the first returned object) with the indices
              that sorts the flattened list of good coordinate
              values.
            - A vector with the sorted coordinate data.
            - A vector with the data sorted by the respective
              coordinates.
            - A vector with the slit illumination profile.

        To construct the empirical 2D illumination profile::

            illum = np.zeros(norm_spec.shape, dtype=float)
            illum[_spat_gpm] = profile[np.argsort(srt)]

        where ``norm_spec`` is the input array, ``_spat_gpm`` is the
        first returned object, ``srt`` is the second returned object,
        and ``profile`` is the last returned object.

    """
    if illum_rej is None and illum_iter > 0:
        msgs.warn('Cannot use iterative rejection to construct the illumination function if the '
                  'rejection is not provided.  Continuing without iteration.')

    _spat_gpm = np.ones(norm_spec.shape, dtype=bool) if spat_gpm is None else np.copy(spat_gpm)
    _spat_gpm, spat_srt, spat_coo_data, spat_flat_data_raw \
            = sorted_flat_data(norm_spec, spat_coo, gpm=_spat_gpm)
    spat_gpm_data_raw = np.ones(spat_flat_data_raw.size, dtype=bool)

    # Assume the density of samples in any given spatial coordinate is
    # roughly the same at all spatial positions. Calculate the fraction
    # of the slit width for the median filter as set by the
    # ``spat_samp`` parameter.
    med_width = int(np.ceil(np.sum(spat_gpm) * spat_samp / slitwidth))

    # Construct the filtered illumination profile (iteratively if requested)
    for i in range(illum_iter+1):
        spat_flat_data = illum_filter(spat_flat_data_raw[spat_gpm_data_raw], med_width)

        if illum_iter == 0 or illum_rej is None:
            # No iterations so we're done (skips debug plot)
            return spat_gpm, spat_srt, spat_coo_data, spat_flat_data_raw, spat_flat_data

        if i == illum_iter:
            # Don't perform the rejection on the last iteration
            break

        # Iteration does not keep previous rejections. NOTE: Rejections
        # at either end of the data array would cause the interpolation
        # below to fault, which is why I set bound_error to False. This
        # may be a problem though because I set the fill value to 0...
        interp = scipy.interpolate.interp1d(spat_coo_data[spat_gpm_data_raw], spat_flat_data,
                                      bounds_error=False, fill_value=0.0, assume_sorted=True)
        resid = spat_flat_data_raw - interp(spat_coo_data)
        sigma = np.std(resid)
        spat_gpm_data_raw = np.absolute(resid) < illum_rej*sigma

    # TODO: Provide a report?

    if debug:
        plt.clf()
        ax = plt.gca()
        ax.scatter(spat_coo_data[spat_gpm_data_raw], spat_flat_data_raw[spat_gpm_data_raw],
                   marker='.', lw=0, s=10, color='k', zorder=1, label='used data')
        ax.scatter(spat_coo_data[np.invert(spat_gpm_data_raw)],
                   spat_flat_data_raw[np.invert(spat_gpm_data_raw)],
                   marker='.', lw=0, s=10, color='C3', zorder=2, label='rejected data')
        ax.plot(spat_coo_data[spat_gpm_data_raw], spat_flat_data, color='C2', zorder=3,
                label='filtered profile')
        ax.legend()
        ax.set_title('Optimized slit illumination profile')
        ax.set_xlabel('Spatial coordinate')
        ax.set_ylabel('Spectrally collapsed, normalized flux')
        plt.show()

    # Include the rejected data in the full image good-pixel mask
    _spat_gpm[_spat_gpm] = spat_gpm_data_raw[np.argsort(spat_srt)]
    # Recreate the illumination profile data
    _spat_gpm, spat_srt, spat_coo_data, spat_flat_data_raw \
            = sorted_flat_data(norm_spec, spat_coo, gpm=_spat_gpm)
    return _spat_gpm, spat_srt, spat_coo_data, spat_flat_data_raw, \
                illum_filter(spat_flat_data_raw, med_width)


def illum_profile_spectral_poly(rawimg, waveimg, slitmask, slitmask_trim, model, slit_illum_ref_idx=0,
                                gpmask=None, thismask=None, nbins=20, debug=False):
    """
    Use a polynomial fit to control points along the spectral direction to determine the relative
    spectral illumination of all slits. Currently, this routine is only used for image slicer IFUs.

    Parameters
    ----------
    rawimg : `numpy.ndarray`_
        Image data that will be used to estimate the spectral relative sensitivity
    waveimg : `numpy.ndarray`_
        Wavelength image
    slitmask : `numpy.ndarray`_
        A 2D int mask, the same shape as rawimg, indicating which pixels are on a slit. A -1 value
        indicates not on a slit, while any pixels on a slit should have the value of the slit spatial ID
        number.
    slitmask_trim : `numpy.ndarray`_
        Same as slitmask, but the slit edges are trimmed.
    model : `numpy.ndarray`_
        A model of the rawimg data.
    slit_illum_ref_idx : :obj:`int`
        Index of slit that is used as the reference.
    gpmask : `numpy.ndarray`_, optional
        Boolean good pixel mask (True = Good)
    thismask : `numpy.ndarray`_, optional
        A boolean mask (True = good) that indicates all pixels where the scaleImg should be constructed.
        If None, the slitmask that is generated by this routine will be used.
    nbins : :obj:`int`
        Number of bins in the spectral direction to sample the relative spectral sensitivity
    debug : :obj:`bool`
        If True, some plots will be output to test if the fitting is working correctly.

    Returns
    -------
    scale_model: `numpy.ndarray`_
        An image containing the appropriate scaling
    """
    msgs.info(f"Performing relative spectral sensitivity correction (reference slit = {slit_illum_ref_idx})")
    # Generate the mask
    _thismask = thismask if (thismask is not None) else (slitmask > 0)
    gpm = gpmask if (gpmask is not None) else np.ones_like(rawimg, dtype=bool)
    # Extract the list of  spatial IDs from the slitmask
    slitmask_spatid = np.unique(slitmask)
    slitmask_spatid = np.sort(slitmask_spatid[slitmask_spatid > 0])
    # Initialise the scale image that will be returned
    scaleImg = np.ones_like(rawimg)
    # Divide the slit into several bins and calculate the median of each bin
    for sl, spatid in enumerate(slitmask_spatid):
        # Prepare the masks, edges, and fitting variables
        this_slit = (slitmask == spatid)
        this_slit_trim = (slitmask_trim == spatid)
        this_slit_mask = gpm & this_slit_trim
        this_wave = waveimg[this_slit_mask]
        wavedg = np.linspace(np.min(this_wave), np.max(this_wave), nbins + 1)
        wavcen = 0.5 * (wavedg[1:] + wavedg[:-1])
        scale_all = rawimg[this_slit_mask] * utils.inverse(model[this_slit_mask])
        scale_bin = np.zeros(nbins)
        scale_err = np.zeros(nbins)
        for bb in range(nbins):
            cond = (this_wave >= wavedg[bb]) & (this_wave <= wavedg[bb + 1])
            scale_bin[bb] = np.median(scale_all[cond])
            scale_err[bb] = 1.4826 * np.median(np.abs(scale_all[cond] - scale_bin[bb]))
        wgd = np.where(scale_err > 0)
        coeff = np.polyfit(wavcen[wgd], scale_bin[wgd], w=1/scale_err[wgd], deg=2)
        scaleImg[this_slit] *= np.polyval(coeff, waveimg[this_slit])
        if debug:
            mod = np.polyval(coeff, wavcen[wgd])
            plt.errorbar(wavcen[wgd], scale_bin[wgd], yerr=scale_err[wgd], fmt='o')
            plt.plot(wavcen[wgd], mod, 'r-')
            plt.show()
        if sl == slit_illum_ref_idx:
            scaleImg[_thismask] *= utils.inverse(np.polyval(coeff, waveimg[_thismask]))
    minv, maxv = np.min(scaleImg[_thismask]), np.max(scaleImg[_thismask])
    msgs.info("Minimum/Maximum scales = {0:.5f}, {1:.5f}".format(minv, maxv))
    return scaleImg


def smooth_scale(arr, wave_ref=None, polydeg=None, sn_smooth_npix=None):
    """
    Smooth the relative sensitivity array using a polynomial fit or a boxcar filter.

    Parameters
    ----------
    arr : `numpy.ndarray`_
        Array containing the relative sensitivity
    wave_ref : `numpy.ndarray`_, optional
        Wavelength array corresponding to the relative sensitivity array. Only used if polydeg is not None.
    polydeg : :obj:`int`, optional
        Degree of the polynomial fit to the relative sensitivity array. If None, a boxcar filter will be used.
    sn_smooth_npix : :obj:`int`, optional
        Number of pixels to use for the boxcar filter. Only used if polydeg is None.

    Returns
    -------
    arr_smooth : `numpy.ndarray`
        Smoothed relative sensitivity array
    """
    # Do some checks on the input
    if polydeg is not None and wave_ref is None:
        msgs.error("Must provide a wavelength array if polydeg is not None")
    if polydeg is None and sn_smooth_npix is None:
        msgs.error("Must provide either polydeg or sn_smooth_npix")
    # Smooth the relative sensitivity array
    if polydeg is not None:
        gd = (arr != 0)
        wave_norm = (wave_ref - wave_ref[0]) / (wave_ref[1] - wave_ref[0])
        coeff = np.polyfit(wave_norm[gd], arr[gd], polydeg)
        ref_relscale = np.polyval(coeff, wave_norm)
    else:
        ref_relscale = coadd.smooth_weights(arr, (arr != 0), sn_smooth_npix)
    # Return the smoothed relative sensitivity array
    return ref_relscale


# TODO:: See pypeit/deprecated/flat.py for a spline version. The following polynomial version is faster, but
#        the spline version is more versatile.
def poly_map(rawimg, rawivar, waveimg, slitmask, slitmask_trim, modelimg, deg=3,
             slit_illum_ref_idx=0, gpmask=None, thismask=None, debug=False):
    """
    Use a polynomial fit to control points along the spectral direction to construct a map between modelimg and
    rawimg. Currently, this routine is only used for image slicer IFUs.

    This problem needs to be recast into a chi-squared problem, that can take advantage of polynomial fitting.
    Refer to the following for variable assignments:
    https://numpy.org/doc/stable/reference/generated/numpy.polyfit.html
    where:
    y = science/model
    w = model/science_error
    resid = w*(y-f(x))
    Then, iterate over this. The reason to iterate is that there should be a mapping between science and
    science_error. Once we have an estimate of f(x), we can do the following:
    spl = spline(science, science_error)
    new_science_error = spl(model*f(x))
    Now iterate a few times until new_science_error (and the fit) is converged.

    Parameters
    ----------
    rawimg : `numpy.ndarray`_
        Image data that will be used to estimate the spectral relative sensitivity
    rawivar : `numpy.ndarray`_
        Inverse variance image of rawimg
    waveimg : `numpy.ndarray`_
        Wavelength image
    slitmask : `numpy.ndarray`_
        A 2D int mask, the same shape as rawimg, indicating which pixels are on a slit. A -1 value
        indicates not on a slit, while any pixels on a slit should have the value of the slit spatial ID
        number.
    slitmask_trim : `numpy.ndarray`_
        Same as slitmask, but the slit edges are trimmed.
    model : `numpy.ndarray`_
        A model of the rawimg data.
    slit_illum_ref_idx : :obj:`int`
        Index of slit that is used as the reference.
    gpmask : `numpy.ndarray`_, optional
        Boolean good pixel mask (True = Good)
    thismask : `numpy.ndarray`_, optional
        A boolean mask (True = good) that indicates all pixels where the scaleImg should be constructed.
        If None, the slitmask that is generated by this routine will be used.
    debug : :obj:`bool`
        If True, some plots will be output to test if the fitting is working correctly.

    Returns
    -------
    modelmap : `numpy.ndarray`_
        A 2D image with the same shape as rawimg, that contains the modelimg mapped to the rawimg
    relscale : `numpy.ndarray`_
        A 2D image with the same shape as rawimg, that contains the relative spectral sensitivity
    """
    # Some variables to consider putting as function arguments
    numiter = 4

    # Start by calculating a ratio of the raming and the modelimg
    nspec = rawimg.shape[0]
    _ratio = rawimg * utils.inverse(modelimg)
    _ratio_ivar = rawivar * modelimg**2
    _fit_wghts = modelimg * np.sqrt(rawivar)

    # Generate the mask
    _thismask = thismask if (thismask is not None) else (slitmask > 0)
    gpm = gpmask if (gpmask is not None) else np.ones_like(rawimg, dtype=bool)
    # Extract the list of  spatial IDs from the slitmask
    slitmask_spatid = np.unique(slitmask)
    slitmask_spatid = np.sort(slitmask_spatid[slitmask_spatid > 0])

    # Create a spline between the raw data and the error
    flxsrt = np.argsort(np.ravel(rawimg))
    spl = scipy.interpolate.interp1d(np.ravel(rawimg)[flxsrt], np.ravel(rawivar)[flxsrt], kind='linear',
                                     bounds_error=False, fill_value=0.0, assume_sorted=True)
    modelmap = np.ones_like(rawimg)
    relscale = np.ones_like(rawimg)
    msgs.info("Generating a polynomial map between the model and the raw data")
    for sl, spatid in enumerate(slitmask_spatid):
        # Prepare the masks, edges, and fitting variables
        this_slit = (slitmask == spatid)
        this_slit_trim = (slitmask_trim == spatid)
        this_slit_mask = gpm & this_slit_trim
        this_wave = waveimg[this_slit_mask]
        wmin, wmax = np.min(this_wave), np.max(this_wave)
        this_wghts = _fit_wghts[this_slit_mask]
        for ii in range(numiter):
            # Generate the map between model and data
            coeff = np.polyfit((this_wave-wmin)/(wmax-wmin), _ratio[this_slit_mask], deg, w=this_wghts)
            # Construct the mapping, and use this to make a model of the rawdata
            this_modmap = np.polyval(coeff, (this_wave-wmin)/(wmax-wmin))
            this_modflx = modelimg[this_slit_mask] * this_modmap
            # Update the fit weights
            this_wghts = modelimg[this_slit_mask] * np.sqrt(spl(this_modflx))
        # Produce the final model for this slit
        modelmap[this_slit] *= np.polyval(coeff, (waveimg[this_slit]-wmin)/(wmax-wmin))
        relscale[this_slit] *= np.polyval(coeff, (waveimg[this_slit]-wmin)/(wmax-wmin))
        # Check if this is the reference slit, and if so, set the scale relative to the reference slit.
        if sl == slit_illum_ref_idx:
            for slidx, spatid in enumerate(slitmask_spatid):
                # Get the slit pixels
                norm_slit = (slitmask == spatid)
                relscale[norm_slit] /= np.polyval(coeff, (waveimg[norm_slit]-wmin)/(wmax-wmin))
    # Return the modelmap and the relative scale
    return modelmap, relscale


# TODO: See pypeit/deprecated/flat.py for the previous version. We need
# to continue to vet this algorithm to make sure there are no
# unforeseen corner cases that cause errors.
def tweak_slit_edges(left, right, spat_coo, norm_flat, thresh=0.93, maxfrac=0.1, debug=False):
    r"""
    Adjust slit edges based on the normalized slit illumination profile.

    Args:
        left (`numpy.ndarray`_):
            Array with the left slit edge for a single slit. Shape is
            :math:`(N_{\rm spec},)`.
        right (`numpy.ndarray`_):
            Array with the right slit edge for a single slit. Shape
            is :math:`(N_{\rm spec},)`.
        spat_coo (`numpy.ndarray`_):
            Spatial pixel coordinates in fractions of the slit width
            at each spectral row for the provided normalized flat
            data. Coordinates are relative to the left edge (with the
            left edge at 0.). Shape is :math:`(N_{\rm flat},)`.
            Function assumes the coordinate array is sorted.
        norm_flat (`numpy.ndarray`_)
            Normalized flat data that provide the slit illumination
            profile. Shape is :math:`(N_{\rm flat},)`.
        thresh (:obj:`float`, optional):
            Threshold of the normalized flat profile at which to
            place the two slit edges.
        maxfrac (:obj:`float`, optional):
            The maximum fraction of the slit width that the slit edge
            can be adjusted by this algorithm. If ``maxfrac = 0.1``,
            this means the maximum change in the slit width (either
            narrowing or broadening) is 20% (i.e., 10% for either
            edge).
        debug (:obj:`bool`, optional):
            Show flow interrupting plots that show illumination
            profile in the case of a failure and the placement of the
            tweaked edge for each side of the slit regardless.

    Returns:
        tuple: Returns six objects:

            - The threshold used to set the left edge
            - The fraction of the slit that the left edge is shifted to
              the right
            - The adjusted left edge
            - The threshold used to set the right edge
            - The fraction of the slit that the right edge is shifted to
              the left
            - The adjusted right edge

    """
    # Check input
    nspec = len(left)
    if len(right) != nspec:
        msgs.error('Input left and right traces must have the same length!')

    # Median slit width
    slitwidth = np.median(right - left)

    # Setup the masked array for finding the continuous left and right
    # regions
    masked_flat = np.ma.MaskedArray(norm_flat)

    # ------------------------------------------------------------------
    # Adjust the left edge

    # Get the maximum to the left of the center
    # TODO: Set a parameter for this
    ileft = (spat_coo > 0.1) & (spat_coo < 0.4)
    if not np.any(ileft):
        msgs.error('No coordinates toward the left of the slit center.  Slit boundaries are '
                   'likely in error, and you probably have a bad (very short) slit.  Slit center '
                   'at center row is {0:.1f}.'.format((left[nspec//2] + right[nspec//2])/2))
    left_thresh = thresh * np.amax(norm_flat[ileft])

    # Find the data that are less than the provided threshold and
    # within the limits set by the offset
    masked_flat[(spat_coo >= maxfrac) | (norm_flat >= left_thresh)] = np.ma.masked

    # To tweak, there must be at least one pixel that meet the above
    # criteria
    left_shift = 0.
    new_left = np.copy(left)
    if not np.all(masked_flat.mask):
        # Find the last index of the first contiguous region
        contiguous_region = np.ma.flatnotmasked_contiguous(masked_flat)[0]
        if contiguous_region.stop is None:
            if debug:
                plt.scatter(spat_coo[masked_flat.mask], norm_flat[masked_flat.mask], marker='.',
                            s=10, color='C3', lw=0)
                plt.scatter(spat_coo[np.invert(masked_flat.mask)],
                            norm_flat[np.invert(masked_flat.mask)], marker='.', s=10, color='k',
                            lw=0)
                plt.show()
            msgs.error('Tweak left edge has failed!  Bad continuous region.')
        i = contiguous_region.stop-1
        if i >= 0 and norm_flat[i-1] > norm_flat[i]:
            msgs.warn('When adjusting left edge, found noisy illumination profile structure.')
        if debug:
            plt.scatter(spat_coo[masked_flat.mask], norm_flat[masked_flat.mask], marker='.', s=10,
                        color='C3', lw=0)
            plt.scatter(spat_coo[np.invert(masked_flat.mask)],
                        norm_flat[np.invert(masked_flat.mask)], marker='.', s=10, color='k', lw=0)
            plt.scatter(spat_coo[i], norm_flat[i], marker='o', facecolor='none', s=50, color='C1')
            plt.show()
        if norm_flat[i+1] < left_thresh:
            msgs.warn('Left slit boundary tweak limited by maximum allowed shift: {:.1f}%'.format(
                        100*maxfrac))
            left_shift = maxfrac
        else:
            left_shift = linear_interpolate(norm_flat[i], spat_coo[i], norm_flat[i+1],
                                           spat_coo[i+1], left_thresh)
        msgs.info('Tweaking left slit boundary by {0:.1f}%'.format(100*left_shift) +
                  ' % ({0:.2f} pixels)'.format(left_shift*slitwidth))
        new_left += left_shift * slitwidth

    # ------------------------------------------------------------------
    # Adjust the right edge

    # Get the maximum to the right of the center
    # TODO: Set a parameter for this
    iright = (spat_coo > 0.6) & (spat_coo < 0.9)
    if not np.any(iright):
        msgs.error('No coordinates toward the right of the slit center.  Slit boundaries are '
                   'likely in error, and you probably have a bad (very short) slit.  Slit center '
                   'at center row is {0:.1f}.'.format((left[nspec//2] + right[nspec//2])/2))
    right_thresh = thresh * np.amax(norm_flat[iright])

    # Find the data that are less than the provided threshold and
    # within the limits set by the offset
    masked_flat.mask = np.ma.nomask
    masked_flat[(spat_coo <= 1 - maxfrac) | (norm_flat >= right_thresh)] = np.ma.masked

    # To tweak, there must be at least one pixel that meets the above
    # criteria
    right_shift = 0.
    new_right = np.copy(right)
    if not np.all(masked_flat.mask):
        # Find the first index of the last contiguous region
        contiguous_region = np.ma.flatnotmasked_contiguous(masked_flat)[-1]
        if contiguous_region.start is None:
            if debug:
                plt.scatter(spat_coo[masked_flat.mask], norm_flat[masked_flat.mask], marker='.',
                            s=10, color='C3', lw=0)
                plt.scatter(spat_coo[np.invert(masked_flat.mask)],
                            norm_flat[np.invert(masked_flat.mask)], marker='.', s=10, color='k',
                            lw=0)
                plt.show()
            msgs.error('Tweak right edge has failed!  Bad continuous region.')
        i = contiguous_region.start
        if i < norm_flat.size-1 and norm_flat[i+1] > norm_flat[i]:
            msgs.warn('When adjusting right edge, found noisy illumination profile structure.')
        if debug:
            plt.scatter(spat_coo[masked_flat.mask], norm_flat[masked_flat.mask], marker='.', s=10,
                        color='C3', lw=0)
            plt.scatter(spat_coo[np.invert(masked_flat.mask)],
                        norm_flat[np.invert(masked_flat.mask)], marker='.', s=10, color='k', lw=0)
            plt.scatter(spat_coo[i], norm_flat[i], marker='o', facecolor='none', s=50, color='C1')
            plt.show()
        if norm_flat[i-1] < right_thresh:
            msgs.warn('Right slit boundary tweak limited by maximum allowed shift: {:.1f}%'.format(
                        100*maxfrac))
            right_shift = maxfrac
        else:
            right_shift = 1-linear_interpolate(norm_flat[i-1], spat_coo[i-1], norm_flat[i],
                                               spat_coo[i], right_thresh)
        msgs.info('Tweaking right slit boundary by {0:.1f}%'.format(100*right_shift) +
                  ' % ({0:.2f} pixels)'.format(right_shift*slitwidth))
        new_right -= right_shift * slitwidth

    return left_thresh, left_shift, new_left, right_thresh, right_shift, new_right

#def flatfield(sciframe, flatframe, bpm=None, illum_flat=None, snframe=None, varframe=None):
def flatfield(sciframe, flatframe, varframe=None):
    r"""
    Field flatten the input image.

    This is a simple scaling of the provided science frame by the inverse of the
    flat frame.

    Args:
        sciframe (`numpy.ndarray`_):
            The science frame to flat-field correct
        flatframe (`numpy.ndarray`_):
            The flat-field image to use for the correction.  Shape must match
            ``sciframe``.
        varframe (`numpy.ndarray`_, optional):
            The variance in the science frame (``sciframe``).  If provided, the
            flat-fielding operation is propagated to the variance, and the
            result is returned.  Shape must match ``sciframe``.

    Returns:
        :obj:`tuple`: A tuple of two or three `numpy.ndarray`_ objects.  The
        first two are the rescaled science frame and a boolean bad-pixel mask
        indicating where the flat frame was not positive.  If a variance image
        is provided, the third object is the propagated variance in the rescaled
        science frame.
    """
    if flatframe.shape != sciframe.shape:
        msgs.error('Shape of flat frame does not match science frame.')
    if varframe is not None and varframe.shape != sciframe.shape:
        msgs.error('Shape of variance frame does not match science frame.')

    # New image
    retframe = np.zeros_like(sciframe)
    gpm = (flatframe > 0.0) & np.isfinite(flatframe)
    retframe[gpm] = sciframe[gpm]/flatframe[gpm]
    if varframe is None:
        return retframe, np.logical_not(gpm)

    # Propagate the variance
    retvar = np.zeros_like(sciframe)
    retvar[gpm] = varframe[gpm]/flatframe[gpm]**2
    return retframe, np.logical_not(gpm), retvar


