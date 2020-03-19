"""
Core module for methods related to flat fielding.

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""
import inspect
import copy
import os

import numpy as np
from scipy import interpolate, ndimage
from matplotlib import pyplot as plt

from IPython import embed

from pypeit import msgs
from pypeit.core import parse
from pypeit.core import pixels
from pypeit.core import tracewave
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
    return ndimage.filters.gaussian_filter1d(spat_flat_data, np.fmax(med_width/20.0, 0.5),
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
    using :func:`illum_profile` and then rejecting deviant residuals.
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
        interp = interpolate.interp1d(spat_coo_data[spat_gpm_data_raw], spat_flat_data,
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

# TODO: See pypeit/deprecated/flat.py for the previous version. We need
# to continue to vet this algorithm to make sure there are no
# unforeseen corner cases that cause errors.
def tweak_slit_edges(left, right, spat_coo, norm_flat, thresh=0.93, maxfrac=0.1):
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
            msgs.error('Tweak left edge has failed!  Bad continuous region.')
        i = contiguous_region.stop-1
        if i >= 0 and norm_flat[i-1] > norm_flat[i]:
            msgs.error('Tweak left edge has failed!  Unexpected illumination profile structure.')
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
            msgs.error('Tweak right edge has failed!  Bad continuous region.')
        i = contiguous_region.start
        if i < norm_flat.size-1 and norm_flat[i+1] > norm_flat[i]:
            msgs.error('Tweak left edge has failed!  Unexpected illumination profile structure.')
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


# Deprecated.  Replaced by tweak_slit_edges above.
#def tweak_slit_edges_brute_force(slit_left_in, slit_righ_in, ximg_fit, normimg, tweak_slits_thresh,
#                                 tweak_slits_maxfrac):
#    """
#    DOC THIS!
#    """
#
#
#    tweak_left = False
#    tweak_righ = False
#    slit_left_out = np.copy(slit_left_in)
#    slit_righ_out = np.copy(slit_righ_in)
#    # How many pixels wide is the slit at each Y?
#    slitwidth = np.median(slit_righ_in - slit_left_in)
#    # Determine the maximum at the left and right end of the slit
#    ileft = (ximg_fit > 0.1) & (ximg_fit < 0.4)
#    irigh = (ximg_fit > 0.6) & (ximg_fit < 0.9)
##    if (not np.any(ileft)) or (not np.any(irigh)):
##        msgs.error('Cannot tweak slits because much of the slit is masked. You probably have a bad slit')
##        tweak_dict = {'xleft': 0.0, 'xrigh': 0.0,
##                      'norm_max_left': 0.0, 'norm_max_righ': 0.0,
##                      'tweak_left': tweak_left, 'tweak_righ': tweak_righ}
##        return slit_left_out, slit_righ_out, tweak_dict
#
#    #xleft = ximg_fit[ileft]
#    #xrigh = ximg_fit[irigh]
#    norm_max_left = normimg[ileft].max()
#    norm_max_righ = normimg[irigh].max()
#
#    msgs.info('Tweaking slit boundaries using slit illumination function')
#    step = 0.001
#    # march out from middle to find left edge
#    msgs.info('Left threshold = {:5.3f}'.format(tweak_slits_thresh * norm_max_left) +
#              ' --  or {:5.3f}'.format(
#                  100.0 * tweak_slits_thresh) + ' % of left side max of illumination function = {:5.3f}'.format(norm_max_left))
#    for xleft in np.arange(0.5, ximg_fit.min(), -step):
#        norm_now = np.interp(xleft, ximg_fit, normimg)
#        if (norm_now < tweak_slits_thresh * norm_max_left) & (xleft < tweak_slits_maxfrac):
#            slit_left_out += xleft * slitwidth
#            tweak_left = True
#            msgs.info('Tweaking left slit boundary by {:5.3f}'.format(100 * xleft) +
#                      ' %, or {:7.3f}'.format(xleft * slitwidth) + ' pixels')
#            if np.abs(xleft - tweak_slits_maxfrac) < 0.01:
#                msgs.warn(
#                    'Left slit boundary tweak limited by maximum changed allowed by tweak_slits_maxfracn={:5.3f}'.format(
#                        100.0 * tweak_slits_maxfrac) + ' %')
#            break
#    msgs.info('Right threshold = {:5.3f}'.format(tweak_slits_thresh * norm_max_righ) +
#              ' --  or {:5.3f}'.format(
#                  100.0 * tweak_slits_thresh) + ' % of right side max of illumination function = {:5.3f}'.format(norm_max_righ))
#    # march out from middle  to find right edge
#    for xrigh in np.arange(0.5, ximg_fit.max(), step):
#        norm_now = np.interp(xrigh, ximg_fit, normimg)
#        if (norm_now < tweak_slits_thresh * norm_max_righ) & ((1.0 - xrigh) < tweak_slits_maxfrac):
#            slit_righ_out -= (1.0 - xrigh) * slitwidth
#            tweak_righ = True
#            msgs.info('Tweaking right slit boundary by {:5.3f}'.format(100 * (1.0 - xrigh)) +
#                      ' %, or {:7.3f}'.format((1.0 - xrigh) * slitwidth) + ' pixels')
#            if np.abs((1.0 - xrigh) - tweak_slits_maxfrac) < 0.01:
#                msgs.warn(
#                    'Right slit boundary tweak limited by maximum changed allowed by tweak_slits_maxfracn={:5.3f}'.format(
#                        100.0 * tweak_slits_maxfrac) + ' %')
#            break
#
#    tweak_dict = {'xleft': xleft, 'xrigh': xrigh,
#                  'norm_max_left': norm_max_left, 'norm_max_righ': norm_max_righ,
#                  'tweak_left': tweak_left, 'tweak_righ': tweak_righ}
#
#    return slit_left_out, slit_righ_out, tweak_dict

# Deprecated.  Replaced by FlatField.fit()
#def fit_flat(flat, tilts_dict, tslits_dict_in, slit, inmask = None,
#             spec_samp_fine = 1.2, spec_samp_coarse = 50.0, spat_samp = 5.0, npoly = None, trim_edg = (3.0,3.0), pad =5.0,
#             tweak_slits = True, tweak_slits_thresh = 0.93, tweak_slits_maxfrac = 0.10, nonlinear_counts =1e10, debug = False):
#
#
#    """ Compute pixelflat and illumination flat from a flat field image.
#
#    Parameters
#    ----------
#    flat :  float ndarray, shape (nspec, nspat)
#        Flat field image in units of electrons.
#
#
#    tilts_dict: dict
#          Dictionary containing wavelength tilts image and other information indicating how wavelengths move across the slit
#
#    tslits_dict: dict
#          Dictionary with information on the slit boundaries
#    slit: int
#          Slit currently being considered
#    inmask: boolean ndarray, shape (nspec, nspat), default inmask = None, optional
#      Input mask for pixels not to be included in sky subtraction fits. True = Good (not masked), False = Bad (masked)
#
#    spec_samp_fine: float, default = 1.2, optional
#      bspline break point spacing in units of pixels for spectral fit to flat field blaze function.
#
#    spec_samp_coarse: float, default = 50.0, optional
#      bspline break point spacing in units of pixels for 2-d bspline-polynomial fit to flat field image residuals.
#      This should be a large number unless you are trying to fit a sky flat with lots of features.
#
#    spat_samp: float, default = 5.0, optional
#      Spatial sampling for spatial slit illumination function. This is the width of the median filter in pixels used to
#      determine the slit illumination function, and thus sets the minimum scale on which the illumination function will
#      have features.
#
#    trim_edg: tuple of floats  (left_edge, right_edge), default (3,3), optional
#      indicates how many pixels to trim from left and right slit edges for creating the edgemask, which is used to mask
#      the edges from the initial (fine) spectroscopic fit to the blaze function.
#
#    pad: int, default = 5, optional
#      Padding window used to create expanded slitmask images used for re-determining slit boundaries. Tilts are also
#      computed using this expanded slitmask in cases the slit boundaries need to be moved outward.
#
#    npoly: int, default = None, optional
#      Order of polynomial for 2-d bspline-polynomial fit to flat field image residuals. The code determines the order of
#      these polynomials to each slit automatically depending on the slit width, which is why the default is None.
#      Do not attempt to set this paramter unless you know what you are doing.
#
#
#    tweak_slits: bool, default = True, optional
#      Slit edges will be tweaked such the left and right bounadaries intersect the location where the illumination
#      function falls below tweak_slits_thresh (see below) of its maximum value near the center (moving out from the center)
#
#    tweak_slits_thresh: float, default = 0.93, optional
#      If tweak_slits is True, this sets the illumination function threshold used to tweak the slits
#
#    tweak_slits_maxfrac: float, default = 0.10, optional
#      Maximum fractinoal amount (of slit width) allowed for each trimming the left and right slit boundaries, i.e. the
#      default is 10% which means slits would shrink by at most 20% (10% on each side)
#
#    debug: bool, default = False, optional
#      Show plots useful for debugging. This will block further execution of the code until the plot windows are closed.
#
#    Returns
#    -------
#    pixeflat:   ndarray with same shape as flat
#      Pixelflat gives pixel-to-pixel variations of detector response. Values are centered about unity.
#
#    illumflat:  ndarray with same shape as flat
#      Illumination flat gives variations of the slit illumination function across the spatial direction of the detect.
#      Values are centered about unity. The slit illumination function is computed by dividing out the spectral response and
#      collapsing out the spectral direction.
#
#    flat_model:  ndarray with same shape as flat
#      Full 2-d model image of the input flat image in units of electrons.  The pixelflat is defined to be flat/flat_model.
#
#    tilts: ndarray with same shape as flat
#      Tilts image fit for this slit evaluated using the new slit boundaries
#
#    thismask_out: ndarray with same shape as flat, bool
#       Boolean mask indicating which pixels are on the slit now with the new slit boundaries
#
#    slit_left_out: ndarray with shape (nspec,)
#       Tweaked left slit bounadries
#
#    slit_righ_out: ndarray with shape (nspec,)
#       Tweaked right slit bounadries
#
#    Notes
#    -----
#
#    Revision History
#        - 11-Mar-2005  First version written by Scott Burles.
#        - 2005-2018    Improved by J. F. Hennawi and J. X. Prochaska
#        - 3-Sep-2018 Ported to python by J. F. Hennawi and significantly improved
#    """
#
#    shape = flat.shape
#    nspec = shape[0]
#    nspat = shape[1]
#
#    # Get the thismask_in and input slit bounadries from the tslits_dict
#    slit_left_in = tslits_dict_in['slit_left'][:,slit]
#    slit_righ_in = tslits_dict_in['slit_righ'][:,slit]
#
#    # ... This is created using the padding in the slits dict
#    thismask_in = pixels.tslits2mask(tslits_dict_in) == slit
#
#    # Check for saturation of the flat. If there are not enough pixels do not attempt a fit
#    good_frac = np.sum(thismask_in & (flat < nonlinear_counts))/np.sum(thismask_in)
#    if good_frac < 0.5:
#        msgs.warn(msgs.newline() + 'Only {:4.2f}'.format(100*good_frac) + '% of the pixels on this slit are not saturated.' +
#                  msgs.newline() + 'Consider raising nonlinear_counts={:5.3f}'.format(nonlinear_counts) +
#                  msgs.newline() + 'Not attempting to flat field slit# {:d}'.format(slit))
#        return np.ones_like(flat), np.ones_like(flat), np.zeros_like(flat), tilts_dict['tilts'], thismask_in, slit_left_in, slit_righ_in
#
#    # Approximate number of pixels sampling each spatial pixel for this (original) slit.
#    npercol = np.fmax(np.floor(np.sum(thismask_in)/nspec),1.0)
#    # Demand at least 10 pixels per row (on average) per degree of the polynomial
#    if npoly is None:
#        npoly_in = 7
#        npoly  = np.clip(npoly_in, 1, np.ceil(npercol/10.).astype(int))
#        #npoly = np.fmax(np.fmin(npoly_in, (np.ceil(npercol/10.)).astype(int)),1)
#
#    # ... Here, thismask_in is only used to set the size of the image
#    ximg_in, edgmask_in = pixels.ximg_and_edgemask(slit_left_in, slit_righ_in, thismask_in, trim_edg=trim_edg)
#    # Create a fractional position image ximg that encompasses the whole image, rather than just the thismask_in slit pixels
#    spat_img = np.outer(np.ones(nspec), np.arange(nspat)) # spatial position everywhere along image
#    slit_left_img = np.outer(slit_left_in, np.ones(nspat))   # left slit boundary replicated spatially
#    slitwidth_img = np.outer(slit_righ_in - slit_left_in, np.ones(nspat)) # slit width replicated spatially
#    ximg = (spat_img - slit_left_img)/slitwidth_img
#
#    # Create a wider slitmask image with shift pixels padded on each side
#    # ... This is created using the padding provided to the function.
#    # This is always 5 because pad was never passed to fit_flat from
#    # FlatField.run
#    slitmask_pad = pixels.tslits2mask(tslits_dict_in, pad = pad)
#    thismask = (slitmask_pad == slit) # mask enclosing the wider slit bounadries
#    # Create a tilts image using this padded thismask, rather than using the original thismask_in slit pixels
#    tilts = tracewave.fit2tilts(shape, tilts_dict['coeffs'], tilts_dict['func2d'])
#    piximg = tilts * (nspec-1)
#    pixvec = np.arange(nspec)
#
#    # ... inmask was always passed as the bpm or unity from FlatField so
#    # this was never called
#    if inmask is None:
#        inmask = np.copy(thismask)
#
#    # Fit the spectral direction of the blaze. We do this in the log
#    log_flat = np.log(np.fmax(flat, 1.0))
#    inmask_log = ((flat > 1.0) & inmask)
#    log_ivar = inmask_log.astype(float)/0.5**2 # set errors to just be 0.5 in the log
#
#    # Flat field pixels for fitting spectral direction. Restrict to original slit pixels
#    # ... This selects good pixels in the trimmed slit
#    fit_spec = thismask_in & inmask & np.invert(edgmask_in) #& (flat < nonlinear_counts)
#    nfit_spec = np.sum(fit_spec)
#    spec_frac = nfit_spec/np.sum(thismask_in)
#    msgs.info('Spectral fit of flatfield for {:}'.format(nfit_spec) + ' pixels')
#    if spec_frac < 0.5:
#        msgs.warn('Spectral flatfield fit is to only {:4.2f}'.format(100*spec_frac) + '% of the pixels on this slit.' +
#                  msgs.newline() + '          Something appears to be wrong here')
#
#    isrt_spec = np.argsort(piximg[fit_spec])
#    pix_fit = piximg[fit_spec][isrt_spec]
#    log_flat_fit = log_flat[fit_spec][isrt_spec]
#    log_ivar_fit = log_ivar[fit_spec][isrt_spec]
#    inmask_log_fit = inmask_log[fit_spec][isrt_spec]
#    logrej = 0.5 # rejectino threshold for spectral fit in log(image)
#
#    # ToDo Figure out how to deal with the fits going crazy at the edges of the chip in spec direction
#    spec_set_fine, outmask_spec, specfit, _, exit_status = \
#        utils.bspline_profile(pix_fit, log_flat_fit, log_ivar_fit,np.ones_like(pix_fit), inmask = inmask_log_fit,
#        nord = 4, upper=logrej, lower=logrej,
#        kwargs_bspline = {'bkspace':spec_samp_fine},kwargs_reject={'groupbadpix':True, 'maxrej': 5})
#
#
#    # Debugging/checking spectral fit
#    if debug:
#        goodbk = spec_set_fine.mask
#        specfit_bkpt, _ = spec_set_fine.value(spec_set_fine.breakpoints[goodbk])
#        was_fit_and_masked = (outmask_spec == False)
#        plt.clf()
#        ax = plt.gca()
#        ax.plot(pix_fit,log_flat_fit, color='k', marker='o', markersize=0.4, mfc='k', fillstyle='full',
#                linestyle='None', label = 'all pixels')
#        ax.plot(pix_fit[was_fit_and_masked],log_flat_fit[was_fit_and_masked], color='red', marker='+',
#                markersize=1.5, mfc='red', fillstyle='full', linestyle='None', label='masked')
#        ax.plot(pix_fit, specfit, color='cornflowerblue', label = 'fit to blaze')
#        ax.plot(spec_set_fine.breakpoints[goodbk], specfit_bkpt, color='lawngreen', marker='o', markersize=2.0,
#                mfc='lawngreen', fillstyle='full', linestyle='None', label='bspline breakpoints')
#        ax.set_ylim((0.99*specfit.min(),1.01*specfit.max()))
#        plt.legend()
#        plt.xlabel('Spectral Pixel')
#        plt.ylabel('log(flat counts)')
#        plt.title('Spectral Fit for slit={:d}'.format(slit))
#        plt.show()
#
#    # Evaluate and save
#    spec_model = np.ones_like(flat)
#    spec_model[thismask], _ = np.exp(spec_set_fine.value(piximg[thismask]))
#    norm_spec = np.ones_like(flat)
#    norm_spec[thismask] = flat[thismask]/np.fmax(spec_model[thismask], 1.0)
#
#    # Flat field pixels for fitting spatial direction
#    # Determine maximum counts in median filtered flat spectrum. Only fit pixels > 0.1 of this maximum
#    specfit_interp = interpolate.interp1d(pix_fit, specfit, kind='linear', bounds_error=False,
#                                          fill_value=-np.inf)
#    log_specfit = specfit_interp(pixvec)
#    specvec = np.exp(log_specfit)
#    spec_sm = utils.fast_running_median(specvec,np.fmax(np.ceil(0.10*nspec).astype(int),10))
#    spec_sm_max = spec_sm.max()
#    fit_spat = thismask & inmask &  (spec_model > 1.0) & (spec_model > 0.1*spec_sm_max) & \
#               (norm_spec > 0.0) & (norm_spec < 1.7)  #& (flat < nonlinear_counts)
#    nfit_spat = np.sum(fit_spat)
#    spat_frac = nfit_spat/np.sum(thismask)
#    # ... Below should be nfit_spat
#    msgs.info('Spatial fit to flatfield for {:}'.format(nfit_spec) + ' pixels')
#    if spat_frac < 0.5:
#        msgs.warn('Spatial flatfield fit is to only {:4.2f}'.format(100*spat_frac) + '% of the pixels on this slit.' +
#                  msgs.newline() + '              Something apperas to be wrong here')
#
#
#    isrt_spat = np.argsort(ximg[fit_spat])
#    ximg_fit = ximg[fit_spat][isrt_spat]
#    norm_spec_fit = norm_spec[fit_spat][isrt_spat]
#    #norm_spec_ivar = np.ones_like(norm_spec_fit)/(spat_illum_thresh**2)
#    nfit_spat = np.sum(fit_spat)
#
#    slitwidth = np.median(slit_righ_in - slit_left_in) # How many pixels wide is the slit at each Y?
#    ximg_resln = spat_samp/slitwidth
#
#    med_width = (np.ceil(nfit_spat*ximg_resln)).astype(int)
#    normimg_raw = utils.fast_running_median(norm_spec_fit,med_width)
#    sig_res = np.fmax(med_width/20.0,0.5)
#    normimg = ndimage.filters.gaussian_filter1d(normimg_raw,sig_res, mode='nearest')
#
#    # mask regions where illumination function takes on extreme values
#    if np.any(np.invert(np.isfinite(normimg))):
#        msgs.error('Inifinities in slit illumination function computation normimg')
#
#    # Determine the breakpoint spacing from the sampling of the ximg
#    ximg_samp = np.median(ximg_fit - np.roll(ximg_fit,1))
#    ximg_1pix = 1.0/slitwidth
#    # Use breakpoints at a spacing of a 1/10th of a pixel, but do not allow a bsp smaller than the typical sampling
#    ximg_bsp  = np.fmax(ximg_1pix/10.0, ximg_samp*1.2)
#    bsp_set = pydl.bspline(ximg_fit,nord=4, bkspace=ximg_bsp)
#    fullbkpt = bsp_set.breakpoints
#
#    spat_set, outmask_spat, spatfit, _, exit_status = \
#        utils.bspline_profile(ximg_fit, normimg, np.ones_like(normimg),np.ones_like(normimg),
#        nord=4,upper=5.0, lower=5.0,fullbkpt = fullbkpt)
#
#    # Evaluate and save
#    illumflat = np.ones_like(flat)
#    illumflat[thismask], _ = spat_set.value(ximg[thismask])
#    norm_spec_spat = np.ones_like(flat)
#    norm_spec_spat[thismask] = flat[thismask]/np.fmax(spec_model[thismask], 1.0)/np.fmax(illumflat[thismask],0.01)
#
#    if tweak_slits:
#        slit_left_out, slit_righ_out, tweak_dict = tweak_slit_edges_brute_force(
#            slit_left_in, slit_righ_in, ximg_fit, normimg, tweak_slits_thresh, tweak_slits_maxfrac)
#        # Recreate all the quantities we need based on the tweaked slits
#        tslits_dict_out = copy.deepcopy(tslits_dict_in)
#        tslits_dict_out['slit_left'][:,slit] = slit_left_out
#        tslits_dict_out['slit_righ'][:,slit] = slit_righ_out
#        # ... This uses the padding in the dict.
#        slitmask_out = pixels.tslits2mask(tslits_dict_out)
#        thismask_out = (slitmask_out == slit)
#        ximg_out, edgmask_out = pixels.ximg_and_edgemask(slit_left_out, slit_righ_out, thismask_out, trim_edg=trim_edg)
#        # Note that nothing changes with the tilts, since these were already extrapolated across the whole image.
#    else:
#        # Generate the edgemask using the original slit boundaries and thismask_in
#        slit_left_out = np.copy(slit_left_in)
#        slit_righ_out = np.copy(slit_righ_in)
#        thismask_out = thismask_in
#        ximg_out = ximg_in
#
#    # Add an approximate pixel axis at the top
#    if debug:
#        plt.clf()
#        ax = plt.gca()
#        ax.plot(ximg_fit, norm_spec_fit, color='k', marker='o', markersize=0.4, mfc='k', fillstyle='full',linestyle='None',
#                label = 'all pixels')
#        #ax.plot(ximg_fit[~imed], norm_spec_fit[~imed], color='darkred', marker='+',markersize=4.0, mfc='red',
#        #        fillstyle='full', linestyle='None', label = 'masked')
#        #ax.plot(ximg_fit[imed], normfit[imed], color='orange', label = 'median spatial profile')
#        ax.plot(ximg_fit, spatfit, color='cornflowerblue', label = 'final slit illumination function')
#        ymin = np.fmax(0.8 * spatfit.min(), 0.5)
#        ymax = 1.2*spatfit.max()
#        ax.set_ylim((np.fmax(0.8 * spatfit.min(), 0.5), 1.2 * spatfit.max()))
#        ax.set_xlim(ximg_fit.min(), ximg_fit.max())
#        plt.vlines(0.0, ymin, ymax, color='lightgreen', linestyle=':', linewidth=2.0, label='original left edge',zorder=8)
#        plt.vlines(1.0,ymin,ymax, color='red',linestyle=':', linewidth = 2.0, label='original right edge',zorder=9)
#        if tweak_slits:
#            if tweak_dict['tweak_left']:
#                label = 'threshold = {:5.2f}'.format(tweak_slits_thresh) + ' % of max of left illumprofile'
#
#                plt.hlines(tweak_slits_thresh*tweak_dict['norm_max_left'], ximg_fit.min(), 0.5, color='lightgreen',
#                           linewidth=3.0,label=label, zorder=10)
#                plt.vlines(tweak_dict['xleft'],ymin,ymax, color='lightgreen',linestyle='--', linewidth = 3.0, label='tweaked left edge',zorder=11)
#            if tweak_dict['tweak_righ']:
#                label = 'threshold = {:5.2f}'.format(tweak_slits_thresh) + ' % of max of right illumprofile'
#
#                plt.hlines(tweak_slits_thresh * tweak_dict['norm_max_righ'], 0.5, ximg_fit.max(), color='red', linewidth=3.0,
#                           label=label, zorder=10)
#                plt.vlines(tweak_dict['xrigh'],ymin,ymax, color='red',linestyle='--', linewidth = 3.0, label='tweaked right edge',zorder=20)
#        plt.legend()
#        plt.xlabel('Normalized Slit Position')
#        plt.ylabel('Normflat Spatial Profile')
#        plt.title('Illumination Function Fit for slit={:d}'.format(slit))
#        plt.show()
#
#    msgs.info('Performing illumination + scattembedered light flat field fit')
#
#    # Flat field pixels for fitting spectral direction
#    isrt_spec = np.argsort(piximg[thismask_out])
#    pix_twod = piximg[thismask_out][isrt_spec]
#    ximg_twod = ximg_out[thismask_out][isrt_spec]
#    norm_twod = norm_spec_spat[thismask_out][isrt_spec]
#
#    fitmask = inmask[thismask_out][isrt_spec] & (np.abs(norm_twod - 1.0) < 0.30)
#    # Here we ignore the formal photon counting errors and simply assume that a typical error per pixel.
#    # This guess is somewhat aribtrary. We then set the rejection threshold with sigrej_illum
#    var_value = 0.01
#    norm_twod_ivar = fitmask.astype(float)/(var_value**2)
#    sigrej_illum = 4.0
#
#    poly_basis = pydl.fpoly(2.0*ximg_twod - 1.0, npoly).T
#
#    # Perform the full 2d fit now
#    twod_set, outmask_twod, twodfit, _ , exit_status = \
#        utils.bspline_profile(pix_twod, norm_twod, norm_twod_ivar,poly_basis,inmask = fitmask, nord = 4,
#        upper=sigrej_illum, lower=sigrej_illum,
#        kwargs_bspline = {'bkspace':spec_samp_coarse},kwargs_reject={'groupbadpix':True, 'maxrej': 10})
#
#    if debug:
#        resid = (norm_twod  - twodfit)
#        badpix = np.invert(outmask_twod) & fitmask
#        goodpix = outmask_twod & fitmask
#        plt.clf()
#        ax = plt.gca()
#        ax.plot(pix_twod[goodpix], resid[goodpix], color='k', marker='o', markersize=0.2, mfc='k', fillstyle='full',linestyle='None',
#                label = 'good points')
#        ax.plot(pix_twod[badpix],resid[badpix], color='red', marker='+',markersize=0.5, mfc='red', fillstyle='full', linestyle='None', label='masked')
#        plt.hlines(sigrej_illum*var_value,pix_twod.min(),pix_twod.max(), color='lawngreen',linestyle='--',
#                   label='rejection thresholds',zorder=10,linewidth=2.0)
#        plt.hlines(-sigrej_illum*var_value,pix_twod.min(),pix_twod.max(), color='lawngreen',linestyle='--',
#                   zorder=10,linewidth=2.0)
#        ax.set_ylim((-0.05,0.05))
#        ax.set_xlim((pix_twod.min(), pix_twod.max()))
#        plt.legend()
#        plt.xlabel('Spectral Pixel')
#        plt.ylabel('Residuals from pixelflat 2-d fit')
#        plt.title('Spectral Residuals for slit={:d}'.format(slit))
#        plt.show()
#
#        plt.clf()
#        ax = plt.gca()
#        ax.plot(ximg_twod[goodpix], resid[goodpix], color='k', marker='o', markersize=0.2, mfc='k', fillstyle='full',
#                linestyle='None',
#                label='good points')
#        ax.plot(ximg_twod[badpix], resid[badpix], color='red', marker='+', markersize=0.5, mfc='red', fillstyle='full',
#                linestyle='None', label='masked')
#        plt.hlines(sigrej_illum * var_value, ximg_twod.min(), ximg_twod.max(), color='lawngreen', linestyle='--',
#                   label='rejection thresholds', zorder=10,linewidth=2.0)
#        plt.hlines(-sigrej_illum * var_value, ximg_twod.min(), ximg_twod.max(), color='lawngreen', linestyle='--',
#                   zorder=10,linewidth=2.0)
#        ax.set_ylim((-0.05, 0.05))
#        ax.set_xlim(-0.02, 1.02)
#        plt.legend()
#        plt.xlabel('Normalized Slit Position')
#        plt.ylabel('Residuals from pixelflat 2-d fit')
#        plt.title('Spatial Residuals for slit={:d}'.format(slit))
#        plt.show()
#
#    # Evaluate and save
#    twod_model = np.ones_like(flat)
#    twod_this = np.zeros_like(twodfit)
#    twod_this[isrt_spec] = twodfit
#    twod_model[thismask_out] = twod_this
#
#    # Compute all the final output images output
#    pixelflat = np.ones_like(flat)
#    flat_model = np.ones_like(flat)
#    flat_model[thismask_out] = twod_model[thismask_out]*np.fmax(illumflat[thismask_out],0.05)*np.fmax(spec_model[thismask_out],1.0)
#    pixelflat[thismask_out] = flat[thismask_out]/flat_model[thismask_out]
#
#    # ToDo Add some code here to treat the edges and places where fits go bad?
#    # Set the pixelflat to 1.0 wherever the flat was nonlinear
#    pixelflat[flat >= nonlinear_counts] = 1.0
#    # Do not apply pixelflat field corrections that are greater than 100% to avoid creating edge effects, etc.
#    # TODO Should we do the same for the illumflat??
#    #pixelflat = np.fmax(np.fmin(pixelflat, 2.0), 0.5)
#    pixelflat = np.clip(pixelflat, 0.5, 2.0)
#
#    return pixelflat, illumflat, flat_model, tilts, thismask_out, slit_left_out, slit_righ_out

# TODO: How much of the rest of this is used?

def flatfield(sciframe, flatframe, bpix, illum_flat=None, snframe=None, varframe=None):
    """ Flat field the input image

    .. todo::
        - Is bpix required?

    Parameters
    ----------
    sciframe : 2d image
    flatframe : 2d image
    illum_flat : 2d image, optional
      slit profile image
    snframe : 2d image, optional
    det : int
      Detector index
    varframe : ndarray
      variance image

    Returns
    -------
    flat-field image
    and updated sigma array if snframe is input
    or updated variance array if varframe is input

    """
    if (varframe is not None) & (snframe is not None):
        msgs.error("Cannot set both varframe and snframe")

    # Fold in the slit profile
    final_flat = flatframe.copy()
    if illum_flat is not None:
        if np.any(illum_flat != 1.0):
            msgs.info('Applying illumination flat')
            final_flat *= illum_flat  # Previous code was modifying flatframe!

    # New image
    retframe = np.zeros_like(sciframe)
    w = np.where(final_flat > 0.0)
    retframe[w] = sciframe[w]/final_flat[w]
    if w[0].size != final_flat.size:
        ww = np.where(final_flat <= 0.0)
        bpix[ww] = 1.0
    # Variance?
    if varframe is not None:
        # This is risky -- Be sure your flat is well behaved!!
        retvar = np.zeros_like(sciframe)
        retvar[w] = varframe[w]/final_flat[w]**2
        return retframe, retvar
    # Error image
    if snframe is None:
        return retframe
    else:
        errframe = np.zeros_like(sciframe)
        wnz = np.where(snframe>0.0)
        errframe[wnz] = retframe[wnz]/snframe[wnz]
        return retframe, errframe


