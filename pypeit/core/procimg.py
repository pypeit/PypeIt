""" Module for image processing core methods

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
from IPython import embed

import numpy as np
import scipy.ndimage
import scipy.optimize
import scipy.signal

from astropy.convolution import convolve, Box2DKernel
from astropy.timeseries import LombScargle
from astropy import stats

from pypeit import msgs
from pypeit import utils
from pypeit.core import parse


# NOTE: This is slower than utils.rebin_evlist by a factor of ~2, but avoids an
# error when the shape of arr is not an integer multiple of boxcar.
def boxcar_average(arr, boxcar):
    """
    Boxcar average an array.

    The dimensionality of the boxcar matches the dimensionality of the provided
    array.  Averages are only performed over the valid region of the array; see
    the return description.  In contrast with a boxcar smoothing algorithm, the
    step between each boxcar average is the size of the box, not one pixel.

    Args:
        arr (`numpy.ndarray`_):
            Array to average.  Can be any shape and dimensionality.
        boxcar (:obj:`int`, :obj:`tuple`):
            Integer number of pixels to average.  If a single integer,
            all axes are averaged with the same size box.  If a
            :obj:`tuple`, the integer is defined separately for each
            array axis; length of tuple must match the number of array
            dimensions.

    Returns:
        `numpy.ndarray`_: The averaged array.  If boxcar is a single
        integer, the returned array shape is::
            
            tuple([s//boxcar for s in arr.shape])

        A similar operation gives the shape when boxcar has elements
        defined for each array dimension.  If the input array is not an
        integer number of boxcar pixels along a given dimension, the
        remainder of the array elements along that dimension are ignored
        (i.e., pixels within the modulus of the array shape and boxcar
        of the end of the array dimension are ignored).
    """
    # Check and configure the input
    _boxcar = (boxcar,)*arr.ndim if isinstance(boxcar, int) else boxcar
    if not isinstance(_boxcar, tuple):
        raise TypeError('Input `boxcar` must be an integer or a tuple.')
    if len(_boxcar) != arr.ndim:
        raise ValueError('Must provide an integer or tuple with one number per array dimension.')

    # Perform the boxcar average over each axis and return the result
    _arr = arr.copy()
    for axis, box in zip(range(arr.ndim), _boxcar):
        _arr = np.add.reduceat(_arr, np.arange(0, _arr.shape[axis], box), axis=axis)/box
    return _arr


# NOTE: This is faster than utils.subsample by a factor of 2-3.
def boxcar_replicate(arr, boxcar):
    """
    Boxcar replicate an array.

    Args:
        arr (`numpy.ndarray`_):
            Array to replicate.  Can be any shape and dimensionality.
        boxcar (:obj:`int`, :obj:`tuple`):
            Integer number of times to replicate each pixel. If a
            single integer, all axes are replicated the same number
            of times. If a :obj:`tuple`, the integer is defined
            separately for each array axis; length of tuple must
            match the number of array dimensions.

    Returns:
        `numpy.ndarray`_: The block-replicated array.
    """
    # Check and configure the input
    _boxcar = (boxcar,)*arr.ndim if isinstance(boxcar, int) else boxcar
    if not isinstance(_boxcar, tuple):
        raise TypeError('Input `boxcar` must be an integer or a tuple.')
    if len(_boxcar) != arr.ndim:
        raise ValueError('Must provide an integer or tuple with one number per array dimension.')

    # Perform the boxcar replication over each axis and return the result
    _arr = arr.copy()
    for axis, box in zip(range(arr.ndim), _boxcar):
        _arr = np.repeat(_arr, box, axis=axis)
    return _arr


def lacosmic(sciframe, saturation=None, nonlinear=1., bpm=None, varframe=None, maxiter=1, grow=1.5,
             remove_compact_obj=True, sigclip=5.0, sigfrac=0.3, objlim=5.0, rm_false_pos=True):
    r"""
    Identify cosmic rays using the L.A.Cosmic algorithm.

    See Peter van Dokkum's `L.A.Cosmic`_ website and `van Dokkum (2001, PASP,
    113, 1420)`_.

    This routine is mostly courtesy of Malte Tewes with some updates/alterations
    by the ``PypeIt`` developers.

    Args:
        sciframe (`numpy.ndarray`_):
            Science frame to process.
        saturation (:obj:`float`, `numpy.ndarray`_, optional):
            The saturation level of the detector in units that match the
            provided frame.  This is used to flag pixels to ignore.  Can be a
            single value or an array; if the latter, the shape must match
            ``sciframe``.  If None, no pixels are flagged as being saturated.
        nonlinear (:obj:`float`, `numpy.ndarray`_, optional):
            The fraction of the saturation level at which the detector response
            becomes non-linear.  Can be a single value or an array; if the
            latter, the shape must match ``sciframe``.
        bpm (`numpy.ndarray`_, optional):
            Bad-pixel mask for the input science frame.  If None, all pixels are
            considered available to be flagged as cosmic rays.  Shape must match
            ``sciframe``.
        varframe (`numpy.ndarray`_, optional):
            The variance in the science frame to process.  If None, the variance
            is estimated by the absolute value of a :math:`5\times 5` boxcar
            median filter of the image.
        maxiter (:obj:`int`, optional):
            Maximum number of detection iterations.
        grow (:obj:`float`, optional):
            The radius (in pixels) used to grow the region masked around
            detected cosmic rays.  See :func:`grow_mask`.  If :math:`\leq 0`,
            the cosmic-ray mask regions are not grown.
        remove_compact_obj (:obj:`bool`, optional):
            Remove likely compact objects from the set of detected cosmic rays.
            This is performed by default here and in the original L.A.Cosmic
            algorithm.
        sigclip (:obj:`float`, optional):
            Threshold for identifying a cosmic ray
        sigfrac (:obj:`float`, optional):
            Fraction of ``sigclip`` used to define the lower threshold used for
            pixels neighboring identified cosmic-rays.
        objlim (:obj:`float`, optional):
            Contrast limit between a cosmic-ray and an underlying object
        rm_false_pos (:obj:`bool`, optional):
            Apply algorithm to detect and remove false positives.  This is not a
            traditional component of the L.A.Cosmic algorithm.

    Returns:
        `numpy.ndarray`_: Boolean array flagging pixels with detected cosmic
        rays; True means the pixel has a cosmic ray.
    """
    # Check input
    if saturation is not None and isinstance(saturation, np.ndarray) \
            and saturation.shape != sciframe.shape:
        msgs.error('Detector pixel saturation array has incorrect shape.')
    if isinstance(nonlinear, np.ndarray) and nonlinear.shape != sciframe.shape:
        msgs.error('Detector nonlinear pixel scale array has incorrect shape.')
    if bpm is not None and bpm.shape != sciframe.shape:
        msgs.error('Bad-pixel mask must match shape of science frame.')
    if varframe is not None and varframe.shape != sciframe.shape:
        msgs.error('Variance frame must match shape of science frame.')

    msgs.info("Detecting cosmic rays with the L.A.Cosmic algorithm")

    # Setup
    # NOTE: We only need a copy of the image if we're performing more than one
    # iteration
    _sciframe = sciframe.copy() if maxiter > 1 else sciframe
    crmask = np.zeros(sciframe.shape, dtype=bool)
    sigcliplow = sigclip * sigfrac

    # Pixels flagged as bad for other reasons are excluded from the list of
    # returned cosmic rays
    _bpm = None if bpm is None else bpm.copy()
    # Include pixel saturation in bad pixel mask if provided
    # Determine if there are saturated pixels
    if saturation is not None:
        satpix = sciframe >= saturation*nonlinear
        if np.any(satpix):
            _bpm = satpix if _bpm is None else _bpm & satpix

    # TODO: Should we be executing boxcar_fill here before the first iteration
    # to remove bad pixels?

    # Initialize the noise model
    if varframe is not None:
        noise = np.sqrt(varframe)
        # NOTE: Inverting the error avoids division by 0 errors
        _inv_err = utils.inverse(noise)

    # Define the kernels
    # - Laplacian kernel
    laplkernel = np.array([[0.0, -1.0, 0.0],
                           [-1.0, 4.0, -1.0],
                           [0.0, -1.0, 0.0]])
    # - Growth kernel
    growkernel = np.ones((3,3), dtype=bool)
    for i in range(maxiter):

        if varframe is None:
            msgs.info("Updating the noise model")
            m5 = scipy.ndimage.median_filter(_sciframe, size=5, mode='mirror')
            noise = np.sqrt(np.absolute(m5))
            # NOTE: Inverting the error avoids division by 0 errors
            _inv_err = utils.inverse(noise)

        # Use the Laplacian transform to construct the image 2nd derivative and
        # get its S/N.  NOTE: the division by 2 in the S/N calculation is from
        # the 2x2 subsampling.  astropy.convolution.convolve gives the same
        # result as scipy.signal.convolve2d, but is nearly a factor of 2 faster.
        msgs.info("Convolving image with Laplacian kernel")
        deriv = convolve(boxcar_replicate(_sciframe, 2), laplkernel, normalize_kernel=False,
                         boundary='extend')
        s = utils.rebin_evlist(np.clip(deriv, 0, None), _sciframe.shape) * _inv_err / 2.0 

        # Remove the large structures
        sp = s - scipy.ndimage.median_filter(s, size=5, mode='mirror')

        # Candidate cosmic rays
        cosmics = sp > sigclip
        ncr = np.sum(cosmics)
        msgs.info(f'Found {ncr} candidate cosmic-ray pixels')

        if _bpm is not None:
            # Remove known bad pixels
            cosmics &= np.logical_not(_bpm)
            ncr = np.sum(cosmics)
            msgs.info(f'Reduced to {ncr} candidates after excluding known bad pixels.')

        if remove_compact_obj:
            # Build the fine structure image
            m3 = scipy.ndimage.median_filter(_sciframe, size=3, mode='mirror')
            m37 = scipy.ndimage.median_filter(m3, size=7, mode='mirror')
            # TODO: How does clip treat NaNs?
            f = np.clip((m3 - m37) * _inv_err, 0.01, None)
            # Require cosmics to have significant contrast
            cosmics &= sp/f > objlim
            ncr = np.sum(cosmics)
            msgs.info(f'Reduced to {ncr} candidates after excluding compact objects.')

        # What follows is a special treatment for neighbors, with more relaxed
        # constraints.
        msgs.info("Finding neighboring pixels affected by cosmic rays")
        # We grow these cosmics a first time to determine the immediate
        # neighborhod, keeping those that also meet the S/N requirement
        cosmics = scipy.ndimage.binary_dilation(cosmics, structure=growkernel)
        cosmics &= sp > sigclip

        # Now we repeat this procedure, but lower the detection limit to sigmalimlow :
        cosmics = scipy.ndimage.binary_dilation(cosmics, structure=growkernel)
        cosmics &= sp > sigcliplow
        ncr = np.sum(cosmics)
        msgs.info(f'Changed to {ncr} candidates after evaluating neighboring pixels.')

        if _bpm is not None:
            # Remove known bad pixels
            cosmics &= np.logical_not(_bpm)
            ncr = np.sum(cosmics)
            msgs.info(f'Reduced to {ncr} candidates after excluding known bad pixels.')

        # Determine how many new cosmics were found
        nnew = np.sum(np.logical_not(crmask) & cosmics)
        crmask |= cosmics
        msgs.info(f'Iteration {i+1}: {np.sum(crmask)} pixels identified as cosmic rays '
                  f'({nnew} are new)')
        if nnew == 0 or i == maxiter - 1:
            # TODO: Warn the user if the maximum number of iterations was
            # reached (and maxiter!=1)?
            break

        # Prepare for the next iteration
        msgs.info('Preparing for next iteration')
        _sciframe = boxcar_fill(_sciframe, 5, bpm=crmask if _bpm is None else crmask | _bpm)

    if not rm_false_pos:
        return grow_mask(crmask, grow) if grow > 0 else crmask

    # Additional algorithms (not traditionally implemented by LA cosmic) to
    # remove some false positives.
    #msgs.work("The following algorithm would be better on the rectified, tilts-corrected image")
    filt  = scipy.ndimage.sobel(sciframe, axis=1, mode='constant')
    _inv_mad = utils.inverse(np.sqrt(np.abs(sciframe))) # Avoid divisions by 0
    filty = scipy.ndimage.sobel(filt * _inv_mad, axis=0, mode='constant')
    # TODO: Can we skip this now that we're not dividing by 0?
    filty[np.isnan(filty)] = 0.0

    sigimg = cr_screen(filty)

    sigsmth = scipy.ndimage.gaussian_filter(sigimg, 1.5)
    sigsmth[np.isnan(sigsmth)] = 0.0

    crmask &= sigsmth > sigclip
    msgs.info(f'{np.sum(crmask)} pixels identified as cosmic rays after removing false positives')
    return grow_mask(crmask, grow) if grow > 0 else crmask


def boxcar_fill(img, width, bpm=None, maxiter=None, fill_value=np.nan):
    """
    Use convolution with a boxcar kernel to iteratively fill masked regions of
    an image.

    .. warning::
        
        Depending the size of the kernel, the masked regions, and the maximum
        number of iterations, the procedure may not be able to fill all masked
        pixels.  Any left-over masked pixels are returned with the provided
        ``fill_value``.

    Args:
        img (`numpy.ndarray`_, `numpy.ma.MaskedArray`_):
            Image to fill.  Must be 2D.  If an unmasked array, function will
            simply return ``img`` if ``mask`` is not provided.
        width (:obj:`int`):
            Width of the square box to use for the convolution kernel.  Must be
            odd.
        bpm (`numpy.ndarray`_, optional):
            Image bad-pixel mask.  If ``img`` is a masked array, this mask is
            combined with the ``img`` mask attribute.
        maxiter (:obj:`int`, optional):
            The maximum number of fill iterations.  If None, iterations continue
            until all masked pixels are filled.
        fill_value (:obj:`float`, optional):
            If the number of allowed iterations are insufficient to fill the
            array, replace masked pixels with this value.

    Returns:
        `numpy.ndarray`_: Array with the same shape as ``img`` with the masked
        pixels filled by the average of the unmasked pixels in the boxcar
        kernel.
    """
    # TODO: A 2D median filter that accounts for masked pixels may be better,
    # but I couldn't find a canned algorithm.  Could also imagine using
    # different kernels (e.g., a Gaussian).  See: 
    #   https://docs.astropy.org/en/stable/convolution/index.html

    # Check input
    if isinstance(img, np.ndarray) and bpm is None:
        return img
    # Handle masked array input
    if isinstance(img, np.ma.MaskedArray):
        _msk = np.ma.getmaskarray(img).copy()
        if bpm is not None:
            _msk |= bpm
        _img = img.data.copy()
    else:
        _img = img.copy()
        _msk = bpm.copy()

    # Construct the convolution kernel
    fillkernel = Box2DKernel(width)

    # Replace masked pixels with np.nan so that astropy.convolution.convolve
    # ignores them.
    _img[_msk] = np.nan

    # Iteratively fill the masked pixels
    filliter = 0
    while np.any(_msk) and (maxiter is None or filliter < maxiter):
        # TODO: Suppress the astropy warning when the returned array still has
        # NaNs in it.
        _img[_msk] = convolve(_img, fillkernel)[_msk]
        _msk = np.isnan(_img)
        filliter += 1
    if np.any(_msk):
        # Fill any left-over masked pixels
        _img[_msk] = fill_value
    return _img


def cr_screen(a, mask_value=0.0, spatial_axis=1):
    r"""
    Calculate the significance of pixel deviations from the median along
    the spatial direction.

    No type checking is performed of the input array; however, the
    function assumes floating point values.

    Args:
        a (`numpy.ndarray`_): Input 2D array
        mask_value (float): (**Optional**) Values to ignore during the
            calculation of the median.  Default is 0.0.
        spatial_axis (int): (**Optional**) Axis along which to calculate
            the median.  Default is 1.

    Returns:
        `numpy.ndarray`_: Returns a map of :math:`|\Delta_{i,j}|/\sigma_j`,
        where :math:`\Delta_{i,j}` is the difference between the pixel
        value and the median along axis :math:`i` and :math:`\sigma_j`
        is robustly determined using the median absolute deviation,
        :math:`sigma_j = 1.4826 MAD`.
    """
    # Check input
    if len(a.shape) != 2:
        msgs.error('Input array must be two-dimensional.')
    if spatial_axis not in [0,1]:
        msgs.error('Spatial axis must be 0 or 1.')

    # Mask the pixels equal to mask value: should use np.isclose()
    _a = np.ma.MaskedArray(a, mask=(a==mask_value))
    # Get the median along the spatial axis
    meda = np.ma.median(_a, axis=spatial_axis)
    # Get a robust measure of the standard deviation using the median
    # absolute deviation; 1.4826 factor is the ratio of sigma/MAD
    d = np.absolute(_a - meda[:,None])
    mada = 1.4826*np.ma.median(d, axis=spatial_axis)
    # Return the ratio of the difference to the standard deviation
    return np.ma.divide(d, mada[:,None]).filled(mask_value)


# NOTE: This is a factor of a few times faster than the previous version.  The
# speed improvement is better for smaller images.  For larger images (e.g.,
# 2048x2048), the improvement is about a factor of 3.
def grow_mask(mask, radius):
    """
    Grow pixels flagged as True in a boolean mask by the provided radius.

    This is largely a convience wrapper for `scipy.ndimage.binary_dilation`_.

    Args:
        mask (`numpy.ndarray`_):
            Boolean mask to process.  Pixels flagged as True are expanded into
            circles with the provided radius.
        radius (scalar-like):
            Radius in pixels to grow the mask.

    Returns:
        `numpy.ndarray`_: The boolean mask grown with the masked region grown by
        the provided radius.
    """
    # Prep for the dilation structure
    size = int(radius*2+1)
    if size % 2 == 0:
        size += 1
    x, y = np.meshgrid(np.arange(size) - size//2, np.arange(size) - size//2)
    # Dilate the mask
    return scipy.ndimage.binary_dilation(mask, structure=x**2 + y**2 <= radius**2)


def gain_frame(amp_img, gain):
    """
    Generate an image with the gain for each pixel.

    Args:
        amp_img (`numpy.ndarray`_):
            Integer array that identifies which (1-indexed) amplifier
            was used to read each pixel.
        gain (array-like):
            List of amplifier gain values in e-/ADU.  Must be that the gain for
            amplifier 1 is provided by `gain[0]`, etc.

    Returns:
        `numpy.ndarray`_: Image with the gain for each pixel.
    """
    # TODO: Remove this or actually do it.
    # msgs.warn("Should probably be measuring the gain across the amplifier boundary")
    # Build and return the gain image
    gain_img = np.zeros_like(amp_img, dtype=float)
    for i,_gain in enumerate(gain):
        gain_img[amp_img == i+1] = _gain
    return gain_img


def rn2_frame(datasec_img, ronoise, units='e-', gain=None, digitization=False):
    r"""
    Construct a readnoise variance image.

    Provided the detector readnoise and gain for each amplifier, this constructs
    an image with the combination of the readnoise and digitization (or
    quantization) noise expected for a single detector readout.  Digitization
    noise is a fixed :math:`\sqrt{1/12}` ADU [1]_ [2]_, derived as the second
    moment of a uniform distribution between values of -1/2 to 1/2 (i.e., the
    variance associated with converting a number of electrons into an ADU
    integer quantized by the gain).  The digitization noise is typically much
    smaller than the readnoise, unless the gain is very large, and, depending on
    how it was measured, the digitization noise is most often incorporated in
    the documented readnoise of the given instrument.  To include the
    digitization noise in the variance, you must provide ``gain`` and set
    ``digitization=True``.

    The variance calculation in electrons is :math:`V = {\rm RN}^2 +
    \gamma^2/12`, when including the digitization noise, and simply :math:`V =
    {\rm RN}^2` otherwise; where RN is the readnoise and :math:`\gamma` is the
    gain in e-/ADU.  In the rare case one would need the units in ADU, the
    returned variance is :math:`V/\gamma^2`.
    
    .. [1] `Newberry (1991, PASP, 103, 122) <https://ui.adsabs.harvard.edu/abs/1991PASP..103..122N/abstract>`_
    .. [2] `Merline & Howell (1995, ExA, 6, 163) <https://ui.adsabs.harvard.edu/abs/1995ExA.....6..163M/abstract>`_

    Args:
        datasec_img (`numpy.ndarray`_):
            An integer array indicating the 1-indexed amplifier used to read
            each pixel in the main data section of the detector.  Values of 0
            are ignored.  Amplifier numbers are expected sequential and match
            the number of readnoise and gain values provided.  The shape of this
            image dictates the shape of the output readnoise variance image.
        ronoise (:obj:`float`, array-like):
            The value of the readnoise for each amplifier in electrons (e-).  If
            there is only one amplifier, this can be provided as a single float.
        units (:obj:`str`, optional):
            Units for the output variance.  Options are ``'e-'`` for variance in
            square electrons (counts) or ``'ADU'`` for square ADU.
        gain (:obj:`float`, array-like, optional):
            The value of the gain for each amplifier in e-/ADU.  If
            ``digitization`` is False, this is ignored.
        digitization (:obj:`bool`, optional):
            Include digitization error in the calculation.  If True, ``gain``
            *must* be provided.

    Returns:
        `numpy.ndarray`_: The image variance resulting from reading the detector
        in the selected units for each pixel.  The shape is the same as
        ``datasec_img``.  Pixels where ``datasec_img`` is 0 are set to 0.
    """
    # Check units
    if units not in ['e-', 'ADU']:
        msgs.error(f"Unknown units: {units}.  Must be 'e-' or 'ADU'.")
    if gain is None and (digitization or units == 'ADU'):
        msgs.error('If including digitization error or return units in ADU, must provide gain.')

    # Determine the number of amplifiers from the datasec image
    _datasec_img = datasec_img.astype(int)
    numamplifiers = np.amax(_datasec_img)
    if numamplifiers == 0:
        msgs.error('Amplifier identification image (datasec_img) does not have any values larger '
                   'than 0!  The image should indicate the 1-indexed integer of the amplifier '
                   'used to read each pixel.')

    # Check the number of RN values
    _ronoise = np.atleast_1d(ronoise) if isinstance(ronoise, (list, np.ndarray)) \
                        else np.array([ronoise])
    if len(_ronoise) != numamplifiers:
        msgs.error('Must provide a read-noise for each amplifier.')

    # Get the amplifier indices
    indx = np.logical_not(_datasec_img == 0)
    amp = _datasec_img[indx] - 1

    # Instantiate the output image.  Any pixels without an assigned amplifier
    # are given a noise of 0.
    var = np.zeros(_datasec_img.shape, dtype=float)
    var[indx] = (_ronoise**2)[amp]

    if not digitization and units == 'e-':
        return var

    # Check the number of gain values
    _gain = np.atleast_1d(gain) if isinstance(gain, (list, np.ndarray)) else np.array([gain])
    if len(_gain) != numamplifiers:
        msgs.error('Must provide a gain for each amplifier.')

    if digitization:
        # Add in the digitization error
        var[indx] += (_gain**2/12)[amp]

    if units == 'ADU':
        # Convert to ADUs
        var[indx] /= (_gain**2)[amp]

    return var


def rect_slice_with_mask(image, mask, mask_val=1):
    """
    Generate rectangular slices from a mask image.

    Args:
        image (`numpy.ndarray`_):
            Image to mask
        mask (`numpy.ndarray`_):
            Mask image
        mask_val (:obj:`int`, optional):
            Value to mask on

    Returns:
        :obj:`tuple`: The image at mask values and a 2-tuple with the
        :obj:`slice` objects that select the masked data.
    """
    pix = np.where(mask == mask_val)
    slices = (slice(np.min(pix[0]), np.max(pix[0])+1), slice(np.min(pix[1]), np.max(pix[1])+1))
    return image[slices], slices


def subtract_overscan(rawframe, datasec_img, oscansec_img, method='savgol', params=[5,65],
                      var=None):
    """
    Subtract overscan.

    Args:
        rawframe (`numpy.ndarray`_):
            Frame from which to subtract overscan.  Must be 2d.
        datasec_img (`numpy.ndarray`_):
            An array the same shape as ``rawframe`` that identifies the pixels
            associated with the data on each amplifier; 0 for no data, 1 for
            amplifier 1, 2 for amplifier 2, etc.
        oscansec_img (:obj:`numpy.ndarray`):
            An array the same shape as ``rawframe`` that identifies the pixels
            associated with the overscan region on each amplifier; 0 for no
            data, 1 for amplifier 1, 2 for amplifier 2, etc.
        method (:obj:`str`, optional):
            The method used to fit the overscan region.  Options are
            polynomial, savgol, median.
        params (:obj:`list`, optional):
            Parameters for the overscan subtraction.  For ``method=polynomial``,
            set ``params`` to the order, number of pixels, number of repeats;
            for ``method=savgol``, set ``params`` to the order and window size;
            for ``method=median``, ``params`` are ignored.
        var (`numpy.ndarray`_, optional):
            Variance in the raw frame.  If provided, must have the same shape as
            ``rawframe`` and used to estimate the error in the overscan
            subtraction.  The estimated error is the standard error in the
            median for the pixels included in the overscan correction.  This
            estimate is also used for the ``'savgol'`` method as an upper limit.
            If None, no variance in the overscan subtraction is calculated, and
            the 2nd object in the returned tuple is None.

    Returns:
        :obj:`tuple`: The input frame with the overscan region subtracted and an
        estimate of the variance in the overscan subtraction; both have the same
        shape as the input ``rawframe``.  If ``var`` is no provided, the 2nd
        returned object is None.
    """
    # Check input
    if method.lower() not in ['polynomial', 'savgol', 'median']:
        msgs.error(f'Unrecognized overscan subtraction method: {method}')
    if rawframe.ndim != 2:
        msgs.error('Input raw frame must be 2D.')
    if datasec_img.shape != rawframe.shape:
        msgs.error('Datasec image must have the same shape as the raw frame.')
    if oscansec_img.shape != rawframe.shape:
        msgs.error('Overscan section image must have the same shape as the raw frame.')
    if var is not None and var.shape != rawframe.shape:
        msgs.error('Variance image must have the same shape as the raw frame.')

    # Copy the data so that the subtraction is not done in place
    no_overscan = rawframe.copy()
    _var = None if var is None else np.zeros(var.shape, dtype=float)

    # Amplifiers
    amps = np.unique(datasec_img[datasec_img > 0]).tolist()

    # Perform the overscan subtraction for each amplifier
    for amp in amps:
        # Pull out the overscan data
        if np.sum(oscansec_img == amp) == 0:
            msgs.error(f'No overscan region for amplifier {amp+1}!')
        overscan, os_slice = rect_slice_with_mask(rawframe, oscansec_img, amp)
        if var is not None:
            osvar = var[os_slice]
        # Pull out the real data
        if np.sum(datasec_img == amp) == 0:
            msgs.error(f'No data region for amplifier {amp+1}!')
        data, data_slice = rect_slice_with_mask(rawframe, datasec_img, amp)

        # Shape along at least one axis must match
        if not np.any([dd == do for dd, do in zip(data.shape, overscan.shape)]):
            msgs.error('Overscan sections do not match amplifier sections for '
                       'amplifier {0}'.format(amp))
        compress_axis = 1 if data.shape[0] == overscan.shape[0] else 0

        # Fit/Model the overscan region
        osfit = np.median(overscan) if method.lower() == 'median' \
                    else np.median(overscan, axis=compress_axis)
        if var is not None:
            # pi/2 coefficient yields asymptotic variance in the median relative
            # to the error in the mean
            osvar = np.pi/2*(np.sum(osvar)/osvar.size**2 if method.lower() == 'median' 
                             else np.sum(osvar, axis=compress_axis)/osvar.shape[compress_axis]**2)
        if method.lower() == 'polynomial':
            # TODO: Use np.polynomial.polynomial.polyfit instead?
            c = np.polyfit(np.arange(osfit.size), osfit, params[0])
            ossub = np.polyval(c, np.arange(osfit.size))
        elif method.lower() == 'savgol':
            ossub = scipy.signal.savgol_filter(osfit, params[1], params[0])
        elif method.lower() == 'median':
            # Subtract scalar and continue
            no_overscan[data_slice] -= osfit
            if var is not None:
                _var[data_slice] = osvar
            continue

        # Subtract along the appropriate axis
        no_overscan[data_slice] -= (ossub[:, None] if compress_axis == 1 else ossub[None, :])
        if var is not None:
            _var[data_slice] = (osvar[:,None] if compress_axis == 1 else osvar[None,:])

    return no_overscan, _var


def subtract_pattern(rawframe, datasec_img, oscansec_img, frequency=None, axis=1): #, debug=False):
    """
    Subtract a sinusoidal pattern from the input rawframe. The algorithm
    calculates the frequency of the signal, generates a model, and subtracts
    this signal from the data. This sinusoidal pattern noise was first
    identified in KCWI, but the source of this pattern noise is not currently
    known. The pattern model is generated with a three step process:

        #. Given a first guess at the frequency, calculate how frequency depends
           on pixel row (slight linear dependence)

        #. Using the frequency model, calculate the amplitude of the signal
           (usually a constant for all pixel rows)

        #. Using the model of the frequency and the amplitude, calculate the
           phase of the signal for each pixel row.  Note that the phase is
           different for each pixel row.

    A complete detector model is generated for the pattern noise using the
    frequency+amplitude+phase, and an estimate of the improved effective read
    noise is provided.

    Args:
        rawframe (`numpy.ndarray`_):
            Frame from which to subtract overscan
        numamplifiers (:obj:`int`):
            Number of amplifiers for this detector.
        datasec_img (`numpy.ndarray`_):
            An array the same shape as rawframe that identifies
            the pixels associated with the data on each amplifier.
            0 for not data, 1 for amplifier 1, 2 for amplifier 2, etc.
        oscansec_img (`numpy.ndarray`_):
            An array the same shape as rawframe that identifies
            the pixels associated with the overscan region on each
            amplifier.
            0 for not data, 1 for amplifier 1, 2 for amplifier 2, etc.
        frequency (:obj:`float`, :obj:`list`, optional):
            The frequency (or list of frequencies - one for each amplifier)
            of the sinusoidal pattern. If None, the frequency of each amplifier
            will be determined from the overscan region.
        axis (:obj:`int`, optional):
            Which axis should the pattern subtraction be applied?

    Returns:
        `numpy.ndarray`_: The input frame with the pattern subtracted
    """
    msgs.info("Analyzing detector pattern")

    # Copy the data so that the subtraction is not done in place
    frame_orig = rawframe.copy()
    outframe = rawframe.copy()
    tmp_oscan = oscansec_img.copy()
    tmp_data = datasec_img.copy()
    if axis == 0:
        frame_orig = rawframe.copy().T
        outframe = rawframe.copy().T
        tmp_oscan = oscansec_img.copy().T
        tmp_data = datasec_img.copy().T

    # Amplifiers
    amps = np.sort(np.unique(tmp_data[tmp_data > 0])).tolist()

    # Estimate the frequency in each amplifier (then average over all amps)
    if frequency is None:
        frq = np.zeros(len(amps))
        for aa, amp in enumerate(amps):
            pixs = np.where(tmp_oscan == amp)
            #pixs = np.where((tmp_oscan == amp) | (tmp_data ==  amp))
            cmin, cmax = np.min(pixs[0]), np.max(pixs[0])
            rmin, rmax = np.min(pixs[1]), np.max(pixs[1])
            frame = frame_orig[cmin:cmax, rmin:rmax].astype(np.float64)
            frq[aa] = pattern_frequency(frame)
        frequency = np.mean(frq)

    # Perform the overscan subtraction for each amplifier
    for aa, amp in enumerate(amps):
        # Get the frequency to use for this amplifier
        if isinstance(frequency, list):
            # if it's a list, then use a different frequency for each amplifier
            use_fr = frequency[aa]
        else:
            # float
            use_fr = frequency

        # Extract overscan
        overscan, os_slice = rect_slice_with_mask(frame_orig, tmp_oscan, amp)
        # Extract overscan+data
        oscandata, osd_slice = rect_slice_with_mask(frame_orig, tmp_oscan+tmp_data, amp)
        # Subtract the DC offset
        overscan -= np.median(overscan, axis=1)[:, np.newaxis]

        # STEP 1 - calculate how frequency depends on pixel row
        # Estimate the frequency at each pixel
        all_rows = np.arange(overscan.shape[0])
        all_freq = np.zeros(overscan.shape[0])
        pixels = np.arange(overscan.shape[1])
        for ii in range(overscan.shape[0]):
            sgnl = overscan[ii,:]
            LSfreq, power = LombScargle(pixels, sgnl).autopower(minimum_frequency=use_fr*(1-100/frame_orig.shape[1]), maximum_frequency=use_fr*(1+100/frame_orig.shape[1]), samples_per_peak=10)
            bst = np.argmax(power)
            cc = np.polyfit(LSfreq[bst-2:bst+3],power[bst-2:bst+3],2)
            all_freq[ii] = -0.5*cc[1]/cc[0]
        cc = np.polyfit(all_rows, all_freq, 1)
        frq_mod = np.polyval(cc, all_rows) * (overscan.shape[1]-1)

        # Convert frequency to the size of the overscan region
        msgs.info("Subtracting detector pattern from amplifier {0:d} with frequency = {1:f}".format(amp, use_fr))

        # Get a first guess of the amplitude and phase information
        xdata, step = np.linspace(0.0, 1.0, overscan.shape[1], retstep=True)
        xdata_all = (np.arange(osd_slice[1].start, osd_slice[1].stop) - os_slice[1].start) * step
        tmpamp = np.fft.rfft(overscan, axis=1)
        idx = (np.arange(overscan.shape[0]), np.argmax(np.abs(tmpamp), axis=1))
        # Convert result to amplitude and phase
        amps = (np.abs(tmpamp))[idx] * (2.0 / overscan.shape[1])

        # STEP 2 - Using th emodel frequency, calculate how amplitude depends on pixel row (usually constant)
        # Use the above to as initial guess parameters for a chi-squared minimisation of the amplitudes
        msgs.info("Measuring amplitude-pixel dependence of amplifier {0:d}".format(amp))
        nspec = overscan.shape[0]
        model_pattern = np.zeros_like(oscandata)
        cosfunc = lambda xarr, *p: p[0] * np.cos(2.0 * np.pi * xarr + p[1])
        nsamp = xdata.size // 10
        bins = np.linspace(0.0, 1.0, nsamp)
        cent = 0.5 * (bins[1:] + bins[:-1])
        amps_fit = np.zeros(nspec)
        for ii in range(nspec):
            # Register all values between 0-1
            resid = (frq_mod[ii] * xdata) % 1
            hist, _ = np.histogram(resid, bins=bins, weights=overscan[ii, :])
            norm, _ = np.histogram(resid, bins=bins)
            hist *= utils.inverse(norm)
            # Only use the good pixels
            wgd = norm != 0
            try:
                # Now fit it
                popt, pcov = scipy.optimize.curve_fit(
                    cosfunc, cent[wgd], hist[wgd], p0=[amps[ii], 0.0],
                    bounds=([0, -np.inf],[np.inf, np.inf])
                )
            except ValueError:
                msgs.warn("Input data invalid for pattern subtraction of row {0:d}/{1:d}".format(ii + 1, overscan.shape[0]))
                continue
            except RuntimeError:
                msgs.warn("Pattern subtraction fit failed for row {0:d}/{1:d}".format(ii + 1, overscan.shape[0]))
                continue
            amps_fit[ii] = popt[0]
        # Construct a model of the amplitudes as a fucntion of spectral pixel
        xspec = np.arange(nspec)
        amp_mod = np.polyval(np.polyfit(xspec, amps_fit, 1), xspec)

        # STEP 3 - Using the model frequency and amplitude, calculate the phase of every pixel row
        # Now determine the phase, given a prior on the amplitude and frequency
        msgs.info("Calculating pattern phases of amplifier {0:d}".format(amp))
        cosfunc = lambda xarr, *p: np.cos(2.0 * np.pi * xarr + p[0])
        cosfunc_full = lambda xarr, *p: p[0] * np.cos(2.0 * np.pi * p[1] * xarr + p[2])
        for ii in range(nspec):
            resid = (frq_mod[ii] * xdata) % 1
            hist, _ = np.histogram(resid, bins=bins, weights=overscan[ii, :])
            norm, _ = np.histogram(resid, bins=bins)
            hist *= utils.inverse(norm)
            hist /= amp_mod[ii]  # Normalise so that the amplitude is ~1
            # Only use the good pixels
            wgd = norm != 0
            try:
                # Now fit it
                popt, pcov = scipy.optimize.curve_fit(
                    cosfunc, cent[wgd], hist[wgd], p0=[0.0],
                    bounds=([-np.inf], [np.inf])
                )
            except ValueError:
                msgs.warn("Input data invalid for pattern subtraction of row {0:d}/{1:d}".format(ii + 1, overscan.shape[0]))
                continue
            except RuntimeError:
                msgs.warn("Pattern subtraction fit failed for row {0:d}/{1:d}".format(ii + 1, overscan.shape[0]))
                continue
            # Calculate the model pattern, given the amplitude, frequency and phase information
            model_pattern[ii, :] = cosfunc_full(xdata_all, amp_mod[ii], frq_mod[ii], popt[0])

        # Estimate the improvement of the effective read noise
        tmp = outframe.copy()
        tmp[osd_slice] -= model_pattern
        mod_oscan, _ = rect_slice_with_mask(tmp, tmp_oscan, amp)
        old_ron = stats.sigma_clipped_stats(overscan, sigma=5)[-1]
        new_ron = stats.sigma_clipped_stats(overscan-mod_oscan, sigma=5)[-1]
        msgs.info(f'Effective read noise of amplifier {amp} reduced by a factor of {old_ron/new_ron:.2f}x')

        # Subtract the model pattern from the full datasec
        outframe[osd_slice] -= model_pattern

    # Transpose if the input frame if applied along a different axis
    if axis == 0:
        outframe = outframe.T
    # Return the result
    return outframe


def pattern_frequency(frame, axis=1):
    """
    Using the supplied 2D array, calculate the pattern frequency
    along the specified axis.

    Args:
        frame (`numpy.ndarray`_):
            2D array to measure the pattern frequency
        axis (:obj:`int`, optional):
            Which axis should the pattern frequency be measured?

    Returns:
        :obj:`float`: The frequency of the sinusoidal pattern.
    """
    # For axis=0, transpose
    arr = frame.copy()
    if axis == 0:
        arr = frame.T
    elif axis != 1:
        msgs.error("frame must be a 2D image, and axis must be 0 or 1")

    # Calculate the output image dimensions of the model signal
    # Subtract the DC offset
    arr -= np.median(arr, axis=1)[:, np.newaxis]

    # Compute the Fourier transform to obtain an estimate of the dominant frequency component
    amp = np.fft.rfft(arr, axis=1)
    idx = (np.arange(arr.shape[0]), np.argmax(np.abs(amp), axis=1))

    # Construct the variables of the sinusoidal waveform
    frqs = idx[1]

    min_fr = np.median(frqs-1)/(arr.shape[1]-1)
    max_fr = np.median(frqs+1)/(arr.shape[1]-1)
    all_freq = np.zeros(arr.shape[0])
    pixels = np.arange(arr.shape[1])
    for ii in range(arr.shape[0]):
        sgnl = arr[ii, :]
        LSfreq, power = LombScargle(pixels, sgnl).autopower(minimum_frequency=min_fr, maximum_frequency=max_fr, samples_per_peak=10)
        bst = 2+np.argmax(power[2:-2])  # Ignore edges, and add 2 to get the index of the original array - allows a quadratic function to be fit to the highest power.
        cc = np.polyfit(LSfreq[bst-2:bst+3], power[bst-2:bst+3], 2)
        all_freq[ii] = -0.5*cc[1]/cc[0]
    medfrq = np.median(all_freq)
    return medfrq


# TODO: Provide a replace_pixels method that does this on a pixel by
# pixel basis instead of full columns.
def replace_columns(img, bad_cols, replace_with='mean', copy=False):
    """
    Replace bad image columns.

    Args:
        img (`numpy.ndarray`_):
            A 2D array with image values to replace.
        bad_cols (`numpy.ndarray`_):
            Boolean array selecting bad columns in `img`.  Must have the
            correct shape.
        replace_with (:obj:`str`, optional):
            Method to use for the replacements.  Can be 'mean' (see
            :func:`replace_column_mean`) or 'linear' (see
            :func:`replace_column_linear`).
        copy (:obj:`bool`, optional):
            Copy `img` to a new array before making any
            modifications.  Otherwise, `img` is modified in-place.

    Returns:
        `numpy.ndarray`_: The modified image, which is either a new
        array or points to the in-place modification of `img` according
        to the value of `copy`.
    """
    # Check
    if img.ndim != 2:
        msgs.error('Images must be 2D!')
    if bad_cols.size != img.shape[1]:
        msgs.error('Bad column array has incorrect length!')
    if np.all(bad_cols):
        msgs.error('All columns are bad!')

    _img = img.copy() if copy else img

    if np.sum(bad_cols) == 0:
        # No bad columns
        return _img

    # Find the starting/ending indices of adjacent bad columns
    borders = np.zeros(img.shape[1], dtype=int)
    borders[bad_cols] = 1
    borders = borders - np.roll(borders,1)
    if borders[0] == -1:
        borders[0] = 0

    # Get edge indices and deal with edge cases
    lindx = borders == 1
    ledges = np.where(lindx)[0] if np.any(lindx) else [0]
    rindx = borders == -1
    redges = np.where(rindx)[0] if np.any(rindx) else [img.shape[1]]
    if ledges[0] > redges[0]:
        ledges = np.append([0], ledges)
    if ledges[-1] > redges[-1]:
        redges = np.append(redges, [img.shape[1]])
    # If this is tripped, there's a coding error
    assert len(ledges) == len(redges), 'Problem in edge setup'

    # Replace the image values
    if replace_with == 'mean':
        for l,r in zip(ledges, redges):
            replace_column_mean(_img, l, r)
    elif replace_with == 'linear':
        for l,r in zip(ledges, redges):
            replace_column_linear(_img, l, r)
    else:
        msgs.error('Unknown replace_columns method.  Must be mean or linear.')
    return _img


def replace_column_mean(img, left, right):
    """
    Replace the column values between left and right indices for all
    rows by the mean of the columns just outside the region.

    Columns at the end of the image with no left or right reference
    column (`left==0` or `right==img.shape[1]`) are just replaced by the
    closest valid column.

    Args:
        img (`numpy.ndarray`_):
            Image with values to both use and replace.
        left (:obj:`int`):
            Inclusive starting column index.
        right (:obj:`int`):
            Exclusive ending column index.
    """
    if left == 0:
        img[:,left:right] = img[:,right][:,None]
        return
    if right == img.shape[1]:
        img[:,left:] = img[:,left-1][:,None]
        return
    img[:,left:right] = 0.5*(img[:,left-1]+img[:,right])[:,None]


def replace_column_linear(img, left, right):
    """
    Replace the column values between left and right indices for all
    rows by a linear interpolation between the columns just outside the
    region.

    If possible, extrapolation is used for columns at the end of the
    image with no left or right reference column (`left==0` or
    `right==img.shape[1]`) using the two most adjacent columns.
    Otherwise, this function calls :func:`replace_column_mean`.

    Args:
        img (`numpy.ndarray`_):
            Image with values to both use and replace.
        left (:obj:`int`):
            Inclusive starting column index.
        right (:obj:`int`):
            Exclusive ending column index.
    """
    if left == 0 and right > img.shape[1]-2 or right == img.shape[1] and left < 2:
        # No extrapolation available so revert to mean
        return replace_column_mean(img, left, right)
    if left == 0:
        # Extrapolate down
        img[:,:right] = (img[:,right+1]-img[:,right])[:,None]*np.arange(right)[None,:] \
                            + img[:,right][:,None]
        return
    if right == img.shape[1]:
        # Extrapolate up
        img[:,left:] = (img[:,left-1]-img[:,left-2])[:,None]*np.arange(right-left)[None,:] \
                            + img[:,left-2][:,None]
        return
    # Interpolate
    img[:,left:right] = np.divide(img[:,right]-img[:,left-1],right-left+1)[:,None] \
                            * (np.arange(right-left)+1)[None,:] + img[:,left-1][:,None]


def trim_frame(frame, mask):
    """
    Trim the masked regions from a frame.

    Args:
        frame (`numpy.ndarray`_):
            Image to be trimmed
        mask (`numpy.ndarray`_):
            Boolean image set to True for values that should be trimmed
            and False for values to be returned in the output trimmed
            image.

    Return:
        `numpy.ndarray`_: Trimmed image

    Raises:
        PypitError:
            Error raised if the trimmed image includes masked values
            because the shape of the valid region is odd.
    """
    # TODO: Should check for this failure mode earlier
    if np.any(mask[np.logical_not(np.all(mask,axis=1)),:][:,np.logical_not(np.all(mask,axis=0))]):
        msgs.error('Data section is oddly shaped.  Trimming does not exclude all '
                   'pixels outside the data sections.')
    return frame[np.logical_not(np.all(mask,axis=1)),:][:,np.logical_not(np.all(mask,axis=0))]


def base_variance(rn_var, darkcurr=None, exptime=None, proc_var=None, count_scale=None):
    r"""
    Calculate the "base-level" variance in a processed image driven by the
    detector properties and the additive noise from the image processing steps.

    The full variance model (see :func:`variance_model`), :math:`V`, is:

    .. math::

        V = s^2\ \left[ {\rm max}(0, C) + D t_{\rm exp} / 3600 +
                V_{\rm rn} + V_{\rm proc} \right] 
                + \epsilon^2 {\rm max}(0, c)^2

    where:

        - :math:`c=s\ C` are the rescaled observed sky + object counts,
        - :math:`C` is the observed number of sky + object counts,
        - :math:`s=s\prime / N_{\rm frames}` is a scale factor derived
          from the (inverse of the) flat-field frames plus the number
          of frames contributing to the object counts (see ``count_scale``),
        - :math:`D` is the dark current in electrons per **hour** (see
          ``darkcurr``),
        - :math:`t_{\rm exp}` is the effective exposure time in seconds (see
          ``exptime``),
        - :math:`V_{\rm rn}` is the detector readnoise variance (i.e.,
          read-noise squared; see ``rn_var``),
        - :math:`V_{\rm proc}` is added variance from image processing (e.g.,
          bias subtraction; see ``proc_var``), and
        - :math:`\epsilon` is an added error term that imposes a maximum
          signal-to-noise on the observed counts.

    This function consolidates terms that do not change with the forward
    modeling of the sky + object counts.  That is, this function calculates

    .. math::

        V_{\rm base} = s^2\ \left[ D t_{\rm exp} / 3600 + V_{\rm rn} + V_{\rm
        proc} \right]

    such that the first equation can be re-written as

    .. math::

        V = s {\rm max}(0,c) + V_{\rm base} + \epsilon^2 {\rm max}(0, c)^2.

    .. warning::

        - If :math:`s` (``count_scale``) is provided, the variance will be 0
          wherever :math:`s \leq 0`.

        - Note that dark current is typically given in electrons per second *per
          pixel*.  If on-chip binning was used for the detector readout, each
          binned pixel will have accummulated the expected dark-current (in
          e-/s/pixel) multiplied by the number of binned pixels.  Beware the
          units of ``darkcurr``, both in that it is dark-current per *hour* and
          that it is the dark-current expected in the *binned* pixel.  For
          example, see the calling function
          :func:`pypeit.images.rawimage.RawImage.build_ivar`.
        
        - The input arrays can have any dimensionality (i.e., they can be single
          2D images or a 3D array containing multiple 2D images); however, the
          exposure time must be a scalar applied to all array values.

    Args:
        rn_var (`numpy.ndarray`_):
            A 2D array with the readnoise variance (i.e., readnoise squared)
            from the instrument detector; see :func:`rn2_frame`.  This should
            include digitization noise and any difference in the readnoise
            across the detector due to the use of multiple amplifiers.
            Readnoise should be in e-, meaning this is in elections squared.
        darkcurr (:obj:`float`, `numpy.ndarray`_, optional):
            Dark current in electrons per **hour** (as is the convention for the
            :class:`~pypeit.images.detector_container.DetectorContainer` object)
            if the exposure time is provided, otherwise in electrons.  Note that
            this is the dark-current in each read pixel, meaning you likely need
            to multiply the quoted detector dark-current by the number of pixels
            in a bin (e.g., 4 for 2x2 binning) for binned data.  If None, set to
            0.  If a single float, assumed to be constant across the full image.
            If an array, the shape must match ``rn_var``.
        exptime (:obj:`float`, optional):
            Exposure time in seconds.  If None, dark current *must* be
            in electrons.
        proc_var (:obj:`float`, `numpy.ndarray`_, optional):
            Additional variance terms to include that are due to the image
            processing steps (e.g., bias subtraction).  If None, set to 0.  If a
            single float, assumed to be constant across the full image.  If an
            array, the shape must match ``rn_var``.
        count_scale (:obj:`float`, `numpy.ndarray`_, optional):
            A scale factor that *has already been applied* to the provided
            counts. It accounts for the number of frames contributing to
            the provided counts, and the relative throughput factors that
            can be measured from flat-field frames. For example, if the image
            has been flat-field corrected, this is the inverse of the flat-field counts.
            If None, set to 1. If a single float, assumed to be constant across the full image.
            If an array, the shape must match ``rn_var``.  The variance will be 0
            wherever :math:`s \leq 0`, modulo the provided ``noise_floor``.

    Returns:
        `numpy.ndarray`_: Base-level variance image computed via the equation
        above with the same shape as ``rn_var``.
    """
    # Check input
    if count_scale is not None and isinstance(count_scale, np.ndarray) \
            and count_scale.shape != rn_var.shape:
        msgs.error('Count scale and readnoise variance have different shape.')
    if proc_var is not None and isinstance(proc_var, np.ndarray) \
            and proc_var.shape != rn_var.shape:
        msgs.error('Processing variance and readnoise variance have different shape.')
    if darkcurr is not None and isinstance(darkcurr, np.ndarray) \
            and darkcurr.shape != rn_var.shape:
        msgs.error('Dark image and readnoise variance have different shape.')

    # Build the variance
    #   - First term is the read-noise
    var = rn_var.copy()
    #   - Add the processing noise
    if proc_var is not None:
        var += proc_var
    #   - Add the dark current
    if darkcurr is not None:
        var += darkcurr if exptime is None else darkcurr * exptime / 3600
    #   - Include the rescaling
    if count_scale is not None:
        _count_scale = count_scale.copy() if isinstance(count_scale, np.ndarray) \
                        else np.full(var.shape, count_scale, dtype=float)
        var *= _count_scale**2
    # Done
    return var


def variance_model(base, counts=None, count_scale=None, noise_floor=None):
    r"""
    Calculate the expected variance in an image.

    The full variance model (see :func:`variance_model`), :math:`V`, is:

    .. math::

        V = s^2\ \left[ {\rm max}(0, C) + D t_{\rm exp} / 3600 +
                V_{\rm rn} + V_{\rm proc} \right] 
                + \epsilon^2 {\rm max}(0, c)^2

    where:

        - :math:`c=s\ C` are the rescaled observed sky + object counts (see
          ``counts``),
        - :math:`C` is the observed number of sky + object counts,
        - :math:`s=s\prime / N_{\rm frames}` is a scale factor derived
          from the (inverse of the) flat-field frames plus the number
          of frames contributing to the object counts (see ``count_scale``),
        - :math:`D` is the dark current in electrons per **hour**,
        - :math:`t_{\rm exp}` is the effective exposure time in seconds,
        - :math:`V_{\rm rn}` is the detector readnoise variance (i.e.,
          read-noise squared),
        - :math:`V_{\rm proc}` is added variance from image processing (e.g.,
          bias subtraction), and
        - :math:`\epsilon` is an added error term that imposes a maximum
          signal-to-noise on the observed counts (see ``noise_floor``).

    The function :func:`base_variance` consolidates all terms
    that do not change with the forward
    modeling of the sky + object counts into a single "base-level" variance

    .. math::

        V_{\rm base} = s^2\ \left[ D t_{\rm exp} / 3600 + V_{\rm rn} + V_{\rm
        proc} \right]

    such that the first equation can be re-written as

    .. math::

        V = s {\rm max}(0,c) + V_{\rm base} + \epsilon^2 {\rm max}(0, c)^2,

    which is the quantity returned by this function.

    We emphasize that this is a *model* for the per-pixel image variance.  In
    real data, the as-observed pixel values are used to estimate the Poisson
    error in the observed counts.  Because of the variance in the image, this
    systematically overestimates the variance toward low counts (:math:`\lesssim
    2 \sigma_{\rm rn}`), with a bias of approximately :math:`1.4/\sigma_{\rm
    rn}` for :math:`C=0` (i.e., about 20% for a readnoise of 2 e-) and less than
    10% for :math:`C=1`.

    .. warning::

        - If :math:`s` (``count_scale``) is provided, the variance will be 0
          wherever :math:`s \leq 0`, modulo the provided ``noise_floor``.

        - The input arrays can have any dimensionality (i.e., they can be single
          2D images or a 3D array containing multiple 2D images); however,
          currently the noise floor must be a single scalar applied to all array
          values.

    Args:
        base (`numpy.ndarray`_):
            The "base-level" variance in the data set by the detector properties
            and the image processing steps.  See :func:`base_variance`;
            :math:`V_{\rm base}` in the equations above.
        counts (`numpy.ndarray`_, optional):
            An array with the number of source-plus-sky counts, possibly
            rescaled by a relative throughput; see :math:`c` in the equations
            above.  Because this is used to calculate the noise floor, this
            *must* be provided if ``noise_floor`` is not None.  Shape must match
            ``base``.  Shape must match ``base``.
        count_scale (:obj:`float`, `numpy.ndarray`_, optional):
            A scale factor that *has already been applied* to the provided
            counts; see :math:`s` in the equations above.  It accounts for
            the number of frames contributing to  the provided counts, and
            the relative throughput factors that  can be measured from flat-field frames.
            For example, if the image has been flat-field corrected, this is the inverse
            of the flat-field counts.  If None, no scaling is expected, meaning
            ``counts`` are exactly the observed detector counts.  If a single
            float, assumed to be constant across the full image.  If an array,
            the shape must match ``base``.  The variance will be 0 wherever
            :math:`s \leq 0`, modulo the provided ``noise_floor``.
        noise_floor (:obj:`float`, optional):
            A fraction of the counts to add to the variance, which has the
            effect of ensuring that the S/N is never greater than
            ``1/noise_floor``; see :math:`epsilon` in the equations above.  If
            None, no noise floor is added.  If not None, ``counts`` *must* be
            provided.

    Returns:
        `numpy.ndarray`_: Variance image computed via the equation above with
        the same shape as ``base``.
    """
    # Check input
    if noise_floor is not None and noise_floor > 0. and counts is None:
        msgs.error('To impose a noise floor, must provide counts.')
    if counts is not None and counts.shape != base.shape:
        msgs.error('Counts image and base-level variance have different shape.')
    if count_scale is not None and isinstance(count_scale, np.ndarray) \
            and count_scale.shape != base.shape:
        msgs.error('Count scale and base-level variance have different shape.')

    # Clip the counts
    _counts = None if counts is None else np.clip(counts, 0, None)

    # Build the variance
    #   - Start with the base-level variance
    var = base.copy()
    #   - Add the sky + object counts
    if counts is not None:
        var += _counts if count_scale is None else count_scale * _counts
        #   - Add the noise floor
        if noise_floor is not None and noise_floor > 0.:
            var += (noise_floor * _counts)**2
    # Done
    return var


