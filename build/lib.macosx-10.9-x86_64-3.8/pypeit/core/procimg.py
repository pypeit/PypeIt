""" Module for image processing core methods

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import numpy as np
from scipy import signal, ndimage
from scipy.optimize import curve_fit
from IPython import embed

from pypeit import msgs
from pypeit import utils
from pypeit.core import parse


def lacosmic(sciframe, saturation, nonlinear, varframe=None, maxiter=1, grow=1.5,
             remove_compact_obj=True, sigclip=5.0, sigfrac=0.3, objlim=5.0):
    """
    Identify cosmic rays using the L.A.Cosmic algorithm
    U{http://www.astro.yale.edu/dokkum/lacosmic/}
    (article : U{http://arxiv.org/abs/astro-ph/0108003})
    This routine is mostly courtesy of Malte Tewes

    Args:
        sciframe:
        saturation:
        nonlinear:
        varframe:
        maxiter:
        grow:
        remove_compact_obj:
        sigclip (float):
            Threshold for identifying a CR
        sigfrac:
        objlim:

    Returns:
        ndarray: mask of cosmic rays (0=no CR, 1=CR)

    """
    msgs.info("Detecting cosmic rays with the L.A.Cosmic algorithm")
#    msgs.work("Include these parameters in the settings files to be adjusted by the user")
    # Set the settings
    scicopy = sciframe.copy()
    crmask = np.cast['bool'](np.zeros(sciframe.shape))
    sigcliplow = sigclip * sigfrac

    # Determine if there are saturated pixels
    satpix = np.zeros_like(sciframe)
#    satlev = settings_det['saturation']*settings_det['nonlinear']
    satlev = saturation*nonlinear
    wsat = np.where(sciframe >= satlev)
    if wsat[0].size == 0: satpix = None
    else:
        satpix[wsat] = 1.0
        satpix = np.cast['bool'](satpix)

    # Define the kernels
    laplkernel = np.array([[0.0, -1.0, 0.0], [-1.0, 4.0, -1.0], [0.0, -1.0, 0.0]])  # Laplacian kernal
    growkernel = np.ones((3,3))
    for i in range(1, maxiter+1):
        msgs.info("Convolving image with Laplacian kernel")
        # Subsample, convolve, clip negative values, and rebin to original size
        subsam = utils.subsample(scicopy)
        conved = signal.convolve2d(subsam, laplkernel, mode="same", boundary="symm")
        cliped = conved.clip(min=0.0)
        lplus = utils.rebin_evlist(cliped, np.array(cliped.shape)/2.0)

        msgs.info("Creating noise model")
        # Build a custom noise map, and compare  this to the laplacian
        m5 = ndimage.filters.median_filter(scicopy, size=5, mode='mirror')
        if varframe is None:
            noise = np.sqrt(np.abs(m5))
        else:
            noise = np.sqrt(varframe)
        msgs.info("Calculating Laplacian signal to noise ratio")

        # Laplacian S/N
        s = lplus / (2.0 * noise)  # Note that the 2.0 is from the 2x2 subsampling

        # Remove the large structures
        sp = s - ndimage.filters.median_filter(s, size=5, mode='mirror')

        msgs.info("Selecting candidate cosmic rays")
        # Candidate cosmic rays (this will include HII regions)
        candidates = sp > sigclip
        nbcandidates = np.sum(candidates)

        msgs.info("{0:5d} candidate pixels".format(nbcandidates))

        # At this stage we use the saturated stars to mask the candidates, if available :
        if satpix is not None:
            msgs.info("Masking saturated pixels")
            candidates = np.logical_and(np.logical_not(satpix), candidates)
            nbcandidates = np.sum(candidates)

            msgs.info("{0:5d} candidate pixels not part of saturated stars".format(nbcandidates))

        msgs.info("Building fine structure image")

        # We build the fine structure image :
        m3 = ndimage.filters.median_filter(scicopy, size=3, mode='mirror')
        m37 = ndimage.filters.median_filter(m3, size=7, mode='mirror')
        f = m3 - m37
        f /= noise
        f = f.clip(min=0.01)

        msgs.info("Removing suspected compact bright objects")

        # Now we have our better selection of cosmics :

        if remove_compact_obj:
            cosmics = np.logical_and(candidates, sp/f > objlim)
        else:
            cosmics = candidates
        nbcosmics = np.sum(cosmics)

        msgs.info("{0:5d} remaining candidate pixels".format(nbcosmics))

        # What follows is a special treatment for neighbors, with more relaxed constains.

        msgs.info("Finding neighboring pixels affected by cosmic rays")

        # We grow these cosmics a first time to determine the immediate neighborhod  :
        growcosmics = np.cast['bool'](signal.convolve2d(np.cast['float32'](cosmics), growkernel, mode="same", boundary="symm"))

        # From this grown set, we keep those that have sp > sigmalim
        # so obviously not requiring sp/f > objlim, otherwise it would be pointless
        growcosmics = np.logical_and(sp > sigclip, growcosmics)

        # Now we repeat this procedure, but lower the detection limit to sigmalimlow :

        finalsel = np.cast['bool'](signal.convolve2d(np.cast['float32'](growcosmics), growkernel, mode="same", boundary="symm"))
        finalsel = np.logical_and(sp > sigcliplow, finalsel)

        # Unmask saturated pixels:
        if satpix is not None:
            msgs.info("Masking saturated stars")
            finalsel = np.logical_and(np.logical_not(satpix), finalsel)

        ncrp = np.sum(finalsel)

        msgs.info("{0:5d} pixels detected as cosmics".format(ncrp))

        # We find how many cosmics are not yet known :
        newmask = np.logical_and(np.logical_not(crmask), finalsel)
        nnew = np.sum(newmask)

        # We update the mask with the cosmics we have found :
        crmask = np.logical_or(crmask, finalsel)

        msgs.info("Iteration {0:d} -- {1:d} pixels identified as cosmic rays ({2:d} new)".format(i, ncrp, nnew))
        if ncrp == 0: break
    # Additional algorithms (not traditionally implemented by LA cosmic) to remove some false positives.
    msgs.work("The following algorithm would be better on the rectified, tilts-corrected image")
    filt  = ndimage.sobel(sciframe, axis=1, mode='constant')
    filty = ndimage.sobel(filt/np.sqrt(np.abs(sciframe)), axis=0, mode='constant')
    filty[np.where(np.isnan(filty))]=0.0

    sigimg = cr_screen(filty)

    sigsmth = ndimage.filters.gaussian_filter(sigimg,1.5)
    sigsmth[np.where(np.isnan(sigsmth))]=0.0
    sigmask = np.cast['bool'](np.zeros(sciframe.shape))
    sigmask[np.where(sigsmth>sigclip)] = True
    crmask = np.logical_and(crmask, sigmask)
    msgs.info("Growing cosmic ray mask by 1 pixel")
    crmask = grow_masked(crmask.astype(np.float), grow, 1.0)

    return crmask.astype(bool)


def cr_screen(a, mask_value=0.0, spatial_axis=1):
    r"""
    Calculate the significance of pixel deviations from the median along
    the spatial direction.

    No type checking is performed of the input array; however, the
    function assumes floating point values.

    Args:
        a (numpy.ndarray): Input 2D array
        mask_value (float): (**Optional**) Values to ignore during the
            calculation of the median.  Default is 0.0.
        spatial_axis (int): (**Optional**) Axis along which to calculate
            the median.  Default is 1.

    Returns:
        numpy.ndarray: Returns a map of :math:`|\Delta_{i,j}|/\sigma_j`,
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


def grow_masked(img, grow, growval):

    if not np.any(img == growval):
        return img

    _img = img.copy()
    sz_x, sz_y = img.shape
    d = int(1+grow)
    rsqr = grow*grow

    # Grow any masked values by the specified amount
    for x in range(sz_x):
        for y in range(sz_y):
            if img[x,y] != growval:
                continue

            mnx = 0 if x-d < 0 else x-d
            mxx = x+d+1 if x+d+1 < sz_x else sz_x
            mny = 0 if y-d < 0 else y-d
            mxy = y+d+1 if y+d+1 < sz_y else sz_y

            for i in range(mnx,mxx):
                for j in range(mny, mxy):
                    if (i-x)*(i-x)+(j-y)*(j-y) <= rsqr:
                        _img[i,j] = growval
    return _img


def gain_frame(amp_img, gain):
    """
    Generate an image with the gain for each pixel.

    Args:
        amp_img (`numpy.ndarray`_):
            Integer array that identifies which (1-indexed) amplifier
            was used to read each pixel.
        gain (:obj:`list`):
            List of amplifier gain values.  Must be that the gain for
            amplifier 1 is provided by `gain[0]`, etc.

    Returns:
        `numpy.ndarray`_: Image with the gain for each pixel.
    """
    # TODO: Remove this or actually do it.
    # msgs.warn("Should probably be measuring the gain across the amplifier boundary")

    # Build the gain image
    gain_img = np.zeros_like(amp_img, dtype=float)
    for i,_gain in enumerate(gain):
        gain_img[amp_img == i+1] = _gain

    # Return the image, trimming if requested
    return gain_img


def rn_frame(datasec_img, gain, ronoise):
    """ Generate a RN image

    Parameters
    ----------
    datasec_img : ndarray
        Index image (2D array) assigning each pixel in the data section to an amplifier.
        A value of 0 means not a data section
    gain : ndarray, list
        A list of the gains for each amplifier
    ronoise : ndarray, list
        A list of read noise values for each amplifier. If any element of the array is 0.0,
        the read noise will be determined from the overscan region

    Returns
    -------
    rn_img : ndarray
        Read noise *variance* image (i.e. RN**2)
    """
    # Determine the number of amplifiers from the datasec image
    numamplifiers = np.amax(datasec_img)

    # Check the input types
    _gain = np.asarray(gain) if isinstance(gain, (list, np.ndarray)) else np.array([gain])
    _ronoise = np.asarray(ronoise) if isinstance(ronoise, (list, np.ndarray)) \
                        else np.array([ronoise])
    if len(_gain) != numamplifiers:
        raise ValueError('Must provide a gain for each amplifier.')
    if len(_ronoise) != numamplifiers:
        raise ValueError('Must provide a read-noise for each amplifier.')
    if np.any(datasec_img > numamplifiers):
        raise ValueError('Pixel amplifier IDs do not match number of amplifiers.')

    # ToDO We should not be using numpy masked arrays!!!

    # Get the amplifier indices
    indx = datasec_img.astype(int) == 0
    amp = np.ma.MaskedArray(datasec_img.astype(int) - 1, mask=indx).filled(0)

    # Return the read-noise image.  Any pixels without an assigned
    # amplifier are given a noise of 0.
    return np.ma.MaskedArray(np.square(_ronoise[amp]) + np.square(0.5*_gain[amp]),
                             mask=indx).filled(0.0)


def rect_slice_with_mask(image, mask, mask_val=1):
    """
    Generate rectangular slices from a mask image

    Args:
        image (np.ndarray): Image to mask
        mask (np.ndarray): Mask image
        mask_val (int,optiona): Value to mask on

    Returns:
        :obj:`tuple`: The image at mask values and a 2-tuple with the
        :obj:`slice` objects that select the masked data.
    """
    pix = np.where(mask == mask_val)
    slices = (slice(np.min(pix[0]), np.max(pix[0])+1), slice(np.min(pix[1]), np.max(pix[1])+1))
    return image[slices], slices


def subtract_overscan(rawframe, datasec_img, oscansec_img,
                          method='savgol', params=[5, 65]):
    """
    Subtract overscan

    Args:
        rawframe (:obj:`numpy.ndarray`):
            Frame from which to subtract overscan
        numamplifiers (int):
            Number of amplifiers for this detector.
        datasec_img (:obj:`numpy.ndarray`):
            An array the same shape as rawframe that identifies
            the pixels associated with the data on each amplifier.
            0 for not data, 1 for amplifier 1, 2 for amplifier 2, etc.
        oscansec_img (:obj:`numpy.ndarray`):
            An array the same shape as rawframe that identifies
            the pixels associated with the overscan region on each
            amplifier.
            0 for not data, 1 for amplifier 1, 2 for amplifier 2, etc.
        method (:obj:`str`, optional):
            The method used to fit the overscan region.  Options are
            polynomial, savgol, median.
        params (:obj:`list`, optional):
            Parameters for the overscan subtraction.  For
            method=polynomial, set params = order, number of pixels,
            number of repeats ; for method=savgol, set params = order,
            window size ; for method=median, params are ignored.

    Returns:
        :obj:`numpy.ndarray`: The input frame with the overscan region
        subtracted
    """
    # Copy the data so that the subtraction is not done in place
    no_overscan = rawframe.copy()

    # Amplifiers
    amps = np.unique(datasec_img[datasec_img > 0]).tolist()

    # Perform the overscan subtraction for each amplifier
    for amp in amps:
        # Pull out the overscan data
        overscan, _ = rect_slice_with_mask(rawframe, oscansec_img, amp)
        # Pull out the real data
        data, data_slice = rect_slice_with_mask(rawframe, datasec_img, amp)

        # Shape along at least one axis must match
        data_shape = data.shape
        if not np.any([dd == do for dd, do in zip(data_shape, overscan.shape)]):
            msgs.error('Overscan sections do not match amplifier sections for'
                       'amplifier {0}'.format(amp))
        compress_axis = 1 if data_shape[0] == overscan.shape[0] else 0

        # Fit/Model the overscan region
        osfit = np.median(overscan) if method.lower() == 'median' \
            else np.median(overscan, axis=compress_axis)
        if method.lower() == 'polynomial':
            # TODO: Use np.polynomial.polynomial.polyfit instead?
            c = np.polyfit(np.arange(osfit.size), osfit, params[0])
            ossub = np.polyval(c, np.arange(osfit.size))
        elif method.lower() == 'savgol':
            ossub = signal.savgol_filter(osfit, params[1], params[0])
        elif method.lower() == 'median':
            # Subtract scalar and continue
            no_overscan[data_slice] -= osfit
            continue
        else:
            raise ValueError('Unrecognized overscan subtraction method: {0}'.format(method))

        # Subtract along the appropriate axis
        no_overscan[data_slice] -= (ossub[:, None] if compress_axis == 1 else ossub[None, :])

    return no_overscan


def subtract_pattern(rawframe, datasec_img, oscansec_img, frequency=None, axis=1, debug=False):
    """
    Subtract a sinusoidal pattern from the input rawframe. The algorithm
    calculates the frequency of the signal, generates a model, and subtracts
    this signal from the data. This sinusoidal pattern noise was first
    identified in KCWI, but the source of this pattern noise is not
    currently known.

    Args:
        rawframe (`numpy.ndarray`_):
            Frame from which to subtract overscan
        numamplifiers (int):
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
        frequency (float, list, optional):
            The frequency (or list of frequencies - one for each amplifier)
            of the sinusoidal pattern. If None, the frequency of each amplifier
            will be determined from the overscan region.
        axis (int):
            Which axis should the pattern subtraction be applied?
        debug (bool):
            Debug the code (True means yes)

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

        # Convert frequency to the size of the overscan region
        msgs.info("Subtracting detector pattern with frequency = {0:f}".format(use_fr))
        use_fr *= (overscan.shape[1]-1)

        # Get a first guess of the amplitude and phase information
        amp = np.fft.rfft(overscan, axis=1)
        idx = (np.arange(overscan.shape[0]), np.argmax(np.abs(amp), axis=1))
        # Convert result to amplitude and phase
        amps = (np.abs(amp))[idx] * (2.0 / overscan.shape[1])
        phss = np.arctan2(amp.imag, amp.real)[idx]

        # Use the above to as initial guess parameters in chi-squared minimisation
        cosfunc = lambda xarr, *p: p[0] * np.cos(2.0 * np.pi * p[1] * xarr + p[2])
        xdata, step = np.linspace(0.0, 1.0, overscan.shape[1], retstep=True)
        xdata_all = (np.arange(osd_slice[1].start, osd_slice[1].stop) - os_slice[1].start) * step
        model_pattern = np.zeros_like(oscandata)
        val = np.zeros(overscan.shape[0])
        # Get the best estimate of the amplitude
        for ii in range(overscan.shape[0]):
            try:
                popt, pcov = curve_fit(cosfunc, xdata, overscan[ii, :], p0=[amps[ii], use_fr, phss[ii]],
                                       bounds=([-np.inf, use_fr * 0.99999999, -np.inf], [+np.inf, use_fr * 1.00000001, +np.inf]))
            except ValueError:
                msgs.warn("Input data invalid for pattern subtraction of row {0:d}/{1:d}".format(ii + 1, overscan.shape[0]))
                continue
            except RuntimeError:
                msgs.warn("Pattern subtraction fit failed for row {0:d}/{1:d}".format(ii + 1, overscan.shape[0]))
                continue
            val[ii] = popt[0]
            model_pattern[ii, :] = cosfunc(xdata_all, *popt)
        use_amp = np.median(val)
        # Get the best estimate of the phase, and generate a model
        for ii in range(overscan.shape[0]):
            try:
                popt, pcov = curve_fit(cosfunc, xdata, overscan[ii, :], p0=[use_amp, use_fr, phss[ii]],
                                       bounds=([use_amp * 0.99999999, use_fr * 0.99999999, -np.inf],
                                               [use_amp * 1.00000001, use_fr * 1.00000001, +np.inf]))
            except ValueError:
                msgs.warn("Input data invalid for pattern subtraction of row {0:d}/{1:d}".format(ii + 1, overscan.shape[0]))
                continue
            except RuntimeError:
                msgs.warn("Pattern subtraction fit failed for row {0:d}/{1:d}".format(ii + 1, overscan.shape[0]))
                continue
            model_pattern[ii, :] = cosfunc(xdata_all, *popt)
        outframe[osd_slice] -= model_pattern

    debug = False
    if debug:
        embed()
        import astropy.io.fits as fits
        hdu = fits.PrimaryHDU(rawframe)
        hdu.writeto("tst_raw.fits", overwrite=True)
        hdu = fits.PrimaryHDU(outframe)
        hdu.writeto("tst_sub.fits", overwrite=True)
        hdu = fits.PrimaryHDU(rawframe - outframe)
        hdu.writeto("tst_mod.fits", overwrite=True)

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
        frame (:obj:`numpy.ndarray`):
            2D array to measure the pattern frequency
        axis (int, optional):
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
    # Find significant deviations and ignore those rows
    mad = 1.4826*np.median(np.abs(arr))
    ww = np.where(arr > 10*mad)
    # Create a mask of these rows
    msk = np.sort(np.unique(ww[0]))

    # Compute the Fourier transform to obtain an estimate of the dominant frequency component
    amp = np.fft.rfft(arr, axis=1)
    idx = (np.arange(arr.shape[0]), np.argmax(np.abs(amp), axis=1))

    # Construct the variables of the sinusoidal waveform
    amps = (np.abs(amp))[idx] * (2.0 / arr.shape[1])
    phss = np.arctan2(amp.imag, amp.real)[idx]
    frqs = idx[1]

    # Use the above to as initial guess parameters in chi-squared minimisation
    cosfunc = lambda xarr, *p: p[0] * np.cos(2.0 * np.pi * p[1] * xarr + p[2])
    xdata = np.linspace(0.0, 1.0, arr.shape[1])
    # Calculate the amplitude distribution
    amp_dist = np.zeros(arr.shape[0])
    frq_dist = np.zeros(arr.shape[0])
    # Loop over all rows to new independent values that can be averaged
    for ii in range(arr.shape[0]):
        if ii in msk:
            continue
        try:
            popt, pcov = curve_fit(cosfunc, xdata, arr[ii, :], p0=[amps[ii], frqs[ii], phss[ii]],
                                   bounds=([-np.inf, frqs[ii]-1, -np.inf],
                                           [+np.inf, frqs[ii]+1, +np.inf]))
        except ValueError:
            msgs.warn("Input data invalid for pattern frequency fit of row {0:d}/{1:d}".format(ii+1, arr.shape[0]))
            continue
        except RuntimeError:
            msgs.warn("Pattern frequency fit failed for row {0:d}/{1:d}".format(ii+1, arr.shape[0]))
            continue
        amp_dist[ii] = popt[0]
        frq_dist[ii] = popt[1]
    ww = np.where(amp_dist > 0.0)
    use_amp = np.median(amp_dist[ww])
    use_frq = np.median(frq_dist[ww])
    # Calculate the frequency distribution with a prior on the amplitude
    frq_dist = np.zeros(arr.shape[0])
    for ii in range(arr.shape[0]):
        if ii in msk:
            continue
        try:
            popt, pcov = curve_fit(cosfunc, xdata, arr[ii, :], p0=[use_amp, use_frq, phss[ii]],
                                   bounds=([use_amp * 0.99999999, use_frq-1, -np.inf],
                                           [use_amp * 1.00000001, use_frq+1, +np.inf]))
        except ValueError:
            msgs.warn("Input data invalid for patern frequency fit of row {0:d}/{1:d}".format(ii+1, arr.shape[0]))
            continue
        except RuntimeError:
            msgs.warn("Pattern frequency fit failed for row {0:d}/{1:d}".format(ii+1, arr.shape[0]))
            continue
        frq_dist[ii] = popt[1]
    # Ignore masked values, and return the best estimate of the frequency
    ww = np.where(frq_dist > 0.0)
    medfrq = np.median(frq_dist[ww])
    return medfrq/(arr.shape[1]-1)


'''
def subtract_overscan(rawframe, numamplifiers, datasec, oscansec, method='savgol', params=[5,65]):
    """
    Subtract overscan

    TODO: Describe the method.

    Args:
        rawframe (:obj:`numpy.ndarray`):
            Frame from which to subtract overscan
        numamplifiers (int):
            Number of amplifiers for this detector.
        datasec (list):
            List of tuples, one per amplifier, with the slice along each
            dimension that selects the data section of the provided raw
            frame.  See :func:`parse.sec2slice` to convert a string
            section (as read from a file header) into a list of slices.
        oscansec (list):
            List of tuples, one per amplifier, with the slice along each
            dimension that selects the overscan section of the provided
            raw frame.  See :func:`parse.sec2slice` to convert a string
            section (as read from a file header) into a list of slices.
        method (:obj:`str`, optional):
            The method used to fit the overscan region.  Options are
            polynomial, savgol, median.
        params (:obj:`list`, optional):
            Parameters for the overscan subtraction.  For
            method=polynomial, set params = order, number of pixels,
            number of repeats ; for method=savgol, set params = order,
            window size ; for method=median, params are ignored.

    Returns:
        :obj:`numpy.ndarray`: The input frame with the overscan region
        subtracted
    """
    # Check input
    if len(datasec) != numamplifiers or len(oscansec) != numamplifiers:
        msgs.error('Number of amplifiers does not match provided image sections.')

    # If the input image sections are strings, convert them
    if not isinstance(datasec[0], tuple):
        msgs.error('Data section must be a tuple of slice objects.')
    if not isinstance(oscansec[0], tuple):
        msgs.error('Overscan section must be a tuple of slice objects.')

    # Check that there are no overlapping data sections
    testframe = np.zeros_like(rawframe, dtype=int)
    for i in range(numamplifiers):
        testframe[datasec[i]] += 1
    if np.any(testframe > 1):
        raise ValueError('Image has overlapping data sections!')

    # Copy the data so that the subtraction is not done in place
    nobias = rawframe.copy()

    # Perform the bias subtraction for each amplifier
    for i in range(numamplifiers):
        # Pull out the overscan data
        overscan = rawframe[oscansec[i]]

        # Shape along at least one axis must match
        data_shape = rawframe[datasec[i]].shape
        if not np.any([ dd == do for dd, do in zip(data_shape, overscan.shape)]):
            msgs.error('Overscan sections do not match amplifier sections for'
                       'amplifier {0}'.format(i+1))
        compress_axis = 1 if data_shape[0] == overscan.shape[0] else 0

        # Fit/Model the overscan region
        osfit = np.median(overscan) if method.lower() == 'median' \
                        else np.median(overscan, axis=compress_axis)
        if method.lower() == 'polynomial':
            # TODO: Use np.polynomial.polynomial.polyfit instead?
            c = np.polyfit(np.arange(osfit.size), osfit, params[0])
            ossub = np.polyval(c, np.arange(osfit.size))
        elif method.lower() == 'savgol':
            ossub = signal.savgol_filter(osfit, params[1], params[0])
        elif method.lower() == 'median':
            # Subtract scalar and continue
            nobias[datasec[i]] -= osfit
            continue
        else:
            raise ValueError('Unrecognized overscan subtraction method: {0}'.format(method))

        # Subtract along the appropriate axis
        nobias[datasec[i]] -= (ossub[:,None] if compress_axis == 1 else ossub[None,:])

    return nobias
'''


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


def old_replace_columns(img, bad_cols, replace_with='mean'):
    """ Replace bad columns with values from the neighbors

    Parameters
    ----------
    img : ndarray
    bad_cols: ndarray (bool, 1D, shape[1] of img)
      True = bad column
      False = ok column
    replace_with : str, optional
      Option for replacement
       mean -- Use the mean of the closest left/right columns

    Returns
    -------
    img2 : ndarray
      Copy of the input image with the bad columns replaced
    """
    # Prep
    img2 = img.copy()
    # Find the starting/ends of the bad column sets
    tmp = np.zeros(img.shape[1], dtype=int)
    tmp[bad_cols] = 1
    tmp2 = tmp - np.roll(tmp,1)
    # Deal with first column
    if bad_cols[0]:
        tmp2[0]=1
    # Deal with last column
    if bad_cols[-1]:
        tmp2[-1]=-1
    ledges = np.where(tmp2 == 1)[0]
    redges = np.where(tmp2 == -1)[0]
    # Last column?
    if tmp2[-1] == 1:
        redges = np.concatenate([redges, np.array([bad_cols.size-1])])
    # Loop on em
    for kk, ledge in enumerate(ledges):
        lval = img[:,redges[kk]+1] if ledge == 0 else img[:,ledge-1]
        rval = img[:, redges[kk]]
        # First columns?
        # Replace
        if replace_with == 'mean':
            mval = (lval+rval)/2.
            for ii in range(ledge, redges[kk]+1):
                img2[:,ii] = mval
        else:
            msgs.error("Bad option to replace_columns")
    # Return
    return img2


def trim_frame(frame, mask):
    """
    Trim the masked regions from a frame.

    Args:
        frame (:obj:`numpy.ndarray`):
            Image to be trimmed
        mask (:obj:`numpy.ndarray`):
            Boolean image set to True for values that should be trimmed
            and False for values to be returned in the output trimmed
            image.

    Return:
        :obj:`numpy.ndarray`:
            Trimmed image

    Raises:
        PypitError:
            Error raised if the trimmed image includes masked values
            because the shape of the valid region is odd.
    """
    # TODO: Should check for this failure mode earlier
    if np.any(mask[np.invert(np.all(mask,axis=1)),:][:,np.invert(np.all(mask,axis=0))]):
        msgs.error('Data section is oddly shaped.  Trimming does not exclude all '
                   'pixels outside the data sections.')
    return frame[np.invert(np.all(mask,axis=1)),:][:,np.invert(np.all(mask,axis=0))]


def variance_frame(datasec_img, sciframe, gain, ronoise, darkcurr=None,
                   exptime=None, skyframe=None, objframe=None, adderr=0.01, rnoise=None):
    """
    Calculate the variance image including detector noise.

    .. todo::
        This needs particular attention because exptime and darkcurr;
        used to be dnoise

    Args:
        datasec_img (:obj:`numpy.ndarray`):
            Image that identifies which amplifier (1-indexed) was used
            to read each pixel.  Anything less than 1 is ignored.
        sciframe (:obj:`numpy.ndarray`):
            Science frame with counts in ?
        gain (:obj:`float`, array-like):
            Gain for each amplifier
        ronoise (:obj:`float`, array-like):
            Read-noise for each amplifier
        darkcurr (:obj:`float`, optional):
            Dark current in electrons per second if the exposure time is
            provided, otherwise in electrons.  If None, set to 0.
        exptime (:obj:`float`, optional):
            Exposure time in **hours**.  If None, darkcurrent *must* be
            in electrons.
        skyframe (`numpy.ndarray`_, optional):
            Sky image.
        objframe (`numpy.ndarray`_, optional):
            Model of object counts
        adderr (:obj:`float`, optional):
            Error floor. The quantity adderr**2*sciframe**2 is added in
            qudarature to the variance to ensure that the S/N is never >
            1/adderr, effectively setting a floor on the noise or a
            ceiling on the S/N.
        rnoise (:obj:`numpy.ndarray`, optional):
            Read noise image.  If not provided, it will be generated

    Returns:
        `numpy.ndarray`_: Variance image
    """

    # ToDO JFH: I would just add the darkcurrent here into the effective read noise image
    # The effective read noise (variance image)
    if rnoise is None:
        rnoise = rn_frame(datasec_img, gain, ronoise)

    # No sky frame provided
    if skyframe is None:
        _darkcurr = 0 if darkcurr is None else darkcurr
        if exptime is not None:
            _darkcurr *= exptime/3600.
        var = np.abs(sciframe - np.sqrt(2.0)*np.sqrt(rnoise)) + rnoise + _darkcurr
        var = var + adderr**2*(np.abs(sciframe))**2
        return var

    # TODO: There's some complicated logic here.  Why is objframe
    # needed?  Can't a users just use objframe in place of sciframe and
    # get the same behavior?  Why is darkcurr (what was dnoise) used
    # with sciframe and not objframe?

    # ToDO JFH: shouldn't dark current be added here as well??
    _objframe = np.zeros_like(skyframe) if objframe is None else objframe
    var = np.abs(skyframe + _objframe - np.sqrt(2.0)*np.sqrt(rnoise)) + rnoise
    var = var + adderr ** 2 * (np.abs(sciframe)) ** 2
    embed(header='this appears to be broken!')
    return




