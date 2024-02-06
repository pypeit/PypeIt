""" Module for sky subtraction

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import numpy as np

from scipy.optimize import least_squares
from scipy import signal, interpolate, ndimage
from IPython import embed

from pypeit import msgs, utils


def pad_frame(_frame, detpad=300):
    """
    Clean the edges of the input frame and then pad the frame to avoid edge effects.

    Parameters
    ----------
    _frame : `numpy.ndarray`_
        Frame to be padded
    detpad : int
        Number of pixels to pad the frame on each side

    Returns
    -------
    _frame_pad : `numpy.ndarray`_
        Padded frame
    """
    _frame[0, :] = np.median(_frame[0:10, :], axis=0)
    _frame[-1, :] = np.median(_frame[-10:, :], axis=0)
    return np.pad(_frame, detpad, mode='edge')  # Model should be generated on padded data


def scattered_light_model_pad(param, img, detpad=300):
    """
    Construct a scattered light model for the input image, with the model parameters
    defined by param. This function is used to generate a model of the scattered light,
    based on a set of model parameters that have first been optimized using scattered_light().
    The model is generated on a padded version of the input image, and then trimmed to
    match the input image size.

    Parameters
    ----------
    param : `numpy.ndarray`_
        Model parameters that determine the scattered light based on the input img.
        Here is a list of the individual parameter meanings:

        * param[0] = Gaussian kernel width in the spectral direction
        * param[1] = Gaussian kernel width in the spatial direction
        * param[2] = Lorentzian kernel width in the spectral direction
        * param[3] = Lorentzian kernel width in the spatial direction
        * param[4] = Pixel shift of the scattered light in the spectral direction
        * param[5] = Pixel shift of the scattered light in the spatial direction
        * param[6] = Zoom factor of the scattered light (~1)
        * param[7] = constant offset for scattered light (independent of img)
        * param[8] = Kernel angle
        * param[9] = Relative importance of Gaussian vs Lorentzian.
                        0 < value < 1 means Lorentzian is weighted more
                        value > 1 means Gaussian is weighted more.
        * param[10:] = Polynomial scaling coefficients
    img : `numpy.ndarray`_
        Image used to generate the scattered light model
    detpad : int
        Number of pixels to pad the frame on each side

    Returns
    -------
    model : `numpy.ndarray`_
        Model of the scattered light for the input
    """
    _frame_pad = pad_frame(img, detpad=detpad)
    return scattered_light_model(param, _frame_pad)[detpad:-detpad, detpad:-detpad]


def scattered_light_model(param, img):
    """ Model used to calculate the scattered light.

    The current model to generate the scattered light is a shifted, scaled, and blurred version of the
    data recorded on the detector. We allow for a shift in the spatial and spectral direction, and the
    blurring kernel can have different widths in the spectral and spatial directions. Finally, the model
    includes a polynomial correction in the spectral direction to account for a change in the spectral
    efficiency of the scattered light relative to the detector image.

    Parameters
    ----------
    param : `numpy.ndarray`_
        Model parameters that determine the scattered light based on the input img.
        Here is a list of the individual parameter meanings:

        * param[0] = Gaussian kernel width in the spectral direction
        * param[1] = Gaussian kernel width in the spatial direction
        * param[2] = Lorentzian kernel width in the spectral direction
        * param[3] = Lorentzian kernel width in the spatial direction
        * param[4] = Pixel shift of the scattered light in the spectral direction
        * param[5] = Pixel shift of the scattered light in the spatial direction
        * param[6] = Zoom factor of the scattered light (~1)
        * param[7] = constant offset for scattered light (independent of img)
        * param[8] = Kernel angle
        * param[9] = Relative importance of Gaussian vs Lorentzian.
                        0 < value < 1 means Lorentzian is weighted more
                        value > 1 means Gaussian is weighted more.
        * param[10:] = Polynomial scaling coefficients
    img : `numpy.ndarray`_
        Raw image that you want to compute the scattered light model.
        shape is (nspec, nspat)
        Model used to calculate the scattered light. This function is used to
        generate a model of the scattered light, based on a set of model parameters
        that have been optimized using self.scattered_light().

    Returns
    -------
    model : `numpy.ndarray`_
        Model of the scattered light for the input
    """
    # Extract the parameters into more conveniently named variables
    sigmx_g, sigmy_g, sigmx_l, sigmy_l = param[0], param[1], param[2], param[3]
    shft_spec, shft_spat, zoom_spec, zoom_spat = param[4], param[5], param[6], param[7]
    constant, kern_angle, kern_scale = param[8], param[9], param[10]
    polyterms_spat = param[11:13]
    polyterms_spec = param[13:]

    # Make a grid of coordinates
    specvec, spatvec = np.arange(img.shape[0]), np.arange(img.shape[1])
    spat, spec = np.meshgrid(spatvec/(spatvec.size-1), specvec/(specvec.size - 1))
    # Generate the polynomial efficiency scaling in the spatial direction
    polyscale = spat*(polyterms_spat[0] + polyterms_spat[1]*spec)  # linear term and a spectral cross-term
    # Now include the spectral direction
    for pp in range(polyterms_spec.size):
        polyscale += polyterms_spec[pp] * spec**pp

    # Generate a 2D smoothing kernel, composed of a 2D Gaussian and a 2D Lorentzian
    sigmx, sigmy = max(sigmx_g, sigmx_l), max(sigmy_g, sigmy_l),
    xkern, ykern = np.meshgrid(np.arange(int(10 * sigmx)) - 5 * sigmx,
                               np.arange(int(10 * sigmy)) - 5 * sigmy)
    # Rotate the kernel
    xkernrot = (xkern * np.cos(kern_angle) - ykern * np.sin(kern_angle))
    ykernrot = (xkern * np.sin(kern_angle) + ykern * np.cos(kern_angle))
    # Create and normalise the Gaussian kernel
    kernel_gaussian = np.exp(-((xkernrot/sigmx_g)**2 + (ykernrot/sigmy_g)**2))
    kernel_gaussian /= np.sum(kernel_gaussian)
    # Create and normalise the Lorenztian kernel
    kernel_lorentzian = 1 / ((xkernrot/sigmx_l) ** 2 + (ykernrot/sigmy_l) ** 2 + 1)
    kernel_lorentzian /= np.sum(kernel_lorentzian)
    # Add the individual kernels into a single kernel. Arbitrarily scale the Gaussian kernel (either is fine).
    # The point of this is to make it so either a lorentzian, gaussian, or something in-between can be
    # used as the kernel, making it more flexible for different spectrographs.
    kernel = kernel_lorentzian + kern_scale * kernel_gaussian
    kernel /= np.sum(kernel)

    scale_img = polyscale * signal.oaconvolve(img, kernel, mode='same')
    spl = interpolate.RectBivariateSpline(specvec, spatvec, scale_img, kx=1, ky=1)
    return constant + spl(zoom_spec * (specvec + shft_spec), zoom_spat * (spatvec + shft_spat))


def scattlight_resid(param, wpix, img):
    """ Residual function used to optimize the model parameters

    Parameters
    ----------
    param : `numpy.ndarray`_
        1D array of model parameters to use for the fitting function.
    wpix : tuple
        A tuple containing the x,y coordinates of the pixels in img
        to be used for computing the residual.
    img : `numpy.ndarray`_
        Data image to be used to compute the residual. Shape is (nspec, nspat)

    Returns
    -------
    resid : `numpy.ndarray`_
        A 1D vector of the residuals
    """
    model = scattered_light_model(param, img)
    return img[wpix] - model[wpix]


def scattered_light(frame, bpm, offslitmask, x0, bounds, detpad=300, debug=False):
    """ Calculate a model of the scattered light of the input frame.

    Parameters
    ----------
    frame : `numpy.ndarray`_
        Raw 2D data frame (nspec, nspat) to be used to compute the scattered light.
    bpm : `numpy.ndarray`_
        2D boolean array indicating the bad pixels (True=bad), same shape as frame
    offslitmask : `numpy.ndarray`_
        A boolean mask indicating the pixels that are on/off the slit (True = off the slit), same shape as frame
    x0 : `numpy.ndarray`_
        A 1D array containing the best-fitting model parameters
    bounds : :obj:`tuple`
        A tuple of two elements, containing two `numpy.ndarray`_ of the same length as x0. These
        two arrays contain the lower (first element of the tuple) and upper (second element of the tuple)
        bounds to consider on the scattered light model parameters.
    debug : :obj:`bool`, optional
        If True, debug the final model fit that's been output

    Returns
    -------
    scatt_img : `numpy.ndarray`_
        A 2D image of the scattered light determined from the input frame.
        Alternatively, if a constant value is used, a constant floating point
        value can be returned as well.
    modelpar : `numpy.ndarray`_
        A 1D array containing the best-fitting model parameters
    success : :obj:`bool`
        True if the fit was successful, False otherwise
    """
    # Convert the BPM to a GPM for convenience
    gpm = np.logical_not(bpm)

    # Replace bad pixels with the nearest good pixel
    _frame = utils.replace_bad(frame, bpm)

    # Pad the edges of the data
    _frame_pad = pad_frame(_frame, detpad)
    offslitmask_pad = np.pad(offslitmask * gpm, detpad, mode='constant', constant_values=0)  # but don't include padded data in the fit
    # Grab the pixels to be included in the fit
    wpix = np.where(offslitmask_pad)

    # Compute the best-fitting model parameters
    msgs.info("Performing a least-squares fit to the scattered light")
    res_lsq = least_squares(scattlight_resid, x0, bounds=bounds, args=(wpix, _frame_pad),
                            verbose=2, ftol=1.0E-4)

    # Store if this is a successful fit
    success = res_lsq.success
    if success:
        msgs.info("Generating best-fitting scattered light model")
        scatt_img = scattered_light_model(res_lsq.x, _frame_pad)[detpad:-detpad, detpad:-detpad]
    else:
        msgs.warn("Scattered light model fitting failed")
        scatt_img = np.zeros_like(frame)

    if debug:
        # Do some checks on the results
        embed()
        scatt_img_alt = scattered_light_model(x0, _frame_pad)[detpad:-detpad, detpad:-detpad]
        from matplotlib import pyplot as plt
        vmin, vmax = 0, np.max(scatt_img)#40
        plt.imshow(frame - scatt_img, vmin=-vmax/2, vmax=vmax/2)
        plt.show()
        print(res_lsq.x)

        plt.subplot(221)
        plt.imshow(_frame, vmin=vmin, vmax=vmax)
        plt.subplot(222)
        plt.imshow(scatt_img, vmin=vmin, vmax=vmax)
        plt.subplot(223)
        plt.imshow(frame - scatt_img, vmin=-vmax/2, vmax=vmax/2)
        plt.subplot(224)
        plt.imshow(_frame, vmin=-vmax/2, vmax=vmax/2)
        # plt.subplot(235)
        # plt.imshow(scatt_img_alt, vmin=vmin, vmax=vmax)
        # plt.subplot(236)
        # plt.imshow(_frame - scatt_img_alt, vmin=vmin, vmax=vmax)
        plt.show()
    return scatt_img, res_lsq.x, success


def mask_slit_regions(offslitmask, centrace, mask_regions=None):
    """ Exclude some user-specified inter-slit regions for the fine correction to the scattered light determination

    Parameters
    ----------
    offslitmask : `numpy.ndarray`_
        A boolean mask indicating the pixels that are on/off the slit (True = off the slit)
    centrace : `numpy.ndarray`_
        A 2D array, shape is (nspec, nslit), containing the central trace of each slit.
    mask_regions : :obj:`int`, :obj:`list`, optional
        A list of regions that should be excluded in the fine correction to the scattered light. A zero
        value indicates that all pixels left of the first slit will be masked. A value of one indicates
        that all pixels between the first and second slit will be masked, and so forth. A list of these
        integers can be supplied. If None, no inter-slit regions will be excluded. If an integer is
        specified, just one region will be excluded.

    Returns
    -------
    good_mask : `numpy.ndarray`_
        A 2D boolean array indicating the pixels that are suitable to use for the fine scattered light correction.
    """
    # Check if there are regions to be masked
    if mask_regions is None:
        msgs.warn("There are no inter-slit regions specified that need to be masked")
        return offslitmask
    elif isinstance(mask_regions, int):
        # Convert this to a list
        _mask_regions = [mask_regions]
    else:
        _mask_regions = mask_regions

    # Grab the dimensions
    nspec, nspat = offslitmask.shape
    nslit = centrace.shape[1]

    # Setup the pixel coordinates
    spat = np.arange(nspat)

    # Loop through all user-specified inter-slit regions to mask
    bad_mask = np.zeros_like(offslitmask)
    for ii in _mask_regions:
        if ii == 0:  # All pixels to the left of the first slit
            mask_pix = spat[None, :] < centrace[:, ii, None]
        elif ii == nslit:  # All pixels to the right of the last slit
            mask_pix = spat[None, :] > centrace[:, ii-1, None]
        else:  # Everything else in between
            mask_pix = (spat[None, :] > centrace[:, ii-1, None]) \
                   & (spat[None, :] < centrace[:, ii, None])
        # Exclude these pixels
        bad_mask[mask_pix] = True
    # Return the mask of good inter-slit pixels
    return offslitmask & np.logical_not(bad_mask)


def fine_correction(frame, bpm, offslitmask, method='median', polyord=2, debug=False):
    """ Calculate a fine correction to the residual scattered light of the input frame.

    Parameters
    ----------
    frame : `numpy.ndarray`_
        Raw 2D data frame (nspec, nspat) to be used to compute the fine correction of the scattered light.
        This frame should be the raw frame, minus the first estimate of the scattered light
        that has been derived from the :func:`scattered_light_model` function.
    bpm : `numpy.ndarray`_
        2D boolean array indicating the bad pixels (True=bad), same shape as frame
    offslitmask : `numpy.ndarray`_
        A boolean mask indicating the pixels that are on/off the slit (True = off the slit), same shape as frame
    method : :obj:`str`, optional
        Method to use to determine the fine correction to the scattered light. Options are:
            - 'median': Use the median of the off-slit pixels to determine the scattered light
            - 'poly': Use a polynomial fit to the off-slit pixels to determine the scattered light. If this
                        option is chosen, the polynomial order must be specified. See `polyord` below.
    polyord : :obj:`int`, optional
        Polynomial order to use for fitting the residual scattered light in the spatial direction.
    debug : :obj:`bool`, optional
        If True, debug the final model fit that's been output

    Returns
    -------
    scatt_img : `numpy.ndarray`_
        A 2D image (nspec, nspat) of the fine correction to the scattered light determined from the input frame.
    """
    if method not in ['median', 'poly']:
        msgs.error("Unrecognized method to determine the fine correction to the scattered light: {:s}".format(method))
    msgs.info("Performing a fine correction to the scattered light using the {:s} method".format(method))
    nspec, nspat = frame.shape
    if method == 'median':
        # Use the median of the off-slit pixels to determine the scattered light
        mask = bpm | np.logical_not(offslitmask)
        frame_msk = np.ma.array(frame, mask=mask)
        scatt_light_fine = np.repeat(np.ma.median(frame_msk, axis=1).data[:, np.newaxis], nspat, axis=1)
    elif method == 'poly':
        # Convert the BPM to a GPM for convenience
        gpm = np.logical_not(bpm)

        # Define some useful variables
        xspat = np.linspace(0, 1, nspat)
        model = np.zeros_like(frame)

        # Loop over the residual scattered light in the spectral direction and perform
        # a low order polynomial fit to the scattered light in the spatial direction.
        for yy in range(nspec):
            ext = frame[yy, :]
            gd = np.where(offslitmask[yy, :] & gpm[yy, :])
            coeff = np.polyfit(xspat[gd], ext[gd], polyord)
            model[yy, :] = np.polyval(coeff, xspat)
        # Median filter in the spectral direction to smooth out irregularities in the fine correction
        model_med = ndimage.median_filter(model, size=(50, 1))  # Median filter to get rid of CRs
        scatt_light_fine = ndimage.gaussian_filter(model_med, sigma=10)  # Gaussian filter to smooth median filter
    if debug:
        from matplotlib import pyplot as plt
        vmin, vmax = -np.max(scatt_light_fine), np.max(scatt_light_fine)
        plt.subplot(121)
        plt.imshow(frame, vmin=vmin, vmax=vmax)
        plt.subplot(122)
        plt.imshow(frame-scatt_light_fine, vmin=vmin, vmax=vmax)
        plt.show()
    # Return the fine correction model of the scattered light
    return scatt_light_fine
