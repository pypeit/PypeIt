""" Module for sky subtraction

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import numpy as np

from scipy.optimize import least_squares
from scipy import signal, interpolate

from IPython import embed

from pypeit import msgs


def scattered_light_model(param, img, kernel='gaussian'):
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
        The first 5 parameters are the same for all spectrographs. The sixth and all
        subsequent parameters in this list are polynomial corrections to the spectral
        direction. Here are the individual parameters:

        * param[0] = Kernel width in the spectral direction
        * param[1] = Kernel width in the spatial direction
        * param[2] = Pixel shift of the scattered light in the spectral direction
        * param[3] = Pixel shift of the scattered light in the spatial direction
        * param[4] = Zoom factor of the scattered light (~1)
        * param[5:] = Polynomial scaling coefficients in the spectral direction
    img : `numpy.ndarray`_
        Raw image that you want to compute the scattered light model.
        shape is (nspec, nspat)
        Model used to calculate the scattered light. This function is used to
        generate a model of the scattered light, based on a set of model parameters
        that have been optimized using self.scattered_light().
    kernel : :obj:`str`_, optional
        The shape of the kernel to use. The allowed values are 'gaussian' and 'lorentzian',
        but the default value ('gaussian') is currently believed to provide the best fit.

    Returns
    -------
    model : `numpy.ndarray`_
        Model of the scattered light for the input
    """
    # Extract the parameters into more conveniently named variables
    sigmx, sigmy, shft_spec, shft_spat, zoom = param[0], param[1], param[2], param[3], param[4]
    polyterms = param[5:]
    # Generate a 2D smoothing kernel
    if kernel == 'gaussian':
        # Gaussian
        subkrnx = np.exp(-0.5 * (np.arange(int(6 * sigmx)) - 3 * sigmx) ** 2 / sigmx ** 2)
        subkrny = np.exp(-0.5 * (np.arange(int(6 * sigmy)) - 3 * sigmy) ** 2 / sigmy ** 2)
    elif kernel == 'lorentzian':
        # Lorentzian
        subkrnx = sigmx / ((np.arange(int(10 * sigmx)) - 5 * sigmx) ** 2 + sigmx ** 2)
        subkrny = sigmy / ((np.arange(int(10 * sigmy)) - 5 * sigmy) ** 2 + sigmy ** 2)
    else:
        msgs.error(f"Unknown kernel: {kernel}")
    kernel = np.outer(subkrnx, subkrny)
    kernel /= np.sum(kernel)
    # Make a grid of coordinates
    specvec, spatvec = np.arange(img.shape[0]), np.arange(img.shape[1])
    spat, spec = np.meshgrid(spatvec, specvec / (specvec.size - 1))
    # Generate the polynomial efficiency scaling in the spectral direction
    polyscale = np.zeros_like(img)
    for pp in range(polyterms.size):
        polyscale += polyterms[pp] * spec ** pp
    # Convolve the input image (note: most of the time is spent here)
    # oaconvolve is the fastest option when the kernel is much smaller dimensions than the image
    # scale_img = polyscale * signal.fftconvolve(img, kernel, mode='same')
    scale_img = polyscale * signal.oaconvolve(img, kernel, mode='same')
    spl = interpolate.RectBivariateSpline(specvec, spatvec, scale_img, kx=1, ky=1)
    return spl(zoom * (specvec + shft_spec), zoom * (spatvec + shft_spat))


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
        Raw 2D data frame to be used to compute the scattered light.
    bpm : `numpy.ndarray`_
        2D boolean array indicating the bad pixels (True=bad)
    offslitmask : `numpy.ndarray`_
        A boolean mask indicating the pixels that are on/off the slit (True = off the slit)
    x0 : `numpy.ndarray`_
        A 1D array containing the best-fitting model parameters
    bounds : :obj:`tuple`_
        A tuple of two elements, containing two `np.ndarray`_ of the same length as x0. These
        two arrays contain the lower (first element of the tuple) and upper (second element of the tuple)
        bounds to consider on the scattered light model parameters.
    detpad : :obj:`int`_, optional
        Number of pixels to pad to each of the detector edges to reduce edge effects.
    debug : :obj:`bool`_, optional
        If True, debug the final model fit that's been output

    Returns
    -------
    scatt_img : `numpy.ndarray`_
        A 2D image of the scattered light determined from the input frame.
        Alternatively, if a constant value is used, a constant floating point
        value can be returned as well.
    modelpar : `numpy.ndarray`_
        A 1D array containing the best-fitting model parameters
    success : :obj:`bool`_
        True if the fit was successful, False otherwise
    """

    # Grab a copy of the input frame, and do some pre-processing on it
    _frame = frame.copy()

    # First pad the edges to minimize edge effects
    # Do a median filter near the edges
    _frame[0, :] = np.median(_frame[0:10, :], axis=0)
    _frame[-1, :] = np.median(_frame[-10:, :], axis=0)
    img = np.pad(_frame, detpad, mode='edge')  # Model should be generated on padded data
    offslitmask_pad = np.pad(offslitmask * np.logical_not(bpm), detpad, mode='constant',
                             constant_values=0)  # but don't include padded data in the fit
    # Grab the pixels to be included in the fit
    wpix = np.where(offslitmask_pad)

    # Compute the best-fitting model parameters
    msgs.info("Computing best-fitting model parameters of the scattered light")
    res_lsq = least_squares(scattlight_resid, x0, bounds=bounds, args=(wpix, img), verbose=2, ftol=1.0E-4)

    # Store if this is a successful fit
    success = res_lsq.success
    if success:
        msgs.info("Generating best-fitting scattered light model")
        scatt_img = scattered_light_model(res_lsq.x, img)[detpad:-detpad, detpad:-detpad]
    else:
        msgs.warn("Scattered light model fitting failed")
        scatt_img = np.zeros_like(_frame)
    if debug:
        # Do some checks on the results
        embed()
        from matplotlib import pyplot as plt
        scatt_img_alt = scattered_light_model(x0, img)[detpad:-detpad, detpad:-detpad]
        vmin, vmax = 0, np.max(scatt_img_alt)
        plt.subplot(231)
        plt.imshow(_frame, vmin=vmin, vmax=vmax)
        plt.subplot(232)
        plt.imshow(scatt_img, vmin=vmin, vmax=vmax)
        plt.subplot(233)
        plt.imshow(_frame - scatt_img, vmin=-vmax/2, vmax=vmax/2)
        plt.subplot(234)
        plt.imshow(_frame, vmin=vmin, vmax=vmax)
        plt.subplot(235)
        plt.imshow(scatt_img_alt, vmin=vmin, vmax=vmax)
        plt.subplot(236)
        plt.imshow(_frame - scatt_img_alt, vmin=vmin, vmax=vmax)
        plt.show()
    return scatt_img, res_lsq.x, success
