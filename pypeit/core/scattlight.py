""" Module for sky subtraction

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import numpy as np

from scipy.optimize import least_squares
from scipy import signal, interpolate

from IPython import embed

from pypeit import msgs


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
        The first 5 parameters are the same for all spectrographs. The sixth and all
        subsequent parameters in this list are polynomial corrections to the spectral
        direction. Here are the individual parameters:

        * param[0] = Gaussian kernel width in the spectral direction
        * param[1] = Gaussian kernel width in the spatial direction
        * param[2] = Lorentzian kernel width in the spectral direction
        * param[3] = Lorentzian kernel width in the spatial direction
        * param[4] = Pixel shift of the scattered light in the spectral direction
        * param[5] = Pixel shift of the scattered light in the spatial direction
        * param[6] = Zoom factor of the scattered light (~1)
        * param[7] = Kernel angle
        * param[8] = Relative importance of Gaussian vs Lorentzian.
                     0 < value < 1 means Lorentzian is weighted more
                     value > 1 means Gaussian is weighted more.
        * param[9:] = Polynomial scaling coefficients
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
    shft_spec, shft_spat, zoom = param[4], param[5], param[6]
    kern_angle, kern_scale = param[7], param[8]
    polyterms = param[9:]

    # Make a grid of coordinates
    specvec, spatvec = np.arange(img.shape[0]), np.arange(img.shape[1])
    spat, spec = np.meshgrid(spatvec/(spatvec.size-1), specvec/(specvec.size - 1))
    # Generate the polynomial efficiency scaling in the spectral direction
    polyscale = polyterms[0] + polyterms[1]*spec + polyterms[2]*spat + polyterms[3]*spec*spat
    # polyscale = np.zeros_like(img)
    # for pp in range(polyterms_spec.size):
    #     polyscale += polyterms_spec[pp] * spec ** pp
    # for pp in range(polyterms_spat.size):
    #     polyscale += polyterms_spat[pp] * spat ** pp

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

    # Convolve the input image (note: most of the time is spent here)
    # oaconvolve is the fastest option when the kernel is much smaller dimensions than the image
    scale_img = polyscale * signal.fftconvolve(img, kernel, mode='same')
    # scale_img = polyscale * signal.oaconvolve(img, kernel, mode='same')
    spl = interpolate.RectBivariateSpline(specvec, spatvec, scale_img, kx=1, ky=1)
    return spl(zoom * (specvec + shft_spec), zoom * (spatvec + shft_spat))


def evaluate(param, frame, detpad=300):
    """ Evaluate the scattered light model, allowing for some image padding to minimise edge effects

    Parameters
    ----------
    param : `numpy.ndarray`_
        Model parameters that determine the scattered light based on the input img.
        See `scattered_light_model()`_
    img : `numpy.ndarray`_
        Raw image that you want to compute the scattered light model.
        shape is (nspec, nspat)

    Returns
    -------
    model : `numpy.ndarray`_
        Model of the scattered light for the input
    """
    _frame = frame.copy()
    _frame[0, :] = np.median(_frame[0:10, :], axis=0)
    _frame[-1, :] = np.median(_frame[-10:, :], axis=0)
    img = np.pad(_frame, detpad, mode='edge')  # Model should be generated on padded data
    return scattered_light_model(param, img)[detpad:-detpad, detpad:-detpad]


def scattlight_resid(param, wpix, img, modimg):
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
    model = scattered_light_model(param, modimg)
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
    _frame = frame.copy()
    _frame[0, :] = np.median(_frame[0:10, :], axis=0)
    _frame[-1, :] = np.median(_frame[-10:, :], axis=0)
    img = np.pad(_frame, detpad, mode='edge')  # Model should be generated on padded data
    offslitmask_pad = np.pad(offslitmask * np.logical_not(bpm), detpad, mode='constant',
                             constant_values=0)  # but don't include padded data in the fit
    # Grab the pixels to be included in the fit
    wpix = np.where(offslitmask_pad)

    # Iterate on the scattered light
    scatt_img = 0.0
    niter = 1
    for ii in range(niter):
        # Grab a copy of the input frame, and do some pre-processing on it
        _frame = frame.copy()-scatt_img

        # First pad the edges to minimize edge effects
        # Do a median filter near the edges
        _frame[0, :] = np.median(_frame[0:10, :], axis=0)
        _frame[-1, :] = np.median(_frame[-10:, :], axis=0)
        modimg = np.pad(_frame, detpad, mode='edge')  # Model should be generated on padded data

        # Compute the best-fitting model parameters
        msgs.info(f"ITERATION {ii+1}: Computing best-fitting model parameters of the scattered light")
        res_lsq = least_squares(scattlight_resid, x0, bounds=bounds, args=(wpix, img, modimg),
                                verbose=2, ftol=1.0E-4)

        # Store if this is a successful fit
        success = res_lsq.success
        if success:
            msgs.info("Generating best-fitting scattered light model")
            scatt_img = scattered_light_model(res_lsq.x, modimg)[detpad:-detpad, detpad:-detpad]
        else:
            msgs.warn("Scattered light model fitting failed")
            scatt_img = np.zeros_like(_frame)
            break
    if debug:
        # Do some checks on the results
        embed()
        scatt_img_alt = scattered_light_model(x0, img)[detpad:-detpad, detpad:-detpad]
        from matplotlib import pyplot as plt
        vmin, vmax = 0, 40.0#np.max(scatt_img_alt)
        plt.imshow(frame - scatt_img, vmin=-vmax/2, vmax=vmax/2)
        plt.show()
        print(res_lsq.x)

        plt.subplot(231)
        plt.imshow(_frame, vmin=vmin, vmax=vmax)
        plt.subplot(232)
        plt.imshow(scatt_img, vmin=vmin, vmax=vmax)
        plt.subplot(233)
        plt.imshow(frame - scatt_img, vmin=-vmax/2, vmax=vmax/2)
        plt.subplot(234)
        plt.imshow(_frame, vmin=vmin, vmax=vmax)
        plt.subplot(235)
        plt.imshow(scatt_img_alt, vmin=vmin, vmax=vmax)
        plt.subplot(236)
        plt.imshow(_frame - scatt_img_alt, vmin=vmin, vmax=vmax)
        plt.show()
    return scatt_img, res_lsq.x, success
