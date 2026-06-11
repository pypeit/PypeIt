import numpy as np
import matplotlib.pyplot as plt
import astropy.convolution

from pypeit import utils


def conv2res(wavelength, flux, resolution, central_wl='midpt',
             debug=False):
    """Convolve an imput spectrum to a specific resolution. This is only
    approximate. It takes a fix FWHM for the entire spectrum given by:
    fwhm = wl_cent / resolution

    Parameters
    ----------
    wavelength : `numpy.ndarray`_
        wavelength
    flux : `numpy.ndarray`_
        flux
    resolution : float
        resolution of the spectrograph
    central_wl 
        if 'midpt' the central pixel of wavelength is used, otherwise
        the central_wl will be used.
    debug : bool
        If True will show debug plots

    Returns
    -------
    flux_convolved : `numpy.ndarray`_
        Resulting flux after convolution
    px_sigma : float
        Size of the sigma in pixels at central_wl
    px_bin : float
        Size of one pixel at central_wl
    """

    if central_wl == 'midpt':
        wl_cent = np.median(wavelength)
    else:
        wl_cent = float(central_wl)
    wl_sigma =  wl_cent / resolution / 2.355
    wl_bin = np.abs((wavelength - \
                     np.roll(wavelength,1))[np.where( 
                         np.abs(wavelength-wl_cent) == np.min(np.abs(wavelength-wl_cent)) )])
    px_bin = wl_bin[0]
    px_sigma = wl_sigma / px_bin

    gauss_kernel = astropy.convolution.Gaussian1DKernel(px_sigma)

    flux_convolved = astropy.convolution.convolve(flux, gauss_kernel)

    if debug:
        utils.pyplot_rcparams()
        plt.figure()
        plt.plot(wavelength, flux,
                 color='navy', linestyle='-', alpha=0.8,
                 label=r'Original')
        plt.plot(wavelength, flux_convolved, 
                 color='crimson', linestyle='-', alpha=0.8,
                 label=r'Convolved')
        plt.legend()
        plt.xlabel(r'Wavelength')
        plt.ylabel(r'Flux')
        plt.title(r'Spectrum Convolved at R = {}'.format(resolution))
        plt.show(block=True)
        plt.close()
        utils.pyplot_rcparams_default()

    return flux_convolved, px_sigma, px_bin
