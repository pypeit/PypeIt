# Module for creating models of arc lines.
from __future__ import absolute_import, division, print_function

import inspect
import os

from pkg_resources import resource_filename
from astropy.io import fits
from astropy.convolution import convolve, Gaussian1DKernel
from scipy import interpolate

import numpy as np
import matplotlib.pyplot as plt

from pypeit import msgs
from pypeit.core import arc



PLANCK  = 6.62606885e-27 # erg*s
C_LIGHT = 2.99792458e+10 # cm/s
K_BOLTZ = 1.38065040e-16 # erg/K

RADIAN_PER_ARCSEC = 1. / 3600. * 3.14159 / 180.

def transparency(wavelength, debug=False):
    """ Interpolate atmospheric transmission model in the IR over
    a given wavelength range.

    Parameters
    ----------
    wavelength : np.array
        wavelength vector in microns

    Returns
    -------
    transparency : np.array
        Transmission of the sky over the considered wavelength rage.
        1. means fully transparent and 0. fully opaque
    """
    msgs.info("Reading in the atmospheric transmission model")
    skisim_dir = resource_filename('pypeit', 'data/skisim/')
    transparency = np.loadtxt(skisim_dir+'atm_transmission_secz1.5_1.6mm.dat')
    wave_mod = transparency[:,0]
    tran_mod = transparency[:,1]

    # Limit model between 0.8 and np.max(wavelength) microns
    filt_wave_mod = (wave_mod>0.8) & (wave_mod<np.max(wavelength))
    wave_mod = wave_mod[filt_wave_mod]
    tran_mod = tran_mod[filt_wave_mod]
    
    
    # Interpolate over input wavelengths
    interp_tran = interpolate.interp1d(wave_mod, tran_mod,
                                       kind='cubic',
                                       fill_value='extrapolate')
    transmission = interp_tran(wavelength)
    transmission[wavelength<0.9] = 1.

    # Clean for spourious values due to interpolation
    transmission[transmission<0.] = 0.
    transmission[transmission>1.] = 1.

    if debug:
        prettyplot()
        msgs.info("Plot of the sky transmission template")
        plt.figure()
        plt.plot(wave_mod, tran_mod,
                 color='navy', linestyle='-', alpha=0.8,
                 label=r'Original')
        plt.plot(wavelength, transmission, 
                 color='crimson', linestyle='-', alpha=0.8,
                 label=r'Resampled')
        plt.legend()
        plt.xlabel(r'Wavelength [microns]')
        plt.ylabel(r'Transmission')
        plt.title(r' IR Transmission Spectra ')
        plt.show()

    # Returns
    return transmission

