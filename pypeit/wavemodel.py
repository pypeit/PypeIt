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
from pypeit import utils


PLANCK  = 6.62606885e-27 # erg*s
C_LIGHT = 2.99792458e+10 # cm/s
K_BOLTZ = 1.38065040e-16 # erg/K

RADIAN_PER_ARCSEC = 1. / 3600. * 3.14159 / 180.

def transparency(wavelength, debug=False):
    """ Interpolate the atmospheric transmission model in the IR over
    a given wavelength (in microns) range.

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
        utils.pyplot_rcparams()
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

def blackbody(wavelength, T_BB=250., debug=False):
    """ Given wavelength [in microns] and Temperature in Kelvin
    it returns the black body emission.

    Parameters
    ----------
    wavelength : np.array
        wavelength vector in microns
    T_BB : float 
        black body temperature in Kelvin. Default is set to:
        T_BB = 250.

    Returns
    -------
    blackbody : np.array
        spectral radiance of the black body in cgs units:
        B_lambda = 2.*h*c^2/lambda^5.*(1./(exp(h*c/(lambda*k_b*T_BB))-1.)
    blackbody_counts : np.array
        Same as above but in flux density
    """

    msgs.info("Creating BB spectrum at T={}K".format(T_BB))
    lam = wavelength / 1e4 # convert wave in cm.
    
    blackbody_pol = 2.*PLANCK*np.power(C_LIGHT,2) / np.power(lam,5)
    blackbody_exp = np.exp(PLANCK*C_LIGHT/(lam*K_BOLTZ*T_BB)) - 1.
    blackbody = blackbody_pol / blackbody_exp

    blackbody_counts = blackbody / (PLANCK * C_LIGHT / lam) * 1e-4 \
                 * np.power(RADIAN_PER_ARCSEC, 2.)

    if debug:
        utils.pyplot_rcparams()
        msgs.info("Plot of the blackbody spectrum.")
        plt.figure()
        plt.plot(wavelength, blackbody,
                 color='navy', linestyle='-', alpha=0.8,
                 label=r'T_BB={}'.format(T_BB))
        plt.legend()
        plt.xlabel(r"Wavelength [micron]")
        plt.ylabel(r"Spectral Radiance")
        plt.title(r"Planck's law")
        plt.show()

    return blackbody, blackbody_counts

def oh_lines():
    """ Reads in the Rousselot (2000) OH line list"

    Returns
    -------
    wavelength, amplitute : np.arrays
        Wavelength [in microns] and amplitude of the OH lines.
    """

    msgs.info("Reading in the Rousselot (2000) OH line list")
    skisim_dir = resource_filename('pypeit', 'data/skisim/')
    oh = np.loadtxt(skisim_dir+"rousselot2000.dat", usecols=(0, 1))
    return oh[:,0]/10000., oh[:,1] # wave converted to microns

def h2o_lines():
    """ Reads in the H2O atmospheric spectrum"

    Returns
    -------
    wavelength, flux : np.arrays
        Wavelength [in microns] and flux of the H2O atmospheric 
        spectrum.
    """

    msgs.info("Reading in the water atmsopheric spectrum")
    skisim_dir = resource_filename('pypeit', 'data/skisim/')
    h2o = np.loadtxt(skisim_dir+"HITRAN.txt", usecols=(0, 1))
    h2o_wv = 1./ h2o[:,0] * 1e4 # microns
    h2o_rad = h2o[:,1] * 5e11

    return h2o_wv, h2o_rad

def addlines2spec(wavelength, wl_line, fl_line, resolution,
                  scale_spec=1., debug=False):
    """ Create a spectrum with a set of (gaussian) emission lines.
    
    Parameters
    ----------
    wavelength : np.array
        wavelength vector of the input spectrum
    wl_line, fl_line : np.arrays
        wavelength and flux of each individual line
    resolution : np.float
        resolution of the spectrograph. In other words, the lines
        will have a FWHM equal to:
        fwhm_line = wl_line / resolution
    scale_spec : np.float
        rescale all the  normalization of the final spectrum.
        Default scale_spec=1.

    Returns
    -------
    line_spec : np.array
        Spectrum with lines

    """
    line_spec = np.zeros_like(wavelength)
    wl_line_min, wl_line_max = np.min(wavelength), np.max(wavelength)
    good_lines = (wl_line>wl_line_min) & (wl_line<wl_line_max)
    wl_line_good = wl_line[good_lines]
    fl_line_good = fl_line[good_lines]

    # define sigma of the gaussians
    sigma = wl_line_good / resolution / 2.355

    msgs.info("Creating line spectrum")
    for ii in np.arange(len(wl_line_good)):
        # ToDo EMA: This was ported from XIDL, but at the moment does not 
        # have correct normalization for the line flux.
        line_spec = line_spec + \
                    ( fl_line_good[ii] * scale_spec * \
                     np.exp(-np.power((wl_line_good[ii]-wavelength),2.)/(2.*np.power(sigma[ii],2.)) ) )

    if debug:
        utils.pyplot_rcparams()
        msgs.info("Plot of the line spectrum.")
        plt.figure()
        plt.plot(wavelength, line_spec,
                 color='navy', linestyle='-', alpha=0.8,
                 label=r'Spectrum with lines included')
        plt.legend()
        plt.xlabel(r'Wavelength')
        plt.ylabel(r'Flux')
        plt.show()

    return line_spec

def nearIR_modelsky(resolution, waveminmax=(0.8,2.6), dlam=40.0,
                    flgd=True, outfile='SKY_MODEL.fits', T_BB=250.,
                    SCL_BB=1., SCL_OH=1., SCL_H2O=1.,
                    WAVE_WATER=2.3, debug=False):
    """ Generate a model sky in the near-IR. This includes a continuum model
    to match to gemini broadband level, a black body at T_BB, OH lines, and 
    H2O lines (but only at lambda>WAVE_WATER). Everythins is smoothed at the
    given resolution.

    Parameters
    ----------
    resolution : np.float
        resolution of the spectrograph. The OH and H2O lines will have a 
        FWHM equal to:
        fwhm_line = wl_line / resolution
    waveminmax : tuple
        wavelength range in microns to be covered by the model.
        Default is: (0.8, 2.6)
    dlam : 
        bin to be used to create the wavelength grid of the model.
        If flgd='True' it is a bin in velocity (km/s). If flgd='False'
        it is a bin in linear space (microns).
        Default is: 40.0 (with flgd='True')
    flgd : boolean
        if flgd='True' (default) wavelengths are created with 
        equal steps in log space. If 'False', wavelengths will be
        created wit equal steps in linear space.
    outfile : str
        name of the fits file where the model sky spectrum will be stored.
        default is: 'SKY_MODEL.fits'
    T_BB : float 
        black body temperature in Kelvin. Default is set to:
        T_BB = 250.
    SCL_BB : float
        scale factor for modelling the sky black body emssion.
        Default: SCL_BB=1.
    SCL_OH : float
        scale factor for modelling the OH emssion.
        Default: SCL_OH=1.
    SCL_H2O : float
        scale factor for modelling the H2O emssion.
        Default: SCL_H2O=1.
    WAVE_WATER : float
        wavelength (in microns) at which the H2O are inclued.
        Default: WAVE_WATER = 2.3

    Returns
    -------
    wave, sky_model : np.arrays
        wavelength (in Ang.) and flux of the final model of the sky.
    """

    # Create the wavelength array:
    wv_min = waveminmax[0]
    wv_max = waveminmax[1]
    if flgd :
        msgs.info("Creating wavelength vector in velocity space.")
        velpix = dlam # km/s
        loglam = np.log10(1.0 + velpix/299792.458)
        wave = np.power(10.,np.arange(np.log10(wv_min), np.log10(wv_max), loglam))
    else :
        msgs.info("Creating wavelength vector in linear space.")
        wave = np.arange(wv_min, wv_max, dlam)

    # Calculate transparency
    # trans = transparency(wave, debug=False)

    # Empirical match to gemini broadband continuum level
    logy = - 0.55 - 0.55 * (wave-1.0)
    y = np.power(10.,logy)

    msgs.info("Add in a blackbody for the atmosphere.")
    bb, bb_counts = blackbody(wave, T_BB=T_BB, debug=False)
    bb_counts = bb_counts

    msgs.info("Add in OH lines")
    oh_wv, oh_fx = oh_lines()
    oh_wv = oh_wv  # convert to microns
    # produces better wavelength solutions with 1.0 threshold
    msgs.info("Selecting stronger OH lines")
    filt_oh = oh_fx > 1.
    oh_wv, oh_fx = oh_wv[filt_oh], oh_fx[filt_oh]
    ohspec = addlines2spec(wave, oh_wv, oh_fx, resolution=resolution,
                           scale_spec=((resolution/1000.)/160.),
                           debug=False)

    if wv_max > WAVE_WATER :
        msgs.info("Add in H2O lines")
        h2o_wv, h2o_rad = h2o_lines()
        filt_h2o = (h2o_wv>wv_min-0.1) & (h2o_wv<wv_max+0.1)
        h2o_wv  = h2o_wv[filt_h2o]
        h2o_rad = h2o_rad[filt_h2o]
        # calculate sigma at the mean wavelenght of the H2O spectrum
        filt_h2o_med = h2o_wv>WAVE_WATER
        mn_wv = np.mean(h2o_wv[filt_h2o_med])
        # Convolve to the instrument resolution.  This is only
        # approximate.
        smooth_fx, dwv, h2o_dwv = conv2res(h2o_wv, h2o_rad,
                                           resolution,
                                           central_wl = mn_wv,
                                           debug=False)
        # Interpolate over input wavelengths
        interp_h2o = interpolate.interp1d(h2o_wv, smooth_fx,
                                          kind='cubic', 
                                          fill_value='extrapolate')
        h2ospec = interp_h2o(wave)
        # Zero out below WAVE_WATER microns (reconsider)
        h2ospec[wave<WAVE_WATER] = 0.
        h2ospec[wave>np.max(h2o_wv)] = 0.

    sky_model = y+bb_counts*SCL_BB+ohspec*SCL_OH+h2ospec*SCL_H2O


    msgs.info("Saving the sky model in: SKY_MODEL.fits")
    hdu = fits.PrimaryHDU(np.array(sky_model))
    header = hdu.header
    if flgd :
        header['CRVAL1'] = np.log10(wv_min)
        header['CDELT1'] = loglam
        header['DC-FLAG'] = 1
    else :
        header['CRVAL1'] = wv_min
        header['CDELT1'] = dlam
        header['DC-FLAG'] = 0
    hdu.writeto('SKY_MODEL.fits', overwrite = True)

    if debug:
        utils.pyplot_rcparams()
        msgs.info("Plot of the sky emission at R={}".format(resolution))
        plt.figure()
        plt.plot(wave, sky_model,
                 color='black', linestyle='-', alpha=0.8,
                 label=r'Sky Model')
        plt.plot(wave, y,
                 color='darkorange', linestyle='-', alpha=0.6,
                 label=r'Continuum')
        plt.plot(wave, bb_counts*SCL_BB, 
                 color='green', linestyle='-', alpha=0.6,
                 label=r'Black Body at T={}K'.format(T_BB))
        plt.plot(wave, ohspec*SCL_OH,
                 color='darkviolet', linestyle='-', alpha=0.6,
                 label=r'OH')
        plt.plot(wave, h2ospec*SCL_H2O,
                 color='dodgerblue', linestyle='-', alpha=0.6,
                 label=r'H2O')
        plt.legend()
        plt.xlabel(r'Wavelength [microns]')
        plt.ylabel(r'Emission')
        plt.title(r' Sky Emission Spectrum at R={}'.format(resolution))
        plt.show()

    return np.array(wave*10000.), np.array(sky_model)
