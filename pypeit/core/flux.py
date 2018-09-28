""" Module for fluxing routines
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import glob
import numpy as np
import scipy

from pkg_resources import resource_filename

from astropy import units
from astropy import constants
from astropy import coordinates
from astropy.table import Table, Column
from astropy.io import ascii
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from matplotlib import pyplot as plt

from matplotlib import pyplot as plt

try:
    from linetools.spectra.xspectrum1d import XSpectrum1D
except ImportError:
    pass

from pypeit.core import pydl
from pypeit import msgs
from pypeit import utils
from pypeit import debugger
from pypeit.core import qa

TINY = 1e-15
MAGFUNC_MAX = 25.0
MAGFUNC_MIN = -25.0
SN2_MAX = (20.0) ** 2

def apply_sensfunc(spec_obj, sensfunc, airmass, exptime, 
                   spectrograph, MAX_EXTRAP=0.05):
    """ Apply the sensitivity function to the data
    We also correct for extinction.

    Parameters
    ----------
    spec_obj : dict
      SpecObj
    sensfunc : dict
      Sens Function dict
    airmass : float
      Airmass
    exptime : float
      Exposure time in seconds
    spectrograph : dict
      Instrument specific dict
      Used for extinction correction
    MAX_EXTRAP : float, optional [0.05]
      Fractional amount to extrapolate sensitivity function
    """

    # ToDo Is MAX_EXTRAP necessary?
    # Loop on extraction modes
    for extract_type in ['boxcar', 'optimal']:
        extract = getattr(spec_obj, extract_type)
        if len(extract) == 0:
            continue
        msgs.info("Fluxing {:s} extraction for:".format(extract_type) + msgs.newline() +
                  "{}".format(spec_obj))
        wave = np.copy(np.array(extract['WAVE']))
        magfit, _ = sensfunc['mag_set'].value(wave)
        sensfit = np.power(10.0, 0.4 * np.maximum(np.minimum(magfit, MAGFUNC_MAX), MAGFUNC_MIN))
        
        msgs.warn("Extinction correction applyed only if the spectra covers <10000Ang.")
        # Apply Extinction if optical bands
        if np.max(wave) < 10000.:
            msgs.info("Applying extinction correction")
            extinct = load_extinction_data(spectrograph.telescope['longitude'],
                                           spectrograph.telescope['latitude'])
            ext_corr = extinction_correction(wave* units.AA, airmass, extinct)
            sensfit = sensfit * ext_corr
        else:
            msgs.info("Extinction correction not applied")
        extract['FLAM'] = extract['COUNTS'] * sensfit / exptime
        extract['FLAM_SIG'] = (sensfit / exptime) / (np.sqrt(extract['COUNTS_IVAR']))
        extract['FLAM_IVAR'] = extract['COUNTS_IVAR'] / (sensfit / exptime) **2 

def generate_sensfunc(wave, counts, counts_ivar, airmass, exptime, spectrograph, telluric=False, star_type=None,
                      star_mag=None, RA=None, DEC=None, BALM_MASK_WID=5., nresln=None):
    """ Function to generate the sensitivity function.
    This can work in different regimes:
    - If telluric=False and RA=None and Dec=None
      the code creates a sintetic standard star spectrum using the Kurucz models,
      and from this it generates a sens func using nresln=20.0 and masking out
      telluric regions.
    - If telluric=False and RA and Dec are assigned
      the standard star spectrum is extracted from the archive, and a sens func 
      is generated using nresln=20.0 and masking out telluric regions.
    - If telluric=True
      the code creates a sintetic standard star spectrum using the Kurucz models,
      the sens func is created setting nresln=1.5 it contains the correction for
      telluric lines.

    Parameters:
    ----------
    wave : array
      Wavelength of the star with units
    counts : array
      Flux (in counts) of the star
    counts_ivar : array
      Inverse variance of the star
    airmass : float
      Airmass
    exptime : float
      Exposure time in seconds
    spectrograph : dict
      Instrument specific dict
      Used for extinction correction
    telluric : bool
      if True performs a telluric correction
    star_type : str
      Spectral type of the telluric star (used if telluric=True)
    star_mag : float
      Apparent magnitude of telluric star (used if telluric=True)
    RA : float
      deg, RA of the telluric star
      if assigned, the standard star spectrum will be extracted from
      the archive
    DEC : float
      deg, DEC of the telluric star
      if assigned, the standard star spectrum will be extracted from
      the archive
    BALM_MASK_WID : float
      Mask parameter for Balmer absorption. A region equal to
      BALM_MASK_WID*resln is masked wher resln is the estimate
      for the spectral resolution.
    nresln : float
      Number of resolution elements for break-point placement.
      If assigned, overwrites the settings imposed by the code.

    Returns:
    -------
    sens_dict : dict
      sensitivity function described by a dict
    """

    # Create copy of the arrays to avoid modification and convert to 
    # electrons / s
    wave_star = wave.copy()
    flux_star = counts.copy() / exptime
    ivar_star = counts_ivar.copy() * exptime ** 2

    # ToDo
    # This should be changed. At the moment the extinction correction procedure
    # requires the spectra to be in the optical. For the NIR is probably enough
    # to extend the tables to longer wavelength setting the extinction to 0.0mag.
    msgs.warn("Extinction correction applyed only if the spectra covers <10000Ang.")
    # Apply Extinction if optical bands
    if np.max(wave_star) < 10000. * units.AA:
        msgs.info("Applying extinction correction")
        extinct = load_extinction_data(spectrograph.telescope['longitude'],
                                       spectrograph.telescope['latitude'])
        ext_corr = extinction_correction(wave_star, airmass, extinct)
        # Correct for extinction
        flux_star = flux_star * ext_corr
        ivar_star = ivar_star / ext_corr ** 2
    else:
        msgs.info("Extinction correction not applied")

    # Create star model
    if (RA is not None) and (DEC is not None):
        # Pull star spectral model from archive
        msgs.info("Get standard model")
        # Grab closest standard within a tolerance
        std_dict = find_standard_file((RA, DEC))
        # Load standard
        load_standard_file(std_dict)
        # Interpolate onto observed wavelengths
        std_xspec = XSpectrum1D.from_tuple((std_dict['wave'], std_dict['flux']))
        xspec = std_xspec.rebin(wave_star)  # Conserves flambda
        flux_true = xspec.flux.value
    else:
        # Create star spectral model
        msgs.info("Creating standard model")
        # Create star model
        star_loglam, star_flux, std_dict = telluric_sed(star_mag, star_type)
        star_lam = 10 ** star_loglam
        # Generate a dict matching the output of find_standard_file
        std_dict = dict(file='KuruczTelluricModel', name=star_type, fmt=1,
                        ra=None, dec=None)
        std_dict['wave'] = star_lam * units.AA
        std_dict['flux'] = 1e17 * star_flux * units.erg / units.s / units.cm ** 2 / units.AA
        # ToDO If the Kuruck model is used, rebin create weird features
        # I using scipy interpolate to avoid this
        flux_true = scipy.interpolate.interp1d(std_dict['wave'], std_dict['flux'], fill_value='extrapolate')(wave_star)
    
    """
    plt.figure(1)
    plt.plot(std_dict['wave'],std_dict['flux'],label='Orig')
    plt.scatter(wave_star,flux_true,s=5,c='red', label='rebin scipy')
    plt.plot(wave_star,flux_true2,label='rebin X')

    plt.xlim(np.min(wave_star.value),np.max(wave_star.value))
    plt.ylim(np.min(flux_true),np.max(flux_true))
    plt.legend()
    plt.show()
    """
    
    if np.min(flux_true) == 0.:
        msgs.warn('Your spectrum extends beyond calibrated standard star.')

    # Set nresln
    if nresln is None:
        if telluric:
            nresln = 1.5
            msgs.info("Set nresln to 1.5")
        else:
            nresln = 20.0
            msgs.info("Set nresln to 20.0")

    # ToDo
    # Compute an effective resolution for the standard. This could be improved
    # to setup an array of breakpoints based on the resolution. At the
    # moment we are using only one number
    msgs.work("Should pull resolution from arc line analysis")
    msgs.work("At the moment the resolution is taken as the PixelScale")
    msgs.work("This needs to be changed!")
    std_pix = np.median(np.abs(wave_star - np.roll(wave_star, 1)))
    std_res = std_pix
    resln = std_res
    if (nresln * std_res) < std_pix:
        msgs.warn("Bspline breakpoints spacing shoud be larger than 1pixel")
        msgs.warn("Changing input nresln to fix this")
        nresln = std_res / std_pix

    # Mask bad pixels, edges, and Balmer, Paschen, Brackett, and Pfund lines
    # Mask (True = good pixels)
    msgs.info("Masking spectral regions:")
    msk_star = np.ones_like(flux_star).astype(bool)

    # Mask bad pixels
    msgs.info(" Masking bad pixels")
    msk_star[ivar_star <= 0.] = False
    msk_star[flux_star <= 0.] = False

    # Mask edges
    msgs.info(" Masking edges")
    msk_star[:1] = False
    msk_star[-1:] = False

    # Mask Balmer
    msgs.info(" Masking Balmer")
    lines_balm = np.array([3836.4, 3969.6, 3890.1, 4102.8, 4102.8, 4341.6, 4862.7, 5407.0,
                           6564.6, 8224.8, 8239.2]) * units.AA
    for line_balm in lines_balm:
        ibalm = np.abs(wave_star - line_balm) <= BALM_MASK_WID * resln
        msk_star[ibalm] = False

    # Mask Paschen
    msgs.info(" Masking Paschen")
    # air wavelengths from:
    # https://www.subarutelescope.org/Science/Resources/lines/hi.html
    lines_pasc = np.array([8203.6, 9229.0, 9546.0, 10049.4, 10938.1,
                           12818.1, 18751.0]) * units.AA
    for line_pasc in lines_pasc:
        ipasc = np.abs(wave_star - line_pasc) <= BALM_MASK_WID * resln
        msk_star[ipasc] = False

    # Mask Brackett
    msgs.info(" Masking Brackett")
    # air wavelengths from:
    # https://www.subarutelescope.org/Science/Resources/lines/hi.html
    lines_brac = np.array([14584.0, 18174.0, 19446.0, 21655.0,
                           26252.0, 40512.0]) * units.AA
    for line_brac in lines_brac:
        ibrac = np.abs(wave_star - line_brac) <= BALM_MASK_WID * resln
        msk_star[ibrac] = False

    # Mask Pfund
    msgs.info(" Masking Pfund")
    # air wavelengths from:
    # https://www.subarutelescope.org/Science/Resources/lines/hi.html
    lines_pfund = np.array([22788.0, 32961.0, 37395.0, 46525.0,
                            74578.0]) * units.AA
    for line_pfund in lines_pfund:
        ipfund = np.abs(wave_star - line_pfund) <= BALM_MASK_WID * resln
        msk_star[ipfund] = False

    # Mask Atm. cutoff
    msgs.info(" Masking Below the atmospheric cutoff")
    atms_cutoff = wave_star <= 3000.0 * units.AA
    msk_star[atms_cutoff] = False

    if ~telluric:
        # Mask telluric absorption
        msgs.info("Masking Telluric")
        tell = np.any([((wave >= 7580.00 * units.AA) & (wave <= 7750.00 * units.AA)),
                       ((wave >= 7160.00 * units.AA) & (wave <= 7340.00 * units.AA)),
                       ((wave >= 6860.00 * units.AA) & (wave <= 6930.00 * units.AA)),
                       ((wave >= 9310.00 * units.AA) & (wave <= 9665.00 * units.AA)),
                       ((wave >= 11120.0 * units.AA) & (wave <= 11615.0 * units.AA)),
                       ((wave >= 12610.0 * units.AA) & (wave <= 12720.0 * units.AA)),
                       ((wave >= 13160.0 * units.AA) & (wave <= 15065.0 * units.AA)),
                       ((wave >= 15700.0 * units.AA) & (wave <= 15770.0 * units.AA)),
                       ((wave >= 16000.0 * units.AA) & (wave <= 16100.0 * units.AA)),
                       ((wave >= 16420.0 * units.AA) & (wave <= 16580.0 * units.AA)),
                       ((wave >= 17310.0 * units.AA) & (wave <= 20775.0 * units.AA)),
                       (wave >= 22680.0 * units.AA)], axis=0)
        msk_star[tell] = False

    # Apply mask
    ivar_star[~msk_star] = 0.0

    # Fit in magnitudes
    kwargs_bspline = {'bkspace': resln.value * nresln}
    kwargs_reject = {'maxrej': 5}
    mag_set = bspline_magfit(wave_star.value, flux_star, ivar_star, flux_true, inmask=msk_star,
                             kwargs_bspline=kwargs_bspline, kwargs_reject=kwargs_reject)

    # Creating the dict
    msgs.work("Is min, max and wave_min, wave_max a duplicate?")
    sens_dict = dict(bspline=mag_set, func='bspline', min=None, max=None, std=std_dict)

    # Add in wavemin,wavemax
    sens_dict['wave_min'] = np.min(wave_star)
    sens_dict['wave_max'] = np.max(wave_star)
    sens_dict['wave'] = wave_star
    sens_dict['msk_star'] = msk_star
    sens_dict['mag_set'] = mag_set

    """
    # Write the sens_dict to a json file
    import json
    msgs.info("Writing sens_dict into .json file")
    with open('sens_dict.json', 'w') as fp:
        json.dump(sens_dict, fp, sort_keys=True, indent=4)
    """

    return sens_dict


def bspline_magfit(wave, flux, ivar, flux_std, inmask=None, maxiter=35, upper=2, lower=2,
                   kwargs_bspline={}, kwargs_reject={}, debug=False):
    """
    Perform a bspline fit to the flux ratio of standard to
    observed counts. Used to generate a sensitivity function.

    Parameters
    ----------
    wave : ndarray
      wavelength as observed
    flux : ndarray
      counts/s as observed
    ivar : ndarray
      inverse variance
    flux_std : Quantity array
      standard star true flux (erg/s/cm^2/A)
    inmask : ndarray
      bspline mask
    maxiter : integer
      maximum number of iterations for bspline_iterfit
    upper : integer
      number of sigma for rejection in bspline_iterfit
    lower : integer
      number of sigma for rejection in bspline_iterfit
    kwargs_bspline : dict, optional
      keywords for bspline_iterfit
    kwargs_reject : dict, optional
      keywords for bspline_iterfit
    debug : bool
      if True shows some dubugging plots

    Returns
    -------
    bset_log1
    """
    # Create copy of the arrays to avoid modification
    wave_obs = wave.copy()
    flux_obs = flux.copy()
    ivar_obs = ivar.copy()

    # preparing arrays to run in bspline_iterfit

    if np.all(~np.isfinite(ivar_obs)):
        msgs.warn("NaN are present in the inverse variance")

    # Preparing arrays to run in bspline_iterfit

    if np.all(~np.isfinite(ivar_obs)):
        msgs.warn("NaN are present in the inverse variance")

    # Removing outliers

    # Calculate log of flux_obs setting a floor at TINY
    logflux_obs = 2.5 * np.log10(np.maximum(flux_obs, TINY))
    # Set a fix value fro the variance of logflux
    logivar_obs = np.ones_like(logflux_obs) * (10.0 ** 2)

    # Calculate log of flux_std model setting a floor at TINY
    logflux_std = 2.5 * np.log10(np.maximum(flux_std, TINY))

    # Calculate ratio setting a floor at MAGFUNC_MIN and a ceiling at
    # MAGFUNC_MAX
    magfunc = logflux_std - logflux_obs
    magfunc = np.maximum(np.minimum(magfunc, MAGFUNC_MAX), MAGFUNC_MIN)
    magfunc_mask = (magfunc < 0.99 * MAGFUNC_MAX) & (magfunc > 0.99 * MAGFUNC_MIN)

    # Mask outliners
    if inmask is None:
        masktot = (ivar_obs > 0.0) & np.isfinite(logflux_obs) & np.isfinite(ivar_obs) & \
                  np.isfinite(logflux_std) & magfunc_mask
    else:
        masktot = inmask & (ivar_obs > 0.0) & np.isfinite(logflux_obs) & np.isfinite(ivar_obs) & \
                  np.isfinite(logflux_std) & magfunc_mask
    logivar_obs[~masktot] = 0.


    # Calculate sensfunc
    sensfunc = 10.0 ** (0.4 * magfunc)

    msgs.info("Initialize bspline for flux calibration")

    init_bspline = pydl.bspline(wave_obs, bkspace=kwargs_bspline['bkspace'])
    fullbkpt = init_bspline.breakpoints

    # TESTING turning off masking for now
    # remove masked regions from breakpoints
    msk_obs = np.ones_like(wave_obs).astype(bool)
    msk_obs[~masktot] = False
    msk_bkpt = scipy.interpolate.interp1d(wave_obs, msk_obs, kind='nearest', fill_value='extrapolate')(fullbkpt)
    init_breakpoints = fullbkpt[msk_bkpt > 0.999]

    # init_breakpoints = fullbkpt

    #  First round of the fit:
    msgs.info("Bspline fit: step 1")
    bset1, bmask = pydl.iterfit(wave_obs, magfunc, invvar=logivar_obs, inmask=masktot, upper=upper, lower=lower,
                                fullbkpt=init_breakpoints, maxiter=maxiter, kwargs_bspline=kwargs_bspline,
                                kwargs_reject=kwargs_reject)
    logfit1, _ = bset1.value(wave_obs)
    logfit_bkpt, _ = bset1.value(init_breakpoints)

    if debug:
        # Check for calibration
        plt.figure(1)
        plt.plot(wave_obs, magfunc, drawstyle='steps-mid', color='black', label='magfunc')
        plt.plot(wave_obs, logfit1, color='cornflowerblue', label='logfit1')
        plt.plot(wave_obs[~masktot], magfunc[~masktot], '+', color='red', markersize=5.0, label='masked magfunc')
        plt.plot(wave_obs[~masktot], logfit1[~masktot], '+', color='red', markersize=5.0, label='masked logfit1')
        plt.plot(init_breakpoints, logfit_bkpt, '.', color='green', markersize=4.0, label='breakpoints')
        plt.plot(init_breakpoints, np.interp(init_breakpoints, wave_obs, magfunc), '.', color='green', markersize=4.0,
                 label='breakpoints')
        plt.plot(wave_obs, 1.0 / np.sqrt(logivar_obs), color='orange', label='sigma')
        plt.legend()
        plt.xlabel('Wavelength [ang]')
        plt.ylim(0.0, 1.2 * MAGFUNC_MAX)
        plt.title('1st Bspline fit')
        plt.show()

    modelfit1 = np.power(10.0, 0.4 * np.maximum(np.minimum(logfit1, MAGFUNC_MAX), MAGFUNC_MIN))
    residual = sensfunc / (modelfit1 + (modelfit1 == 0)) - 1.
    # new_mask = masktot & (sensfunc > 0)

    # residual_ivar = (modelfit1 * flux_obs / (sensfunc + (sensfunc == 0.0))) ** 2 * ivar_obs
    residual_ivar = np.ones_like(residual) / (0.1 ** 2)
    residual_ivar = residual_ivar * masktot

    (mean, med, stddev) = sigma_clipped_stats(residual[masktot], sigma_lower=3.0, sigma_upper=3.0)

    if np.median(stddev > 0.01):
        #  Second round of the fit:
        msgs.info("Bspline fit: step 2")
        #  Now do one more fit to the ratio of data/model - 1.
        bset_residual, bmask2 = pydl.iterfit(wave_obs, residual, invvar=residual_ivar, inmask=masktot, upper=upper,
                                             lower=lower, maxiter=maxiter, fullbkpt=bset1.breakpoints,
                                             kwargs_bspline=kwargs_bspline, kwargs_reject=kwargs_reject)
        bset_log1 = bset1.copy()
        bset_log1.coeff = bset_log1.coeff + bset_residual.coeff
        if debug:
            # Check for calibration
            resid_fit, _ = bset_residual.value(wave_obs)
            logfit2, _ = bset_log1.value(wave_obs)
            logfit2_bkpt, _ = bset_log1.value(bset1.breakpoints)
            plt.figure(1)
            plt.plot(wave_obs, residual, drawstyle='steps-mid', color='black', label='residual')
            plt.plot(wave_obs, resid_fit, color='cornflowerblue', label='resid_fit')
            plt.plot(wave_obs[~masktot], residual[~masktot], '+', color='red', markersize=5.0, label='masked residual')
            plt.plot(wave_obs[~masktot], resid_fit[~masktot], '+', color='red', markersize=5.0, label='masked resid_fit')
            plt.plot(init_breakpoints, logfit2_bkpt, '.', color='green', markersize=4.0, label='breakpoints')
            plt.plot(wave_obs, 1.0 / np.sqrt(residual_ivar), color='orange', label='sigma')
            plt.legend()
            plt.xlabel('Wavelength [ang]')
            plt.ylim(-0.1, 0.1)
            plt.title('2nd Bspline fit')
            plt.show()
    else:
        bset_log1 = bset1.copy()

    # ToDo JFH I think we should move towards writing this out as a vector in a fits table
    # rather than the b-spline.

    # Create sensitivity function
    newlogfit, _ = bset_log1.value(wave_obs)
    sensfit = np.power(10.0, 0.4 * np.maximum(np.minimum(newlogfit, MAGFUNC_MAX), MAGFUNC_MIN))

    sensfit[~magfunc_mask] = 0.0

    if debug:

        # Check for calibration
        plt.figure(1)
        plt.plot(wave_obs, sensfunc, drawstyle='steps-mid', color='black', label='sensfunc')
        plt.plot(wave_obs, sensfit, color='cornflowerblue', label='sensfunc fit')
        plt.plot(wave_obs[~masktot], sensfunc[~masktot], '+', color='red', markersize=5.0, label='masked sensfunc')
        plt.plot(wave_obs[~masktot], sensfit[~masktot], '+', color='red', markersize=5.0, label='masked sensfuncfit')
        plt.legend()
        plt.xlabel('Wavelength [ang]')
        plt.ylim(0.0, 100.0)
        plt.show()

    # Check quality of the fit
    absdev = np.median(np.abs(sensfit / modelfit1 - 1))
    msgs.info('Difference between fits is {:g}'.format(absdev))

    # Check for residual of the fit
    if debug:

        # scale = np.power(10.0, 0.4 * sensfit)
        flux_cal = flux_obs * sensfit
        ivar_cal = ivar_obs / sensfit ** 2.

        plt.rcdefaults()
        plt.rcParams['font.family']= 'times new roman'
        plt.figure(figsize=(11, 8.5))
        plt.clf()
        plt.plot(wave_obs,flux_cal, label='Calibrated Spectrum')
        plt.plot(wave_obs,flux_std, label='Model')
        plt.plot(wave_obs,np.sqrt(1/ivar_cal))
        plt.legend()
        plt.xlabel('Wavelength [ang]')
        plt.ylabel('Flux [erg/s/cm2/Ang.]')
        plt.ylim(0,np.median(flux_std)*2.5)
        plt.show()
        plt.close()




    # QA
    msgs.work("Add QA for sensitivity function")
    qa_bspline_magfit(wave_obs, bset_log1, magfunc, masktot)
    
    
    """
    bspline_magfit_new_qa(wave_obs, magfunc, logfit1,
                          newlogfit, bset1.breakpoints,
                          outfile=None, title=None)
    """

    return bset_log1

def qa_bspline_magfit(wave, bset, magfunc, mask):
    plt.close("all")
    plt.rcParams['savefig.dpi'] = 600
    plt.rcParams['xtick.top'] = True
    plt.rcParams['ytick.right'] = True
    plt.rcParams['xtick.minor.visible'] = True
    plt.rcParams['ytick.minor.visible'] = True
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['xtick.major.size'] = 6
    plt.rcParams['ytick.major.size'] = 6
    plt.rcParams['xtick.minor.size'] = 3
    plt.rcParams['ytick.minor.size'] = 3
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['xtick.minor.width'] = 1
    plt.rcParams['ytick.minor.width'] = 1
    plt.rcParams['axes.linewidth'] = 1
    plt.rcParams['lines.linewidth'] = 2
    plt.rcParams['legend.frameon'] = False
    plt.rcParams['legend.handletextpad'] = 1.0
    final_fit, _ = bset.value(wave)
    final_fit_bkpt, _ = bset.value(bset.breakpoints)

    plt.figure(1)
    plt.plot(bset.breakpoints, final_fit_bkpt, '.', color='green', markersize=4.0, label='breakpoints')
    plt.plot(wave, magfunc, drawstyle='steps-mid', color='black', label='magfunc')
    plt.plot(wave, final_fit, color='cornflowerblue', label='bspline fit')
    plt.plot(wave[~mask], magfunc[~mask], '+', color='red', markersize=5.0, label='masked points')
    plt.legend()
    plt.xlabel('Wavelength [ang]')
    plt.title('Final Result of the Bspline fit')

    plt.show()
    return

def extinction_correction(wave, airmass, extinct):
    """
    Derive extinction correction
    Based on algorithm in LowRedux (long_extinct)

    Parameters
    ----------
    wave : Quantity array
      Wavelengths for interpolation. Should be sorted
    airmass : float
      Airmass
    extinct : Table
      Table of extinction values

    Returns
    -------
    flux_corr : ndarray
      Flux corrections at the input wavelengths
    """
    # Checks
    if airmass < 1.:
        msgs.error("Bad airmass value in extinction_correction")
    # Interpolate
    f_mag_ext = scipy.interpolate.interp1d(extinct['wave'],
                                           extinct['mag_ext'], bounds_error=False, fill_value=0.)
    mag_ext = f_mag_ext(wave.to('AA').value)

    # Deal with outside wavelengths
    gdv = np.where(mag_ext > 0.)[0]
    if len(gdv) == 0:
        msgs.error(
            "None of the input wavelengths are in the extinction correction range. Presumably something was wrong.")
    if gdv[0] != 0:  # Low wavelengths
        mag_ext[0:gdv[0]] = mag_ext[gdv[0]]
        msgs.warn("Extrapolating at low wavelengths using last valid value")
    if gdv[-1] != (mag_ext.size - 1):  # High wavelengths
        mag_ext[gdv[-1] + 1:] = mag_ext[gdv[-1]]
        msgs.warn("Extrapolating at high wavelengths using last valid value")
    # Evaluate
    flux_corr = 10.0 ** (0.4 * mag_ext * airmass)
    # Return
    return flux_corr


def find_standard_file(radec, toler=20. * units.arcmin, check=False):
    """
    Find a match for the input file to one of the archived
    standard star files (hopefully).  Priority is by order of search.

    Parameters
    ----------
    radec : tuple
      ra, dec in string format ('05:06:36.6','52:52:01.0')
    toler : Angle
      Tolerance on matching archived standards to input
    check : bool
      If True, the routine will only check to see if a
      standard star exists within the input ra, dec, and toler range.

    Returns
    -------
    sdict : dict
      'file': str -- Filename
      'fmt': int -- Format flag
           1=Calspec style FITS binary table
      'name': str -- Star name
      'ra': str -- RA(2000)
      'dec': str -- DEC(2000)
    """
    # Priority
    std_sets = [load_calspec]
    std_file_fmt = [1]  # 1=Calspec style FITS binary table

    # SkyCoord
    obj_coord = coordinates.SkyCoord(radec[0], radec[1], unit=(units.hourangle, units.deg))
    # Loop on standard sets
    closest = dict(sep=999 * units.deg)
    for qq, sset in enumerate(std_sets):
        # Stars
        path, star_tbl = sset()
        star_coords = coordinates.SkyCoord(star_tbl['RA_2000'], star_tbl['DEC_2000'],
                                           unit=(units.hourangle, units.deg))
        # Match
        idx, d2d, d3d = coordinates.match_coordinates_sky(obj_coord, star_coords, nthneighbor=1)
        if d2d < toler:
            if check:
                return True
            else:
                # Generate a dict
                std_dict = dict(file=path + star_tbl[int(idx)]['File'], name=star_tbl[int(idx)]['Name'],
                                fmt=std_file_fmt[qq], ra=star_tbl[int(idx)]['RA_2000'],
                                dec=star_tbl[int(idx)]['DEC_2000'])
                # Return
                msgs.info("Using standard star {:s}".format(std_dict['name']))
                return std_dict
        else:  # Save closest, if it is
            imind2d = np.argmin(d2d)
            mind2d = d2d[imind2d]
            if mind2d < closest['sep']:
                closest['sep'] = mind2d
                closest.update(dict(name=star_tbl[int(idx)]['Name'],ra=star_tbl[int(idx)]['RA_2000'],
                                    dec=star_tbl[int(idx)]['DEC_2000']))
    # Standard star not found
    if check: return False
    msgs.warn("No standard star was found within a tolerance of {:g}".format(toler))
    msgs.info("Closest standard was {:s} at separation {:g}".format(closest['name'], closest['sep'].to('arcmin')))
    msgs.warn("Flux calibration will not be performed")
    return None


def load_calspec():
    """
    Load the list of calspec standards

    Parameters
    ----------

    Returns
    -------
    calspec_path : str
      Path from pypeitdir to calspec standard star files
    calspec_stds : Table
      astropy Table of the calspec standard stars (file, Name, RA, DEC)
    """
    # Read
    calspec_path = '/data/standards/calspec/'
    calspec_file = resource_filename('pypeit', calspec_path + 'calspec_info.txt')
    calspec_stds = Table.read(calspec_file, comment='#', format='ascii')
    # Return
    return calspec_path, calspec_stds


def load_extinction_data(longitude, latitude, toler=5. * units.deg):
    """
    Find the best extinction file to use, based on longitude and latitude
    Loads it and returns a Table

    Parameters
    ----------
    toler : Angle, optional
      Tolerance for matching detector to site (5 deg)

    Returns
    -------
    ext_file : Table
      astropy Table containing the 'wavelength', 'extinct' data for AM=1.
    """
    # Mosaic coord
    mosaic_coord = coordinates.SkyCoord(longitude, latitude, frame='gcrs', unit=units.deg)
    # Read list
    extinct_path = resource_filename('pypeit', '/data/extinction/')
    extinct_summ = extinct_path + 'README'
    extinct_files = Table.read(extinct_summ, comment='#', format='ascii')
    # Coords
    ext_coord = coordinates.SkyCoord(extinct_files['Lon'], extinct_files['Lat'], frame='gcrs',
                                     unit=units.deg)
    # Match
    idx, d2d, d3d = coordinates.match_coordinates_sky(mosaic_coord, ext_coord, nthneighbor=1)
    if d2d < toler:
        extinct_file = extinct_files[int(idx)]['File']
        msgs.info("Using {:s} for extinction corrections.".format(extinct_file))
    else:
        msgs.warn("No file found for extinction corrections.  Applying none")
        msgs.warn("You should generate a site-specific file")
        return None
    # Read
    extinct = Table.read(extinct_path + extinct_file, comment='#', format='ascii',
                         names=('iwave', 'mag_ext'))
    wave = Column(np.array(extinct['iwave']) * units.AA, name='wave')
    extinct.add_column(wave)
    # Return
    return extinct[['wave', 'mag_ext']]


def load_standard_file(std_dict):
    """Load standard star data

    Parameters
    ----------
    std_dict : dict
      Info on standard star indcluding filename in 'file'
      May be compressed

    Returns
    -------
    std_wave : Quantity array
      Wavelengths of standard star array
    std_flux : Quantity array
      Flux of standard star
    """
    root = resource_filename('pypeit', std_dict['file'] + '*')
    fil = glob.glob(root)
    if len(fil) == 0:
        msgs.error("No standard star file: {:s}".format(fil))
    else:
        fil = fil[0]
        msgs.info("Loading standard star file: {:s}".format(fil))
        msgs.info("Fluxes are flambda, normalized to 1e-17")

    if std_dict['fmt'] == 1:
        std_spec = fits.open(fil)[1].data
        # Load
        std_dict['wave'] = std_spec['WAVELENGTH'] * units.AA
        std_dict['flux'] = 1e17 * std_spec['FLUX'] * units.erg / units.s / units.cm ** 2 / units.AA
    else:
        msgs.error("Bad Standard Star Format")
    return


def find_standard(specobjs):
    """Take the median boxcar and then the max object as the standard

    Parameters
    ----------
    specobjs : list

    Returns
    -------

    """

    # Repackage as necessary (some backwards compatibility)
    all_specobj = utils.unravel_specobjs(specobjs)

    # Do it
    medfx = []
    for indx, spobj in enumerate(all_specobj):
        if spobj is None:
            medfx.append(0.)
        else:
            medfx.append(np.median(spobj.boxcar['COUNTS']))
    try:
        mxix = np.argmax(np.array(medfx))
    except:
        debugger.set_trace()
    msgs.info("Putative standard star {} has a median boxcar count of {}".format(all_specobj[mxix], np.max(medfx)))

    # Return
    return mxix


def telluric_params(sptype):
    """Compute physical parameters for a given stellar type.
    This is used by telluric_sed(V, sptype) to create the model spectrum.

    Parameters:
    ----------
    sptype: str
      Spectral type of telluric star

    Returns:
    ----------
    tell_param: dict
      Star parameters
    """

    # log(g) of the Sun
    logg_sol = np.log10(6.67259e-8) + np.log10(1.989e33) - 2.0 * np.log10(6.96e10)

    # Load Schmidt-Kaler (1982) table
    sk82_file = resource_filename('pypeit', 'data/standards/kurucz93/schmidt-kaler_table.txt')
    sk82_tab = ascii.read(sk82_file, names=('Sp', 'logTeff', 'Teff', '(B-V)_0', 'M_V', 'B.C.', 'M_bol', 'L/L_sol'))

    # Match input type
    mti = np.where(sptype == sk82_tab['Sp'])[0]
    if len(mti) != 1:
        raise ValueError('Not ready to interpolate yet.')

    # Calculate final quantities
    # Relation between radius, temp, and bolometric luminosity
    logR = 0.2 * (42.26 - sk82_tab['M_bol'][mti[0]] - 10.0 * sk82_tab['logTeff'][mti[0]])

    # Mass-bolometric luminosity relation from schimdt-kaler p28 valid for M_bol < 7.5
    logM = 0.46 - 0.10 * sk82_tab['M_bol'][mti[0]]
    logg = logM - 2.0 * logR + logg_sol
    M_V = sk82_tab['M_V'][mti[0]]
    tell_param = dict(logR=logR, logM=logM, logg=logg, M_V=M_V,
                      T=sk82_tab['Teff'][mti[0]])

    # Return
    return tell_param


def telluric_sed(V, sptype):
    """Parse Kurucz SED given T and g
    Also convert absolute/apparent magnitudes

    Parameters:
    ----------
    V: float
      Apparent magnitude of the telluric star
    sptype: str
      Spectral type of the telluric star

    Returns:
    ----------
    loglam: ndarray
      log wavelengths
    flux: ndarray
      SED f_lambda (cgs units, I think, probably per Ang)
    """

    # Grab telluric star parameters
    tell_param = telluric_params(sptype)

    # Flux factor (absolute/apparent V mag)

    # Define constants
    parsec = constants.pc.cgs  # 3.086e18
    R_sol = constants.R_sun.cgs  # 6.96e10

    # Distance modulus
    logd = 0.2 * (V - tell_param['M_V']) + 1.0
    D = parsec * 10. ** logd
    R = R_sol * 10. ** tell_param['logR']

    # Factor converts the kurucz surface flux densities to flux observed on Earth
    flux_factor = (R / D.value) ** 2

    # Grab closest T in Kurucz SEDs
    T1 = 3000. + np.arange(28) * 250
    T2 = 10000. + np.arange(6) * 500
    T3 = 13000. + np.arange(22) * 1000
    T4 = 35000. + np.arange(7) * 2500
    Tk = np.concatenate([T1, T2, T3, T4])
    indT = np.argmin(np.abs(Tk - tell_param['T']))

    # Grab closest g in Kurucz SEDs
    loggk = np.arange(11) * 0.5
    indg = np.argmin(np.abs(loggk - tell_param['logg']))

    # Grab Kurucz filename
    std_file = resource_filename('pypeit', '/data/standards/kurucz93/kp00/kp00_{:d}.fits.gz'.format(int(Tk[indT])))
    std = Table.read(std_file)

    # Grab specific spectrum
    loglam = np.array(np.log10(std['WAVELENGTH']))
    gdict = {0: 'g00', 1: 'g05', 2: 'g10', 3: 'g15', 4: 'g20',
             5: 'g25', 6: 'g30', 7: 'g35', 8: 'g40', 9: 'g45',
             10: 'g50'}
    flux = std[gdict[indg]]

    # Generate the standard star dict
    std_dict = dict(stellar_type=sptype, Vmag=V)

    # Return
    return loglam, flux.data * flux_factor, std_dict
