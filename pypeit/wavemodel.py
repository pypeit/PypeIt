"""
Module to create models of arc lines.
"""
import os

import astropy
import re
import scipy

import numpy as np
import matplotlib.pyplot as plt


from astropy.io import fits
from astropy.convolution import convolve, Gaussian1DKernel
from astropy.table import Table
from astropy import units

from pypeit import msgs
from pypeit.core import arc
from pypeit import utils
from pypeit.core.wave import airtovac
from pypeit import io
from pypeit import data

from IPython import embed

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

    # Define constants in cgs
    PLANCK  = astropy.constants.h.cgs.value    # erg*s
    C_LIGHT = astropy.constants.c.cgs.value    # cm/s
    K_BOLTZ = astropy.constants.k_B.cgs.value  # erg/K
    RADIAN_PER_ARCSEC = 1./3600.*np.pi/180.

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
        msgs.info("Close the Figure to continue.")
        plt.show(block=True)
        plt.close()
        utils.pyplot_rcparams_default()

    return blackbody, blackbody_counts


def addlines2spec(wavelength, wl_line, fl_line, resolution,
                  scale_spec=1., debug=False):
    """ Create a spectrum with a set of (gaussian) emission lines.
    
    Parameters
    ----------
    wavelength : np.array
        wavelength vector of the input spectrum
    wl_line, fl_line : np.arrays
        wavelength and flux of each individual line
    resolution : float
        resolution of the spectrograph. In other words, the lines
        will have a FWHM equal to:
        fwhm_line = wl_line / resolution
    scale_spec : float
        rescale all the  normalization of the final spectrum.
        Default scale_spec=1.
    debug : boolean
        If True will show debug plots

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
        line_spec += scale_spec*fl_line_good[ii]*\
                     np.exp(-np.power((wl_line_good[ii]-wavelength),2.)/(2.*np.power(sigma[ii],2.)))

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
        msgs.info("Close the Figure to continue.")
        plt.show(block=True)
        plt.close()
        utils.pyplot_rcparams_default()

    return line_spec


def oh_lines():
    """ Reads in the Rousselot (2000) OH line list"

    Returns
    -------
    wavelength, amplitude : np.arrays
        Wavelength [in microns] and amplitude of the OH lines.
    """

    msgs.info("Reading in the Rousselot (2000) OH line list")
    oh = np.loadtxt(data.get_skisim_filepath('rousselot2000.dat'),
                    usecols=(0, 1))
    return oh[:,0]/10000., oh[:,1] # wave converted to microns


def transparency(wavelength, debug=False):
    """ Interpolate the atmospheric transmission model in the IR over
    a given wavelength (in microns) range.

    Parameters
    ----------
    wavelength : np.array
        wavelength vector in microns
    debug : boolean
        If True will show debug plots

    Returns
    -------
    transparency : np.array
        Transmission of the sky over the considered wavelength rage.
        1. means fully transparent and 0. fully opaque
    """

    msgs.info("Reading in the atmospheric transmission model")
    transparency = np.loadtxt(data.get_skisim_filepath('atm_transmission_secz1.5_1.6mm.dat'))
    wave_mod = transparency[:,0]
    tran_mod = transparency[:,1]

    # Limit model between 0.8 and np.max(wavelength) microns
    filt_wave_mod = (wave_mod>0.8) & (wave_mod<np.max(wavelength))
    wave_mod = wave_mod[filt_wave_mod]
    tran_mod = tran_mod[filt_wave_mod]

    # Interpolate over input wavelengths
    interp_tran = scipy.interpolate.interp1d(wave_mod, tran_mod,
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
        msgs.info("Close the Figure to continue.")
        plt.show(block=True)
        plt.close()
        utils.pyplot_rcparams_default()

    # Returns
    return transmission


def h2o_lines():
    """ Reads in the H2O atmospheric spectrum"

    Returns
    -------
    wavelength, flux : np.arrays
        Wavelength [in microns] and flux of the H2O atmospheric 
        spectrum.
    """

    msgs.info("Reading in the water atmsopheric spectrum")
    h2o = np.loadtxt(data.get_skisim_filepath('HITRAN.txt'),
                     usecols=(0, 1))
    h2o_wv = 1./ h2o[:,0] * 1e4 # microns
    h2o_rad = h2o[:,1] * 5e11 # added to match XIDL

    return h2o_wv, h2o_rad


def thar_lines():
    """ Reads in the H2O atmospheric spectrum"
    Detailed information are here: http://astronomy.swin.edu.au/~mmurphy/thar/index.html

    Returns
    -------
    wavelength, flux : np.arrays
        Wavelength [in angstrom] and flux of the ThAr lamp 
        spectrum.
    """

    msgs.info("Reading in the ThAr spectrum")
    thar = data.load_thar_spec()
    
    # create pixel array
    thar_pix = np.arange(thar[0].header['CRPIX1'],len(thar[0].data[0,:])+1)
    # convert pixels to wavelength in Angstrom
    thar_wv = thar[0].header['UP_WLSRT']*10**((thar_pix-thar[0].header['CRPIX1'])*thar[0].header['CD1_1'])
    # read in spectrum
    thar_spec = thar[0].data[0,:]

    return thar_wv, thar_spec


def nearIR_modelsky(resolution, waveminmax=(0.8,2.6), dlam=40.0,
                    flgd=True, nirsky_outfile=None, T_BB=250.,
                    SCL_BB=1., SCL_OH=1., SCL_H2O=10.,
                    WAVE_WATER=2.3, debug=False):
    """ Generate a model sky in the near-IR. This includes a continuum model
    to match to gemini broadband level, a black body at T_BB, OH lines, and 
    H2O lines (but only at lambda>WAVE_WATER). Everythins is smoothed at the
    given resolution.

    Parameters
    ----------
    resolution : float
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
    nirsky_outfile : str
        name of the fits file where the model sky spectrum will be stored.
        default is 'None' (i.e., no file will be written).
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
        Default: SCL_H2O=10.
    WAVE_WATER : float
        wavelength (in microns) at which the H2O are inclued.
        Default: WAVE_WATER = 2.3
    debug : boolean
        If True will show debug plots

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
    bb, bb_counts = blackbody(wave, T_BB=T_BB, debug=debug)
    bb_counts = bb_counts

    msgs.info("Add in OH lines")
    oh_wv, oh_fx = oh_lines()
    # produces better wavelength solutions with 1.0 threshold
    msgs.info("Selecting stronger OH lines")
    filt_oh = oh_fx > 1.
    oh_wv, oh_fx = oh_wv[filt_oh], oh_fx[filt_oh]
    # scale_spec was added to match the XIDL code
    ohspec = addlines2spec(wave, oh_wv, oh_fx, resolution=resolution,
                           scale_spec=((resolution/1000.)/40.),
                           debug=debug)

    if wv_max > WAVE_WATER :
        msgs.info("Add in H2O lines")
        h2o_wv, h2o_rad = h2o_lines()
        filt_h2o = (h2o_wv>wv_min-0.1) & (h2o_wv<wv_max+0.1)
        h2o_wv  = h2o_wv[filt_h2o]
        h2o_rad = h2o_rad[filt_h2o]
        # calculate sigma at the mean wavelenght of the H2O spectrum
        filt_h2o_med = h2o_wv>WAVE_WATER
        mn_wv = np.mean(h2o_wv[filt_h2o_med])
        # Convolve to the instrument resolution. This is only
        # approximate.
        smooth_fx, dwv, h2o_dwv = conv2res(h2o_wv, h2o_rad,
                                           resolution,
                                           central_wl = mn_wv,
                                           debug=debug)
        # Interpolate over input wavelengths
        interp_h2o = scipy.interpolate.interp1d(h2o_wv, smooth_fx,
                                                kind='cubic', 
                                                fill_value='extrapolate')
        h2ospec = interp_h2o(wave)
        # Zero out below WAVE_WATER microns (reconsider)
        h2ospec[wave<WAVE_WATER] = 0.
        h2ospec[wave>np.max(h2o_wv)] = 0.
    else:
        h2ospec = np.zeros(len(wave),dtype='float')

    sky_model = y+bb_counts*SCL_BB+ohspec*SCL_OH+h2ospec*SCL_H2O

    if nirsky_outfile is not None:
        msgs.info("Saving the sky model in: {}".format(nirsky_outfile))
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
        hdu.writeto(nirsky_outfile, overwrite = True)

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
        plt.title(r'Sky Emission Spectrum at R={}'.format(resolution))
        msgs.info("Close the Figure to continue.")
        plt.show(block=True)
        plt.close()
        utils.pyplot_rcparams_default()

    return np.array(wave*10000.), np.array(sky_model)


def optical_modelThAr(resolution, waveminmax=(3000.,10500.), dlam=40.0,
                      flgd=True, thar_outfile=None, debug=False):
    """ Generate a model of a ThAr lamp in the uvb/optical. This is based on the
    Murphy et al. ThAr spectrum. Detailed information are here:
    http://astronomy.swin.edu.au/~mmurphy/thar/index.html
    Everythins is smoothed at the given resolution.

    Parameters
    ----------
    resolution : float
        resolution of the spectrograph. The ThAr lines will have a 
        FWHM equal to:
        fwhm_line = wl_line / resolution
    waveminmax : tuple
        wavelength range in angstrom to be covered by the model.
        Default is: (3000.,10500.)
    dlam : 
        bin to be used to create the wavelength grid of the model.
        If flgd='True' it is a bin in velocity (km/s). If flgd='False'
        it is a bin in linear space (microns).
        Default is: 40.0 (with flgd='True')
    flgd : boolean
        if flgd='True' (default) wavelengths are created with 
        equal steps in log space. If 'False', wavelengths will be
        created wit equal steps in linear space.
    thar_outfile : str
        name of the fits file where the model sky spectrum will be stored.
        default is 'None' (i.e., no file will be written).
    debug : boolean
        If True will show debug plots

    Returns
    -------
    wave, thar_model : np.arrays
        wavelength (in Ang.) and flux of the final model of the ThAr lamp emission.
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

    msgs.info("Add in ThAr lines")
    th_wv, th_fx = thar_lines()

    # select spectral region
    filt_wl = (th_wv>=wv_min) & (th_wv<=wv_max)
    # calculate sigma at the mean wavelenght of the ThAr spectrum
    mn_wv = np.mean(th_wv[filt_wl])
    # Convolve to the instrument resolution. This is only
    # approximate.
    smooth_fx, dwv, thar_dwv = conv2res(th_wv, th_fx,
                                        resolution,
                                        central_wl = mn_wv,
                                        debug=debug)
    # Interpolate over input wavelengths
    interp_thar = scipy.interpolate.interp1d(th_wv, smooth_fx,
                                             kind='cubic',
                                             fill_value='extrapolate')
    thar_spec = interp_thar(wave)

    # remove negative artifacts
    thar_spec[thar_spec<0.] = 0.
    # Remove regions of the spectrum outside the wavelength covered by the ThAr model
    if wv_min<np.min(th_wv):
        msgs.warn("Model of the ThAr spectrum outside the template coverage.")
        thar_spec[wave<np.min(th_wv)] = 0.
    if wv_max<np.max(th_wv):
        msgs.warn("Model of the ThAr spectrum outside the template coverage.")
        thar_spec[wave>np.max(th_wv)] = 0.

    if thar_outfile is not None:
        msgs.info("Saving the ThAr model in: {}".format(thar_outfile))
        hdu = fits.PrimaryHDU(np.array(thar_spec))
        header = hdu.header
        if flgd :
            header['CRVAL1'] = np.log10(wv_min)
            header['CDELT1'] = loglam
            header['DC-FLAG'] = 1
        else :
            header['CRVAL1'] = wv_min
            header['CDELT1'] = dlam
            header['DC-FLAG'] = 0
        hdu.writeto(thar_outfile, overwrite = True)

    if debug:
        utils.pyplot_rcparams()
        msgs.info("Plot of the Murphy et al. template at R={}".format(resolution))
        plt.figure()
        plt.plot(th_wv, th_fx,
                 color='navy', linestyle='-', alpha=0.3,
                 label=r'Original')
        plt.plot(th_wv, smooth_fx, 
                 color='crimson', linestyle='-', alpha=0.6,
                 label=r'Convolved at R={}'.format(resolution))
        plt.plot(wave, thar_spec, 
                 color='maroon', linestyle='-', alpha=1.0,
                 label=r'Convolved at R={} and resampled'.format(resolution))
        plt.legend()
        plt.xlabel(r'Wavelength [Ang.]')
        plt.ylabel(r'Emission')
        plt.title(r'Murphy et al. ThAr spectrum at R={}'.format(resolution))
        msgs.info("Close the Figure to continue.")
        plt.show(block=True)
        plt.close()
        utils.pyplot_rcparams_default()

    return np.array(wave), np.array(thar_spec)


def conv2res(wavelength, flux, resolution, central_wl='midpt',
             debug=False):
    """Convolve an imput spectrum to a specific resolution. This is only
    approximate. It takes a fix FWHM for the entire spectrum given by:
    fwhm = wl_cent / resolution

    Parameters
    ----------
    wavelength : np.array
        wavelength
    flux : np.array
        flux
    resolution : float
        resolution of the spectrograph
    central_wl 
        if 'midpt' the central pixel of wavelength is used, otherwise
        the central_wl will be used.
    debug : boolean
        If True will show debug plots

    Returns
    -------
    flux_convolved :np.array
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
    wl_bin = np.abs((wavelength - np.roll(wavelength,1))[np.where( np.abs(wavelength-wl_cent) == np.min(np.abs(wavelength-wl_cent)) )])
    msgs.info("The binning of the wavelength array at {} is: {}".format(wl_cent, wl_bin[0]))
    px_bin = wl_bin[0]
    px_sigma = wl_sigma / px_bin

    msgs.info("Covolving with a Gaussian kernel with sigma = {} pixels".format(px_sigma))
    gauss_kernel = Gaussian1DKernel(px_sigma)

    flux_convolved = convolve(flux, gauss_kernel)

    if debug:
        utils.pyplot_rcparams()
        msgs.info("Spectrum Convolved at R = {}".format(resolution))
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
        msgs.info("Close the Figure to continue.")
        plt.show(block=True)
        plt.close()
        utils.pyplot_rcparams_default()

    return flux_convolved, px_sigma, px_bin


def iraf_datareader(database_dir, id_file):
    """Reads in a line identification database created with IRAF
    identify. These are usually locate in a directory called 'database'.
    This read pixel location and wavelength of the lines that
    have been id with IRAF. Note that the first pixel in IRAF
    is '1', while is '0' in python. The pixel location is thus
    shifted of one pixel while reading the database.

    Parameters
    ----------
    database_dir : string
        directory where the id files are located.
    id_file : string
        filename that is going to be read.

    Returns
    -------
    pixel, line_id : np.arrays
        Position of the line in pixel and ID of the line. 
        For IRAF output, these are usually in Ang.
    """

    lines_database = []

    # Open file for reading of text data.
    with open (database_dir+id_file, 'rt') as id_file_iraf:
        for line in id_file_iraf:
            lines_database.append(line.split())
            feat_line = re.search(r'features\t(\d+)', line)
            if feat_line is not None:
                N_lines = int(feat_line.group(1))

    msgs.info("The number of IDs in the IRAF database {} is {}".format(id_file, N_lines))

    pixel = np.zeros(N_lines)
    line_id = np.zeros(N_lines)
    for iii in range(0,N_lines):
        pixel[iii] = lines_database[10:N_lines+10][iii][0]
        line_id[iii] = lines_database[10:N_lines+10][iii][2]

    # Moving from IRAF 1-based to Python 0-based convention.
    pixel = pixel - 1.

    return pixel, line_id


def create_linelist(wavelength, spec, fwhm, sigdetec=2.,
                    cont_samp=10., line_name=None, file_root_name=None,
                    iraf_frmt=False, debug=False, convert_air_to_vac=True):
    """ Create list of lines detected in a spectrum in a PypeIt
    compatible format. The name of the output file is
    file_root_name+'_lines.dat'.

    Parameters
    ----------
    wavelength : np.array
        wavelength
    spec : np.array
        spectrum
    fwhm : float
        fwhm in pixels used for filtering out arc lines that are too
        wide and not considered in fits. Parameter of arc.detect_lines().
    sigdetec : float
        sigma threshold above fluctuations for line detection. Parameter
        of arc.detect_lines(). Default = 2.
    cont_samp : float
        the number of samples across the spectrum used for continuum
        subtraction. Parameter of arc.detect_lines().  Default = 10.
    line_name : str
        name of the lines to listed in the file.
    file_root_name : str
        name of the file where the identified lines will be stored.
        The code automatically add '_lines.dat' at the end of the
        root name.
    iraf_frmt : bool
        if True, the file is written in the IRAF format (i.e. wavelength,
        ion name, amplitude).
    convert_air_to_vac (bool):
        If True, convert the wavelengths of the created linelist from air to vacuum
    """

    msgs.info("Searching for peaks {} sigma above background".format(sigdetec))
    tampl_true, tampl, tcent, twid, centerr, ww, arcnorm, nsig = arc.detect_lines(spec, sigdetect=sigdetec,
                                                                                  fwhm=fwhm, cont_samp=cont_samp,
                                                                                  debug=debug)

    peaks_good = tcent[ww]
    ampl_good = tampl[ww]
    # convert from pixel location to wavelength
    pixvec = np.arange(spec.size)
    wave_peak = scipy.interpolate.interp1d(pixvec, wavelength, bounds_error=False, fill_value='extrapolate')(peaks_good)
    # Convert to vacuum?
    if convert_air_to_vac:
        msgs.info("Converting wavelengths from air to vacuum")
        wave_peak = airtovac(wave_peak * units.AA).value

    npeak = len(wave_peak)
    ion = npeak*[str(line_name)]
    NIST = npeak*[1]
    Instr = npeak*[32]
    Source = npeak*['wavemodel.py']

    if iraf_frmt:
        msgs.info("Printing file in IRAF format: {}".format(file_root_name+'_iraf_lines.dat'))
        ion = np.array(ion)
        id_lines_iraf = np.vstack( (np.round(wave_peak,5), ion, np.round(ampl_good,5)) ).T
        np.savetxt(file_root_name+'_iraf_lines.dat', id_lines_iraf, fmt="%15s %6s %15s", delimiter="  ")
    else:
        msgs.info("Printing file: {}".format(file_root_name+'_lines.dat'))
        dat = Table([wave_peak, ion, NIST, Instr, ampl_good, Source], names=('wave', 'ion','NIST','Instr','amplitude','Source'))
        dat.write(file_root_name+'_lines.dat',format='ascii.fixed_width')


def create_OHlinelist(resolution, waveminmax=(0.8,2.6), dlam=40.0, flgd=True, nirsky_outfile=None,
                      fwhm=None, sigdetec=3., line_name='OH', file_root_name=None, iraf_frmt=False,
                      debug=False):
    """Create a synthetic sky spectrum at a given resolution, extract significant lines, and
    store them in a PypeIt compatibile file. The skymodel is built from nearIR_modelsky and
    includes black body at 250K, OH lines, and H2O lines (but only at lambda>2.3microns).

    Parameters
    ----------
    resolution : float
        resolution of the spectrograph
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
    nirsky_outfile : str
        name of the fits file where the model sky spectrum will be stored.
        default is 'None' (i.e., no file will be written).
    fwhm : float
        fwhm in pixels used for filtering out arc lines that are too
        wide and not considered in fits. Parameter of arc.detect_lines().
        If set to 'None' the fwhm will be derived from the resolution as:
        2. * central_wavelength / resolution
    sigdetec : float
        sigma threshold above fluctuations for line detection. Parameter
        of arc.detect_lines(). Default = 2.
    line_name : str
        name of the lines to listed in the file. Default is 'OH'.
    file_root_name : str
        name of the file where the identified lines will be stored.
        The code automatically add '_lines.dat' at the end of the
        root name.
    iraf_frmt : bool
        if True, the file is written in the IRAF format (i.e. wavelength,
        ion name, amplitude).
    debug : boolean
        If True will show debug plots
    """

    wavelength, spec = nearIR_modelsky(resolution, waveminmax=waveminmax, dlam=dlam,
                                       flgd=flgd, nirsky_outfile=nirsky_outfile, debug=debug)

    if fwhm is None:
        msgs.warn("No min FWHM for the line detection set. Derived from the resolution at the center of the spectrum")
        wl_cent = np.average(wavelength)
        wl_fwhm = wl_cent / resolution
        wl_bin = np.abs((wavelength-np.roll(wavelength,1))[np.where(np.abs(wavelength-wl_cent)==np.min(np.abs(wavelength-wl_cent)))])
        # In order not to exclude all the lines, fwhm is set to 5 times
        # the minimum fwhm of the spectrum
        fwhm = 1.1 * wl_fwhm / wl_bin[0]
        if fwhm < 1.:
             msgs.warn("Lines are unresolved. Setting FWHM=2.pixels")
             fwhm = 2.

    if line_name is None:
        msgs.warn("No line_name as been set. The file will contain XXX as ion")
        line_name = 'XXX'

    if file_root_name is None:
        msgs.warn("No file_root_name as been set. The file will called OH_SKY_lines.dat")
        file_root_name = 'OH_SKY'

    create_linelist(wavelength, spec, fwhm=fwhm, sigdetec=sigdetec, line_name=line_name,
                    file_root_name=file_root_name, iraf_frmt=iraf_frmt, debug=debug, convert_air_to_vac=False)


def create_ThArlinelist(resolution, waveminmax=(3000.,10500.), dlam=40.0, flgd=True, thar_outfile=None,
                        fwhm=None, sigdetec=3., line_name='ThAr', file_root_name=None, iraf_frmt=False,
                        debug=False, convert_air_to_vac=True):
    """Create a syntetic ThAr spectrum at a given resolution, extract significant lines, and
    store them in a PypeIt compatibile file. This is based on the Murphy et al. ThAr spectrum.
    Detailed information are here: http://astronomy.swin.edu.au/~mmurphy/thar/index.html

    Parameters
    ----------
    resolution : float
        resolution of the spectrograph
    waveminmax : tuple
        wavelength range in ang. to be covered by the model.
        Default is: (3000., 10500.)
    dlam : 
        bin to be used to create the wavelength grid of the model.
        If flgd='True' it is a bin in velocity (km/s). If flgd='False'
        it is a bin in linear space (angstrom).
        Default is: 40.0 (with flgd='True')
    flgd : boolean
        if flgd='True' (default) wavelengths are created with 
        equal steps in log space. If 'False', wavelengths will be
        created wit equal steps in linear space.
    thar_outfile : str
        name of the fits file where the model sky spectrum will be stored.
        default is 'None' (i.e., no file will be written).
    fwhm : float
        fwhm in pixels used for filtering out arc lines that are too
        wide and not considered in fits. Parameter of arc.detect_lines().
        If set to 'None' the fwhm will be derived from the resolution as:
        2. * central_wavelength / resolution
    sigdetec : float
        sigma threshold above fluctuations for line detection. Parameter
        of arc.detect_lines(). Default = 2.
    line_name : str
        name of the lines to listed in the file.
    file_root_name : str
        name of the file where the identified lines will be stored.
        The code automatically add '_lines.dat' at the end of the
        root name.
    iraf_frmt : bool
        if True, the file is written in the IRAF format (i.e. wavelength,
        ion name, amplitude).
    debug : boolean
        If True will show debug plots
    convert_air_to_vac (bool):
        If True, convert the wavelengths of the created linelist from air to vacuum
    """

    wavelength, spec = optical_modelThAr(resolution, waveminmax=waveminmax, dlam=dlam,
                                         flgd=flgd, thar_outfile=thar_outfile, debug=debug)
    if fwhm is None:
        msgs.warn("No min FWHM for the line detection set. Derived from the resolution at the center of the spectrum")
        wl_cent = np.average(wavelength)
        wl_fwhm = wl_cent / resolution
        wl_bin = np.abs((wavelength-np.roll(wavelength,1))[np.where(np.abs(wavelength-wl_cent)==np.min(np.abs(wavelength-wl_cent)))])
        # In order not to exclude all the lines, fwhm is set to 5 times
        # the minimum fwhm of the spectrum
        fwhm = 1.1 * wl_fwhm / wl_bin[0]
        if fwhm < 1.:
             msgs.warn("Lines are unresolved. Setting FWHM=2.*pixels")
             fwhm = 2.

    if line_name is None:
        msgs.warn("No line_name as been set. The file will contain XXX as ion")
        line_name = 'XXX'

    if file_root_name is None:
        msgs.warn("No file_root_name as been set. The file will called ThAr_lines.dat")
        file_root_name = 'ThAr'

    create_linelist(wavelength, spec, fwhm=fwhm, sigdetec=sigdetec, line_name=line_name,
                    file_root_name=file_root_name, iraf_frmt=iraf_frmt, debug=debug, convert_air_to_vac=convert_air_to_vac)

