"""
Module to create models of arc lines.

.. include:: ../include/links.rst

"""

import astropy
import re
import scipy

import numpy as np
import matplotlib.pyplot as plt


from astropy.io import fits
from astropy.convolution import convolve, Gaussian1DKernel
from astropy.table import Table
from astropy import units
from astropy import constants

from pypeit.core import arc
from pypeit.core.plot import pyplot_rcparams, pyplot_rcparams_default
from pypeit.core.wave import airtovac

def blackbody(wavelength, T_BB=250., debug=False):
    """ Given wavelength [in microns] and Temperature in Kelvin
    it returns the black body emission.

    Parameters
    ----------
    wavelength : `numpy.ndarray`_
        wavelength vector in microns
    T_BB : float 
        black body temperature in Kelvin. Default is set to:
        T_BB = 250.

    Returns
    -------
    blackbody : `numpy.ndarray`_
        spectral radiance of the black body in cgs units:
        B_lambda = 2.*h*c^2/lambda^5.*(1./(exp(h*c/(lambda*k_b*T_BB))-1.)
    blackbody_counts : `numpy.ndarray`_
        Same as above but in flux density
    """

    # Define constants in cgs
    PLANCK  = constants.h.cgs.value    # erg*s
    C_LIGHT = constants.c.cgs.value    # cm/s
    K_BOLTZ = constants.k_B.cgs.value  # erg/K
    RADIAN_PER_ARCSEC = 1./3600.*np.pi/180.

    log.info("Creating BB spectrum at T={}K".format(T_BB))
    lam = wavelength / 1e4 # convert wave in cm.
    blackbody_pol = 2.*PLANCK*np.power(C_LIGHT,2) / np.power(lam,5)
    blackbody_exp = np.exp(PLANCK*C_LIGHT/(lam*K_BOLTZ*T_BB)) - 1.
    blackbody = blackbody_pol / blackbody_exp
    blackbody_counts = blackbody / (PLANCK * C_LIGHT / lam) * 1e-4 \
                 * np.power(RADIAN_PER_ARCSEC, 2.)

    if debug:
        pyplot_rcparams()
        log.info("Plot of the blackbody spectrum.")
        plt.figure()
        plt.plot(wavelength, blackbody,
                 color='navy', linestyle='-', alpha=0.8,
                 label=r'T_BB={}'.format(T_BB))
        plt.legend()
        plt.xlabel(r"Wavelength [micron]")
        plt.ylabel(r"Spectral Radiance")
        plt.title(r"Planck's law")
        log.info("Close the Figure to continue.")
        plt.show(block=True)
        plt.close()
        pyplot_rcparams_default()

    return blackbody, blackbody_counts

def addlines2spec(wavelength, wl_line, fl_line, resolution,
                  scale_spec=1., debug=False):
    """ Create a spectrum with a set of (gaussian) emission lines.
    
    Parameters
    ----------
    wavelength : `numpy.ndarray`_
        wavelength vector of the input spectrum
    wl_line : `numpy.ndarray`_
        wavelength of each individual line
    fl_line : `numpy.ndarray`_
        flux of each individual line
    resolution : float
        resolution of the spectrograph. In other words, the lines
        will have a FWHM equal to:
        fwhm_line = wl_line / resolution
    scale_spec : float
        rescale all the  normalization of the final spectrum.
        Default scale_spec=1.
    debug : bool
        If True will show debug plots

    Returns
    -------
    line_spec : `numpy.ndarray`_
        Spectrum with lines

    """
    line_spec = np.zeros_like(wavelength)
    wl_line_min, wl_line_max = np.min(wavelength), np.max(wavelength)
    good_lines = (wl_line>wl_line_min) & (wl_line<wl_line_max)
    wl_line_good = wl_line[good_lines]
    fl_line_good = fl_line[good_lines]

    # define sigma of the gaussians
    sigma = wl_line_good / resolution / 2.355

    log.info("Creating line spectrum")
    for ii in np.arange(len(wl_line_good)):
        line_spec += scale_spec*fl_line_good[ii]*\
                     np.exp(-np.power((wl_line_good[ii]-wavelength),2.)/(2.*np.power(sigma[ii],2.)))

    if debug:
        pyplot_rcparams()
        log.info("Plot of the line spectrum.")
        plt.figure()
        plt.plot(wavelength, line_spec,
                 color='navy', linestyle='-', alpha=0.8,
                 label=r'Spectrum with lines included')
        plt.legend()
        plt.xlabel(r'Wavelength')
        plt.ylabel(r'Flux')
        log.info("Close the Figure to continue.")
        plt.show(block=True)
        plt.close()
        pyplot_rcparams_default()

    return line_spec

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
    wl_bin = np.abs((wavelength - np.roll(wavelength,1))[np.where( np.abs(wavelength-wl_cent) == np.min(np.abs(wavelength-wl_cent)) )])
    log.info("The binning of the wavelength array at {} is: {}".format(wl_cent, wl_bin[0]))
    px_bin = wl_bin[0]
    px_sigma = wl_sigma / px_bin

    log.info("Covolving with a Gaussian kernel with sigma = {} pixels".format(px_sigma))
    gauss_kernel = Gaussian1DKernel(px_sigma)

    flux_convolved = convolve(flux, gauss_kernel)

    if debug:
        pyplot_rcparams()
        log.info("Spectrum Convolved at R = {}".format(resolution))
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
        log.info("Close the Figure to continue.")
        plt.show(block=True)
        plt.close()
        pyplot_rcparams_default()

    return flux_convolved, px_sigma, px_bin


def iraf_datareader(database_dir, id_file):
    """
    Reads in a line identification database created with IRAF
    identify. These are usually locate in a directory called 'database'.
    This read pixel location and wavelength of the lines that
    have been id with IRAF. Note that the first pixel in IRAF
    is '1', while is '0' in python. The pixel location is thus
    shifted of one pixel while reading the database.

    Parameters
    ----------
    database_dir : str
        directory where the id files are located.
    id_file : str
        filename that is going to be read.

    Returns
    -------
    pixel : `numpy.ndarray`_
        Position of the line in pixel of the line.  For IRAF output, these are
        usually in Ang.
    line_id : `numpy.ndarray`_
        ID of the line.
    """

    lines_database = []

    # Open file for reading of text data.
    with open (database_dir+id_file, 'rt') as id_file_iraf:
        for line in id_file_iraf:
            lines_database.append(line.split())
            feat_line = re.search(r'features\t(\d+)', line)
            if feat_line is not None:
                N_lines = int(feat_line.group(1))

    log.info("The number of IDs in the IRAF database {} is {}".format(id_file, N_lines))

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
    wavelength : `numpy.ndarray`_
        wavelength
    spec : `numpy.ndarray`_
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
    convert_air_to_vac : bool
        If True, convert the wavelengths of the created linelist from air to vacuum
    """

    log.info("Searching for peaks {} sigma above background".format(sigdetec))
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
        log.info("Converting wavelengths from air to vacuum")
        wave_peak = airtovac(wave_peak * units.AA).value

    npeak = len(wave_peak)
    ion = npeak*[str(line_name)]
    NIST = npeak*[1]
    Instr = npeak*[32]
    Source = npeak*['wavemodel.py']

    if iraf_frmt:
        log.info("Printing file in IRAF format: {}".format(file_root_name+'_iraf_lines.dat'))
        ion = np.array(ion)
        id_lines_iraf = np.vstack( (np.round(wave_peak,5), ion, np.round(ampl_good,5)) ).T
        np.savetxt(file_root_name+'_iraf_lines.dat', id_lines_iraf, fmt="%15s %6s %15s", delimiter="  ")
    else:
        log.info("Printing file: {}".format(file_root_name+'_lines.dat'))
        dat = Table([wave_peak, ion, NIST, Instr, ampl_good, Source], names=('wave', 'ion','NIST','Instr','amplitude','Source'))
        dat.write(file_root_name+'_lines.dat',format='ascii.fixed_width')


