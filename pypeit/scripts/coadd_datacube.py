#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
This script enables the user to convert spec2D FITS files
from IFU instruments into a 3D cube with a defined WCS.
"""

import argparse
from IPython import embed

from astropy import wcs, units
from astropy.io import fits
from astropy.coordinates import AltAz, SkyCoord
import scipy.optimize as opt
from scipy.interpolate import interp1d
from scipy import stats
import numpy as np
import copy

from pypeit import msgs
from pypeit.spectrographs.util import load_spectrograph
from pypeit.core.flux_calib import load_extinction_data, extinction_correction
from pypeit.core.procimg import grow_masked
from pypeit.core.flexure import calculate_image_offset
from pypeit.core import parse
from pypeit.core import fitting
from pypeit import spec2dobj


def parser(options=None):

    parser = argparse.ArgumentParser(description='Read in a spec2D file and convert it to a datacube',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file', type = str, default=None, help='ascii file with list of spec2D files to combine')
    parser.add_argument('--det', default=1, type=int, help="Detector")
    parser.add_argument('-o', '--overwrite', default=False, action='store_true',
                        help='Overwrite any existing files/directories')

    return parser.parse_args() if options is None else parser.parse_args(options)


def dar_fitfunc(radec, coord_ra, coord_dec, datfit, wave, obstime, location, pressure, temperature, rel_humidity):
    """ Generates a fitting function to calculate the offset due to differential atmospheric refraction

    Args:
        radec (tuple):
            A tuple containing two floats representing the shift in ra and dec due to DAR.
        coord_ra (float):
            RA in degrees
        coord_dec (float):
            Dec in degrees
        datfit (`numpy.ndarray`_):
            The RA and DEC that the model needs to match
        wave (float):
            Wavelength to calculate the DAR
        location (astropy EarthLocation):
            observatory location
        pressure (float):
            Outside pressure at `location`
        temperature (float):
            Outside ambient air temperature at `location`
        rel_humidity (float):
            Outside relative humidity at `location`. This should be between 0 to 1.

    Returns:
        chisq (float):
            chi-squared difference between datfit and model
    """
    (diff_ra, diff_dec) = radec
    # Generate the coordinate with atmopheric conditions
    coord_atmo = SkyCoord(coord_ra + diff_ra, coord_dec + diff_dec, unit=(units.deg, units.deg))
    coord_altaz = coord_atmo.transform_to(AltAz(obstime=obstime, location=location, obswl=wave,
                                         pressure=pressure, temperature=temperature,
                                         relative_humidity=rel_humidity))
    # Return chi-squared value
    return np.sum((np.array([coord_altaz.alt.value, coord_altaz.az.value])-datfit)**2)


def dar_correction(wave_arr, coord, obstime, location, pressure, temperature, rel_humidity,
                   wave_ref=None, numgrid = 10):
    """Apply a differental atmospheric refraction correction to the input ra/dec.
    This implementation is based on ERFA, which is called through astropy

    Args:
        wave_arr (`numpy.ndarray`_):
            wavelengths to obtain ra and dec offsets
        coord (astropy SkyCoord):
            ra, dec positions at the centre of the field
        obstime (astropy Time):
            time at the midpoint of observation
        location (astropy EarthLocation):
            observatory location
        pressure (float):
            Outside pressure at `location`
        temperature (float):
            Outside ambient air temperature at `location`
        rel_humidity (float):
            Outside relative humidity at `location`. This should be between 0 to 1.
        wave_ref (float):
            Reference wavelength (The DAR correction will be performed relative to this wavelength)
        numgrid (int):
            Number of grid points to evaluate the DAR correction.

    Returns:
        ra_diff (`numpy.ndarray`_):
            Relative RA shift at each wavelength given by `wave_arr`
        dec_diff (`numpy.ndarray`_):
            Relative DEC shift at each wavelength given by `wave_arr`

    TODO :: There's probably going to be issues when the RA angle is either side of RA=0
    TODO :: Move this routine to the main PypeIt code?
    """
    msgs.info("Performing differential atmospheric refraction correction")

    if wave_ref is None:
        wave_ref = 0.5*(wave_arr.min() + wave_arr.max())

    # First create the reference frame and wavelength grid
    coord_altaz = coord.transform_to(AltAz(obstime=obstime, location=location))
    wave_grid = np.linspace(wave_arr.min(), wave_arr.max(), numgrid) * units.AA
    # Prepare the fit
    ra_grid, dec_grid = np.zeros(numgrid), np.zeros(numgrid)
    datfit = np.array([coord_altaz.alt.value, coord_altaz.az.value])
    # Loop through all wavelengths
    for ww in range(numgrid):
        # Fit the differential
        args = (coord.ra.value, coord.dec.value, datfit, wave_grid[ww], obstime, location, pressure, temperature, rel_humidity)
        #b_popt, b_pcov = opt.curve_fit(dar_fitfunc, tmp, datfit, p0=(0.0, 0.0))
        res_lsq = opt.least_squares(dar_fitfunc, [0.0, 0.0], args=args, xtol=1.0e-6, ftol=None, gtol=None)
        # Store the result
        ra_grid[ww] = res_lsq.x[0]
        dec_grid[ww] = res_lsq.x[1]

    # Generate spline of differentials
    spl_ra = interp1d(wave_grid, ra_grid, kind='cubic')
    spl_dec = interp1d(wave_grid, dec_grid, kind='cubic')

    # Evaluate the differentials at the input wave_arr
    ra_diff = spl_ra(wave_arr) - spl_ra(wave_ref)
    dec_diff = spl_dec(wave_arr) - spl_dec(wave_ref)

    return ra_diff, dec_diff


def generate_whiteLightImgs(all_ra, all_dec, all_wave, all_sci, all_wghts, all_idx,
                            dspat, numfiles=None, all_ivar=None):
    """ Generate a whitelight image of every input frame

    Args:
        all_ra (`numpy.ndarray`_):
            RA values of each pixel from all spec2d files
        all_dec (`numpy.ndarray`_):
            DEC values of each pixel from all spec2d files
        all_wave (`numpy.ndarray`_):
            Wavelength values of each pixel from all spec2d files
        all_sci (`numpy.ndarray`_):
            Counts of each pixel from all spec2d files
        all_wghts (`numpy.ndarray`_):
            Weights attributed to each pixel from all spec2d files
        all_idx (`numpy.ndarray`_):
            An index indicating which spec2d file each pixel originates from
        dspat (float):
            The size of each spaxel on the sky (in degrees)
        numfiles (int, optional):
            Number of spec2d files included. If not provided, it will be calculated from all_idx
        all_ivar (`numpy.ndarray`_, optional):
            Inverse variance of each pixel from all spec2d files. If provided,
            inverse variance images will be calculated and return for each white light image.

    Returns:
        tuple : two 3D arrays will be returned, each of shape [N, M, numfiles],
        where N and M are the spatial dimensions of the combined white light images.
        The first array is a white light image, and the second array is the corresponding
        inverse variance image. If all_ivar is None, this will be an empty array.
    """
    # Determine number of files
    if numfiles is None:
        numfiles = np.unique(all_idx).size

    # Generate coordinates
    cosdec = np.cos(np.mean(all_dec) * np.pi / 180.0)
    numra = int((np.max(all_ra) - np.min(all_ra)) * cosdec / dspat)
    numdec = int((np.max(all_dec) - np.min(all_dec)) / dspat)
    xbins = np.arange(1 + numra) - 1
    ybins = np.arange(1 + numdec) - 1
    spec_bins = np.arange(2) - 1
    bins = (xbins, ybins, spec_bins)

    # Generate a master 2D WCS to register all frames
    coord_min = [np.min(all_ra), np.min(all_dec), np.min(all_wave)]
    coord_dlt = [dspat, dspat, np.max(all_wave) - np.min(all_wave)]
    whitelightWCS = generate_masterWCS(coord_min, coord_dlt)

    whitelight_Imgs = np.zeros((numra, numdec, numfiles))
    whitelight_ivar = np.zeros((numra, numdec, numfiles))
    trim = 3
    for ff in range(numfiles):
        msgs.info("Generating white light image of frame {0:d}/{1:d}".format(ff + 1, numfiles))
        ww = (all_idx == ff)
        # Make the cube
        pix_coord = whitelightWCS.wcs_world2pix(np.vstack((all_ra[ww], all_dec[ww], all_wave[ww] * 1.0E-10)).T, 0)
        wlcube, edges = np.histogramdd(pix_coord, bins=bins, weights=all_sci[ww] * all_wghts[ww])
        norm, edges = np.histogramdd(pix_coord, bins=bins, weights=all_wghts[ww])
        nrmCube = (norm > 0) / (norm + (norm == 0))
        whtlght = (wlcube * nrmCube)[:, :, 0]
        # Create a mask of good pixels (trim the edges)
        msk = grow_masked(whtlght == 0, trim, 1) == 0
        whtlght *= msk
        # Set the masked regions to the minimum value
        minval = np.min(whtlght[msk == 1])
        whtlght[msk == 0] = minval
        # Store the white light image
        whitelight_Imgs[:, :, ff] = whtlght.copy()
        # Now operate on the inverse variance image
        if all_ivar is not None:
            ivar_img, _ = np.histogramdd(pix_coord, bins=bins, weights=all_ivar[ww])
            ivar_img = ivar_img[:, :, 0]
            ivar_img *= msk
            minval = np.min(ivar_img[msk == 1])
            ivar_img[msk == 0] = minval
            whitelight_ivar[:, :, ff] = ivar_img.copy()
    return whitelight_Imgs, whitelight_ivar


def generate_masterWCS(crval, cdelt, equinox=2000.0):
    """ Generate a WCS that will cover all input spec2D files

    Args:
        crval (list):
            3 element list containing the [RA, DEC, WAVELENGTH] of the reference pixel
        cdelt (list):
            3 element list containing the delta values of the [RA, DEC, WAVELENGTH]
        equinox (float):
            Equinox of the WCS

    Returns:
        w (`astropy.WCS`_):
            astropy WCS to be used for the combined cube
    """
    # Create a new WCS object.
    msgs.info("Generating Master WCS")
    w = wcs.WCS(naxis=3)
    w.wcs.equinox = equinox
    w.wcs.name = 'KCWI'
    w.wcs.radesys = 'FK5'
    # Insert the coordinate frame
    w.wcs.cname = ['KCWI RA', 'KCWI DEC', 'KCWI Wavelength']
    w.wcs.cunit = [units.degree, units.degree, units.Angstrom]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN", "AWAV"]
    w.wcs.crval = crval  # RA, DEC, and wavelength zeropoints
    w.wcs.crpix = [0, 0, 0]  # RA, DEC, and wavelength reference pixels
    #w.wcs.cd = np.array([[cdval[0], 0.0, 0.0], [0.0, cdval[1], 0.0], [0.0, 0.0, cdval[2]]])
    w.wcs.cdelt = cdelt
    w.wcs.lonpole = 180.0  # Native longitude of the Celestial pole
    w.wcs.latpole = 0.0  # Native latitude of the Celestial pole
    return w


def calculate_spectral_weights(all_ra, all_dec, all_wave, all_sci, all_ivar, all_wghts, all_idx, dspat, dwv, numfiles=None):
    # Determine number of files
    if numfiles is None:
        numfiles = np.unique(all_idx).size

    # Generate a white light image of *all* data
    msgs.info("Generating global white light image")
    whitelight_Img, ivar = generate_whiteLightImgs(all_ra, all_dec, all_wave, all_sci, all_wghts, np.zeros(all_ra.size),
                                                dspat, numfiles=1)
    whitelight_Img = whitelight_Img[:, :, 0]
    # Find the location of the object with the highest S/N in the combined white light image
    idx_max = np.unravel_index(np.argmax(whitelight_Img), whitelight_Img.shape)
    msgs.info("Highest S/N object located at spaxel (x, y) = {0:d}, {1:d}".format(idx_max[0], idx_max[1]))

    # Generate a master 2D WCS to register all frames
    coord_min = [np.min(all_ra), np.min(all_dec), np.min(all_wave)]
    coord_dlt = [dspat, dspat, dwv]
    whitelightWCS = generate_masterWCS(coord_min, coord_dlt)
    # Make the bin edges to be at +/- 1 pixels around the maximum (i.e. summing 9 pixels total)
    numwav = int((np.max(all_wave) - np.min(all_wave)) / dwv)
    xbins = np.array([idx_max[0]-1, idx_max[0]+2]) - 0.5
    ybins = np.array([idx_max[1]-1, idx_max[1]+2]) - 0.5
    spec_bins = np.arange(1 + numwav) - 0.5
    bins = (xbins, ybins, spec_bins)

    # Extract the spectrum of the highest S/N object
    all_spec = np.zeros((numwav, numfiles))
    all_snr = np.zeros((numwav, numfiles))
    for ff in range(numfiles):
        msgs.info("Extracting spectrum of highest S/N detection from frame {0:d}/{1:d}".format(ff + 1, numfiles))
        ww = (all_idx == ff)
        # Extract the spectrum
        pix_coord = whitelightWCS.wcs_world2pix(np.vstack((all_ra[ww], all_dec[ww], all_wave[ww] * 1.0E-10)).T, 0)
        spec, edges = np.histogramdd(pix_coord, bins=bins, weights=all_sci[ww])
        norm, edges = np.histogramdd(pix_coord, bins=bins)
        ivar, edges = np.histogramdd(pix_coord, bins=bins, weights=all_ivar[ww])
        nrmSpec = (norm > 0) / (norm + (norm == 0))
        spec = (spec*nrmSpec)[0, 0, :]
        all_spec[:, ff] = spec.copy()
        all_snr[:, ff] = spec*np.sqrt(ivar[0, 0, :])

    # Construct the relative weights based on the S/N as a function of wavelength
    # Obtain a wavelength of each pixel
    wcs_res = whitelightWCS.wcs_pix2world(np.vstack((np.zeros(numwav), np.zeros(numwav), np.arange(numwav))).T, 0)
    wave_spec = wcs_res[:, 2] * 1.0E10
    # Identify the reference spectrum (the one with the highest typical S/N ratio)
    medsnr = np.median(all_snr, axis=0)
    ref_idx = np.argmax(medsnr)
    msgs.info("The reference spectrum (frame {0:d}) has a typical S/N = {1:.3f}".format(ref_idx, medsnr[ref_idx]))
    all_wghts = np.ones(all_sci.size)
    for ff in range(numfiles):
        ww = (all_idx == ff)
        # The reference spectrum will have unit weights
        if ff == ref_idx:
            all_wghts[ww] = 1.0
            continue
        # Calculate the relative weights according the S/N
        relsnr = (all_snr[:, ff]/all_snr[:, ref_idx])**2
        # Perform a low order polynomial fit
        wght_fit = fitting.robust_fit(wave_spec, relsnr, 3, function="legendre",
                                      minx=wave_spec.min(), maxx=wave_spec.max(),
                                      lower=5, upper=5)
        # Apply fitting function to all wavelengths
        all_wghts[ww] = wght_fit.eval(all_wave[ww])
    return all_wghts


def main(args):
    # Get a list of files for the combination
    files = open(args.file, 'r').readlines()
    filelist = []
    for fil in files:
        filelist.append(fil.rstrip("\n"))

    # Coadd the files
    coadd_cube(filelist, det=args.det, overwrite=args.overwrite)


def coadd_cube(files, det=1, overwrite=False):
    """ Main routine to coadd spec2D files

    Args:
        files (list):
            List of all spec2D files
        det (int):
            detector
        overwrite (bool):
            Overwrite the output file, if it exists?
    """
    # prep
    numfiles = len(files)
    combine = True if numfiles > 1 else False

    all_ra, all_dec, all_wave = np.array([]), np.array([]), np.array([])
    all_sci, all_ivar, all_idx, all_wghts = np.array([]), np.array([]), np.array([]), np.array([])
    all_wcs = []
    dspat = None  # binning size on the sky (in arcsec)
    ref_scale = None  # This will be used to correct relative scaling among the various input frames
    wave_ref = None
    weights = np.ones(numfiles)  # Weights to use when combining cubes
    for ff, fil in enumerate(files):
        # Load it up
        spec2DObj = spec2dobj.Spec2DObj.from_file(fil, det)

        # Load the spectrograph
        specname = spec2DObj.head0['SPECTROG']
        spec = load_spectrograph(specname)
        detector = spec2DObj.detector

        # Setup for PypeIt imports
        msgs.reset(verbosity=2)

        if ref_scale is None:
            ref_scale = spec2DObj.scaleimg.copy()
        # EXtract the information
        sciimg = (spec2DObj.sciimg-spec2DObj.skymodel) * (ref_scale/spec2DObj.scaleimg)  # Subtract sky and apply relative sky
        ivar = spec2DObj.ivarraw * (ref_scale/spec2DObj.scaleimg)**2
        waveimg = spec2DObj.waveimg
        bpmmask = spec2DObj.bpmmask

        # Grab the slit edges
        slits = spec2DObj.slits

        wave0 = waveimg[waveimg != 0.0].min()
        diff = waveimg[1:, :] - waveimg[:-1, :]
        dwv = float(np.median(diff[diff != 0.0]))
        msgs.info("Using wavelength solution: wave0={0:.3f}, dispersion={1:.3f} Angstrom/pixel".format(wave0, dwv))

        msgs.info("Constructing slit image")
        slitid_img_init = slits.slit_img(pad=0, initial=True, flexure=spec2DObj.sci_spat_flexure)
        onslit_gpm = (slitid_img_init > 0) & (bpmmask == 0)

        # Grab the WCS of this frame
        wcs = spec.get_wcs(spec2DObj.head0, slits, detector.platescale, wave0, dwv)
        all_wcs.append(copy.deepcopy(wcs))

        # Find the largest spatial scale of all images being combined
        # TODO :: probably need to put this in the DetectorContainer
        pxscl = detector.platescale * parse.parse_binning(detector.binning)[1] / 3600.0  # This should be degrees/pixel
        slscl = spec.get_meta_value([spec2DObj.head0], 'slitwid')
        if dspat is None:
            dspat = max(pxscl, slscl)
        elif max(pxscl, slscl) > dspat:
            dspat = max(pxscl, slscl)

        # Generate an RA/DEC image
        msgs.info("Generating RA/DEC image")
        raimg, decimg, minmax = slits.get_radec_image(wcs, initial=True, flexure=spec2DObj.sci_spat_flexure)

        # Perform the DAR correction
        if wave_ref is None:
            wave_ref = 0.5*(np.min(waveimg[onslit_gpm]) + np.max(waveimg[onslit_gpm]))
        # Get DAR parameters
        raval = spec.get_meta_value([spec2DObj.head0], 'ra')
        decval = spec.get_meta_value([spec2DObj.head0], 'dec')
        obstime = spec.get_meta_value([spec2DObj.head0], 'obstime')
        pressure = spec.get_meta_value([spec2DObj.head0], 'pressure')
        temperature = spec.get_meta_value([spec2DObj.head0], 'temperature')
        rel_humidity = spec.get_meta_value([spec2DObj.head0], 'humidity')
        coord = SkyCoord(raval, decval, unit=(units.deg, units.deg))
        location = spec.location  # TODO :: spec.location should probably end up in the TelescopePar (spec.telescope.location)
        ra_corr, dec_corr = dar_correction(waveimg[onslit_gpm], coord, obstime, location,
                                           pressure, temperature, rel_humidity, wave_ref=wave_ref)
        raimg[onslit_gpm] += ra_corr
        decimg[onslit_gpm] += dec_corr

        # Get copies of arrays to be saved
        wave_ext = waveimg[onslit_gpm].copy()
        flux_ext = sciimg[onslit_gpm].copy()
        ivar_ext = ivar[onslit_gpm].copy()

        # Perform extinction correction
        msgs.info("Applying extinction correction")
        longitude = spec.telescope['longitude']
        latitude = spec.telescope['latitude']
        airmass = spec2DObj.head0[spec.meta['airmass']['card']]
        extinct = load_extinction_data(longitude, latitude)
        # extinction_correction requires the wavelength is sorted
        wvsrt = np.argsort(wave_ext)
        ext_corr = extinction_correction(wave_ext[wvsrt] * units.AA, airmass, extinct)
        # Correct for extinction
        flux_sav = flux_ext[wvsrt] * ext_corr
        ivar_sav = ivar_ext[wvsrt] / ext_corr ** 2
        # sort back to the original ordering
        resrt = np.argsort(wvsrt)

        # Calculate the weights relative to the zeroth cube
        if ff != 0:
            weights[ff] = np.median(flux_sav[resrt]*np.sqrt(ivar_sav[resrt]))**2

        # Store the information
        numpix = raimg[onslit_gpm].size
        all_ra = np.append(all_ra, raimg[onslit_gpm].copy())
        all_dec = np.append(all_dec, decimg[onslit_gpm].copy())
        all_wave = np.append(all_wave, wave_ext.copy())
        all_sci = np.append(all_sci, flux_sav[resrt].copy())
        all_ivar = np.append(all_ivar, ivar_sav[resrt].copy())
        all_idx = np.append(all_idx, ff*np.ones(numpix))
        all_wghts = np.append(all_wghts, weights[ff]*np.ones(numpix))

    # Grab cos(dec) for convenience
    cosdec = np.cos(np.mean(all_dec) * np.pi / 180.0)

    # Register spatial offsets between all frames if several frames are being combined
    if combine:
        # Generate white light images
        whitelight_Imgs, _ = generate_whiteLightImgs(all_ra, all_dec, all_wave, all_sci, all_wghts, all_idx,
                                                     dspat, numfiles=numfiles)

        # ref_idx will be the index of the cube with the highest S/N
        ref_idx = np.argmax(weights)
        msgs.info("Calculating the relative spatial translation of each cube (reference cube = {0:d})".format(ref_idx+1))
        # Calculate the image offsets - check the reference is a zero shift
        ra_shift_ref, dec_shift_ref = calculate_image_offset(whitelight_Imgs[:, :, ref_idx], whitelight_Imgs[:, :, ref_idx])
        for ff in range(numfiles):
            # Don't correlate the reference image with itself
            if ff == ref_idx:
                continue
            # Calculate the shift
            ra_shift, dec_shift = calculate_image_offset(whitelight_Imgs[:, :, ff], whitelight_Imgs[:, :, ref_idx])
            # Convert to reference
            ra_shift -= ra_shift_ref
            dec_shift -= dec_shift_ref
            # Convert pixel shift to degress shift
            ra_shift *= dspat/cosdec
            dec_shift *= dspat
            msgs.info("Image shift of cube {0:d}: RA, DEC (arcsec) = {1:+0.3f}, {2:+0.3f}".format(ff+1, ra_shift*3600.0, dec_shift*3600.0))
            # Apply the shift
            all_ra[all_idx == ff] += ra_shift
            all_dec[all_idx == ff] += dec_shift

        # msgs.info("Recomputing white light images")
        # # Generate white light images for each spec2D file
        # whitelight_Imgs, whitelight_ivar = generate_whiteLightImgs(all_ra, all_dec, all_wave, all_sci, all_wghts,
        #                                                            all_idx, dspat, numfiles=numfiles, all_ivar=all_ivar)
        # Calculate the relative spectral weights of all pixels
        all_wghts = calculate_spectral_weights(all_ra, all_dec, all_wave, all_sci, all_ivar, all_wghts, all_idx,
                                               dspat, dwv, numfiles=numfiles)
        # TODO :: Recalculate the white light images to get the relative weights of the images
        # TODO :: SOME NOTES...
        # Probably don't need to recalculate all individual white light images.
        # Here is a better process:
        # (1) generate a global white light image
        # (2) Find spaxel with maximum S/N ratio
        # (3) Histogram each spec2d file, centred around the spaxel of maximum S/N
        #     (i.e. +/-1 pixel in each direction = 9 pixels total), but, perform
        #     the wavelength binning so you get a spectrum of the brightest object.
        # (4) Calculate the S/N of each pixel in the spectrum, relative to a reference spectrum
        # (5) Find a way to apply these weights to all pixels within a given wavelength range.
        #     Could do it by brute force, but this is surely going to be very slow... maybe jit it?

    # Generate a master WCS to register all frames
    coord_min = [np.min(all_ra), np.min(all_dec), np.min(all_wave)]
    coord_dlt = [dspat, dspat, dwv]
    masterwcs = generate_masterWCS(coord_min, coord_dlt)

    # Generate the output binning
    if combine:
        numra = int((np.max(all_ra)-np.min(all_ra)) * cosdec / dspat)
        numdec = int((np.max(all_dec)-np.min(all_dec))/dspat)
        numwav = int((np.max(all_wave)-np.min(all_wave))/dwv)
        xbins = np.arange(1+numra)-0.5
        ybins = np.arange(1+numdec)-0.5
        spec_bins = np.arange(1+numwav)-0.5
    else:
        # TODO :: This is KCWI specific - probably should put this in the spectrograph file, or just delete it.
        msgs.warn("This routine is Keck/KCWI specific")
        slitlength = int(np.round(np.median(slits.get_slitlengths(initial=True, median=True))))
        xbins = np.arange(1 + 24) - 12.0 - 0.5
        ybins = np.linspace(np.min(minmax[:, 0]), np.max(minmax[:, 1]), 1+slitlength) - 0.5
        spec_bins = np.arange(1+int(round((np.max(waveimg)-wave0)/dwv))) - 0.5

    # Make the cube
    msgs.info("Generating datacube")
    # TODO :: Need to include variance weights and make a variance cube
    if combine:
        pix_coord = masterwcs.wcs_world2pix(all_ra, all_dec, all_wave * 1.0E-10, 0)
        hdr = masterwcs.to_header()
    else:
        pix_coord = wcs.wcs_world2pix(np.vstack((all_ra, all_dec, all_wave*1.0E-10)).T, 0)
        hdr = wcs.to_header()

    # Find the NGP coordinates for all input pixels
    bins = (xbins, ybins, spec_bins)
    datacube, edges = np.histogramdd(pix_coord, bins=bins, weights=all_sci*all_wghts)
    norm, edges = np.histogramdd(pix_coord, bins=bins, weights=all_wghts)
    ivarcube, edges = np.histogramdd(pix_coord, bins=bins, weights=all_ivar)
    normCube = (norm > 0) / (norm + (norm == 0))
    varCube = (ivarcube > 0) / (ivarcube + (ivarcube == 0))

    # Save the datacube
    debug = False
    if debug:
        datacube_resid, edges = np.histogramdd(pix_coord, bins=(xbins, ybins, spec_bins), weights=all_sci*np.sqrt(all_ivar))
        norm, edges = np.histogramdd(pix_coord, bins=(xbins, ybins, spec_bins))
        normCube = (norm > 0) / (norm + (norm == 0))
        outfile = "datacube_resid.fits"
        msgs.info("Saving datacube as: {0:s}".format(outfile))
        hdu = fits.PrimaryHDU((datacube_resid*normCube).T, header=masterwcs.to_header())
        hdu.writeto(outfile, overwrite=overwrite)

    outfile = "datacube.fits"
    msgs.info("Saving datacube as: {0:s}".format(outfile))
    primary_hdu = fits.PrimaryHDU(header=spec2DObj.head0)
    sci_hdu = fits.ImageHDU((datacube*normCube).T, name="scicube", header=hdr)
    var_hdu = fits.ImageHDU(varCube.T, name="varcube", header=hdr)
    hdulist = fits.HDUList([primary_hdu, sci_hdu, var_hdu])
    hdulist.writeto(outfile, overwrite=overwrite)
