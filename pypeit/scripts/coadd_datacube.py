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
import numpy as np
import copy

from pypeit import msgs
from pypeit.spectrographs.util import load_spectrograph
from pypeit.core.flux_calib import load_extinction_data, extinction_correction
from pypeit.core.procimg import grow_masked
from pypeit.core.flexure import calculate_image_offset
from pypeit.core import parse
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

    # If several frames are being combined, generate white light images to register the offsets
    if combine:
        numra = int((np.max(all_ra)-np.min(all_ra)) * np.cos(np.mean(all_dec)*np.pi/180.0) / dspat)
        numdec = int((np.max(all_dec)-np.min(all_dec))/dspat)
        xbins = np.arange(1+numra)-1
        ybins = np.arange(1+numdec)-1
        spec_bins = np.arange(2)-1

        # Generate a master 2D WCS to register all frames
        coord_min = [np.min(all_ra), np.min(all_dec), np.min(all_wave)]
        coord_dlt = [dspat, dspat, np.max(all_wave)-np.min(all_wave)]
        whitelightWCS = generate_masterWCS(coord_min, coord_dlt)

        # Register spatial offsets between all frames
        whitelight_Imgs = np.zeros((numra, numdec, numfiles))
        trim = 3
        for ff in range(numfiles):
            msgs.info("Generating white light image of frame {0:d}/{1:d}".format(ff+1, numfiles))
            ww = (all_idx == ff)
            # Make the cube
            pix_coord = whitelightWCS.wcs_world2pix(np.vstack((all_ra[ww], all_dec[ww], all_wave[ww] * 1.0E-10)).T, 0)
            wlcube, edges = np.histogramdd(pix_coord, bins=(xbins, ybins, spec_bins), weights=all_sci[ww] * all_wghts[ww])
            norm, edges = np.histogramdd(pix_coord, bins=(xbins, ybins, spec_bins), weights=all_wghts[ww])
            nrmCube = (norm > 0) / (norm + (norm == 0))
            whtlght = (wlcube * nrmCube)[:, :, 0]
            # Create a mask of good pixels (trim the edges)
            msk = grow_masked(whtlght == 0, trim, 1) == 0
            whtlght *= msk
            # Set the masked regions to the median value
            minval = np.min(whtlght[msk == 1])
            whtlght[msk == 0] = minval
            whitelight_Imgs[:, :, ff] = whtlght.copy()

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

        # Recalculate the white light images to get the relative weights of the images
        weights = np.ones(numfiles)

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
    datacube, edges = np.histogramdd(pix_coord, bins=(xbins, ybins, spec_bins), weights=all_sci*all_wghts)
    norm, edges = np.histogramdd(pix_coord, bins=(xbins, ybins, spec_bins), weights=all_wghts)
    ivarcube, edges = np.histogramdd(pix_coord, bins=(xbins, ybins, spec_bins), weights=all_ivar)
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
