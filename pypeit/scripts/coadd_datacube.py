#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
This script enables the user to view a 2D FITS file
and define the sky background regions interactively.
Run above the Science/ folder.
"""

import argparse
from IPython import embed

from astropy import wcs, units
from astropy.io import fits
from astropy.coordinates import AltAz, SkyCoord
import scipy.optimize as opt
from scipy.interpolate import interp1d
import numpy as np

from pypeit import msgs
from pypeit.spectrographs.util import load_spectrograph
from pypeit.core.parse import get_dnum
from pypeit.core.flux_calib import load_extinction_data, extinction_correction
from pypeit import alignframe
from pypeit import spec2dobj


def parser(options=None):

    parser = argparse.ArgumentParser(description='Read in a spec2D file and convert it to a datacube',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file', type = str, default=None, help='ascii file with list of spec2D files to combine')
    parser.add_argument('--list', default=False, help='List the extensions only?',
                        action='store_true')
    parser.add_argument('--det', default=1, type=int, help="Detector")
    parser.add_argument('-o', '--overwrite', default=False, action='store_true',
                        help='Overwrite any existing files/directories')

    return parser.parse_args() if options is None else parser.parse_args(options)


def dar_fitfunc(radec, coord_ra, coord_dec, datfit, wave, obstime, location, pressure, temperature, rel_humidity):
    (diff_ra, diff_dec) = radec
    # Generate the coordinate with atmopheric conditions
    coord_atmo = SkyCoord(coord_ra + diff_ra, coord_dec + diff_dec, unit=(units.deg, units.deg))
    coord_altaz = coord_atmo.transform_to(AltAz(obstime=obstime, location=location, obswl=wave,
                                         pressure=pressure, temperature=temperature,
                                         relative_humidity=rel_humidity))
    # Return chi-squared value
    return np.sum((np.array([coord_altaz.alt.value, coord_altaz.az.value])-datfit)**2)


def dar_correction(wave_arr, coord, obstime, location, pressure, temperature, rel_humidity):
    """Apply a differental atmospheric refraction correction to the input ra/dec.
    This implementation is based on ERFA, which is called through astropy

    wave_arr (ndarray):
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

    Returns
    -------
    ra_diff (ndarray):
        Relative RA shift at each wavelength given by `wave_arr`
    dec_diff (ndarray):
        Relative DEC shift at each wavelength given by `wave_arr`

    TODO :: There's probably going to be issues when the RA angle is either side of RA=0
    TODO :: Move this routine to the main PypeIt code
    """
    msgs.info("Performing differential atmospheric refraction correction")
    numgrid = 100
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
        res_lsq = opt.least_squares(dar_fitfunc, [0.0, 0.0], args=args)
        # Store the result
        ra_grid[ww] = res_lsq.x[0]
        dec_grid[ww] = res_lsq.x[1]

    # Generate spline of differentials
    spl_ra = interp1d(wave_grid, ra_grid, kind='cubic')
    spl_dec = interp1d(wave_grid, dec_grid, kind='cubic')

    # Evaluate the differentials at the input wave_arr
    ra_diff = spl_ra(wave_arr)
    dec_diff = spl_dec(wave_arr)

    return ra_diff, dec_diff


def generate_masterWCS(crval, cdelt, equinox=2000.0):
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

    # List only?
    if args.list:
        hdu = fits.open(args.file)
        hdu.info()
        return

    # Get a list of files for the combination
    files = open(args.file, 'r').readlines()
    combine = True if len(files) > 1 else False

    all_ra, all_dec, all_wave = np.array([]), np.array([]), np.array([])
    all_sci, all_ivar = np.array([]), np.array([])
    dspat = None  # binning size on the sky (in arcsec)
    ref_scale = None  # This will be used to correct relative scaling among the various input frames
    for fil in files:
        # Load it up
        spec2DObj = spec2dobj.Spec2DObj.from_file(fil.rstrip("\n"), args.det)

        # Load the spectrograph
        specname = spec2DObj.head0['SPECTROG']
        spec = load_spectrograph(specname)

        # Setup for PypeIt imports
        msgs.reset(verbosity=2)

        # Init
        # TODO: get_dnum needs to be deprecated...
        sdet = "{0:s}-".format(get_dnum(args.det, caps=True, prefix=True))

        if ref_scale is None:
            ref_scale = spec2DObj.scaleimg.copy()
        # EXtract the information
        sciimg = (spec2DObj.sciimg-spec2DObj.skymodel) * (ref_scale/spec2DObj.scaleimg)
        ivar = spec2DObj.ivarraw * (ref_scale/spec2DObj.scaleimg)**2
        waveimg = spec2DObj.waveimg
        bpmmask = spec2DObj.bpmmask

        # Grab the slit edges
        slits = spec2DObj.slits

        # Load the master alignments
        msgs.info("Loading alignments")
        # TODO :: Include ALGNMKEY or alignments in Spec2D
        alignfile = fil.split("Science")[0] + "Masters/MasterAlignment_{0:s}_01.fits".format(spec2DObj.head0['FLATMKEY'])
        alignments = alignframe.Alignments.from_file(alignfile)

        wave0 = waveimg[waveimg != 0.0].min()
        diff = waveimg[1:, :] - waveimg[:-1, :]
        dwv = float(np.median(diff[diff != 0.0]))
        msgs.info("Using wavelength solution: wave0={0:.3f}, dispersion={1:.3f} Angstrom/pixel".format(wave0, dwv))

        msgs.info("Constructing slit image")
        slitid_img_init = slits.slit_img(pad=0, initial=True, flexure=spec2DObj.sci_spat_flexure)
        onslit_gpm = (slitid_img_init > 0) & (bpmmask == 0)

        # Grab the WCS of this frame
        wcs = spec.get_wcs(spec2DObj.head0, slits, wave0, dwv)

        # Find the largest spatial scale of all images being combined
        pxscl, slscl = spec.get_scales(spec2DObj.head0)
        if dspat is None:
            dspat = max(pxscl, slscl)
        elif max(pxscl, slscl) > dspat:
            dspat = max(pxscl, slscl)

        # Generate an RA/DEC image
        msgs.info("Generating RA/DEC image")
        raimg, decimg, minmax = spec.get_radec_image(alignments, slits, wcs, flexure=spec2DObj.sci_spat_flexure)

        # Perform and apply the DAR correction
        if False:
            darpar = spec.get_dar_params(spec2DObj.head0)
            ra_corr, dec_corr = dar_correction(waveimg[onslit_gpm], *darpar)
            raimg[onslit_gpm] += ra_corr
            decimg[onslit_gpm] += dec_corr
        # TODO :: THE FOLLOWING IS WRONG!!!  It's a placeholder while I figure out how to put everything onto a masterWCS
        # raimg[onslit_gpm] += (ra_corr-ra_corr[0])
        # decimg[onslit_gpm] += (dec_corr-dec_corr[0])

        # Extinction correction
        if False:
            msgs.info("Applying extinction correction")
            extinct = load_extinction_data(longitude, latitude)
            ext_corr = extinction_correction(wave * units.AA, airmass, extinct)
            # Correct for extinction
            flux_star = flux_star * ext_corr
            ivar_star = ivar_star / ext_corr ** 2

        # TODO :: Need to apply relative spectral scaling so that all spectra are on the same spectral shape

        all_ra = np.append(all_ra, raimg[onslit_gpm].copy())
        all_dec = np.append(all_dec, decimg[onslit_gpm].copy())
        all_wave = np.append(all_wave, waveimg[onslit_gpm].copy())
        all_sci = np.append(all_sci, sciimg[onslit_gpm].copy())
        all_ivar = np.append(all_ivar, ivar[onslit_gpm].copy())

    # Generate a master WCS to register all frames
    coord_min = [np.min(all_ra), np.min(all_dec), np.min(all_wave)]
    coord_dlt = [dspat, dspat, dwv]
    masterwcs = generate_masterWCS(coord_min, coord_dlt)

    # Generate the output binning
    if combine:
        numra = int((np.max(all_ra)-np.min(all_ra)) * np.cos(np.mean(all_dec)*np.pi/180.0) / dspat)
        numdec = int((np.max(all_dec)-np.min(all_dec))/dspat)
        numwav = int((np.max(all_wave)-np.min(all_wave))/dwv)
        xbins = np.arange(1+numra)-0.5
        ybins = np.arange(1+numdec)-0.5
        zbins = np.arange(1+numwav)-0.5
    else:
        slitlength = int(np.round(np.median(slits.get_slitlengths(initial=True, median=True))))
        xbins = np.arange(1 + 24) - 12.0 - 0.5
        ybins = np.linspace(np.min(minmax[:, 0]), np.max(minmax[:, 1]), 1+slitlength) - 0.5
        zbins = np.arange(1+int(round((np.max(waveimg)-wave0)/dwv))) - 0.5

    # Make the cube
    msgs.info("Generating datacube")
    # TODO :: Need to include variance weights and make a variance cube
    if combine:
        pix_coord = masterwcs.wcs_world2pix(all_ra, all_dec, all_wave * 1.0E-10, 0)
        hdr = masterwcs.to_header()
    else:
        pix_coord = wcs.wcs_world2pix(np.vstack((all_ra, all_dec, all_wave*1.0E-10)).T, 0)
        hdr = wcs.to_header()
    #datacube_resid, edges = np.histogramdd(pix_coord, bins=(xbins, ybins, zbins), weights=all_sci*np.sqrt(all_ivar))
    datacube, edges = np.histogramdd(pix_coord, bins=(xbins, ybins, zbins), weights=all_sci*all_ivar)
    norm, edges = np.histogramdd(pix_coord, bins=(xbins, ybins, zbins), weights=all_ivar)
    scaleArr = (norm > 0)/(norm + (norm == 0))

    # Save the datacube
    if False:
        outfile = "datacube_resid.fits"
        msgs.info("Saving datacube as: {0:s}".format(outfile))
        hdu = fits.PrimaryHDU(datacube_resid*norm, header=masterwcs.to_header())
#        hdu = fits.PrimaryHDU(norm, header=masterwcs.to_header())
        hdu.writeto(outfile, overwrite=args.overwrite)

    outfile = "datacube.fits"
    msgs.info("Saving datacube as: {0:s}".format(outfile))
    hdu = fits.PrimaryHDU((datacube*scaleArr).T, header=hdr)
    hdu.writeto(outfile, overwrite=args.overwrite)
