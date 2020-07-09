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

from astropy.io import fits
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

    parser.add_argument('file', type = str, default=None, help='PypeIt file')
    parser.add_argument('--list', default=False, help='List the extensions only?',
                        action='store_true')
    parser.add_argument('--det', default=1, type=int, help="Detector")
    parser.add_argument('-o', '--overwrite', default=False, action='store_true',
                        help='Overwrite any existing files/directories')

    return parser.parse_args() if options is None else parser.parse_args(options)

from astropy.coordinates import AltAz, SkyCoord
import astropy.units as units
import scipy.optimize as opt
from scipy.interpolate import interp1d


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


def main(args):

    # List only?
    if args.list:
        hdu = fits.open(args.file)
        hdu.info()
        return

    # Load it up
    spec2DObj = spec2dobj.Spec2DObj.from_file(args.file, args.det)

    # Load the spectrograph
    specname = spec2DObj.head0['SPECTROG']
    spec = load_spectrograph(specname)

    # Setup for PypeIt imports
    msgs.reset(verbosity=2)

    # Init
    # TODO: get_dnum needs to be deprecated...
    sdet = "{0:s}-".format(get_dnum(args.det, caps=True, prefix=True))

    # EXtract the information
    sciimg = spec2DObj.sciimg-spec2DObj.skymodel
    ivar = spec2DObj.ivarraw#model
    waveimg = spec2DObj.waveimg

    # Grab the slit edges
    slits = spec2DObj.slits

    # Load the master alignments
    msgs.info("Loading alignments")
    # TODO :: Include ALGNMKEY or alignments in Spec2D
    alignfile = "Masters/MasterAlignment_{0:s}_01.fits".format(spec2DObj.head0['FLATMKEY'])
    alignments = alignframe.Alignments.from_file(alignfile)

    wave0 = waveimg[waveimg != 0.0].min()
    diff = waveimg[1:, :] - waveimg[:-1, :]
    dwv = float(np.median(diff[diff != 0.0]))
    msgs.info("Using wavelength solution: wave0={0:.3f}, dispersion={1:.3f} Angstrom/pixel".format(wave0, dwv))

    msgs.info("Constructing slit image")
    slitid_img_init = slits.slit_img(pad=0, initial=True, flexure=spec2DObj.sci_spat_flexure)
    onslit = slitid_img_init > 0

    # Generate a master WCS to register all frames
    masterwcs = spec.get_wcs(spec2DObj.head0, slits, wave0, dwv)

    # Grab the WCS of this frame
    wcs = spec.get_wcs(spec2DObj.head0, slits, wave0, dwv)

    # Generate an RA/DEC image
    raimg, decimg, minmax = spec.get_radec_image(alignments, slits, wcs, flexure=spec2DObj.sci_spat_flexure)

    # Perform and apply the DAR correction
    darpar = spec.get_dar_params(spec2DObj.head0)
    ra_corr, dec_corr = dar_correction(waveimg[onslit], *darpar)
    # TODO :: THE FOLLOWING IS WRONG!!! It should just be:
    # raimg[onslit] += ra_corr
    # decimg[onslit] += dec_corr
    # but we first need to generate a master WCS that takes this effect into account.
    raimg[onslit] += (ra_corr-ra_corr[0])
    decimg[onslit] += (dec_corr-dec_corr[0])

    # Generate the output binning
    slitlength = int(np.round(np.median(slits.get_slitlengths(initial=True, median=True))))
    xbins = np.arange(1 + 24) - 12.0 - 0.5
    ybins = np.linspace(np.min(minmax[:, 0]), np.max(minmax[:, 1]), 1+slitlength) - 0.5
    zbins = np.arange(1+int(round((np.max(waveimg)-wave0)/dwv))) - 0.5

    # Make the cube
    msgs.info("Generating datacube")
    pix_coord = wcs.wcs_world2pix(raimg[onslit], decimg[onslit], waveimg[onslit]*1.0E-10, 0)
    datacube_resid, edges = np.histogramdd(pix_coord, bins=(xbins, ybins, zbins), weights=sciimg[onslit]*np.sqrt(ivar[onslit]))
    datacube, edges = np.histogramdd(pix_coord, bins=(xbins, ybins, zbins), weights=sciimg[onslit])
    norm, edges = np.histogramdd(pix_coord, bins=(xbins, ybins, zbins))
    norm += (norm == 0)

    # Extinction correction
    if False:
        msgs.info("Applying extinction correction")
        extinct = load_extinction_data(longitude, latitude)
        ext_corr = extinction_correction(wave * units.AA, airmass, extinct)
        # Correct for extinction
        flux_star = flux_star * ext_corr
        ivar_star = ivar_star / ext_corr ** 2

#    embed()

    # Save the datacube
    outfile = "datacube_resid.fits"
    msgs.info("Saving datacube as: {0:s}".format(outfile))
    hdu = fits.PrimaryHDU(datacube_resid/norm, header=masterwcs.to_header())
#    hdu = fits.PrimaryHDU(norm, header=masterwcs.to_header())
    hdu.writeto(outfile, overwrite=args.overwrite)

    outfile = "datacube.fits"
    msgs.info("Saving datacube as: {0:s}".format(outfile))
    hdu = fits.PrimaryHDU(datacube/norm, header=masterwcs.to_header())
    hdu.writeto(outfile, overwrite=args.overwrite)
