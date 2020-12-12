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

from astropy import units
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import numpy as np
import copy, os

from pypeit import msgs, par, io, spec2dobj
from pypeit.spectrographs.util import load_spectrograph
from pypeit.core import datacube as dc_utils
from pypeit.core.flux_calib import load_extinction_data, extinction_correction
from pypeit.core.flexure import calculate_image_offset
from pypeit.core import parse

from IPython import embed


def parse_args(options=None, return_parser=False):

    parser = argparse.ArgumentParser(description='Read in an array of spec2D files and convert them into a datacube',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file', type = str, default=None, help='filename.coadd3d file')
    parser.add_argument('--det', default=1, type=int, help="Detector")
    parser.add_argument('-o', '--overwrite', default=False, action='store_true',
                        help='Overwrite any existing files/directories')

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)


def coadd_cube(files, parset, overwrite=False):
    """ Main routine to coadd spec2D files into a 3D datacube

    Args:
        files (list):
            List of all spec2D files
        parset (:class:`pypeit.par.core.PypeItPar`):
            An instance of the parameter set.
        overwrite (bool):
            Overwrite the output file, if it exists?
    """
    # Get the detector number
    det = 1 if parset is None else parset['rdx']['detnum']

    # Load the spectrograph
    spec2DObj = spec2dobj.Spec2DObj.from_file(files[0], det)
    specname = spec2DObj.head0['PYP_SPEC']
    spec = load_spectrograph(specname)

    # Grab the parset, if not provided
    if parset is None: parset = spec.default_pypeit_par()
    cubepar = parset['reduce']['cube']

    # Check the output file
    outfile = cubepar['output_filename'] if ".fits" in cubepar['output_filename'] else cubepar['output_filename']+".fits"
    out_whitelight = outfile.replace(".fits", "_whitelight.fits")
    if os.path.exists(outfile) and not overwrite:
        msgs.error("Output filename already exists:"+msgs.newline()+outfile)
    elif os.path.exists(out_whitelight) and cubepar['save_whitelight'] and not overwrite:
        msgs.error("Output filename already exists:"+msgs.newline()+out_whitelight)
    # Check the reference cube and image exist, if requested
    ref_scale = None  # This will be used to correct relative scaling among the various input frames
    if cubepar['standard_cube'] is not None:
        if not os.path.exists(cubepar['standard_cube']):
            msgs.error("Standard cube does not exist:" + msgs.newline() + cubepar['reference_cube'])
        cube = fits.open(cubepar['standard_cube'])
        ref_scale = cube['REFSCALE'].data
    if cubepar['reference_image'] is not None:
        if not os.path.exists(cubepar['reference_image']):
            msgs.error("Reference cube does not exist:" + msgs.newline() + cubepar['reference_image'])
    if cubepar['flux_calibrate']:
        msgs.error("Flux calibration is not currently implemented" + msgs.newline() +
                   "Please set 'flux_calibrate = False'")

    # prep
    numfiles = len(files)
    combine = cubepar['combine']

    all_ra, all_dec, all_wave = np.array([]), np.array([]), np.array([])
    all_sci, all_ivar, all_idx, all_wghts = np.array([]), np.array([]), np.array([]), np.array([])
    all_wcs = []
    dspat = None if cubepar['spatial_delta'] is None else  cubepar['spatial_delta']/3600.0  # binning size on the sky (/3600 to convert to degrees)
    dwv = cubepar['wave_delta']       # binning size in wavelength direction (in Angstroms)
    wave_ref = None
    whitelight_img = None  # This is the whitelight image based on all input spec2d frames
    weights = np.ones(numfiles)  # Weights to use when combining cubes
    for ff, fil in enumerate(files):
        # Load it up
        spec2DObj = spec2dobj.Spec2DObj.from_file(fil, det)
        detector = spec2DObj.detector

        # Setup for PypeIt imports
        msgs.reset(verbosity=2)

        if ref_scale is None:
            ref_scale = spec2DObj.scaleimg.copy()
        # Extract the information
        sciimg = (spec2DObj.sciimg-spec2DObj.skymodel) * (ref_scale/spec2DObj.scaleimg)  # Subtract sky and apply relative sky
        ivar = spec2DObj.ivarraw / (ref_scale/spec2DObj.scaleimg)**2
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
        ra_corr, dec_corr = dc_utils.dar_correction(waveimg[onslit_gpm], coord, obstime, location,
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

        # Check if a reference whitelight image should be used to register the offsets
        if cubepar["reference_image"] is None:
            # Generate white light images
            whitelight_imgs, _, _ = dc_utils.make_whitelight(all_ra, all_dec, all_wave, all_sci, all_wghts, all_idx,
                                                             dspat)
            # ref_idx will be the index of the cube with the highest S/N
            ref_idx = np.argmax(weights)
            reference_image = whitelight_imgs[:, :, ref_idx].copy()
            msgs.info("Calculating spatial translation of each cube relative to cube #{0:d})".format(ref_idx+1))
        else:
            ref_idx = -1  # Don't use an index
            # Load reference information
            reference_image, whitelight_imgs, wlwcs = \
                dc_utils.make_whitelight_fromref(all_ra, all_dec, all_wave, all_sci, all_wghts, all_idx, dspat,
                                                 cubepar['reference_image'])
            msgs.info("Calculating the spatial translation of each cube relative to user-defined 'reference_image'")
        # Calculate the image offsets - check the reference is a zero shift
        ra_shift_ref, dec_shift_ref = calculate_image_offset(reference_image.copy(), reference_image.copy())
        for ff in range(numfiles):
            # Don't correlate the reference image with itself
            if ff == ref_idx:
                continue
            # Calculate the shift
            ra_shift, dec_shift = calculate_image_offset(whitelight_imgs[:, :, ff], reference_image.copy())
            # Convert to reference
            ra_shift -= ra_shift_ref
            dec_shift -= dec_shift_ref
            # Convert pixel shift to degress shift
            ra_shift *= dspat/cosdec
            dec_shift *= dspat
            msgs.info("Spatial shift of cube #{0:d}: RA, DEC (arcsec) = {1:+0.3f}, {2:+0.3f}".format(ff+1, ra_shift*3600.0, dec_shift*3600.0))
            # Apply the shift
            all_ra[all_idx == ff] += ra_shift
            all_dec[all_idx == ff] += dec_shift

        # Generate a white light image of *all* data
        msgs.info("Generating global white light image")
        if cubepar["reference_image"] is None:
            whitelight_img, _, wlwcs = dc_utils.make_whitelight(all_ra, all_dec, all_wave, all_sci, all_wghts,
                                                                np.zeros(all_ra.size), dspat)
        else:
            _, whitelight_img, wlwcs = \
                dc_utils.make_whitelight_fromref(all_ra, all_dec, all_wave, all_sci, all_wghts, np.zeros(all_ra.size),
                                                 dspat, cubepar['reference_image'])

        # Calculate the relative spectral weights of all pixels
        all_wghts = dc_utils.compute_weights(all_ra, all_dec, all_wave, all_sci, all_ivar, all_idx,
                                             whitelight_img[:, :, 0], dspat, dwv,
                                             relative_weights=cubepar['relative_weights'])
    # Check if a whitelight image should be saved
    if cubepar['save_whitelight']:
        # Check if the white light image still needs to be generated - if so, generate it now
        if whitelight_img is None:
            msgs.info("Generating global white light image")
            if cubepar["reference_image"] is None:
                whitelight_img, _, wlwcs = dc_utils.make_whitelight(all_ra, all_dec, all_wave, all_sci, all_wghts,
                                                                    np.zeros(all_ra.size), dspat)
            else:
                _, whitelight_img, wlwcs = \
                    dc_utils.make_whitelight_fromref(all_ra, all_dec, all_wave, all_sci, all_wghts,
                                                     np.zeros(all_ra.size),
                                                     dspat, cubepar['reference_image'])
        # Prepare and save the fits file
        msgs.info("Saving white light image as: {0:s}".format(out_whitelight))
        img_hdu = fits.PrimaryHDU(whitelight_img.T, header=wlwcs.to_header())
        img_hdu.writeto(out_whitelight, overwrite=overwrite)

    # Setup the cube ranges
    ra_min = cubepar['ra_min'] if cubepar['ra_min'] is not None else np.min(all_ra)
    ra_max = cubepar['ra_max'] if cubepar['ra_max'] is not None else np.max(all_ra)
    dec_min = cubepar['dec_min'] if cubepar['dec_min'] is not None else np.min(all_dec)
    dec_max = cubepar['dec_max'] if cubepar['dec_max'] is not None else np.max(all_dec)
    wav_min = cubepar['wave_min'] if cubepar['wave_min'] is not None else np.min(all_wave)
    wav_max = cubepar['wave_max'] if cubepar['wave_max'] is not None else np.max(all_wave)
    if cubepar['wave_delta'] is not None: dwv = cubepar['wave_delta']
    # Generate a master WCS to register all frames
    coord_min = [ra_min, dec_min, wav_min]
    coord_dlt = [dspat, dspat, dwv]
    masterwcs = dc_utils.generate_masterWCS(coord_min, coord_dlt, name=specname)
    msgs.info(msgs.newline()+"-"*40 +
              msgs.newline() + "Parameters of the WCS:" +
              msgs.newline() + "RA   min, max = {0:f}, {1:f}".format(ra_min, ra_max) +
              msgs.newline() + "DEC  min, max = {0:f}, {1:f}".format(dec_min, dec_max) +
              msgs.newline() + "WAVE min, max = {0:f}, {1:f}".format(wav_min, wav_max) +
              msgs.newline() + "Spaxel size = {0:f}''".format(3600.0*dspat) +
              msgs.newline() + "Wavelength step = {0:f} A".format(dwv) +
              msgs.newline() + "-" * 40)

    # Generate the output binning
    if combine:
        numra = int((ra_max-ra_min) * cosdec / dspat)
        numdec = int((dec_max-dec_min)/dspat)
        numwav = int((wav_max-wav_min)/dwv)
        xbins = np.arange(1+numra)-0.5
        ybins = np.arange(1+numdec)-0.5
        spec_bins = np.arange(1+numwav)-0.5
    else:
        slitlength = int(np.round(np.median(slits.get_slitlengths(initial=True, median=True))))
        numwav = int((np.max(waveimg) - wave0) / dwv)
        xbins, ybins, spec_bins = spec.get_datacube_bins(slitlength, minmax, numwav)

    # Make the cube
    msgs.info("Generating pixel coordinates")
    if combine:
        pix_coord = masterwcs.wcs_world2pix(all_ra, all_dec, all_wave * 1.0E-10, 0)
        hdr = masterwcs.to_header()
    else:
        pix_coord = wcs.wcs_world2pix(np.vstack((all_ra, all_dec, all_wave*1.0E-10)).T, 0)
        hdr = wcs.to_header()

    # Find the NGP coordinates for all input pixels
    msgs.info("Generating data cube")
    bins = (xbins, ybins, spec_bins)
    datacube, edges = np.histogramdd(pix_coord, bins=bins, weights=all_sci*all_wghts)
    norm, edges = np.histogramdd(pix_coord, bins=bins, weights=all_wghts)
    norm_cube = (norm > 0) / (norm + (norm == 0))
    datacube *= norm_cube
    # Create the variance cube, including weights
    msgs.info("Generating variance cube")
    all_var = (all_ivar > 0) / (all_ivar + (all_ivar == 0))
    var_cube, edges = np.histogramdd(pix_coord, bins=bins, weights=all_var * all_wghts**2)
    var_cube *= norm_cube**2

    # Save the datacube
    debug = False
    if debug:
        datacube_resid, edges = np.histogramdd(pix_coord, bins=(xbins, ybins, spec_bins), weights=all_sci*np.sqrt(all_ivar))
        norm, edges = np.histogramdd(pix_coord, bins=(xbins, ybins, spec_bins))
        norm_cube = (norm > 0) / (norm + (norm == 0))
        outfile = "datacube_resid.fits"
        msgs.info("Saving datacube as: {0:s}".format(outfile))
        hdu = fits.PrimaryHDU((datacube_resid*norm_cube).T, header=masterwcs.to_header())
        hdu.writeto(outfile, overwrite=overwrite)

    msgs.info("Saving datacube as: {0:s}".format(outfile))
    final_cube = dc_utils.DataCube(datacube.T, var_cube.T, specname,
                                   refscale=ref_scale, fluxed=cubepar['flux_calibrate'])
    final_cube.to_file(outfile, hdr=hdr, overwrite=overwrite)


def main(args):
    if args.file is None:
        msgs.error('You must input a coadd3d file')
    else:
        spectrograph_name, config_lines, spec2d_files = io.read_spec2d_file(args.file, filetype="coadd3d")
        spectrograph = load_spectrograph(spectrograph_name)

        # Parameters
        spectrograph_def_par = spectrograph.default_pypeit_par()
        parset = par.PypeItPar.from_cfg_lines(cfg_lines=spectrograph_def_par.to_config(),
                                              merge_with=config_lines)
        # If detector was passed as an argument override whatever was in the coadd3d file
        if args.det is not None:
            msgs.info("Restricting to detector={}".format(args.det))
            parset['rdx']['detnum'] = int(args.det)

    # Coadd the files
    coadd_cube(spec2d_files, parset, overwrite=args.overwrite)
