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

import os
import argparse
from IPython import embed

from astropy.io import fits
from linetools import utils as ltu
import numpy as np

from pypeit import msgs
from pypeit.spectrographs.util import load_spectrograph
from pypeit.core.parse import get_dnum
from pypeit.images.imagebitmask import ImageBitMask
from pypeit import alignframe
from pypeit import masterframe
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
    ivar = spec2DObj.ivarmodel
    waveimg = spec2DObj.waveimg

    # Load the wavelength solution and generate a waveimg
    # EXT0001 = 'DET01-SCIIMG'
    # EXT0002 = 'DET01-IVARRAW'
    # EXT0003 = 'DET01-SKYMODEL'
    # EXT0004 = 'DET01-OBJMODEL'
    # EXT0005 = 'DET01-IVARMODEL'
    # EXT0006 = 'DET01-TILTS'
    # EXT0007 = 'DET01-WAVEIMG'
    # EXT0008 = 'DET01-BPMMASK'
    # EXT0009 = 'DET01-SLITS'

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

    # Generate a master WCS to register all frames
    masterwcs = spec.get_wcs(spec2DObj.head0, slits, wave0, dwv)

    pxscl, slscl = spec.get_scales(spec2DObj.head0)

    # Grab the WCS of this frame
    wcs = spec.get_wcs(spec2DObj.head0, slits, wave0, dwv)

    # Generate an RA/DEC image
    raimg, decimg, minmax = spec.get_radec_image(alignments, slits, wcs, flexure=spec2DObj.sci_spat_flexure)

    # Generate the output binning
    slitlength = int(np.round(np.median(slits.get_slitlengths(initial=True, median=True))))
    xbins = np.arange(1 + 24) - 12.0 - 0.5
    ybins = np.linspace(np.min(minmax[:, 0]), np.max(minmax[:, 1]), 1+slitlength) - 0.5
    zbins = np.arange(1+int(round((np.max(waveimg)-wave0)/dwv))) - 0.5

    # Make the cube
    msgs.info("Generating datacube")
    slitid_img_init = slits.slit_img(pad=0, initial=True, flexure=spec2DObj.sci_spat_flexure)
    onslit = slitid_img_init > 0
    pix_coord = wcs.wcs_world2pix(raimg[onslit], decimg[onslit], waveimg[onslit]*1.0E-10, 0)
    datacube_resid, edges = np.histogramdd(pix_coord, bins=(xbins, ybins, zbins), weights=sciimg[onslit]*np.sqrt(ivar[onslit]))
    datacube, edges = np.histogramdd(pix_coord, bins=(xbins, ybins, zbins), weights=sciimg[onslit])

    # Save the datacube
    outfile = "datacube_resid.fits"
    msgs.info("Saving datacube as: {0:s}".format(outfile))
    hdu = fits.PrimaryHDU(datacube_resid, header=masterwcs.to_header())
    hdu.writeto(outfile, overwrite=args.overwrite)

    outfile = "datacube.fits"
    msgs.info("Saving datacube as: {0:s}".format(outfile))
    hdu = fits.PrimaryHDU(datacube, header=masterwcs.to_header())
    hdu.writeto(outfile, overwrite=args.overwrite)
