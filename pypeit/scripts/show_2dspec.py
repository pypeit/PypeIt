#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
This script enables the viewing of a processed FITS file
with extras.  Run above the Science/ folder.
"""
import os

import numpy as np

from IPython import embed

from astropy.io import fits
from astropy.stats import sigma_clipped_stats

from pypeit import msgs
from pypeit import slittrace
from pypeit import specobjs

from pypeit.display import display
from pypeit.core.parse import get_dnum
from pypeit.images.imagebitmask import ImageBitMask
from pypeit import masterframe
from pypeit import spec2dobj


def parse_args(options=None, return_parser=False):
    import argparse
    parser = argparse.ArgumentParser(description='Display sky subtracted, spec2d image in a '
                                                 'Ginga viewer.  Run above the Science/ folder',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file', type = str, default = None, help = 'PYPIT spec2d file')
    parser.add_argument('--list', default=False, help='List the extensions only?',
                        action='store_true')
    parser.add_argument('--det', default=1, type=int, help='Detector number')
    parser.add_argument('--showmask', default=False, help='Overplot masked pixels',
                        action='store_true')
    parser.add_argument('--removetrace', default=False, help="Do not overplot traces in the skysub, "
                                                             "sky_resid and resid channels",
                        action = "store_true")
    parser.add_argument('--embed', default=False, help='Upon completion embed in ipython shell',
                        action='store_true')
    parser.add_argument('--ignore_extract_mask', default=False, help='Ignore the extraction mask',
                        action='store_true')

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)


def show_trace(specobjs, det, viewer, ch):

    if specobjs is None:
        return
    in_det = np.where(specobjs.DET == det)[0]
    for kk in in_det:
        trace = specobjs[kk]['TRACE_SPAT']
        obj_id = specobjs[kk].NAME
        display.show_trace(viewer, ch, trace, obj_id, color='orange') #hdu.name)


def main(args):

    # List only?
    if args.list:
        hdu = fits.open(args.file)
        hdu.info()
        return

    # Load it up -- NOTE WE ALLOW *OLD* VERSIONS TO GO FORTH
    spec2DObj = spec2dobj.Spec2DObj.from_file(args.file, args.det, chk_version=False)

    # Setup for PypeIt imports
    msgs.reset(verbosity=2)

    # Init
    # TODO: get_dnum needs to be deprecated...
    sdet = get_dnum(args.det, prefix=False)


    # Grab the slit edges
    slits = spec2DObj.slits
    if spec2DObj.sci_spat_flexure is not None:
        msgs.info("Offseting slits by {}".format(spec2DObj.sci_spat_flexure))
    all_left, all_right, mask = slits.select_edges(flexure=spec2DObj.sci_spat_flexure)
    # TODO -- This may be too restrictive, i.e. ignore BADFLTCALIB??
    gpm = mask == 0
    left = all_left[:, gpm]
    right = all_right[:, gpm]
    slid_IDs = spec2DObj.slits.slitord_id[gpm]

    bitMask = ImageBitMask()

    # Object traces from spec1d file
    spec1d_file = args.file.replace('spec2d', 'spec1d')
    if args.file[-2:] == 'gz':
        spec1d_file = spec1d_file[:-3]
    if os.path.isfile(spec1d_file):
        sobjs = specobjs.SpecObjs.from_fitsfile(spec1d_file)
    else:
        sobjs = None
        msgs.warn('Could not find spec1d file: {:s}'.format(spec1d_file) + msgs.newline() +
                  '                          No objects were extracted.')

    display.connect_to_ginga(raise_err=True, allow_new=True)

    # Now show each image to a separate channel

    # Show the bitmask?
    mask_in = None
    if args.showmask:
        viewer, ch = display.show_image(spec2DObj.bpmmask, chname="BPM", waveimg=spec2DObj.waveimg, clear=True)
        #bpm, crmask, satmask, minmask, offslitmask, nanmask, ivar0mask, ivarnanmask, extractmask \

    # SCIIMG
    image = spec2DObj.sciimg  # Processed science image
    mean, med, sigma = sigma_clipped_stats(image[spec2DObj.bpmmask == 0], sigma_lower=5.0,
                                           sigma_upper=5.0)
    cut_min = mean - 1.0 * sigma
    cut_max = mean + 4.0 * sigma
    chname_skysub='sciimg-det{:s}'.format(sdet)
    # Clear all channels at the beginning
    viewer, ch = display.show_image(image, chname=chname_skysub, waveimg=spec2DObj.waveimg, clear=True)

    if sobjs is not None:
        show_trace(sobjs, args.det, viewer, ch)
    display.show_slits(viewer, ch, left, right, slit_ids=slid_IDs)

    # SKYSUB
    if args.ignore_extract_mask:
        # TODO -- Is there a cleaner way to do this?
        gpm = (spec2DObj.bpmmask == 0) | (spec2DObj.bpmmask == 2**bitMask.bits['EXTRACT'])
    else:
        gpm = spec2DObj.bpmmask == 0

    image = (spec2DObj.sciimg - spec2DObj.skymodel) * gpm #(spec2DObj.mask == 0)  # sky subtracted image
    mean, med, sigma = sigma_clipped_stats(image[spec2DObj.bpmmask == 0], sigma_lower=5.0,
                                           sigma_upper=5.0)
    cut_min = mean - 1.0 * sigma
    cut_max = mean + 4.0 * sigma
    chname_skysub='skysub-det{:s}'.format(sdet)
    # Clear all channels at the beginning
    # TODO: JFH For some reason Ginga crashes when I try to put cuts in here.
    viewer, ch = display.show_image(image, chname=chname_skysub, waveimg=spec2DObj.waveimg,
                                  bitmask=bitMask, mask=mask_in) #, cuts=(cut_min, cut_max),wcs_match=True)
    if not args.removetrace and sobjs is not None:
            show_trace(sobjs, args.det, viewer, ch)
    display.show_slits(viewer, ch, left, right, slit_ids=slid_IDs)


    # SKRESIDS
    chname_skyresids = 'sky_resid-det{:s}'.format(sdet)
    image = (spec2DObj.sciimg - spec2DObj.skymodel) * np.sqrt(spec2DObj.ivarmodel) * gpm #(spec2DObj.bpmmask == 0)  # sky residual map
    viewer, ch = display.show_image(image, chname_skyresids, waveimg=spec2DObj.waveimg,
                                  cuts=(-5.0, 5.0), bitmask=bitMask, mask=mask_in)
    if not args.removetrace and sobjs is not None:
            show_trace(sobjs, args.det, viewer, ch)
    display.show_slits(viewer, ch, left, right, slit_ids=slid_IDs)

    # RESIDS
    chname_resids = 'resid-det{:s}'.format(sdet)
    # full model residual map
    image = (spec2DObj.sciimg - spec2DObj.skymodel - spec2DObj.objmodel) * np.sqrt(spec2DObj.ivarmodel) * (spec2DObj.bpmmask == 0)
    viewer, ch = display.show_image(image, chname=chname_resids, waveimg=spec2DObj.waveimg,
                                  cuts = (-5.0, 5.0), bitmask=bitMask, mask=mask_in)
    if not args.removetrace and sobjs is not None:
            show_trace(sobjs, args.det, viewer, ch)
    display.show_slits(viewer, ch, left, right, slit_ids=slid_IDs)


    # After displaying all the images sync up the images with WCS_MATCH
    shell = viewer.shell()
    shell.start_global_plugin('WCSMatch')
    shell.call_global_plugin_method('WCSMatch', 'set_reference_channel', [chname_resids], {})

    if args.embed:
        embed()

        # Playing with some mask stuff
        #out = shell.start_operation('TVMask')
        #maskfile = '/Users/joe/python/PypeIt-development-suite/REDUX_OUT/Shane_Kast_blue/600_4310_d55/shane_kast_blue_setup_A/crmask.fits'
        #out = shell.call_local_plugin_method(chname_resids, 'TVMask', 'load_file', [maskfile], {})


