#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
This script enables the viewing of a raw FITS file
"""

def parse_args(options=None, return_parser=False):
    import argparse
    from pypeit.spectrographs import available_spectrographs

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file', type = str, default = None, help = 'FITS file')
    parser.add_argument('spectrograph', type=str,
                        help='A valid spectrograph identifier: {0}'.format(
                             ', '.join(available_spectrographs)))
    parser.add_argument("--list", default=False, help="List the extensions only?", action="store_true")
    parser.add_argument('--exten', type=int, default = 0, help="FITS extension")
    parser.add_argument('--det', type=int, default=1, help="Detector number (ignored for keck_lris, keck_deimos")
    parser.add_argument('--chname', type=str, default='Image', help="Name of Ginga tab")

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(args):

    import subprocess

    from astropy.io import fits

    from pypeit import msgs
    from pypeit.spectrographs import keck_lris
    from pypeit.spectrographs import keck_deimos
    from pypeit.spectrographs import gemini_gmos
    from pypeit.display import display
    from pypeit.spectrographs import mmt_binospec
    from pypeit.spectrographs import mmt_mmirs
    from pypeit.spectrographs import mmt_bluechannel
    from pypeit.spectrographs import util
    from pypeit import msgs
    from pypeit import io

    # List only?
    if args.list:
        hdu = io.fits_open(args.file)
        print(hdu.info())
        return


    # Setup for PYPIT imports
    msgs.reset(verbosity=2)

    # RAW_LRIS??
    if 'keck_lris' in args.spectrograph:
        #
        if args.spectrograph == 'keck_lris_red_orig':
            gen_lris = keck_lris.KeckLRISROrigSpectrograph()
            img = gen_lris.get_rawimage(args.file, 1)[1]
        else:
            gen_lris = keck_lris.KeckLRISRSpectrograph()  # Using LRISr, but this will work for LRISb too
            img = gen_lris.get_rawimage(args.file,  None)[1]
    # RAW_DEIMOS??
    elif args.spectrograph == 'keck_deimos':
        #
        gen_deimos = keck_deimos.KeckDEIMOSSpectrograph()
        img = gen_deimos.get_rawimage(args.file, None)[1]
    # RAW_GEMINI??
    elif 'gemini_gmos' in args.spectrograph:
        # TODO this routine should show the whole mosaic if no detector number is passed in!
        # Need to figure out the number of amps
        spectrograph = util.load_spectrograph(args.spectrograph)
        img = spectrograph.get_rawimage(args.file, args.det)[1]
    # RAW_BinoSpec
    elif args.spectrograph == 'mmt_binospec':
        #
        gen_bino = mmt_binospec.MMTBINOSPECSpectrograph()
        img = gen_bino.get_rawimage(args.file, args.det)[1]
    # RAW_MMIRS
    elif args.spectrograph == 'mmt_mmirs':
        gen_mmirs = mmt_mmirs.MMTMMIRSSpectrograph()
        img = gen_mmirs.get_rawimage(args.file, args.det)[1]
    # RAW MMT blue channel
    elif args.spectrograph == 'mmt_bluechannel':
        gen_bluechan = mmt_bluechannel.MMTBlueChannelSpectrograph()
        img = gen_bluechan.get_rawimage(args.file, args.det)[1]
    else:
        hdu = io.fits_open(args.file)
        img = hdu[args.exten].data
        # Write

    display.show_image(img,chname=args.chname)
