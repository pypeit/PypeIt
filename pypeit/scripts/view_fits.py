#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
This script enables the viewing of a FITS file
"""
import argparse

def parser(options=None):

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file', type = str, default = None, help = 'FITS file')
    parser.add_argument('spectrograph', type=str, help='Name of the spectrograph')
    parser.add_argument("--list", default=False, help="List the extensions only?", action="store_true")
    parser.add_argument('--exten', type=int, default = 0, help="FITS extension")
    parser.add_argument('--det', type=int, default=1, help="Detector number")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args

def main(args):

    import subprocess

    from astropy.io import fits

    from pypeit import msgs
    from pypeit.spectrographs import keck_lris
    from pypeit.spectrographs import keck_deimos
    from pypeit.spectrographs import gemini_gmos
    from pypeit import msgs
    from pypeit import ginga

    # List only?
    if args.list:
        hdu = fits.open(args.file)
        print(hdu.info())
        return


    # Setup for PYPIT imports
    msgs.reset(verbosity=2)

    # RAW_LRIS??
    if 'keck_lris' in args.spectrograph:
        #
        gen_lris = keck_lris.KeckLRISRSpectrograph()  # Using LRISr, but this will work for LRISb too
        img, _, _, _, _ = gen_lris.get_rawimage(args.file,  None)
    # RAW_DEIMOS??
    elif args.spectrograph == 'keck_deimos':
        #
        gen_deimos = keck_deimos.KeckDEIMOSSpectrograph()
        img, _, _, _, _ = gen_deimos.get_rawimage(args.file, None)
    # RAW_GEMINI??
    elif 'gemini_gmos' in args.spectrograph:
        # TODO this routine should show the whole mosaic if no detector number is passed in!
        # Need to figure out the number of amps
        gen_gmos = gemini_gmos.GeminiGMOSSpectrograph()
        img, _, _, _, _ = gen_gmos.get_rawimage(args.file, args.det)
    else:
        hdu = fits.open(args.file)
        img = hdu[args.exten].data
        # Write

    ginga.show_image(img)
