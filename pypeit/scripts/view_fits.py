#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-

"""
This script enables the viewing of a FITS file
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse

def parser(options=None):

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file', type = str, default = None, help = 'FITS file')
    parser.add_argument("--list", default=False, help="List the extensions only?", action="store_true")
    parser.add_argument('--raw_lris', action="store_true")
    parser.add_argument('--raw_deimos', action="store_true")
    parser.add_argument('--raw_gmos', action="store_true")
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
    if args.raw_lris:
        # 
        img, head, _ = keck_lris.read_lris(args.file)
    # RAW_DEIMOS??
    elif args.raw_deimos:
        #
        img, head, _ = keck_deimos.read_deimos(args.file)
    # RAW_GEMINI??
    elif args.raw_gmos:
        # TODO this routine should show the whole mosaic if no detector number is passed in!
        # Need to figure out the number of amps
        img, head, _ = gemini_gmos.read_gmos(args.file, det=args.det)
    else:
        hdu = fits.open(args.file)
        img = hdu[args.exten].data
        # Write

    ginga.show_image(img)