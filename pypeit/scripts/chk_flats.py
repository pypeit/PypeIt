#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
This script displays the flat images
in an RC Ginga window (must be previously launched)
"""
import argparse

def parser(options=None):

    parser = argparse.ArgumentParser(description='Display MasterFlat images in a previously launched RC Ginga viewer',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('master_file', type=str, help='PYPIT MasterFlat file [e.g. MasterFlat_A_1_01.fits]')

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(pargs):

    import time

    from pypeit import ginga
    from pypeit import flatfield

    import subprocess

    # Load up
    flatField = flatfield.FlatField.from_master_file(pargs.master_file)

    try:
        ginga.connect_to_ginga(raise_err=True)
    except ValueError:
        subprocess.Popen(['ginga', '--modules=RC'])
        time.sleep(3)

    # Show RawFlatImage
    viewer, ch = ginga.show_image(flatField.rawflatimg.image, chname='Raw Flat')
    # PixelFlat
    if flatField.mspixelflat is not None:
        viewer, ch = ginga.show_image(flatField.mspixelflat, chname='Pixel Flat')
    # Illumination flat
    if flatField.msillumflat is not None:
        viewer, ch = ginga.show_image(flatField.msillumflat, chname='Illumination Flat')
    # Illumination flat
    if flatField.flat_model is not None:
        viewer, ch = ginga.show_image(flatField.flat_model, chname='Flat Model')

    print("Check your Ginga viewer")


