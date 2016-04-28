#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-

"""
This script pushes a FITS file to ginga
"""

import argparse
import numpy as np
import sys, os
import subprocess
from astropy.io import fits

# Setup for PYPIT imports
this_file = os.path.realpath(__file__)
this_path = this_file[:this_file.rfind('/')]
sys.path.append(os.path.abspath(this_path+'/../src'))
import armsgs as msgs

import ardebug
debug = ardebug.init()
last_updated = "26 November 2015"
version = '0.3'
verbose = False
msgs = msgs.get_logger((None, debug, last_updated, version, verbose))

import arlris

def main() :

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file', type = str, default = None, help = 'FITS file')
    parser.add_argument('--raw_lris', action="store_true")
    parser.add_argument('--exten', type=int, help="FITS extension")

    args = parser.parse_args()

    kludge_fil = 'tmp_ginga.fits'

    # Extension
    if args.exten is not None:
        hdu = fits.open(args.file)
        img = hdu[args.exten].data
        # Write
        msgs.warn('Writing kludge file to {:s}'.format(kludge_fil))
        hdunew = fits.PrimaryHDU(img)
        hdulist = fits.HDUList([hdunew])
        hdulist.writeto(kludge_fil,clobber=True)
        #
        args.file = kludge_fil

    # RAW_LRIS??
    if args.raw_lris:
        # 
        img, head, _ = arlris.read_lris(args.file)
        # Generate hdu
        hdu = fits.PrimaryHDU(img)
        hdulist = fits.HDUList([hdu])
        # Write
        msgs.warn('Writing kludge file to {:s}'.format(kludge_fil))
        hdulist.writeto(kludge_fil,clobber=True)
        args.file = kludge_fil

    # Spawn ginga
    subprocess.call(["ginga", args.file])

    if args.raw_lris:
        msgs.warn('Removing kludge file {:s}'.format(kludge_fil))
        subprocess.call(["rm", args.file])

if __name__ == '__main__':
    main()
