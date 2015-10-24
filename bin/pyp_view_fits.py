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

sys.path.append(os.path.abspath('/Users/xavier/local/Python/PYPIT/src'))
import arlris 
import armsgs as msgs

def main() :

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file', type = str, default = None, help = 'FITS file')
    parser.add_argument('--raw_lris', action="store_true")

    args = parser.parse_args()

    # RAW_LRIS??
    if args.raw_lris:
        # 
        img, head = arlris.read_lris(args.file)
        # Generate hdu
        hdu = fits.PrimaryHDU(img)
        hdulist = fits.HDUList([hdu])
        # Write
        kludge_fil = 'tmp_ginga.fits'
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
