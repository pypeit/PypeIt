#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-

"""
This script pushes a FITS file to ginga
"""

import argparse
from scipy.io.idl import readsav
from linetools.spectra.xspectrum1d import XSpectrum1D

# Setup for PYPIT imports

try:
    from xastropy.xutils import xdebug as debugger
except:
    import pdb as debugger

def main() :

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('lowrdx_sky', type = str, default = None,
                        help = 'LowRedux Sky Spectrum (IDL save file)')
    parser.add_argument('new_file', type = str, default = None, help = 'PYPIT FITS sky spectrum')

    args = parser.parse_args()

    # Read
    lrdx_sky = readsav(args.lowrdx_sky)
    # Generate
    xspec = XSpectrum1D.from_tuple((lrdx_sky['wave_calib'], lrdx_sky['sky_calib']))
    # Write
    xspec.write_to_fits(args.new_file)


if __name__ == '__main__':
    main()
