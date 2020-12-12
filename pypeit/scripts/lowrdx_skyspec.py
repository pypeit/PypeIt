#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
This script generates a sky spectrum from a LowRedux IDL save file
"""

def parse_args(options=None, return_parser=False):
    import argparse

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('lowrdx_sky', type = str, default = None,
                        help = 'LowRedux Sky Spectrum (IDL save file)')
    parser.add_argument('new_file', type = str, default = None, help = 'PYPIT FITS sky spectrum')

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(args):
    from scipy.io.idl import readsav
    from linetools.spectra.xspectrum1d import XSpectrum1D

    # Read
    lrdx_sky = readsav(args.lowrdx_sky)
    # Generate
    xspec = XSpectrum1D.from_tuple((lrdx_sky['wave_calib'], lrdx_sky['sky_calib']))
    # Write
    xspec.write_to_fits(args.new_file)

