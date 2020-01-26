#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
Launch the identify GUI tool.
"""

def parser(options=None):

    import argparse

    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description='Launch PypeIt identify tool, display extracted'
                                                 'MasterArc, and load linelist.'
                                                 'Run above the Masters/ folder',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file', type = str, default=None, help='PypeIt MasterArc file')
    parser.add_argument("--lamps", default=False, help="Comma separated list of calibration lamps (no spaces)",
                        action="store_true")
    parser.add_argument("--wmin", default=3000.0, help="Minimum wavelength range",
                        action="store_true")
    parser.add_argument("--wmax", default=10000.0, help="Maximum wavelength range",
                        action="store_true")

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(args):

    import time
    import os
    import sys
    import numpy as np
    import astropy.io.fits as fits
    from pypeit.spectrographs.util import load_spectrograph
    from pypeit import traceimage, edgetrace, biasframe
    from pypeit.pypeit import PypeIt
    from pypeit.core import parse

    from IPython import embed

    if os.path.exists(args.file):
        msarc = fits.open(args.file)
    else:
        try:
            msarc = fits.open("Masters/{0:s}".format(args.file))
        except FileNotFoundError:
            print("Could not find MasterArc file.")
            sys.exit()
    specname = msarc[0].header['PYP_SPEC']
    spec = load_spectrograph(specname)
    par = spec.default_pypeit_par()
    # Get the lamp list
    if args.lamps is None:
        lamplist = par['calibrations']['wavelengths']['lamps']
        if lamplist is None:
            print("ERROR :: Cannot determine the lamps")
            sys.exit()
    else:
        lamplist = args.lamps.split(",")
