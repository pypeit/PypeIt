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

debug = True
last_updated = "26 November 2015"
version = '0.3'
msgs = msgs.get_logger((None, debug, last_updated, version))

import arlris

try:
    from xastropy.xutils import xdebug as debugger
except:
    import pdb as debugger

def main() :

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('lowrdx_file', type = str, default = None,
                        help = 'LowRedux Pixel Flat FITS file')
    parser.add_argument('new_file', type = str, default = None, help = 'PYPIT FITS file')

    args = parser.parse_args()

    # Assume LRIS for now
    arlris.convert_lowredux_pixflat(args.lowrdx_file, args.new_file)


if __name__ == '__main__':
    main()
