#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-

"""
This script converts a LowRedux pixel flat into a PYPIT ready one
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse

try:
    from xastropy.xutils import xdebug as debugger
except:
    import pdb as debugger

from pypit import arlris

def parser(options=None):

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('lowrdx_file', type = str, default = None,
                        help = 'LowRedux Pixel Flat FITS file')
    parser.add_argument('new_file', type = str, default = None, help = 'PYPIT FITS file')

    args = parser.parse_args()

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args):
    # Assume LRIS for now
    arlris.convert_lowredux_pixflat(args.lowrdx_file, args.new_file)

