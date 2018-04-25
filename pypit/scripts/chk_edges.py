#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-

"""
This script displays the Trace image and the traces
in an RC Ginga window (must be previously launched)
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse

def parser(options=None):

    parser = argparse.ArgumentParser(description='Display MasterTrace image in a previously launched RC Ginga viewer',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file', type = str, default = None, help='PYPIT Master Trace file [e.g. MasterTrace_A_01_aa.fits]')
    parser.add_argument("--chname", default='MTrace', type=str, help="Channel name for image in Ginga")
    #parser.add_argument('--det', default=1, type=int, help="Detector")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(pargs):

    import os

    import pdb as debugger
    import time

    from astropy.io import fits
    from astropy.table import Table

    from pypit import ginga
    from pypit.arspecobj import get_slitid

    import subprocess

    # List only?
    trc_hdu = fits.open(pargs.file)
    head0 = trc_hdu[0].header

    mstrace = trc_hdu[0].data

    try:
        ginga.connect_to_ginga(raise_err=True)
    except ValueError:
        subprocess.Popen(['ginga', '--modules=RC'])
        time.sleep(3)

    # Show Image
    viewer, ch = ginga.show_image(mstrace, chname=pargs.chname)

    lordloc = trc_hdu[5].data  # Should check name
    rordloc = trc_hdu[6].data  # Should check name

    # Get slit ids
    stup = (mstrace.shape, lordloc, rordloc)
    slit_ids = [get_slitid(stup, None, ii)[0] for ii in range(lordloc.shape[1])]
    ginga.show_slits(viewer, ch, lordloc, rordloc, slit_ids)


