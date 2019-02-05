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

    parser.add_argument('root', type=str, default = None, help='PYPIT Master Trace file root [e.g. MasterTrace_A_01_aa]')
    parser.add_argument("--chname", default='MTrace', type=str, help="Channel name for image in Ginga")
    parser.add_argument("--dumb_ids", default=False, action="store_true", help="Slit ID just by order?")
    parser.add_argument("--show", type=str, help="Use the show() method of TraceSlits to show something else")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(pargs):

    import pdb as debugger
    import time

    from pypeit import ginga
    from pypeit import traceslits
    from pypeit.core.trace_slits import get_slitid

    import subprocess

    # Load up
    Tslits = traceslits.TraceSlits.from_master_files(pargs.root)

    try:
        ginga.connect_to_ginga(raise_err=True)
    except ValueError:
        subprocess.Popen(['ginga', '--modules=RC'])
        time.sleep(3)

    if pargs.show is not None:
        Tslits.show(pargs.show)
        print("Check your Ginga viewer")
        return

    # Show Image
    viewer, ch = ginga.show_image(Tslits.mstrace, chname=pargs.chname)

    # Get slit ids
    stup = (Tslits.mstrace.shape, Tslits.lcen, Tslits.rcen)
    if pargs.dumb_ids:
        slit_ids = range(Tslits.lcen.shape[1])
    else:
        slit_ids = [get_slitid(stup[0], stup[1], stup[2], ii)[0] for ii in range(Tslits.lcen.shape[1])]
    ginga.show_slits(viewer, ch, Tslits.lcen, Tslits.rcen, slit_ids, pstep=50)
    print("Check your Ginga viewer")


