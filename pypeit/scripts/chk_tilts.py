#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
This script displays the Arc image from and the traces
in an RC Ginga window (must be previously launched)

.. todo::

    This script is out of date and will not function with the current
    pypeit master branch.

"""
import argparse

def parser(options=None):

    parser = argparse.ArgumentParser(description='Display MasterArc image in a previously launched RC Ginga viewer with tilts',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('option', type = str, default = None, help='Item to show [fweight, model, tilts, final_tilts]')
    parser.add_argument('setup', type = str, default = None, help='setup  -- Run from MF folder (e.g. A_01_aa)')
    parser.add_argument('--slit', type=int, default=None, help='Slit/Order [0,1,2..]')
    #parser.add_argument("--dumb_ids", default=False, action="store_true", help="Slit ID just by order?")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(pargs):

    import pdb as debugger
    import time

    from pypeit import ginga
    from pypeit import wavetilts

    import subprocess

    # Load up
    wTilts = wavetilts.WaveTilts.from_master_files(pargs.setup)

    # Launch ginga if need be
    try:
        ginga.connect_to_ginga(raise_err=True)
    except ValueError:
        subprocess.Popen(['ginga', '--modules=RC'])
        time.sleep(3)

    # Show
    if pargs.slit is not None:
        cname = 'Slit{:03d}'.format(pargs.slit)
    else:
        cname = None
    wTilts.show(pargs.option, slit=pargs.slit, cname=cname)

