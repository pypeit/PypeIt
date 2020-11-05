#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
This script displays the Trace image and the traces
in an RC Ginga window (must be previously launched)
"""

def parse_args(options=None, return_parser=False):
    import argparse

    parser = argparse.ArgumentParser(description='Display MasterEdges image and trace data',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('trace_file', type=str, default = None,
                        help='PypeIt Master Trace file [e.g. MasterEdges_A_01_aa.fits.gz]')
    parser.add_argument('--chname', default='MTrace', type=str,
                        help='Channel name for image in Ginga')
    parser.add_argument('--mpl', default=False, action='store_true',
                        help='Use a matplotlib window instead of ginga to show the trace')
    parser.add_argument('--try_old', default=False, action='store_true',
                        help='Attempt to load old datamodel versions.  A crash may ensue..')

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)

# TODO: JFH I don't see why we are showing the edges and/or what the
# purpose of all this if synced not synced is fore. edgetrace seems to
# crash if the syncing fails. So if we have successfuly run EdgeTrace,
# we create a slittrace object and the slittrace object is the thing we
# should be showing not the edgetrace object. This has the advantage
# that then orders are correctly labeled for Echelle which is not the
# case with the current show method.

def main(pargs):
    from pypeit import edgetrace

    # Load
    edges = edgetrace.EdgeTraceSet.from_file(pargs.trace_file, chk_version=(not pargs.try_old))

    if pargs.mpl:
        edges.show(thin=10, include_img=True, idlabel=True)
    else:
        edges.show(thin=10, in_ginga=True)
    return 0
