#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
This script displays the Trace image and the traces
in an RC Ginga window (must be previously launched)
"""
import argparse

def parser(options=None):
    parser = argparse.ArgumentParser(description='Display MasterEdges image and trace data',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('trace_file', type=str, default = None,
                        help='PypeIt Master Trace file [e.g. MasterEdges_A_01_aa.fits.gz]')
    parser.add_argument('--chname', default='MTrace', type=str,
                        help='Channel name for image in Ginga')
    parser.add_argument('--mpl', default=False, action='store_true',
                        help='Use a matplotlib window instead of ginga to show the trace')

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(pargs):
    from pypeit import edgetrace

    edges = edgetrace.EdgeTraceSet.from_file(pargs.trace_file)
    if pargs.mpl:
        edges.show(thin=10, include_img=True, idlabel=True)
    else:
        edges.show(thin=10, in_ginga=True)
    return 0
