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
    parser = argparse.ArgumentParser(description='Display MasterAlignment image and the trace data',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('master_file', type=str, default = None,
                        help='PypeIt Master Alignment file [e.g. MasterAlignment_A_1_01.fits]')
    parser.add_argument('--chname', default='Alignments', type=str,
                        help='Channel name for image in Ginga')

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(pargs):
    from pypeit import alignframe

    # Load
    alignments = alignframe.Alignments.from_file(pargs.master_file)
    alignments.show()
