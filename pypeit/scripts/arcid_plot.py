#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-

"""
This script generates an ArcID plot from a Master WaveSoln file
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse

def parser(options=None) :

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('wave_soln', type = str, default = None,
                        help = 'MasterWaveSoln file [JSON]')
    parser.add_argument('title', type = str, default = None, help = 'Title for the plot')
    parser.add_argument('outfile', type = str, default = None, help = 'Output PDF file')

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args):
    """
    Parameters
    ----------
    args

    Returns
    -------

    """
    import numpy as np

    try:
        from xastropy.xutils import xdebug as debugger
    except:
        import pdb as debugger

    from linetools.utils import loadjson

    from pypit import arqa
    from pypit import msgs
    msgs.reset(verbosity=2)

    # Read JSON
    fdict = loadjson(args.wave_soln)
    for key in fdict.keys():
        if isinstance(fdict[key], list):
            fdict[key] = np.array(fdict[key])

    # Generate QA
    arqa.arc_fit_qa(None, fdict, outfil=args.outfile, ids_only=True,
                    title=args.title)
    print("Wrote {:s}".format(args.outfile))

