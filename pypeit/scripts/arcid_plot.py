#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-

"""
This script generates an ArcID plot from a Master WaveSoln file

.. todo::

    This script is out of date and will not function with the current
    pypeit master branch.

"""
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

    # TODO: This must be an out-dated script that is never used.
    # Deprecate it?
    from pypeit import arqa
    from pypeit import msgs
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

