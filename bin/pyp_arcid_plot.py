#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-

"""
This script generates an ArcID plot from a Master WaveSoln file
"""

import argparse
import numpy as np
import sys, os
import json

# Setup for PYPIT imports
this_file = os.path.realpath(__file__)
this_path = this_file[:this_file.rfind('/')]
sys.path.append(os.path.abspath(this_path+'/../src'))
import armsgs

import ardebug
debug = ardebug.init()
last_updated = "2 May 2016"
version = '0.1'
msgs = armsgs.get_logger((None, debug, last_updated, version, 1))

import arqa

from linetools.utils import loadjson

try:
    from xastropy.xutils import xdebug as debugger
except:
    import pdb as debugger

def main() :

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('wave_soln', type = str, default = None,
                        help = 'MasterWaveSoln file [JSON]')
    parser.add_argument('title', type = str, default = None, help = 'Title for the plot')
    parser.add_argument('outfile', type = str, default = None, help = 'Output PDF file')

    pargs = parser.parse_args()

    # Read JSON
    fdict = loadjson(pargs.wave_soln)
    for key in fdict.keys():
        if isinstance(fdict[key], list):
            fdict[key] = np.array(fdict[key])

    # Generate QA
    arqa.arc_fit_qa(None, fdict, outfil=pargs.outfile, ids_only=True,
                    title=pargs.title)
    print("Wrote {:s}".format(pargs.outfile))


if __name__ == '__main__':
    main()
