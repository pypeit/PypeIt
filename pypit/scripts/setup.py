#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-


"""
This script generates files to setup a PYPIT run
"""
from __future__ import (print_function, absolute_import, division,
                        unicode_literals)

try:
    from xastropy.xutils import xdebug as debugger
except:
    import pdb as debugger

def parser(options=None):
    import argparse

    parser = argparse.ArgumentParser(description="Script to setup a PYPIT run")
    parser.add_argument("files_root", type=str, help="File root")
    parser.add_argument("spectrograph", type=str, help="Name of spectrograph")
    parser.add_argument("-v", "--verbosity", type=int, default=2, help="(2) Level of verbosity (0-2)")
    parser.add_argument("-d", "--develop", default=False, action='store_true', help="Turn develop debugging on")
    #parser.add_argument("-q", "--quick", default=False, help="Quick reduction", action="store_true")
    #parser.add_argument("-c", "--cpus", default=False, help="Number of CPUs for parallel processing", action="store_true")
    #parser.print_help()

    if options is None:
        pargs = parser.parse_args()
    else:
        pargs = parser.parse_args(options)
    #
    return pargs


def main(args):

    import sys, os
    from pypit import pypit
    from pypit.scripts import run_pypit
    import traceback
    import datetime

    # Generate a dummy .pypit file
    date = str(datetime.date.today().strftime('%Y-%b-%d'))
    root = args.spectrograph+'_'+date
    pyp_file = root+'.pypit'

    with open(pyp_file, 'w') as f:
        f.write("# This is a comment line\n")
        f.write("\n")
        f.write("# Change the default settings\n")
        f.write("run ncpus 1\n")
        f.write("run calcheck True\n")
        f.write("run spectrograph {:s}\n".format(args.spectrograph))
        f.write("output overwrite True\n")
        f.write("output sorted {:s}\n".format(root))
        f.write("\n")
        f.write("# Read in the data\n")
        f.write("data read\n")
        f.write(" {:s}*\n".format(args.files_root))
        f.write("data end\n")
        f.write("\n")
        f.write("spect read\n")
        f.write(" pixelflat number 0\n")
        f.write(" arc number 1\n")
        f.write(" slitflat number 0\n")
        f.write(" bias number 0\n")
        f.write(" standard number 0\n")
        f.write("spect end\n")
    print("Wrote {:s}".format(pyp_file))

    # Run
    args = run_pypit.parser([pyp_file])
    run_pypit.main(args)
