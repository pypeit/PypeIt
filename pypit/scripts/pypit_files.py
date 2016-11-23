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
    parser.add_argument("--extension", default='.fits', action='store_true',
                        help="Extension for data files.  Note any extension for compression (e.g. .gz) is not required.")
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

    import sys
    import glob
    from pypit.scripts import run_pypit
    from pypit import utils as pyp_utils
    import datetime

    # Check for existing setup file
    setup_files = glob.glob('./{:s}*.setup'.format(args.spectrograph))
    if len(setup_files) > 0:
        print("Working directory already includes a .setup file")
        for ifile in setup_files:
            print("Remove: {:s}".format(ifile))
        print("Then you can re-run this script")
        sys.exit()

    # Generate a dummy .pypit file
    date = str(datetime.date.today().strftime('%Y-%b-%d'))
    root = args.spectrograph+'_'+date
    pyp_file = root+'.pypit'

    pyp_utils.make_pypit_file(pyp_file, args.spectrograph,
                              [args.files_root], args.extension)
    print("Wrote {:s}".format(pyp_file))

    # Run
    pinp = [pyp_file]
    if args.develop:
        pinp += ['-d']
    args = run_pypit.parser(pinp)
    run_pypit.main(args)
