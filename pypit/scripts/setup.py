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
    parser.add_argument("files_root", type=str, help="File root or .pypit filename")
    parser.add_argument("spectrograph", type=str, help="Name of spectrograph")
    parser.add_argument("-v", "--verbosity", type=int, default=2, help="(2) Level of verbosity (0-2)")
    parser.add_argument("-d", "--develop", default=False, action='store_true', help="Turn develop debugging on")
    parser.add_argument("--extension", default='.fits',
                        help="Extension for data files.  Note any extension for compression (e.g. .gz) is not required.")
    parser.add_argument("--pypit_file", default=False, action='store_true', help='Input is the .pypit file')
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

    import pdb
    import sys
    import glob
    from pypit.scripts import run_pypit
    from pypit import pyputils
    from pypit.pypit import load_input
    import datetime
    import yaml

    # Check for existing setup file
    setup_files = glob.glob('./{:s}*.setup'.format(args.spectrograph))
    if len(setup_files) > 0:
        print("Working directory already includes a .setup file")
        for ifile in setup_files:
            print("Remove: {:s}".format(ifile))
        print("Then you can re-run this script")
        sys.exit()

    # Generate a dummy .pypit file
    if not args.pypit_file:
        # Name
        date = str(datetime.date.today().strftime('%Y-%b-%d'))
        root = args.spectrograph+'_'+date
        pyp_file = root+'.pypit'
        # Generate
        dfname = "{:s}*{:s}*".format(args.files_root, args.extension)
        pyputils.make_pypit_file(pyp_file, args.spectrograph,
                              [dfname], setup_script=True)
        print("Wrote {:s}".format(pyp_file))
    else:
        pyp_file = args.files_root

    # Run
    pinp = [pyp_file]
    if args.develop:
        pinp += ['-d']
    pargs = run_pypit.parser(pinp)
    run_pypit.main(pargs)

    # #####################
    # Generate custom .pypit files

    # Read group file
    group_file = pyp_file.replace('.pypit', '.group')
    with open(group_file, 'r') as infile:
        group_dict = yaml.load(infile)
    groups = list(group_dict.keys())
    groups.sort()

    # Read master file
    from pypit import pyputils
    msgs = pyputils.get_dummy_logger()
    parlines, datlines, spclines, dfnames = load_input(pyp_file, msgs)
    # Remove run setup from parlines
    for jj,parline in enumerate(parlines):
        if 'run setup' in parline:
            parlines[jj] = 'run setup False\n'

    # Generate .pypit files
    for group in groups:
        root = args.spectrograph+'_setup_'
        pyp_file = root+group+'.pypit'

        pyputils.make_pypit_file(pyp_file, args.spectrograph, dfnames,
                                 parlines=parlines,
                                 spclines=spclines,
                                 calcheck=True)
        print("Wrote {:s}".format(pyp_file))

