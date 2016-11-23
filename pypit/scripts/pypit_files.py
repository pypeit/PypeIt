#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-


"""
This script generates .pypit files for a PYPIT reduction run
"""
from __future__ import (print_function, absolute_import, division,
                        unicode_literals)

try:
    from xastropy.xutils import xdebug as debugger
except:
    import pdb as debugger

def parser(options=None):
    import argparse

    parser = argparse.ArgumentParser(description="Script to generate .pypit files")
    parser.add_argument("master_file", type=str, help="Master .pypit file")
    parser.add_argument("spectrograph", type=str, help="Name of spectrograph")

    if options is None:
        pargs = parser.parse_args()
    else:
        pargs = parser.parse_args(options)
    #
    return pargs


def main(args):

    import sys
    import glob
    from pypit import pyputils

    # Check for existing setup file
    setup_files = glob.glob('./{:s}*.setup'.format(args.spectrograph))
    if len(setup_files) > 1:
        print("Working directory includes multiple .setup files")
        print("Only 1 is allowed.  You should probably start again with pypit_setup")
        print("Then you can re-run this script")
        sys.exit()
    elif len(setup_files) == 1:
        setup_file = setup_files[0]
    else:
        print("Working directory does not include a .setup file")
        print("Generate one first with pypit_setup")

    # Load msgs
    msgs = pyputils.get_dummy_logger()
    from pypit.pypit import load_input

    # Read master file
    parlines, datlines, spclines, dfnames, skip_files = load_input(args.master_file, msgs)
    debugger.set_trace()

    # Generate .pypit files
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
