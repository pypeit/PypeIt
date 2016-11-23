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
    import yaml
    from pypit import pyputils
    msgs = pyputils.get_dummy_logger()
    from pypit.pypit import load_input
    from pypit.arsort import get_setup_file, load_setup

    # Check for existing setup file
    setup_file, nexist = get_setup_file(spectrograph=args.spectrograph)
    if nexist > 1:
        print("Working directory includes multiple .setup files")
        print("Only 1 is allowed.  You should probably start again with pypit_setup")
        print("Then you can re-run this script")
        sys.exit()
    elif nexist == 1:
        pass
    else:
        print("Working directory does not include a .setup file")
        print("Generate one first with pypit_setup")

    # Load msgs

    # Read master file
    parlines, datlines, spclines, dfnames = load_input(args.master_file, msgs)

    # Modify setup line (if present)
    for jj,parline in enumerate(parlines):
        if 'run setup' in parline:
            parlines[jj] = 'run setup False\n'
    # Read setup (may not need the dict)
    setup_dict = load_setup(spectrograph=args.spectrograph)
    setups = setup_dict.keys()
    setups.sort()

    # Read group file
    group_file = setup_file.replace('.setup', '.group')
    with open(group_file, 'r') as infile:
        group_dict = yaml.load(infile)
    groups = group_dict.keys()
    groups.sort()

    # Generate .pypit files
    for group in groups:
        root = args.spectrograph+'_setup_'
        pyp_file = root+group+'.pypit'

        pyputils.make_pypit_file(pyp_file, args.spectrograph, dfnames,
                                 parlines=parlines,
                                 spclines=spclines,
                                 calcheck=True)
        print("Wrote {:s}".format(pyp_file))

