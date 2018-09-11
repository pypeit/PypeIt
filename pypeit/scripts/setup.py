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

import argparse

def parser(options=None):
    parser = argparse.ArgumentParser(description="Script to setup a PypeIt run [v2]")
    parser.add_argument("files_root", type=str, help="File path+root, e.g. /data/Kast/b ")
    parser.add_argument("spectrograph", type=str, help="Name of spectrograph")
    parser.add_argument("-v", "--verbosity", type=int, default=2,
                        help="(2) Level of verbosity (0-2)")
    parser.add_argument("--extension", default='.fits',
                        help='File extension; compression indicators (e.g. .gz) not required.')
    parser.add_argument("--pypeit_file", default=False, action='store_true',
                        help='Input is the .pypeit file')
    parser.add_argument("--setups_path", default=None,
                        help='Path to top-level folder.  Default is current working directory.')
    parser.add_argument("-c", "--custom", default=False, action='store_true',
                        help='Generate custom folders and pypeit files?')
    parser.add_argument('-o', '--overwrite', default=False, action='store_true',
                        help='Overwrite any existing files/directories')
#    parser.add_argument("-q", "--quick", default=False, help="Quick reduction",
#                        action="store_true")
#    parser.add_argument("-c", "--cpus", default=False,
#                        help="Number of CPUs for parallel processing", action="store_true")
#    parser.print_help()

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(args):

    import os
    import pdb as debugger

    from pypeit import msgs
    from pypeit.spectrographs.util import valid_spectrographs
    from pypeit.spectrographs.util import load_spectrograph
    from pypeit import pypeit

    # Check that input spectrograph is supported
    instruments_served = valid_spectrographs()
    if args.spectrograph not in instruments_served:
        raise ValueError('Instrument \'{0}\' unknown to PypeIt.\n'.format(args.spectrograph)
                         + '\tAvailable options are: {0}\n'.format(', '.join(instruments_served))
                         + '\tSelect an available instrument or consult the documentation '
                         + 'on how to add a new instrument.')

    # Instantiate
    spectrograph = load_spectrograph(args.spectrograph)

    # Generate PypeIt
    pypeIt = pypeit.instantiate_me(spectrograph, verbosity=args.verbosity,
                                   overwrite=args.overwrite, setups_path=args.setups_path)

    pypeIt.build_setup_files(args.files_root)

    # Custom?
    if not args.custom:
        return
    else:
        pypeIt.build_custom_pypeitfiles()


