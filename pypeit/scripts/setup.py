#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-


"""
This script generates files to setup a PypeIt run
"""
from __future__ import (print_function, absolute_import, division,
                        unicode_literals)

import argparse
from pypeit.spectrographs.util import valid_spectrographs

def parser(options=None):
    # TODO: Add argument that specifies the log file
    parser = argparse.ArgumentParser(description="Script to setup a PypeIt run [v3]")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-r', '--root', type=str, default=None,
                       help='File path+root, e.g. /data/Kast/b ')
    #group.add_argument('-p', '--pypeit_file', type=str, default=None,
    #                   help='PypeIt file to use')
    parser.add_argument('-s', '--spectrograph', default=None, type=str,
                        help='A valid spectrograph identifier: {0}'.format(
                                ', '.join(valid_spectrographs())))
    parser.add_argument('-e', '--extension', default='.fits',
                        help='File extension; compression indicators (e.g. .gz) not required.')
    parser.add_argument('-d', '--output_path', default=None,
                        help='Path to top-level output directory.  '
                             'Default is the current working directory.')
    parser.add_argument('-o', '--overwrite', default=False, action='store_true',
                        help='Overwrite any existing files/directories')
    parser.add_argument('-c', '--cfg_split', default=None, type=str,
                        help='Generate the PypeIt files and folders by input configuration. [all or A,B or B,D,E or E]')
    parser.add_argument('-b', '--background', default=False, action='store_true',
                        help='Include the background-pair columns for the user to edit')
    parser.add_argument('-v', '--verbosity', type=int, default=2,
                        help='Level of verbosity from 0 to 2; default is 2.')
#    parser.add_argument('-q', '--quick', default=False, help='Quick reduction',
#                        action='store_true')
#    parser.add_argument('-c', '--cpus', default=False,
#                        help='Number of CPUs for parallel processing', action='store_true')
#    parser.print_help()

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(args):

    import os
    import pdb as debugger

    from pypeit import msgs
    from pypeit.pypeitsetup import PypeItSetup

    # Check that the spectrograph is provided if using a file root
    if args.root is not None:
        if args.spectrograph is None:
            raise ValueError('Must provide spectrograph identifier with file root.')
        # Check that input spectrograph is supported
        instruments_served = valid_spectrographs()
        if args.spectrograph not in instruments_served:
            raise ValueError('Instrument \'{0}\' unknown to PypeIt.\n'.format(args.spectrograph)
                             + '\tOptions are: {0}\n'.format(', '.join(instruments_served))
                             + '\tSelect an available instrument or consult the documentation '
                             + 'on how to add a new instrument.')

    # Get the output directory
    output_path = os.getcwd() if args.output_path is None else args.output_path
    sort_dir = os.path.join(output_path, 'setup_files')
   
    # Initialize PypeItSetup based on the arguments
    if args.root is not None:
        ps = PypeItSetup.from_file_root(args.root, args.spectrograph, extension=args.extension,
                                        output_path=sort_dir)
    else:
        # Should never reach here
        raise IOError('Need to set -r !!')

    # Run the setup
    ps.run(setup_only=True, sort_dir=sort_dir, write_bkg_pairs=args.background)

    # Use PypeItMetaData to write the complete PypeIt file
    if args.cfg_split is not None:
        pypeit_file = os.path.join(output_path, '{0}.pypeit'.format(args.spectrograph))
        config_list = args.cfg_split.split(',')
        ps.fitstbl.write_pypeit(pypeit_file, cfg_lines=ps.user_cfg, write_bkg_pairs=args.background,
                                configs=config_list)
    return 0

