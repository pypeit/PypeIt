#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
This script generates files to setup a PypeIt run
"""
import os
from pypeit.spectrographs import available_spectrographs

def parse_args(options=None, return_parser=False):
    import argparse

    # TODO: Add argument that specifies the log file
    parser = argparse.ArgumentParser(description='Parse data files to construct a pypeit file in '
                                                 'preparation for reduction using \'run_pypeit\'',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # TODO: Make root and spectrograph required arguments
    parser.add_argument('-r', '--root', type=str, default=None,
                       help='File path+root, e.g. /data/Kast/b ')
    parser.add_argument('-s', '--spectrograph', default=None, type=str,
                        help='A valid spectrograph identifier: {0}'.format(
                                ', '.join(available_spectrographs)))

    parser.add_argument('-e', '--extension', default='.fits',
                        help='File extension; compression indicators (e.g. .gz) not required.')
    parser.add_argument('-d', '--output_path', default=os.getcwd(),
                        help='Path to top-level output directory.')
    parser.add_argument('-o', '--overwrite', default=False, action='store_true',
                        help='Overwrite any existing files/directories')
    parser.add_argument('-c', '--cfg_split', default=None, type=str,
                        help='Generate the PypeIt files and folders by input configuration.  To '
                             'write all unique configurations identifed, use \'all\', otherwise '
                             'provide the list of configuration letters; e.g., \'A,B\' or '
                             '\'B,D,E\' or \'E\'.')
    parser.add_argument('-b', '--background', default=False, action='store_true',
                        help='Include the background-pair columns for the user to edit')
    parser.add_argument('-v', '--verbosity', type=int, default=2,
                        help='Level of verbosity from 0 to 2.')

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(args):

    from pypeit.pypeitsetup import PypeItSetup

    if args.root is None:
        raise IOError('root is a required argument.  Use the -r, --root command-line option.')
    if args.spectrograph is None:
        raise IOError('spectrograph is a required argument.  Use the -s, --spectrograph '
                      'command-line option.')

    # Check that input spectrograph is supported
    if args.spectrograph not in available_spectrographs:
        raise ValueError('Instrument \'{0}\' unknown to PypeIt.\n'.format(args.spectrograph)
                         + '\tOptions are: {0}\n'.format(', '.join(available_spectrographs))
                         + '\tSelect an available instrument or consult the documentation '
                         + 'on how to add a new instrument.')

    # Get the output directory
    sort_dir = os.path.join(args.output_path, 'setup_files')
   
    # Initialize PypeItSetup based on the arguments
    ps = PypeItSetup.from_file_root(args.root, args.spectrograph, extension=args.extension,
                                    output_path=sort_dir)
    # Run the setup
    ps.run(setup_only=True, sort_dir=sort_dir, write_bkg_pairs=args.background)

    # Use PypeItMetaData to write the complete PypeIt file
    # TODO: Set cfg_split to 'all' by default?
    if args.cfg_split is not None:
        ps.fitstbl.write_pypeit(args.output_path, cfg_lines=ps.user_cfg,
                                write_bkg_pairs=args.background,
                                configs=[item.strip() for item in args.cfg_split.split(',')])

