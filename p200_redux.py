# Automatic Reduction Pipeline for P200 DBSP.

import argparse
import os

import p200_arm_redux
import pypeit.scripts.setup as setup
#import table_edit


def parser(options=None):
    # Define command line arguments.
    parser = argparse.ArgumentParser(description="Automatic Data Reduction Pipeline for P200 DBSP")
    
    # Argument for fully-automatic (i.e. nightly) or with user-checking file typing
    parser.add_argument('-i', '--interactive', type=bool, default=True,
                        help='Interactive file-checking?')

    # Argument for input file directory
    parser.add_argument('-r', '--root', type=str, default=None,
                        help='File path+root, e.g. /data/DBSP_20200127')

    parser.add_argument('-d', '--output_path', default=None,
                        help='Path to top-level output directory.  '
                             'Default is the current working directory.')
    
    # Argument for specifying only red/blue
    
    return parser.parse_args() if options is None else parser.parse_args(options)

def interactive_correction(fitstbl, arm):
    # this needs to actually fix the data's FITS headers
    # deleting entire row from .pypeit file is valid though

    # function for interactively correcting the fits table
    fitstbl.table.sort('filename')
    fitstbl.table.write(f'table_{arm}.dat', format='ascii')
    #table_edit.main(fitstbl)

def main(args):
    pass
    # Build
    options_blue = {
        'root': os.path.join(args.root, 'Blue'),
        'spectrograph': 'p200_dbsp_blue',
        'output_path': args.output_path,
        'extension': '.fits',
        'background': None,
        'cfg_split': 'all'
    }
    options_red = options_blue.copy()
    options_red['spectrograph'] = 'p200_dbsp_red'
    options_red['root'] = os.path.join(args.root, 'Red')

    red_fitstbl, context = p200_arm_redux.setup(options_red)
    # optionally use interactive correction
    if args.interactive:
        interactive_correction(red_fitstbl, 'red')
    pypeit_file_red = p200_arm_redux.write_setup(options_red, context)[0]
    
    blue_fitstbl, context = p200_arm_redux.setup(options_blue)
    if args.interactive:
        interactive_correction(blue_fitstbl, 'blue')
    pypeit_file_blue = p200_arm_redux.write_setup(options_blue, context)[0]

    # Grab pypeit file from write_setup
    options_red['pypeit_file'] = pypeit_file_red
    options_blue['pypeit_file'] = pypeit_file_blue
    #p200_arm_redux.redux(options_red)
    p200_arm_redux.redux(options_blue)


if __name__ == '__main__':
    # args = parser(['-r', '/Users/milan/my_scr2/DBSP_20200716', '-d', '/Users/milan/my_scr2'])
    args = parser(['-r', '/Users/milan/Documents/GitHub/DBSP_20200716', '-d', '/Users/milan/Documents/GitHub/DBSP_20200716'])
    # args = parser(['-r', '/Users/milan/Documents/GitHub/DBSP_20200127'])
    main(args)
