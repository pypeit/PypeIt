"""
This script displays the flat images in an RC Ginga window.
"""

def parse_args(options=None, return_parser=False):

    import argparse

    parser = argparse.ArgumentParser(description='Display MasterFlat images in a previously '
                                                 'launched RC Ginga viewer',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--type", type=str, default='all', help="Which flats to display. Must be one of: pixel, illum, all")
    parser.add_argument('master_file', type=str,
                        help='PypeIt MasterFlat file [e.g. MasterFlat_A_1_01.fits]')
    parser.add_argument('--try_old', default=False, action='store_true',
                        help='Attempt to load old datamodel versions.  A crash may ensue..')

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(args):

    from pypeit import flatfield

    # Load
    flatImages = flatfield.FlatImages.from_file(args.master_file, chk_version=(not args.try_old))
    flatImages.show(args.type)

