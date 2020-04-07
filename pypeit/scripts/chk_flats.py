"""
This script displays the flat images in an RC Ginga window.
"""
import argparse

from astropy.io import fits

from pypeit import flatfield
from pypeit import slittrace
from pypeit import masterframe
from pypeit import msgs
from IPython import embed


def parser(options=None):
    parser = argparse.ArgumentParser(description='Display MasterFlat images in a previously '
                                                 'launched RC Ginga viewer',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('master_file', type=str,
                        help='PypeIt MasterFlat file [e.g. MasterFlat_A_1_01.fits]')
    return parser.parse_args() if options is None else parser.parse_args(options)


def main(pargs):
    # Load
    flatImages = flatfield.FlatImages.from_file(pargs.master_file)
    flatImages.show()
