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
    flatField = flatfield.FlatImages.from_file(pargs.master_file)
    master_key, master_dir = masterframe.grab_key_mdir(pargs.master_file)
    try:
        slit_masterframe_name = masterframe.construct_file_name(slittrace.SlitTraceSet, master_key,
                                                                master_dir=master_dir)
        slits = slittrace.SlitTraceSet.from_file(slit_masterframe_name)
    except:
        msgs.warn('Could not load slits to show with flat-field images. Did you provide the master info??')
        slits = None
    # Show
    # TODO: Add wcs_match as command-line argument?
    flatfield.show_flats(flatField.pixelflat, flatField.illumflat, flatField.procflat, flatField.flat_model,
                         slits=slits)
