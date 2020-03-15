"""
This script displays the flat images in an RC Ginga window.
"""
import argparse
from pypeit import flatfield


def parser(options=None):
    parser = argparse.ArgumentParser(description='Display MasterFlat images in a previously '
                                                 'launched RC Ginga viewer',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('master_file', type=str,
                        help='PYPIT MasterFlat file [e.g. MasterFlat_A_1_01.fits]')
    return parser.parse_args() if options is None else parser.parse_args(options)


def main(pargs):
    flatField = flatfield.FlatField.from_master_file(pargs.master_file)
    # TODO: Add show_slits and wcs_match as command-line arguments?
    flatField.show()

#    import time
#
#    from pypeit import ginga
#    from pypeit import flatfield
#
#    import subprocess
#
#    # Load up
#    flatField = flatfield.FlatField.from_master_file(pargs.master_file)
#
#    try:
#        ginga.connect_to_ginga(raise_err=True)
#    except ValueError:
#        subprocess.Popen(['ginga', '--modules=RC'])
#        time.sleep(3)
#
#    # Show RawFlatImage
#    viewer, ch = ginga.show_image(flatField.rawflatimg.image, chname='Raw Flat')
#    # PixelFlat
#    if flatField.mspixelflat is not None:
#        viewer, ch = ginga.show_image(flatField.mspixelflat, chname='Pixel Flat')
#    # Illumination flat
#    if flatField.msillumflat is not None:
#        viewer, ch = ginga.show_image(flatField.msillumflat, chname='Illumination Flat')
#    # Illumination flat
#    if flatField.flat_model is not None:
#        viewer, ch = ginga.show_image(flatField.flat_model, chname='Flat Model')
#
#    print("Check your Ginga viewer")


