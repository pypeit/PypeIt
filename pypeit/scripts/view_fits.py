#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
This script enables the viewing of a raw FITS file
"""

import subprocess

from astropy.io import fits

from pypeit import msgs
from pypeit.spectrographs import keck_lris
from pypeit.spectrographs import keck_deimos
from pypeit.spectrographs import gemini_gmos
from pypeit.display import display
from pypeit.spectrographs import mmt_binospec
from pypeit.spectrographs import mmt_mmirs
from pypeit.spectrographs import mmt_bluechannel
from pypeit.spectrographs import util
from pypeit import msgs
from pypeit import io
from IPython import embed
from pypeit.images import buildimage

def parse_args(options=None, return_parser=False):
    import argparse
    from pypeit.spectrographs import available_spectrographs

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('spectrograph', type=str,
                        help='A valid spectrograph identifier: {0}'.format(
                             ', '.join(available_spectrographs)))
    parser.add_argument('file', type = str, default = None, help = 'FITS file')
    parser.add_argument("--list", default=False, help="List the extensions only?", action="store_true")
    parser.add_argument("--proc", default=False,
                        help="Process the image (i.e. orient, overscan subtract, multiply by gain) using pypeit.images.buildimage. "
                             "Note det=mosaic will not work with this option", action="store_true")
    parser.add_argument('--exten', type=int, default = None, help="Show a FITS extension in the raw file. Note --proc and --mosaic will not work with this option")
    parser.add_argument('--det', type=str, default=1, help="Detector number. To mosaic keck_deimos or keck_lris images, set equal to mosaic")
    parser.add_argument('--chname', type=str, default='Image', help="Name of Ginga tab")

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(args):

    # List only?
    if args.list:
        hdu = io.fits_open(args.file)
        print(hdu.info())
        return

    # Setup for PYPIT imports
    msgs.reset(verbosity=2)

    if args.proc and args.exten is not None:
        msgs.error('You cannot specify --proc and --exten, since --exten shows the raw image')
    if args.proc and args.det == 'mosaic':
        msgs.error('You cannot specify --proc and --det mosaic, since --mosaic can only display the raw image mosaic')
    if args.exten is not None and args.det == 'mosaic':
        msgs.error('You cannot specify --exten and --det mosaic, since --mosaic displays multiple extensions by definition')


    if args.exten is not None:
        hdu = io.fits_open(args.file)
        img = hdu[args.exten].data
    else:
        spectrograph = util.load_spectrograph(args.spectrograph)
        if args.proc:
            # Use the arc FramePar since this does not do processing
            par = spectrograph.default_pypeit_par()['calibrations']['arcframe']
            par['process']['use_darkimage'] = False
            par['process']['use_biasimage'] = False
            par['process']['mask_cr'] = False
            par['process']['cr_sigrej'] = -1


            img = buildimage.buildimage_fromlist(spectrograph, int(args.det), par, [args.file]).image
        else:
            det = None if args.det == 'mosaic' else int(args.det)
            img = spectrograph.get_rawimage(args.file, det)[1]

    display.connect_to_ginga(raise_err=True, allow_new=True)
    display.show_image(img,chname=args.chname)

def entry_point():
    main(parse_args())


if __name__ == '__main__':
    entry_point()
