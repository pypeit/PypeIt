"""
This script enables the viewing of a raw FITS file

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from IPython import embed

from pypeit.scripts import scriptbase


class ViewFits(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        from pypeit.spectrographs import available_spectrographs
        parser = super().get_parser(description='View FITS files with ginga', width=width)
        parser.add_argument('spectrograph', type=str,
                            help='A valid spectrograph identifier: {0}'.format(
                                 ', '.join(available_spectrographs)))
        parser.add_argument('file', type = str, default = None, help = 'FITS file')
        parser.add_argument("--list", default=False, action="store_true",
                            help="List the extensions only?")
        parser.add_argument("--proc", default=False, action="store_true",
                            help='Process the image (i.e. orient, overscan subtract, multiply by '
                                 'gain) using pypeit.images.buildimage. Note det=mosaic will not '
                                 'work with this option')
        parser.add_argument('--exten', type=int, default=None,
                            help='Show a FITS extension in the raw file. Note --proc and --mosaic '
                                 'will not work with this option.')
        parser.add_argument('--det', type=str, default=1,
                            help='Detector number. To mosaic keck_deimos or keck_lris images, '
                                 'set equal to mosaic.')
        parser.add_argument('--chname', type=str, default='Image', help="Name of Ginga tab")
        return parser

    @staticmethod
    def main(args):

        from pypeit import msgs
        from pypeit.display import display
        from pypeit.spectrographs import util
        from pypeit import io
        from pypeit.images import buildimage

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
            msgs.error('You cannot specify --proc and --det mosaic, since --mosaic can only '
                       'display the raw image mosaic')
        if args.exten is not None and args.det == 'mosaic':
            msgs.error('You cannot specify --exten and --det mosaic, since --mosaic displays '
                       'multiple extensions by definition')


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

                img = buildimage.buildimage_fromlist(spectrograph, int(args.det), par,
                                                     [args.file]).image
            else:
                det = None if args.det == 'mosaic' else int(args.det)
                img = spectrograph.get_rawimage(args.file, det)[1]

        display.connect_to_ginga(raise_err=True, allow_new=True)
        display.show_image(img,chname=args.chname)


