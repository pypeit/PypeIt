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
        parser.add_argument('file', type=str, default=None, help='FITS file')
        parser.add_argument('--list', default=False, action='store_true',
                            help='List the extensions only?')
        parser.add_argument('--proc', default=False, action='store_true',
                            help='Process the image (i.e. orient, overscan subtract, multiply by '
                                 'gain) using pypeit.images.buildimage. Note det=mosaic will not '
                                 'work with this option')
        parser.add_argument('--exten', type=int, default=None,
                            help='Show a FITS extension in the raw file. Note --proc and --mosaic '
                                 'will not work with this option.')
        parser.add_argument('--det', type=str, default=1,
                            help='Detector number. To mosaic keck_deimos or keck_lris images, '
                                 'set equal to mosaic.')
        parser.add_argument('--chname', type=str, default='Image', help='Name of Ginga tab')
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
            bad_read_message = 'Unable to construct image due to a read or image processing ' \
                               'error.  Use case interpreted from command-line inputs requires ' \
                               'a raw image, not an output image product from pypeit.  To show ' \
                               'a pypeit output image, specify the extension using --exten.  ' \
                               'Use --list to show the extension names.'
            if args.proc:
                # Use the biasframe processing parameters because processing
                # these frames is independent of any other frames (ie., does not
                # perform bias subtraction or flat-fielding)
                par = spectrograph.default_pypeit_par()['calibrations']['biasframe']
                try:
                    img = buildimage.buildimage_fromlist(spectrograph, int(args.det), par,
                                                        [args.file]).image
                except Exception as e:
                    msgs.error(bad_read_message 
                               + f'  Original exception -- {type(e).__name__}: {str(e)}')
            else:
                det = None if args.det == 'mosaic' else int(args.det)
                try:
                    img = spectrograph.get_rawimage(args.file, det)[1]
                except Exception as e:
                    msgs.error(bad_read_message 
                               + f'  Original exception -- {type(e).__name__}: {str(e)}')

        display.connect_to_ginga(raise_err=True, allow_new=True)
        display.show_image(img,chname=args.chname)


