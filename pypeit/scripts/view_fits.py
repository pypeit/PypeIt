"""
This script enables the viewing of a raw FITS file

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from IPython import embed

from pypeit.scripts import scriptbase
from pypeit import utils


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
        parser.add_argument('--bkg_file', type=str, default=None, help='FITS file to be subtracted from the image in file.'
                            '--proc must be set in order for this option to work.')

        parser.add_argument('--exten', type=int, default=None,
                            help='Show a FITS extension in the raw file. Note --proc and --mosaic '
                                 'will not work with this option.')
        parser.add_argument('--det', type=str, default='1', nargs='*',
                            help='Detector(s) to show.  If more than one, the list of detectors, i.e. --det 4 8 '
                                 'to show detectors 4 and 8. This combination must be one of the allowed '
                                 'mosaics hard-coded for the selected '
                                 'spectrograph.  Using "mosaic" for gemini_gmos, keck_deimos, or '
                                 'keck_lris will show the mosaic of all detectors.')
        parser.add_argument('--chname', type=str, default='Image', help='Name of Ginga tab')
        parser.add_argument('--embed', default=False, action='store_true',
                            help='Upon completion embed in ipython shell')
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
#        if args.proc and args.det == 'mosaic':
#            msgs.error('You cannot specify --proc and --det mosaic, since --mosaic can only '
#                       'display the raw image mosaic')
        if args.exten is not None and args.det == 'mosaic':
            msgs.error('You cannot specify --exten and --det mosaic, since --mosaic displays '
                       'multiple extensions by definition')

        if args.exten is not None:
            hdu = io.fits_open(args.file)
            img = hdu[args.exten].data
            hdu.close()
        else:
            spectrograph = util.load_spectrograph(args.spectrograph)
            bad_read_message = 'Unable to construct image due to a read or image processing ' \
                               'error.  Use case interpreted from command-line inputs requires ' \
                               'a raw image, not an output image product from pypeit.  To show ' \
                               'a pypeit output image, specify the extension using --exten.  ' \
                               'Use --list to show the extension names.'
            if 'mosaic' in args.det:
                mosaic = True
                _det = spectrograph.default_mosaic 
                if _det is None:
                    msgs.error(f'{args.spectrograph} does not have a known mosaic')
            else:
                try:
                    _det = tuple(int(d) for d in args.det)
                except:
                    msgs.error(f'Could not convert detector input to integer.')
                mosaic = len(_det) > 1
                if not mosaic:
                    _det = _det[0]

            if args.proc:
                # Use the biasframe processing parameters because processing
                # these frames is independent of any other frames (ie., does not
                # perform bias subtraction or flat-fielding)
                par = spectrograph.default_pypeit_par()['calibrations']['biasframe']
                try:
                    Img = buildimage.buildimage_fromlist(spectrograph, _det, par,
                                                         [args.file], mosaic=mosaic)
                except Exception as e:
                    msgs.error(bad_read_message 
                               + f'  Original exception -- {type(e).__name__}: {str(e)}')

                if args.bkg_file is not None:
                    try:
                        bkgImg = buildimage.buildimage_fromlist(spectrograph, _det, par,
                                                                [args.bkg_file], mosaic=mosaic)
                    except Exception as e:
                        msgs.error(bad_read_message
                                   + f'  Original exception -- {type(e).__name__}: {str(e)}')
                    Img = Img.sub(bkgImg)

                img = Img.image

            else:
                try:
                    img = spectrograph.get_rawimage(args.file, _det)[1]
                except Exception as e:
                    msgs.error(bad_read_message 
                               + f'  Original exception -- {type(e).__name__}: {str(e)}')

        display.connect_to_ginga(raise_err=True, allow_new=True)
        display.show_image(img, chname=args.chname)


        if args.embed:
            embed(header=utils.embed_header())

