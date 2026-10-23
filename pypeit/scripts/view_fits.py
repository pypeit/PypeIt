"""
This script enables the viewing of a raw FITS file

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from IPython import embed

from pypeit import log
from pypeit.scripts import scriptbase
from pypeit import utils


class ViewFits(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        from pypeit.spectrographs.util import available_spectrographs
        parser = super().get_parser(description='View FITS files with ginga', width=width)
        parser.add_argument('spectrograph', type=str,
                            help='A valid spectrograph identifier: {0}'.format(
                                 ', '.join(available_spectrographs)))
        parser.add_argument('file', type=str, default=None, help='FITS file.  Either a Raw file or an Intermediate PypeIt file')
        parser.add_argument('--list', default=False, action='store_true',
                            help='List the extensions only?')
        parser.add_argument('--proc', default=False, action='store_true',
                            help='Process the image (i.e. orient, overscan subtract, multiply by '
                                 'gain) using pypeit.images.buildimage.')
        parser.add_argument('--bkg_file', type=str, default=None, help='FITS file to be subtracted from the image in file.'
                            '--proc must be set in order for this option to work.')

        parser.add_argument('--inter', default=False, action='store_true',
                            help='Input file is an Intermediate SciImage file')
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
        parser.add_argument('--showmask', default=False, help='Overplot masked pixels',
                            action='store_true')
        parser.add_argument('--no_clear', dest='clear', default=True, action='store_false',
                            help='If a PypeIt viewer is already open (e.g., from a previous call '
                                 'to this script), show the image in that viewer and return '
                                 'immediately.  Use --chname to avoid replacing an existing '
                                 'channel.')
        parser.add_argument('--embed', default=False, action='store_true',
                            help='Upon completion embed in ipython shell')
        return parser

    @classmethod
    def main(cls, args):

        from pathlib import Path

        from pypeit import PypeItError
        from pypeit import io

        # List only?
        if args.list:
            hdu = io.fits_open(args.file)
            print(hdu.info())
            return

        # Initialize the log
        cls.init_log(args)

        if not Path(args.file).is_file():
            raise PypeItError(f'{args.file} does not exist!')
        args.file = str(Path(args.file).absolute())
        if args.bkg_file is not None:
            args.bkg_file = str(Path(args.bkg_file).absolute())

        # NOTE: The import of the viewer is deferred so that the --list option
        # does not need to import ginga.
        from pypeit.display import viewer

        # Construct the image(s) before opening the viewer, so that any problem
        # reading the file is reported immediately.
        # NOTE: If the file is sent to an open viewer, the image is constructed
        # again by that viewer; constructing it here makes sure that any
        # problem is reported by this script, not in the other viewer.
        images = cls.build_images(args)

        # Show the image in an already open viewer?
        if not args.clear:
            if viewer.send_command('view_fits', **vars(args)):
                log.info(f'Sent {args.file} to the open PypeIt viewer.')
                return
            log.info('No open PypeIt viewer found; opening a new one.')

        v = viewer.get_viewer()
        # Allow other calls to this script (and pypeit_show_2dspec) to add
        # images to this viewer
        v.listen()
        cls.show(args, images=images)

        if args.embed:
            viewer.embed(header=utils.embed_header())

        # Wait for the user to close the window
        v.show(message='close the window to exit')

    @classmethod
    def show(cls, args, images=None):
        """
        Show the image(s) selected by the command-line arguments.

        Args:
            args (`argparse.Namespace`_):
                Parsed command-line arguments.
            images (:obj:`list`, optional):
                The images to show, as returned by :meth:`build_images`.  If
                None, they are constructed from ``args``.
        """
        from pypeit.display import display

        if images is None:
            images = cls.build_images(args)
        for img, chname in images:
            display.show_image(img, chname=chname)

    @classmethod
    def build_images(cls, args):
        """
        Construct the image(s) selected by the command-line arguments.

        Args:
            args (`argparse.Namespace`_):
                Parsed command-line arguments.

        Returns:
            :obj:`list`: List of (image, channel name) tuples with the images
            to show.
        """
        from pypeit import PypeItError
        from pypeit.spectrographs import util
        from pypeit import io
        from pypeit.images import buildimage
        from pypeit.images import pypeitimage

        if args.proc and args.exten is not None:
            raise PypeItError('You cannot specify --proc and --exten, since --exten shows the raw image')
        if args.exten is not None and args.det == 'mosaic':
            raise PypeItError('You cannot specify --exten and --det mosaic, since --mosaic displays '
                       'multiple extensions by definition')

        # Intermediate file?
        if args.inter:
            hdu = io.fits_open(args.file)
            sciImg = pypeitimage.PypeItImage.from_file(args.file)
            img = sciImg.image
        elif args.exten is not None:
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
                    raise PypeItError(f'{args.spectrograph} does not have a known mosaic')
            else:
                try:
                    _det = tuple(int(d) for d in args.det)
                except:
                    raise PypeItError(f'Could not convert detector input to integer.')
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
                    raise PypeItError(bad_read_message 
                               + f'  Original exception -- {type(e).__name__}: {str(e)}')

                if args.bkg_file is not None:
                    try:
                        bkgImg = buildimage.buildimage_fromlist(spectrograph, _det, par,
                                                                [args.bkg_file], mosaic=mosaic)
                    except Exception as e:
                        raise PypeItError(bad_read_message
                                   + f'  Original exception -- {type(e).__name__}: {str(e)}')


                    Img = Img.sub(bkgImg)

                img = Img.image

            else:
                try:
                    img = spectrograph.get_rawimage(args.file, _det)[1]
                except Exception as e:
                    raise PypeItError(bad_read_message 
                               + f'  Original exception -- {type(e).__name__}: {str(e)}')

        images = [(img, args.chname)]

        if args.showmask:
            if not args.proc:
                log.info("You need to use --proc with --showmask to show the mask.  Ignoring your argument")
            else:
                images.append((Img.bpm, "BPM"))
        return images
