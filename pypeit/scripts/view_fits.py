"""
This script enables the viewing of a raw FITS file.

The image is read and displayed inside the ginga viewer by the
``PypeItRemote`` plugin; see :mod:`pypeit.display.remote` and
:func:`pypeit.display.image_data.view_fits_images`.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit.scripts import scriptbase


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
                            help='Do *not* clear all existing channels')
        scriptbase.add_ginga_timeout_argument(parser)
        return parser

    @classmethod
    def main(cls, args):
        from pypeit import io
        from pypeit.display import remote

        # List only?
        if args.list:
            with io.fits_open(args.file) as hdu:
                hdu.info()
            return

        # Initialize the log
        cls.init_log(args)

        # Read and show the image in ginga
        remote.send_command(
            'view_fits', remote.args_to_json(args, path_args=['file', 'bkg_file']),
            timeout=args.timeout
        )
