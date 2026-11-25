"""
Show on a ginga window the archived pixel flat field image.

The images are read and displayed inside the viewer by the ``PypeItRemote``
plugin; see :mod:`pypeit.display.remote`.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit.scripts import scriptbase


class ShowPixFlat(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Show an archived Pixel Flat image in a ginga window.',
                                    width=width)
        parser.add_argument("file", type=str, help="Pixel Flat filename, e.g. pixelflat_keck_lris_blue.fits.gz")
        parser.add_argument('--det', default=None, type=int, nargs='+',
                            help='Detector(s) to show.  If more than one, list the detectors as, e.g. --det 1 2 '
                                 'to show detectors 1 and 2. If not provided, all detectors will be shown.')
        scriptbase.add_ginga_timeout_argument(parser)
        return parser

    @classmethod
    def main(cls, args):
        from pypeit import PypeItError
        from pypeit import dataPaths
        from pypeit.display import remote

        # Initialize the log
        cls.init_log(args)

        # Find the file here, so that a file in the PypeIt cache or one that
        # must be downloaded is reported in the terminal
        file_path = dataPaths.pixelflat.get_file_path(args.file, return_none=True)
        if file_path is None:
            raise PypeItError(f'Provided pixelflat file, {args.file} not found. It is not a '
                              'direct path, a cached file, or a file that can be downloaded '
                              'from a PypeIt repository.')
        args.file = str(file_path)

        # Read and show the pixel flat in ginga
        remote.send_command(
            'show_pixflat', remote.args_to_json(args, path_args=['file']), timeout=args.timeout
        )
