"""
This script displays the Tiltimg image with the traced and fitted arc-line
tilts in a ginga viewer, or in a matplotlib window.

In ginga, the data are read and displayed inside the viewer by the
``PypeItRemote`` plugin; see :mod:`pypeit.display.remote`.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit.scripts import scriptbase


class ChkTilts(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Display Tiltimg image and 2D fitted tilts in Ginga viewer '
                                                'or Matplotlib window. Tiltimg file must be in the same '
                                                'directory as Tilts.',
                                    width=width)

        parser.add_argument('file', type=str, default=None,
                            help='PypeIt Tilts file [e.g. Tilt_A_1_01.fits]')
        parser.add_argument('--mpl', default=False, action='store_true',
                            help='Use a matplotlib window instead of ginga to show the tilts. Faster plotting.')
        parser.add_argument('--show_traces', default=False, action='store_true',
                            help='Show the traced tilts. This slows down the plotting (mostly in Ginga). If not set, '
                                 'only the fitted, masked and rejected in the fit tilts are shown.')
        parser.add_argument('--try_old', default=False, action='store_true',
                            help='Attempt to load old datamodel versions.  A crash may ensue..')
        scriptbase.add_ginga_timeout_argument(parser)
        return parser

    @classmethod
    def main(cls, args):
        from pathlib import Path
        from pypeit import wavetilts
        from pypeit.display import remote

        # Initialize the log
        cls.init_log(args)

        if args.mpl:
            file = Path(args.file).absolute()
            tilts = wavetilts.WaveTilts.from_file(file, chk_version=not args.try_old)
            tilts.show(in_ginga=False, show_traces=args.show_traces, calib_dir=file.parent,
                       chk_version=not args.try_old)
            return

        # Read and show the tilts in ginga
        remote.send_command(
            'chk_tilts', remote.args_to_json(args, path_args=['file']), timeout=args.timeout
        )
