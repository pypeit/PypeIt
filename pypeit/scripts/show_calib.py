"""
This script displays a single-image calibration file (e.g. Bias, Arc, Dark,
Tiltimg or Alignment) in a ginga viewer.

The file is read and displayed inside the ginga viewer by the ``PypeItRemote``
plugin; see :mod:`pypeit.display.remote`.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
from pypeit.scripts import scriptbase


class ShowCalib(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(
            description='Display a calibration image (e.g. Bias, Arc, Dark, Tiltimg or '
                        'Alignment) in a ginga viewer, after clearing the viewer.',
            width=width
        )
        parser.add_argument('file', type=str,
                            help='PypeIt calibration file [e.g. Arc_A_0_DET01.fits]')
        parser.add_argument('--ivar', default=False, action='store_true',
                            help='Add a channel with the inverse variance image')
        parser.add_argument('--mask', default=False, action='store_true',
                            help='Add a channel with the mask bits')
        parser.add_argument('--slits', default=False, action='store_true',
                            help='Overlay the slit edges.  The Slits file is looked for in the '
                                 'calibration directory recorded in the file, and then in the '
                                 'directory of the file.')
        parser.add_argument('--try_old', default=False, action='store_true',
                            help='Attempt to load old datamodel versions.  If the file cannot '
                                 'be parsed, its first image is shown.')
        scriptbase.add_ginga_timeout_argument(parser)
        return parser

    @classmethod
    def main(cls, args):

        from pypeit.display import remote

        # Initialize the log
        cls.init_log(args)

        # Read and show the calibration image in ginga
        remote.send_command(
            'show_calib', remote.args_to_json(args, path_args=['file']), timeout=args.timeout
        )
