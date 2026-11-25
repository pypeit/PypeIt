"""
This script displays the flat images in a ginga viewer.

The images are read and displayed inside the ginga viewer by the
``PypeItRemote`` plugin; see :mod:`pypeit.display.remote`.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit.scripts import scriptbase


class ChkFlats(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Display flat images in Ginga viewer',
                                    width=width)
        parser.add_argument('file', type=str,
                            help='PypeIt Flat file [e.g. Flat_A_1_DET01.fits]')
        parser.add_argument("--type", type=str, default='all', choices=['pixel', 'illum', 'all'],
                            help="Which flats to display.")
        parser.add_argument('--try_old', default=False, action='store_true',
                            help='Attempt to load old datamodel versions.  A crash may ensue..')
        scriptbase.add_ginga_timeout_argument(parser)
        return parser

    @classmethod
    def main(cls, args):

        from pypeit.display import remote

        # Initialize the log
        cls.init_log(args)

        # Read and show the flats in ginga
        remote.send_command(
            'chk_flats', remote.args_to_json(args, path_args=['file']), timeout=args.timeout
        )
