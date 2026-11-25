"""
This script displays the Trace image and the traces in a ginga viewer, or in a
matplotlib window.

In ginga, the data are read and displayed inside the viewer by the
``PypeItRemote`` plugin; see :mod:`pypeit.display.remote`.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit.scripts import scriptbase


class ChkEdges(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Display trace image and edge traces',
                                    width=width)

        parser.add_argument('trace_file', type=str, default=None,
                            help='PypeIt Edges file [e.g. Edges_A_0_DET01.fits.gz]')
        parser.add_argument('--slits_file', type=str, default=None,
                            help='PypeIt Slits file [e.g. Slits_A_1_01.fits]. If this file does '
                                 'not exist or is not provided, PypeIt will attempt to read the '
                                 'default file name (in the calibration directory recorded in '
                                 'the Edges file, or else in the directory of the Edges file).  '
                                 'Ignored if plotting using a matplotlib window instead of '
                                 'ginga.')
        parser.add_argument('--mpl', default=False, action='store_true',
                            help='Use a matplotlib window instead of ginga to show the trace')
        parser.add_argument('--try_old', default=False, action='store_true',
                            help='Attempt to load old datamodel versions.  A crash may ensue..')
        scriptbase.add_ginga_timeout_argument(parser)
        return parser

    @classmethod
    def main(cls, args):

        from pypeit import edgetrace
        from pypeit.display import remote

        # Initialize the log
        cls.init_log(args)

        if args.mpl:
            edges = edgetrace.EdgeTraceSet.from_file(args.trace_file,
                                                     chk_version=not args.try_old)
            edges.show(thin=10, include_img=True, idlabel=True)
            return

        # Read and show the edges in ginga
        remote.send_command(
            'chk_edges', remote.args_to_json(args, path_args=['trace_file', 'slits_file']),
            timeout=args.timeout
        )
