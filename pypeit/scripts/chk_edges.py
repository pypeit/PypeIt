"""
This script displays the Trace image and the traces in an RC Ginga window.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit.scripts import scriptbase


class ChkEdges(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Display trace image and edge traces',
                                    width=width)

        parser.add_argument('trace_file', type=str, default = None,
                            help='PypeIt Edges file [e.g. Edges_A_0_DET01.fits.gz]')
        parser.add_argument('--mpl', default=False, action='store_true',
                            help='Use a matplotlib window instead of ginga to show the trace')
        parser.add_argument('--try_old', default=False, action='store_true',
                            help='Attempt to load old datamodel versions.  A crash may ensue..')
        return parser

    @staticmethod
    def main(pargs):
        from pypeit import edgetrace, slittrace

        # Load
        edges = edgetrace.EdgeTraceSet.from_file(pargs.trace_file, chk_version=(not pargs.try_old))

        if pargs.mpl:
            edges.show(thin=10, include_img=True, idlabel=True)
        else:
            # check if SlitTraceSet file exists(this allows to show the orders id instead of slitids)
            slit_filename = slittrace.SlitTraceSet.construct_file_name(edges.traceimg.calib_key,
                                                                       calib_dir=edges.traceimg.calib_dir)
            slits = slittrace.SlitTraceSet.from_file(slit_filename, chk_version=(not pargs.try_old)) \
                if slit_filename.exists() else None
            edges.show(thin=10, in_ginga=True, slits=slits)
        return 0


