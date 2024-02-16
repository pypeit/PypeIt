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

        parser.add_argument('trace_file', type=str, default=None,
                            help='PypeIt Edges file [e.g. Edges_A_0_DET01.fits.gz]')
        parser.add_argument('--slits_file', type=str, default=None,
                            help='PypeIt Slits file [e.g. Slits_A_1_01.fits]. If this file does '
                                 'not exist or is not provided, PypeIt will attempt to read the '
                                 'default file name (in the Calibrations directory).  Ignored if '
                                 'plotting using a matplotlib window instead of ginga.')
        parser.add_argument('--mpl', default=False, action='store_true',
                            help='Use a matplotlib window instead of ginga to show the trace')
        parser.add_argument('--try_old', default=False, action='store_true',
                            help='Attempt to load old datamodel versions.  A crash may ensue..')
        return parser

    @staticmethod
    def main(args):

        from pathlib import Path

        from pypeit import edgetrace, slittrace, msgs

        chk_version = not args.try_old

        # Load
        edges = edgetrace.EdgeTraceSet.from_file(args.trace_file, chk_version=chk_version)

        if args.mpl:
            edges.show(thin=10, include_img=True, idlabel=True)
            return

        # Set the Slits file name
        slit_filename = args.slits_file
        if slit_filename is not None:
            # File provided by user
            slit_filename = Path(args.slits_file).resolve()
            if not slit_filename.exists():
                # But doesn't exist
                msgs.warn(f'{slit_filename} does not exist!')
                # Set the file name to None so that the code will try to find
                # the default file
                slit_filename = None
        if slit_filename is None:
            slit_filename = slittrace.SlitTraceSet.construct_file_name(
                                edges.traceimg.calib_key, calib_dir=edges.traceimg.calib_dir)
            if not slit_filename.exists():
                msgs.warn(f'{slit_filename} does not exist!')
        # NOTE: At this point, slit_filename *must* be a Path object

        slits = slittrace.SlitTraceSet.from_file(slit_filename, chk_version=chk_version) \
                    if slit_filename.exists() else None
        edges.show(thin=10, in_ginga=True, slits=slits)


