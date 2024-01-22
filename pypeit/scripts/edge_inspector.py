"""
Launch a matplotlib GUI to inspect and interact with the slit edge tracing.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
from IPython import embed

from pypeit.scripts import scriptbase

class EdgeInspector(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        from pathlib import Path
        parser = super().get_parser(description='Interactively inspect/edit slit edge traces',
                                    width=width)
        parser.add_argument('trace_file', type=str, default=None,
                            help='PypeIt Edges file [e.g. Edges_A_0_DET01.fits.gz]')
        parser.add_argument('--try_old', default=False, action='store_true',
                            help='Attempt to load old datamodel versions.  A crash may ensue..')
        return parser

    @staticmethod
    def main(args):

        from pathlib import Path
        from matplotlib import pyplot
        from pypeit import edgetrace
        from pypeit.core.gui import edge_inspector

        chk_version = not args.try_old

        # Set the file name to the full path
        trace_file = Path(args.trace_file).resolve()
        # Load
        edges = edgetrace.EdgeTraceSet.from_file(trace_file, chk_version=chk_version)
        # Inspector object
        pointer = edge_inspector.EdgeInspectorGUI(edges)
        # Run.  Ends when window is closed
        pyplot.show()
        # Close (restore pyplot rc defaults)
        pointer.close()

        # Write the updated edges and slits
        edges.to_file(file_path=trace_file)
        slits_file = trace_file.name.replace('Edges', 'Slits')
        edges.get_slits().to_file(file_path=trace_file.parent / slits_file)


