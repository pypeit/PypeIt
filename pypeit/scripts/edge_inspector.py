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
        return parser

    @staticmethod
    def main(args):

        from matplotlib import pyplot
        from pypeit import edgetrace
        from pypeit.core.gui import edge_inspector

        # Load
        edges = edgetrace.EdgeTraceSet.from_file(args.trace_file)
        # Inspector object
        pointer = edge_inspector.EdgeInspectorGUI(edges)
        # Inspect
        pyplot.show()
        # Close (restore pyplot rc defaults)
        pointer.close()

        # Write the updated edges and slits

        embed()
        exit()
        


