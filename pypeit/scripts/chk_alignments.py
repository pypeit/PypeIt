"""
This script displays the Trace image and the traces in an RC Ginga window.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit.scripts import scriptbase


class ChkAlignments(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        from pypeit.spectrographs import available_spectrographs

        parser = super().get_parser(description='Display MasterAlignment image and the trace data',
                                    width=width)

        parser.add_argument('master_file', type=str, default = None,
                            help='PypeIt Master Alignment file [e.g. MasterAlignment_A_1_01.fits]')
        parser.add_argument('--chname', default='Alignments', type=str,
                            help='Channel name for image in Ginga')
        return parser

    @staticmethod
    def main(pargs):
        from pypeit import alignframe

        # Load
        alignments = alignframe.Alignments.from_file(pargs.master_file)
        # Show
        alignments.show()


