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

        parser = super().get_parser(description='Display Alignment image and the trace data',
                                    width=width)

        parser.add_argument('file', type=str, default = None,
                            help='PypeIt Alignment file [e.g. Alignment_A_1_DET01.fits]')
        parser.add_argument('--chname', default='Alignments', type=str,
                            help='Channel name for image in Ginga')
        parser.add_argument('--try_old', default=False, action='store_true',
                            help='Attempt to load old datamodel versions.  A crash may ensue..')
        return parser

    @staticmethod
    def main(args):
        from pypeit import alignframe

        # Load
        chk_version = not args.try_old
        alignments = alignframe.Alignments.from_file(args.file, chk_version=chk_version)
        # Show
        alignments.show()


