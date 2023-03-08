"""
Parse the calib_id's using the PypeIt file

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit.scripts import scriptbase


class ParseCalibID(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        import argparse

        parser = super().get_parser(description='Parse the pypeit file for the calib IDs',
                                    width=width, formatter=argparse.RawDescriptionHelpFormatter)
        parser.add_argument('pypeit_file', type=str,
                            help='PypeIt reduction file (must have .pypeit extension)')
        return parser

    @staticmethod
    def main(args):

        import os

        from pypeit import pypeit
        from pypeit import msgs

        # Load options from command line
        splitnm = os.path.splitext(args.pypeit_file)
        if splitnm[1] != '.pypeit':
            msgs.error('Input file must have a .pypeit extension!')

        # Instantiate the main pipeline reduction object
        pypeIt = pypeit.PypeIt(args.pypeit_file, verbosity=1,
                               calib_only=True, show=False)

        # Grab the info without running
        calib_dict = pypeIt.calib_all(run=False)
        msgs.info('Data reduction complete')
        msgs.close()

        return calib_dict


