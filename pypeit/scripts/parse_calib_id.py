"""
Parse the calib_id's using the PypeIt file

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit.scripts import scriptbase


class ParseCalibID(scriptbase.ScriptBase):

    @classmethod
    def name(cls):
        """
        Return the name of the executable.
        """
        return 'run_pypeit'

    @classmethod
    def usage(cls):
        """
        Print pypeit usage description.
        """
        import textwrap
        import pypeit
        from pypeit.spectrographs import available_spectrographs

        spclist = ', '.join(available_spectrographs)
        spcl = textwrap.wrap(spclist, width=70)
        descs = '##  '
        descs += 'The PypeIt Parse calib_id script'
        descs += '##  '
        return descs

    @classmethod
    def get_parser(cls, width=None):
        import argparse

        parser = super().get_parser(description=cls.usage(),
                                    width=width, formatter=argparse.RawDescriptionHelpFormatter)
        parser.add_argument('pypeit_file', type=str,
                            help='PypeIt reduction file (must have .pypeit extension)')

    #    parser.add_argument('-c', '--cpus', default=False, action='store_true',
    #                         help='Number of CPUs for parallel processing')

        return parser

    @staticmethod
    def main(args):

        import os

        import numpy as np

        from pypeit import pypeit
        from pypeit import msgs

        # Load options from command line
        splitnm = os.path.splitext(args.pypeit_file)
        if splitnm[1] != '.pypeit':
            msgs.error('Input file must have a .pypeit extension!')

        # Instantiate the main pipeline reduction object
        pypeIt = pypeit.PypeIt(args.pypeit_file, verbosity=1,
                               calib_only=True, show=False)

        pypeIt.calib_all(run=False)
        msgs.info('Data reduction complete')
        # QA HTML
        msgs.info('Generating QA HTML')
        pypeIt.build_qa()
        msgs.close()

        return 0


