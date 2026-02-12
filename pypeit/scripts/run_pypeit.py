"""
Main execution script for ``PypeIt`` reduction pipelines.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit.scripts import scriptbase


class RunPypeIt(scriptbase.ScriptBase):

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
        from pypeit import __version__

        descr = 'PypeIt: The Python Spectroscopic Data Reduction Pipeline\n'
        descr += f'Version {__version__}\n\n'
        import textwrap
        from pypeit.spectrographs import available_spectrographs
        spclist = ', '.join(available_spectrographs)
        spcl = textwrap.wrap(spclist, width=70)
        descr += 'Available spectrographs include:\n'
        for ispcl in spcl:
            descr += f'    {ispcl}\n'
        return descr

    @classmethod
    def get_parser(cls, width=None):
        import argparse

        parser = super().get_parser(description=cls.usage(),
                                    width=width, formatter=argparse.RawDescriptionHelpFormatter,
                                    default_log_file=True)
        parser.add_argument('pypeit_file', type=str,
                            help='PypeIt reduction file (must have .pypeit extension)')

        parser.add_argument('-r', '--redux_path', default=None,
                            help='Path to directory for the reduction.  Only advised for testing')
        parser.add_argument('-m', '--do_not_reuse_calibs', dest='reuse_calibs', default=True,
                            action='store_false',
                            help='Do not load previously generated calibrations, even ones made '
                                 'during the run.')
        parser.add_argument('-s', '--show', default=False, action='store_true',
                            help='Show reduction steps via plots (which will block further '
                                 'execution until clicked on) and outputs to ginga. Requires '
                                 'remote control ginga session via '
                                 '"ginga --modules=RC,SlitWavelength &"')

        parser.add_argument('-o', '--overwrite', default=False, action='store_true',
                            help='Overwrite any existing files/directories')
        parser.add_argument('-c', '--calib_only', default=False, action='store_true',
                            help='Only run on calibrations')

        return parser

    @classmethod
    def main(cls, args):

        from pathlib import Path
        from IPython import embed

        from pypeit import pypeit
        from pypeit import log
        from pypeit import PypeItError

        # Set a default log file based on the name of the pypeit file, not the
        # name of the script
        if args.log_file == 'default':
            _pypeit_file = Path(args.pypeit_file)
            if _pypeit_file.suffix != '.pypeit':
                raise PypeItError('Input file must have a .pypeit extension!')
            args.log_file = _pypeit_file.with_suffix('.log')

        cls.init_log(args)

        # Instantiate the main pipeline reduction object
        pypeIt = pypeit.PypeIt(
            args.pypeit_file, reuse_calibs=args.reuse_calibs, overwrite=args.overwrite,
            redux_path=args.redux_path, calib_only=args.calib_only, show=args.show
        )

        if args.calib_only:
            pypeIt.calib_all()
        else:
            pypeIt.reduce_all()
        log.info('Data reduction complete')

        # QA HTML
        log.info('Generating QA HTML')
        pypeIt.build_qa()

        return 0


