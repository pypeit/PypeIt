"""
Script to check and display the status of a PypeIt reduction
without actually running any processing steps.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit.scripts import scriptbase


class PypeItStatus(scriptbase.ScriptBase):

    @classmethod
    def name(cls):
        """
        Return the name of the executable.
        """
        return 'pypeit_status'

    @classmethod
    def get_parser(cls, width=None):
        import argparse

        parser = super().get_parser(
            description='Check the status of a PypeIt reduction',
            width=width,
            formatter=argparse.RawDescriptionHelpFormatter,
            default_log_file=True)
        parser.add_argument('pypeit_file', type=str,
                            help='PypeIt reduction file (must have .pypeit extension)')
        return parser

    @classmethod
    def main(cls, args):

        from pathlib import Path
        from IPython import embed

        from pypeit import pypeit
        from pypeit.state import science_status
        from pypeit import log
        from pypeit import PypeItError

        # Set a default log file (avoid over-writing the standard one)
        if args.log_file == 'default':
            _pypeit_file = Path(args.pypeit_file)
            if _pypeit_file.suffix != '.pypeit':
                raise PypeItError('Input file must have a .pypeit extension!')
            args.log_file = _pypeit_file.with_suffix('.status.log')

        cls.init_log(args)

        # Instantiate PypeIt in calib_only mode to skip science frame checks
        pypeIt = pypeit.PypeIt(
            args.pypeit_file,
            reuse_calibs=True,
            calib_only=True)

        # Run calibration status check only (no processing)
        pypeIt.calib_all(status_only=True, reload_only=True)

        # Derive the science-frame state from the on-disk products (Science
        # spec2d/spec1d, with Intermediate/ as a fallback); no processing.
        science_status.derive_science_from_disk(
            pypeIt.run_state, pypeIt.par['rdx']['redux_path'],
            fitstbl=pypeIt.fitstbl)

        # Write state to JSON
        #pypeIt.run_state.write()

        # Pretty-print the state (calibrations + science)
        pypeIt.run_state.print_status()

        return 0
