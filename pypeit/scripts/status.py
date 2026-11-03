"""
Script to check and display the status of a PypeIt reduction
without actually running any processing steps.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit.scripts import scriptbase


class PypeItStatus(scriptbase.ScriptBase):
    """
    The ``pypeit_status`` console script: report a reduction's calibration
    and science status without running any processing.
    """

    @classmethod
    def get_parser(cls, width=None):
        """
        Construct the command-line argument parser.

        Parameters
        ----------
        width : :obj:`int`, optional
            Restrict the formatted help width to this many characters.

        Returns
        -------
        `argparse.ArgumentParser`_
            The command-line parser.
        """
        parser = super().get_parser(
            description='Check the status of a PypeIt reduction',
            width=width)
        parser.add_argument('pypeit_file', type=str,
                            help='PypeIt reduction file (must have .pypeit extension)')
        return parser

    @classmethod
    def main(cls, args):
        """
        Check and print a PypeIt reduction's status without processing.

        Loads the reduction in ``calib_only`` mode, derives the calibration
        status from the on-disk ``Calibrations/`` outputs (a status-only read,
        no processing), derives the science-frame status from the on-disk
        ``Science/`` products (see
        :func:`~pypeit.state.science_status.derive_science_from_disk`), and
        prints both to the screen
        (:meth:`~pypeit.state.run_state.RunPypeItState.print_status`).  The
        status is *always* re-derived from the current on-disk products; any
        existing state file is not simply trusted.

        Parameters
        ----------
        args : `argparse.Namespace`_
            The parsed command-line arguments.

        Returns
        -------
        :obj:`int`
            ``0`` on success.
        """
        from pathlib import Path

        from pypeit import pypeit
        from pypeit.state import science_status
        from pypeit import PypeItError

        if Path(args.pypeit_file).suffix != '.pypeit':
            raise PypeItError('Input file must have a .pypeit extension!')

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

        # Pretty-print the state (calibrations + science)
        pypeIt.run_state.print_status()

        return 0
