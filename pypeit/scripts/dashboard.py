"""
Command-line launcher for the PypeIt Dashboard.

The dashboard is launched like every other PypeIt script (e.g.
``run_pypeit``): via this :class:`~pypeit.scripts.scriptbase.ScriptBase`
subclass, registered as the ``pypeit_dashboard`` console entry point.

.. include:: ../include/links.rst
"""

import sys

from pypeit.scripts import scriptbase
# Qt-free helper module (safe to import for, e.g., the auto-generated docs).
from pypeit.dashboard import util


class PypeItDashboard(scriptbase.ScriptBase):
    """
    The ``pypeit_dashboard`` console script.
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
            description='Launch the PypeIt Dashboard: a GUI to monitor and '
                        'inspect a PypeIt reduction.',
            width=width, default_log_file=True)
        # The .pypeit file is a required positional argument, so multiple
        # .pypeit files in the folder are unambiguous.
        parser.add_argument('pypeit_file', type=str,
                            help='PypeIt reduction file (.pypeit) to open.')
        parser.add_argument('--redux_path', type=str, default=None,
                            help='Reduction directory.  Defaults to the '
                                 'directory containing the .pypeit file.')
        return parser

    @classmethod
    def init_log(cls, args):
        """
        Initialize the logger and install the dashboard's exception hook.

        Extends :meth:`~pypeit.scripts.scriptbase.ScriptBase.init_log`:
        after the standard logger setup (which installs the general
        :class:`~pypeit.pkg.logger.PypeItLogger` exception hook), replace
        the hook with :func:`~pypeit.dashboard.util.dashboard_excepthook`,
        which logs the *full* traceback through the logger.  Exceptions
        raised inside Qt slots are otherwise swallowed by the C++ event
        loop, so the dashboard needs every unhandled exception recorded
        loudly, with its complete traceback, in the console and log file.

        Parameters
        ----------
        args : `argparse.Namespace`_
            The parsed command-line arguments.
        """
        super().init_log(args)
        sys.excepthook = util.dashboard_excepthook

    @classmethod
    def main(cls, args):
        """
        Launch the dashboard.

        Parameters
        ----------
        args : `argparse.Namespace`_
            The parsed command-line arguments.

        Returns
        -------
        :obj:`int`
            The Qt application exit code.
        """
        # Import here so that simply importing the script module (e.g. for
        # the auto-generated docs) does not require Qt.
        from pypeit.dashboard.app import launch

        cls.init_log(args)
        return launch(args)
