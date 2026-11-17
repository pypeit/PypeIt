"""
Launch the PypeIt interactive quicklook viewer.

Opens Ginga with the QLView local plugin, which provides a self-contained
file browser, slit-overlay renderer, and one-click quicklook reduction
interface.  Configuration (initial raw/reduced paths, instrument,
reduction backend, etc.) is read from ``~/.pypeit/quicklook.cfg`` at startup.

If a PypeIt viewer is already open, the plugin is opened in that viewer;
otherwise, a new viewer is opened in a separate (background) process.  In
either case, the script returns immediately.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit.scripts import scriptbase


class QLViewer(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(
            description='Launch the PypeIt quicklook viewer (QLView) in Ginga.  If a PypeIt '
                        'viewer is already open (e.g., from a previous call to this script, or '
                        'to pypeit_show_2dspec, pypeit_view_fits, or pypeit_show_1dspec), the '
                        'plugin is opened in that viewer; otherwise, a new viewer is opened.  '
                        'In either case, the script returns immediately.',
            width=width,
        )
        return parser

    @classmethod
    def main(cls, args):
        """Launch Ginga and open the QLView plugin on the QuickLook channel."""
        import subprocess
        import sys

        from pypeit import log
        from pypeit.display import viewer

        cls.init_log(args)

        # Open QLView in an already open viewer?
        if viewer.send_command('qlview'):
            log.info('Opened QLView in the open PypeIt viewer.')
            return

        # Otherwise, open a new viewer in a separate process, which keeps
        # running after this script returns.
        subprocess.Popen([sys.executable, '-c',
                          'from pypeit.scripts.qlview import QLViewer; QLViewer.run_viewer()'],
                         start_new_session=True)
        log.info('Started the PypeIt quicklook viewer in a separate process.')

    @staticmethod
    def run_viewer():
        """
        Open the viewer with QLView in this process and wait for the user to
        close its window.

        This is run by :meth:`main` in a separate process.
        """
        from pypeit.display import viewer
        from pypeit.display.display import launch_qlview

        v = viewer.get_viewer()
        # Allow later calls to this script (and other PypeIt scripts) to use
        # this viewer
        v.listen()
        launch_qlview(v)
        v.show(message='close the window to exit')
