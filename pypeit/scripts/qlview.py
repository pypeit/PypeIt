"""
Launch the PypeIt interactive quicklook viewer.

Opens Ginga with the QLView local plugin, which provides a self-contained
file browser, slit-overlay renderer, and one-click quicklook reduction
interface.  Configuration (initial raw/reduced paths, instrument,
reduction backend, etc.) is read from ``~/.pypeit/quicklook.cfg`` at startup.

By default, the viewer runs in a separate (background) process, and the script
returns immediately.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit.scripts import scriptbase


class QLViewer(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(
            description='Launch the PypeIt quicklook viewer (QLView) in Ginga.  By default, '
                        'the viewer runs in a separate process and this script returns '
                        'immediately.',
            width=width,
        )
        parser.add_argument('--no_clear', default=False, action='store_true',
                            help='If a PypeIt viewer is already open (e.g., from a previous call '
                                 'to this script, or to pypeit_show_2dspec, pypeit_view_fits, '
                                 'or pypeit_show_1dspec), open QLView in that viewer and return '
                                 'immediately.')
        parser.add_argument('--foreground', default=False, action='store_true',
                            help='Run the viewer in this process; the script returns when the '
                                 'viewer window is closed.')
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
        if args.no_clear:
            if viewer.send_command('qlview'):
                log.info('Opened QLView in the open PypeIt viewer.')
                return
            log.info('No open PypeIt viewer found; opening a new one.')

        if not args.foreground:
            # Run the viewer in a separate process, which keeps running after
            # this script returns.
            subprocess.Popen([sys.executable, '-c',
                              'from pypeit.scripts.qlview import QLViewer; '
                              'QLViewer.entry_point()', '--foreground'],
                             start_new_session=True)
            log.info('Started the PypeIt quicklook viewer in a separate process.')
            return

        from pypeit.display.display import launch_qlview

        v = viewer.get_viewer()
        # Allow other PypeIt scripts (with --no_clear) to use this viewer
        v.listen()
        launch_qlview(v)

        # Wait for the user to close the window
        v.show(message='close the window to exit')
