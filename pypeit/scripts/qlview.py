"""
Launch the PypeIt interactive quicklook viewer.

Opens Ginga with the QLView local plugin, which provides a self-contained
file browser, slit-overlay renderer, and one-click quicklook reduction
interface.  Configuration (initial raw/reduced paths, instrument,
reduction backend, etc.) is read from ``~/.quicklook.cfg`` at startup.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit.scripts import scriptbase


class QLViewer(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(
            description='Launch the PypeIt quicklook viewer (QLView) in Ginga.',
            width=width,
        )
        return parser

    @classmethod
    def main(cls, args):
        """Launch Ginga and open the QLView plugin on the QuickLook channel."""
        from pypeit.display.display import launch_qlview

        cls.init_log(args)
        launch_qlview()
