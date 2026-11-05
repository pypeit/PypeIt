"""
The Dashboard status/activity area.

A status bar — distinct from the PypeIt reduction *state* — that reports what
the Dashboard itself is doing and when it is waiting on a task, with busy
indicators for long operations.  The bar is shared by all tabs and is split
into two channels so that build/monitoring status and user-inspection feedback
never overwrite each other:

* **Build** (left) — reports (re)build runs launched from the Dashboard and
  the live monitoring of an active reduction.
* **Inspection** (right) — reports feedback for user-launched viewers
  (inspecting step outputs, QA figures, or input frames).
"""

from qtpy.QtWidgets import QStatusBar, QLabel, QProgressBar


class ActivityBar(QStatusBar):
    """
    The shared Dashboard status/activity bar with a **Build** and an
    **Inspection** channel; see the module docstring for the channel
    semantics.

    Parameters
    ----------
    parent : QWidget, optional
        The parent widget.
    """

    def __init__(self, parent=None):
        super().__init__(parent=parent)
        # Build channel (left): (re)build runs + live monitoring.
        build_title = QLabel('Build:')
        build_title.setStyleSheet('font-weight: bold;')
        self.addWidget(build_title)
        self._build_msg = QLabel('Idle')
        self.addWidget(self._build_msg, stretch=1)
        self._build_busy = self._busy_bar()
        self.addWidget(self._build_busy)

        # Inspection channel (right, permanent): user viewer launches.  Added
        # left-to-right as permanent widgets, so they sit at the right edge.
        inspect_title = QLabel('Inspection:')
        inspect_title.setStyleSheet('font-weight: bold;')
        self.addPermanentWidget(inspect_title)
        self._inspect_msg = QLabel('—')
        self.addPermanentWidget(self._inspect_msg)
        self._inspect_busy = self._busy_bar()
        self.addPermanentWidget(self._inspect_busy)

    @staticmethod
    def _busy_bar():
        """
        Build an indeterminate busy indicator (hidden until shown).

        Returns
        -------
        QProgressBar
            A range-(0,0) spinner, initially hidden.
        """
        bar = QProgressBar()
        bar.setRange(0, 0)               # range (0,0) → indeterminate spinner
        bar.setMaximumWidth(120)
        bar.setVisible(False)
        return bar

    # -- Build channel ---------------------------------------------------

    def set_build(self, message, busy=False):
        """
        Show a message on the **Build** channel ((re)build / monitoring).

        Parameters
        ----------
        message : str
            The message to display.
        busy : bool, optional
            Show the build busy indicator.
        """
        self._build_msg.setText(message)
        self._build_busy.setVisible(busy)

    # -- Inspection channel ----------------------------------------------

    def set_inspection(self, message, busy=False):
        """
        Show a message on the **Inspection** channel (user viewer launches).

        Parameters
        ----------
        message : str
            The message to display.
        busy : bool, optional
            Show the inspection busy indicator.
        """
        self._inspect_msg.setText(message)
        self._inspect_busy.setVisible(busy)

    # -- shared ----------------------------------------------------------

    def idle(self):
        """
        Reset both channels to their idle state.
        """
        self.set_build('Idle', busy=False)
        self.set_inspection('—', busy=False)

    def build_message(self):
        """
        Return the current Build-channel message (used by tests).

        Returns
        -------
        str
            The Build-channel message.
        """
        return self._build_msg.text()

    def inspection_message(self):
        """
        Return the current Inspection-channel message (used by tests).

        Returns
        -------
        str
            The Inspection-channel message.
        """
        return self._inspect_msg.text()
