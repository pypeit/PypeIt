"""
The shared header banner for the PypeIt Dashboard (design requirement R6).

The banner shows the reduction's context — the ``.pypeit`` file name, the
spectrograph (``PYP_SPEC``), the setup/configuration ID, the pipeline
(MultiSlit / Echelle / IFU), and the reduction directory — with the PypeIt
logo in the top-right corner.  It is shared by both the Status and
Calibrations views.
"""

from qtpy.QtWidgets import (QWidget, QHBoxLayout, QVBoxLayout, QGridLayout,
                            QLabel, QFrame, QSizePolicy)
from qtpy.QtCore import Qt
from qtpy.QtGui import QPixmap

from pypeit import __version__, dataPaths

# The logo asset shipped with the package (pypeit/data/images).
_LOGO_PATH = dataPaths.images.path / 'pypeit_logo.png'

# The header fields, in display order: (attribute on HeaderInfo, label text).
_FIELDS = (('pypeit_file', 'PypeIt file'),
           ('spectrograph', 'Spectrograph'),
           ('setup', 'Setup'),
           ('path', 'Pipeline'),
           ('redux_dir', 'Reduction dir'))

# Fields whose values can be long (a filename or a full path).  These use an
# eliding label so they never force a large minimum width on the window
# (otherwise the long monospace path prevents the window from shrinking
# horizontally).
_LONG_FIELDS = ('pypeit_file', 'redux_dir')


class ElidingLabel(QLabel):
    """
    A :obj:`QLabel` that elides its (long) text to fit the available width,
    and reports a small minimum width so it never blocks the window from
    being resized narrower.

    Args:
        mode (``Qt.TextElideMode``, optional):
            How to elide; defaults to eliding the middle (keeps the start
            and end of a path visible).
        parent (:obj:`QWidget`, optional):
            The parent widget.
    """

    def __init__(self, mode=Qt.ElideMiddle, parent=None):
        super().__init__(parent=parent)
        self._full_text = ''
        self._elide_mode = mode
        # Ignored width: the label provides no width hint to the layout, so
        # it can shrink freely; it elides its text to whatever width it gets.
        self.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Preferred)
        self.setMinimumWidth(0)

    def setText(self, text):
        """
        Store the full text and display an elided version.

        Args:
            text (:obj:`str`): The full (un-elided) text.

        Returns:
            None.
        """
        self._full_text = '' if text is None else str(text)
        self._update_elided()

    def fullText(self):
        """
        Return the full, un-elided text.

        Returns:
            str: The stored full text.
        """
        return self._full_text

    def resizeEvent(self, event):
        """
        Re-elide the text whenever the label is resized.

        Args:
            event (``QResizeEvent``): The resize event.

        Returns:
            None.
        """
        super().resizeEvent(event)
        self._update_elided()

    def _update_elided(self):
        """
        Set the displayed text to the full text elided to the current width.

        Returns:
            None.
        """
        metrics = self.fontMetrics()
        elided = metrics.elidedText(self._full_text, self._elide_mode,
                                    max(self.width(), 0))
        super().setText(elided)


class HeaderBanner(QFrame):
    """
    The shared context banner + PypeIt logo (R6).

    Args:
        header_info (:class:`pypeit.dashboard.model.HeaderInfo`):
            The reduction metadata to display.
        parent (:obj:`QWidget`, optional):
            The parent widget.
    """

    def __init__(self, header_info, parent=None):
        super().__init__(parent=parent)
        self.setFrameShape(QFrame.StyledPanel)
        self.setObjectName('dashboardHeader')

        # value_labels maps each field attribute -> its value QLabel, so the
        # banner can be refreshed in place and tests can assert the fields.
        self.value_labels = {}

        layout = QHBoxLayout(self)

        # Left: the field grid (label : value), one row per R6 field.
        grid = QGridLayout()
        for row, (attr, text) in enumerate(_FIELDS):
            name_label = QLabel(f'{text}:')
            name_label.setStyleSheet('font-weight: bold;')
            # Long values (filename, path) use an eliding label so they do
            # not force a large minimum window width; short values use a
            # plain label.
            if attr in _LONG_FIELDS:
                value_label = ElidingLabel()
                # Monospace for file/dir values, matching the design sketches.
                value_label.setStyleSheet('font-family: monospace;')
            else:
                value_label = QLabel('')
            value_label.setTextInteractionFlags(
                Qt.TextSelectableByMouse)
            grid.addWidget(name_label, row, 0, Qt.AlignRight | Qt.AlignTop)
            grid.addWidget(value_label, row, 1, Qt.AlignLeft | Qt.AlignTop)
            self.value_labels[attr] = value_label
        grid.setColumnStretch(1, 1)
        layout.addLayout(grid, stretch=1)

        # Right: the PypeIt logo (top-right corner, R6) with the PypeIt
        # version number underneath it.
        right = QVBoxLayout()
        self._logo_label = QLabel()
        self._logo_label.setAlignment(Qt.AlignTop | Qt.AlignRight)
        self._logo_label.setSizePolicy(QSizePolicy.Fixed,
                                       QSizePolicy.Preferred)
        self._set_logo()
        right.addWidget(self._logo_label, alignment=Qt.AlignRight)
        version_label = QLabel(f'PypeIt v{__version__}')
        version_label.setAlignment(Qt.AlignRight | Qt.AlignTop)
        version_label.setStyleSheet('color: grey; font-size: 10px;')
        version_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        right.addWidget(version_label, alignment=Qt.AlignRight)
        right.addStretch(1)
        layout.addLayout(right)

        # Populate the field values from the supplied metadata.
        self.set_header_info(header_info)

    def _set_logo(self):
        """
        Load and scale the PypeIt logo into the corner label, if present.

        Returns:
            None.
        """
        if not _LOGO_PATH.is_file():
            # Non-fatal: render without a logo rather than failing startup.
            return
        pixmap = QPixmap(str(_LOGO_PATH))
        if pixmap.isNull():
            return
        # Constrain the logo height so it stays clear of the banner text.
        self._logo_label.setPixmap(
            pixmap.scaledToHeight(72, Qt.SmoothTransformation))

    def set_header_info(self, header_info):
        """
        Update the displayed field values from a :class:`HeaderInfo`.

        Args:
            header_info (:class:`pypeit.dashboard.model.HeaderInfo`):
                The reduction metadata to display.

        Returns:
            None.
        """
        for attr, _ in _FIELDS:
            value = getattr(header_info, attr, None)
            self.value_labels[attr].setText('' if value is None else
                                             str(value))
