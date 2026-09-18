"""
Shared helper functions for the PypeIt Dashboard view widgets.

Small Qt utilities used by several of the view modules (the Status,
Calibrations, and Science tabs), collected here so that each is defined only
once.
"""

from qtpy.QtGui import QColor, QPalette


def text_on(hexcolor):
    """
    Pick a readable text color (black or white) for a background color.

    Parameters
    ----------
    hexcolor : str
        The background color as a hex string, e.g. ``'#2E7D32'``.

    Returns
    -------
    str
        ``'#000000'`` (black) for light backgrounds or ``'#FFFFFF'`` (white)
        for dark backgrounds.
    """
    color = QColor(hexcolor)
    luminance = (0.299 * color.redF() + 0.587 * color.greenF()
                 + 0.114 * color.blueF())
    return '#000000' if luminance > 0.6 else '#FFFFFF'


def detect_theme(widget):
    """
    Detect whether a widget is rendered with a light or a dark Qt theme.

    Used by the views to select between the light and dark color sets defined
    by :mod:`pypeit.dashboard.palette`.

    Parameters
    ----------
    widget : QWidget
        The widget whose palette is inspected.

    Returns
    -------
    str
        ``'dark'`` or ``'light'``.
    """
    color = widget.palette().color(QPalette.Window)
    luminance = (0.299 * color.redF() + 0.587 * color.greenF()
                 + 0.114 * color.blueF())
    return 'dark' if luminance < 0.5 else 'light'


def clear_layout(layout):
    """
    Remove and delete every item from a Qt layout.

    Used when a view re-renders itself in place.  Each widget is reparented
    out of the visible widget tree immediately because ``deleteLater`` is
    deferred to the event loop, which would otherwise leave the old widgets
    overlapping a freshly rebuilt panel.

    Parameters
    ----------
    layout : QLayout
        The layout to empty.
    """
    while layout.count() > 0:
        item = layout.takeAt(0)
        widget = item.widget()
        if widget is not None:
            widget.setParent(None)
            widget.deleteLater()
        else:
            child = item.layout()
            if child is not None:
                clear_layout(child)


def fmt_value(value):
    """
    Format a metric value for display in a table cell or label.

    Parameters
    ----------
    value : object
        The value to format, or ``None`` if the value is missing.

    Returns
    -------
    str
        A short display string: ``'yes'``/``'no'`` for booleans, a
        4-significant-figure string for floats, ``str(value)`` otherwise, and
        an em-dash for a missing (``None``) value.
    """
    if value is None:
        return '—'
    if isinstance(value, bool):
        return 'yes' if value else 'no'
    if isinstance(value, float):
        return f'{value:.4g}'
    return str(value)
