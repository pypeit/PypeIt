"""
The Dashboard-wide status palette (Qt-free, plain data).

This module maps a calibration step's ``(required, status, in_pipeline)``
to a color + glyph + label.  It returns **plain data** (hex strings, a
glyph, a label) so it stays Qt-free and unit-testable; the *view* converts
the hex to a ``QColor``.  Both a light and a dark variant are provided,
tuned to keep equivalent contrast; the view picks one by theme.  Every
status pairs a glyph with its color so status is never communicated by
color alone.
"""

from dataclasses import dataclass

# Status-category keys.
SUCCESS = 'success'
RUNNING = 'running'
FAIL = 'fail'
REQUIRED_UNDONE = 'required_undone'
OPTIONAL = 'optional'
NOT_USED = 'not_used'
# Per-slit-only: a slit intentionally skipped (flats SKIPFLATCALIB). Distinct
# from "fail" (a real problem) and "optional" (a not-required step).
SKIP = 'skip'

# Glyph + human label per category (theme-independent), paired with color so
# status never depends on color alone.
GLYPHS = {
    SUCCESS: ('✓', 'success'),          # check mark
    RUNNING: ('⏳', 'running'),           # hourglass
    FAIL: ('✗', 'fail'),                 # ballot X
    REQUIRED_UNDONE: ('○', 'required'),  # open circle
    OPTIONAL: ('–', 'optional'),         # en dash
    NOT_USED: ('–', 'n/a'),              # en dash
    SKIP: ('⊘', 'skipped'),              # circled slash
}

# Light-theme hex colors.
LIGHT_COLORS = {
    SUCCESS: '#2E7D32',
    RUNNING: '#EF6C00',
    FAIL: '#C62828',
    REQUIRED_UNDONE: '#FFFFFF',   # needs an outline to read on light bg
    OPTIONAL: '#9E9E9E',
    NOT_USED: '#BDBDBD',
    SKIP: '#607D8B',              # blue-grey: intentionally skipped
}

# Dark-theme hex colors: brighter fills to keep equivalent contrast on a
# dark background (first pass; tunable).  White/greys still read on dark.
DARK_COLORS = {
    SUCCESS: '#66BB6A',
    RUNNING: '#FFA726',
    FAIL: '#EF5350',
    REQUIRED_UNDONE: '#FFFFFF',
    OPTIONAL: '#9E9E9E',
    NOT_USED: '#616161',
    SKIP: '#90A4AE',             # lighter blue-grey for dark theme
}

# Statuses (from ``pypeit.state``) that count as "generated successfully".
_SUCCESS_STATUS = ('success', 'complete')

# Distinct, bright control colors — deliberately **not** status colors (so a
# button is never mistaken for a "success"/"fail" status).  The (Re)Build
# *action* is blue; "Inspect output" is teal (a different bright hue, so the
# two primary actions are easy to tell apart).  The magenta selection ring
# stays "selected"; the status palette stays "status".
ACTION_COLORS = {'light': '#1565C0', 'dark': '#42A5F5'}
INSPECT_COLORS = {'light': '#00838F', 'dark': '#26C6DA'}

# Neutral **row-selection** fill for the dashboard tables.  A soft
# blue-grey, applied explicitly so a selected row never inherits the desktop
# theme's Highlight color (red on some systems → reads as "failed").
SELECTION_COLORS = {'light': '#CFD8DC', 'dark': '#37474F'}


def selection_style(theme='light'):
    """
    Return a Qt stylesheet giving a table a **neutral** selected-row fill,
    so a selected frame reads as "selected", never as a failure.  Applied
    to the selectable dashboard tables in place of the desktop theme's
    (possibly red) Highlight color.

    Parameters
    ----------
    theme : :obj:`str`, optional
        ``'light'`` (default) or ``'dark'``.

    Returns
    -------
    :obj:`str`
        A ``QTableWidget``/``QListWidget`` ``::item:selected`` stylesheet.
    """
    color = SELECTION_COLORS['dark'] if theme == 'dark' \
        else SELECTION_COLORS['light']
    text = '#FFFFFF' if theme == 'dark' else '#000000'
    return (f'QTableWidget::item:selected, QListWidget::item:selected '
            f'{{ background-color: {color}; color: {text}; }}')


def action_color(theme='light'):
    """
    Return the (Re)Build action-control color for a theme.

    Parameters
    ----------
    theme : :obj:`str`, optional
        ``'light'`` (default) or ``'dark'``.

    Returns
    -------
    :obj:`str`
        The action hex color (a blue, distinct from any status color).
    """
    return ACTION_COLORS['dark'] if theme == 'dark' else ACTION_COLORS['light']


def inspect_color(theme='light'):
    """
    Return the "Inspect output" control color for a theme.

    Parameters
    ----------
    theme : :obj:`str`, optional
        ``'light'`` (default) or ``'dark'``.

    Returns
    -------
    :obj:`str`
        The inspect hex color (a teal, distinct from the action blue and
        any status color).
    """
    return INSPECT_COLORS['dark'] if theme == 'dark' \
        else INSPECT_COLORS['light']


@dataclass
class StepStyle:
    """
    The visual style for one calibration step's state.

    Attributes
    ----------
    category : :obj:`str`
        One of the palette category keys (e.g. ``success``, ``fail``).
    color : :obj:`str`
        Hex color string (e.g. ``#2E7D32``) for the chosen theme.
    glyph : :obj:`str`
        A single-character status glyph (paired with color so status never
        depends on color alone).
    label : :obj:`str`
        A short text label for the status.
    """
    category: str
    color: str
    glyph: str
    label: str


def classify(required, status, in_pipeline):
    """
    Map a step's ``(required, status, in_pipeline)`` to a palette category.

    Parameters
    ----------
    required : :obj:`bool`, optional
        Whether the step is required.  May be ``None`` when unknown
        (treated as not required).
    status : :obj:`str`, optional
        The step status from ``pypeit.state`` (``success``, ``complete``,
        ``running``, ``fail``, ``undone``), or a not-present sentinel
        (e.g. ``absent``/``None``).
    in_pipeline : :obj:`bool`
        Whether the step is part of the active spectrograph's
        ``default_steps()``.

    Returns
    -------
    :obj:`str`
        The palette category key.
    """
    # A step the spectrograph never runs: dimmed, regardless of status.
    if not in_pipeline:
        return NOT_USED
    if status in _SUCCESS_STATUS:
        return SUCCESS
    if status == RUNNING:
        return RUNNING
    if status == FAIL:
        return FAIL
    # Not yet generated (undone / absent / None): required vs optional.
    return REQUIRED_UNDONE if required else OPTIONAL


# Severity order (most → least) for summarizing a set of steps into one
# "worst" category, e.g. coloring a configuration-overview navigator cell.
# optional/not_used rank *below* success so they never worsen a cell that is
# otherwise successful.
_SEVERITY = [FAIL, RUNNING, REQUIRED_UNDONE, SUCCESS, OPTIONAL, NOT_USED]


def worst_category(categories):
    """
    Return the most severe palette category in ``categories``.

    Used to color a navigator cell by the worst status among its steps
    (precedence ``fail > running > required_undone > success``; ``optional``
    and ``not_used`` never worsen an otherwise-successful cell).

    Parameters
    ----------
    categories : iterable
        Palette category keys (e.g. from :func:`classify`).

    Returns
    -------
    :obj:`str`
        The most severe category present, or :data:`NOT_USED` if the input
        is empty.
    """
    present = set(categories)
    for category in _SEVERITY:
        if category in present:
            return category
    return NOT_USED


def step_style(required, status, in_pipeline, theme='light'):
    """
    Return the :class:`StepStyle` for a step's state and theme.

    Parameters
    ----------
    required : :obj:`bool`, optional
        Whether the step is required (``None`` treated as not required).
    status : :obj:`str`, optional
        The step status (see :func:`classify`).
    in_pipeline : :obj:`bool`
        Whether the step is in the spectrograph's ``default_steps()``.
    theme : :obj:`str`, optional
        ``'light'`` (default) or ``'dark'``.

    Returns
    -------
    :class:`StepStyle`
        The color + glyph + label for the step.
    """
    category = classify(required, status, in_pipeline)
    colors = DARK_COLORS if theme == 'dark' else LIGHT_COLORS
    glyph, label = GLYPHS[category]
    return StepStyle(category=category, color=colors[category],
                     glyph=glyph, label=label)


# Per-slit status → palette category (used by the Calibrations view's
# per-slit/order drill-down). Unlike step_style, a slit has only a status —
# no required/in_pipeline — and may be 'skip' (flats SKIPFLATCALIB).
_SLIT_CATEGORY = {
    'success': SUCCESS,
    'complete': SUCCESS,
    'running': RUNNING,
    'fail': FAIL,
    'skip': SKIP,
    'undone': REQUIRED_UNDONE,
}


def slit_style(status, theme='light'):
    """
    Return the :class:`StepStyle` for one per-slit/order status.

    Parameters
    ----------
    status : :obj:`str`
        The per-slit status (``success``/``complete``/``running``/
        ``fail``/``skip``/``undone``).
    theme : :obj:`str`, optional
        ``'light'`` (default) or ``'dark'``.

    Returns
    -------
    :class:`StepStyle`
        The color + glyph + label for the slit.
    """
    category = _SLIT_CATEGORY.get(status, REQUIRED_UNDONE)
    colors = DARK_COLORS if theme == 'dark' else LIGHT_COLORS
    glyph, label = GLYPHS[category]
    return StepStyle(category=category, color=colors[category],
                     glyph=glyph, label=label)
