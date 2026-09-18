"""
Utilities for the PypeIt Dashboard GUI.

The central concern here is *loud failures*: by default an exception raised
inside a Qt slot is swallowed by the C++ event loop, so the application
silently misbehaves and the traceback is lost.  These helpers ensure such
failures are logged with full tracebacks.

The exception hook defined here is installed by the ``pypeit_dashboard``
script (see :meth:`pypeit.scripts.dashboard.PypeItDashboard.init_log`),
which overrides the general hook installed by
:class:`~pypeit.pkg.logger.PypeItLogger`.
"""

import functools
import traceback

from pypeit import log


def dashboard_excepthook(exc_type, exc_value, exc_tb):
    """
    Exception hook that logs the full traceback through PypeIt's logger.

    Intended to be assigned to ``sys.excepthook`` by the dashboard script's
    ``init_log``.  The traceback is emitted via the logger (console stream
    and any log file), rather than printed raw to stderr, so a failure
    inside the Qt event loop always leaves a record in the reduction log.

    Parameters
    ----------
    exc_type : type
        The exception class.
    exc_value : BaseException
        The exception instance.
    exc_tb : traceback
        The associated traceback object.
    """
    tb_text = ''.join(
        traceback.format_exception(exc_type, exc_value, exc_tb))
    log.error(f'Unhandled exception in the PypeIt Dashboard:\n{tb_text}')


def safe_slot(func):
    """
    Decorate a Qt slot so any exception it raises is logged with a full
    traceback and then re-raised, instead of being silently swallowed by
    the Qt event loop.

    Parameters
    ----------
    func : callable
        The slot method to wrap.

    Returns
    -------
    callable
        The wrapped slot.
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception:
            # Deliberate log-and-re-raise guard: record the full traceback,
            # then re-raise so the global excepthook (and the developer)
            # still see the failure.
            log.error(f'Exception in slot {func.__name__}:\n'
                      f'{traceback.format_exc()}')
            raise
    return wrapper
