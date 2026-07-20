"""
Utilities for the PypeIt Dashboard GUI.

The central concern here is *loud failures*: by default an exception raised
inside a Qt slot is swallowed by the C++ event loop, so the application
silently misbehaves and the traceback is lost.  These helpers ensure such
failures are logged with full tracebacks (per the Debugging plan in
``pypeit_dashboard_coding.md``).
"""

import sys
import functools
import traceback

from pypeit import log

# Module-level flag so the excepthook is installed at most once.
_excepthook_installed = False


def _dashboard_excepthook(exc_type, exc_value, exc_tb):
    """
    Global exception hook that logs the full traceback through PypeIt's
    logger before deferring to the default hook.

    Args:
        exc_type (type): The exception class.
        exc_value (BaseException): The exception instance.
        exc_tb (traceback): The associated traceback object.

    Returns:
        None.
    """
    # Format and log the traceback so a failure leaves a record even when
    # the Qt event loop would otherwise discard it.
    tb_text = ''.join(
        traceback.format_exception(exc_type, exc_value, exc_tb))
    log.error(f'Unhandled exception in the PypeIt Dashboard:\n{tb_text}')
    # Still call the default hook so the message reaches stderr as usual.
    sys.__excepthook__(exc_type, exc_value, exc_tb)


def install_excepthook():
    """
    Install the dashboard's global ``sys.excepthook`` (idempotent).

    Returns:
        None.
    """
    global _excepthook_installed
    if _excepthook_installed:
        return
    sys.excepthook = _dashboard_excepthook
    _excepthook_installed = True


def safe_slot(func):
    """
    Decorate a Qt slot so any exception it raises is logged with a full
    traceback and then re-raised, instead of being silently swallowed by
    the Qt event loop.

    Args:
        func (callable):
            The slot method to wrap.

    Returns:
        callable: The wrapped slot.
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception:
            # Log the full traceback, then re-raise so the global
            # excepthook (and the developer) still see the failure.
            log.error(f'Exception in slot {func.__name__}:\n'
                      f'{traceback.format_exc()}')
            raise
    return wrapper
