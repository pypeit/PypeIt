"""
The single-run lock for the PypeIt Dashboard.

The lock is scoped to **one reduction** — the reduction directory the
dashboard was opened on — not to the file system or the machine as a whole:
at most one PypeIt run (``run_pypeit`` or ``pypeit_run_to_calibstep``)
should be active *for that reduction* at a time, and the Dashboard refuses
to launch a second.  Reductions in other directories (or other dashboards
watching them) are unaffected.

The lock engages in two situations:

* a run **the Dashboard launched** is still in progress (driven by the
  :class:`~pypeit.dashboard.launcher.Launcher` via
  :meth:`RunLock.set_dashboard_running`); and
* a run started **outside** the Dashboard is detected by watching the
  reduction's log file, ``<redux_dir>/<pypeit_root>.log`` (see
  :attr:`~pypeit.dashboard.model.DashboardModel.log_path`) — PypeIt writes
  to it continuously while running, so a log whose modification time is
  recent means a run is in progress.

The mtime test (:meth:`RunLock._is_recent`) is a pure, unit-testable
function; the polling :class:`~qtpy.QtCore.QTimer` only calls it.  The lock
emits :attr:`RunLock.lockChanged` when it transitions, which the views
connect to so they can enable/disable their launch controls.
"""

import time
from pathlib import Path

from qtpy.QtCore import QObject, QTimer, Signal


class RunLock(QObject):
    """
    Single-run lock for one reduction directory, with log-mtime detection
    of externally launched runs and a live-monitoring change signal.

    The lock is per reduction: it watches the one
    ``<redux_dir>/<pypeit_root>.log`` and ``<redux_dir>/<pypeit_root>_state.json``
    of the reduction the dashboard was opened on, so it neither sees nor
    blocks runs in other directories.

    The one polling timer serves both roles: it detects whether a run is
    **active** (the ``.log`` mtime is recent) and, while active, whether the
    reduction **state file** has changed (its mtime advanced) — emitting
    :attr:`stateChanged` so the views can refresh live.

    Parameters
    ----------
    log_path : :obj:`str`, :obj:`pathlib.Path`, optional
        The reduction's ``<pypeit_root>.log`` file to watch.  If ``None``,
        only Dashboard-launched runs lock (no external detection).
    state_path : :obj:`str`, :obj:`pathlib.Path`, optional
        The reduction's ``<pypeit_root>_state.json`` file to watch for
        live updates.  If ``None``, :attr:`stateChanged` never fires.
    parent : :obj:`QObject`, optional
        The Qt parent.

    Attributes
    ----------
    lockChanged : Signal
        Emitted with the new locked state (:obj:`bool`) on each transition.
    stateChanged : Signal
        Emitted (no args) when the state file's mtime advances **while a
        run is active** (live monitoring).
    """

    #: Emitted with the new locked state (bool) whenever it changes.
    lockChanged = Signal(bool)

    #: Emitted when the state file changes while a run is active.
    stateChanged = Signal()

    #: A ``.log`` modified within this many seconds means a run is active
    #: (PypeIt writes to it continuously; a quiet log means it finished or
    #: stalled).
    ACTIVE_WINDOW_S = 10.0

    #: Polling cadence (ms) for the external-run + state-change check (one
    #: timer for both).
    POLL_MS = 2500

    def __init__(self, log_path=None, state_path=None, parent=None):
        super().__init__(parent=parent)
        self._log_path = Path(log_path) if log_path is not None else None
        self._state_path = Path(state_path) if state_path is not None else None
        # Locked-by-our-own-run and locked-by-an-external-run, tracked
        # separately so either can hold the lock.
        self._dashboard_running = False
        self._external_active = False
        # Last-seen state-file mtime, for live-monitoring change detection.
        self._last_state_mtime = self._state_mtime()
        self._timer = QTimer(self)
        self._timer.setInterval(self.POLL_MS)
        # Re-check the .log/state mtimes on every polling tick.
        self._timer.timeout.connect(self.poll)

    # -- lifecycle -------------------------------------------------------

    def start(self):
        """
        Start polling the ``.log`` for external runs (and check once now).
        """
        self._timer.start()
        self.poll()

    def stop(self):
        """
        Stop polling.
        """
        self._timer.stop()

    # -- state -----------------------------------------------------------

    def is_locked(self):
        """
        Whether a run is active (Dashboard-launched or external).

        Returns
        -------
        :obj:`bool`
            True if the Dashboard should refuse to launch a run.
        """
        return self._dashboard_running or self._external_active

    def set_dashboard_running(self, running):
        """
        Record whether a run the Dashboard launched is in progress, emitting
        :attr:`lockChanged` if the locked state changes.

        Parameters
        ----------
        running : :obj:`bool`
            True while our launched run is active.
        """
        running = bool(running)
        if running == self._dashboard_running:
            return
        was = self.is_locked()
        self._dashboard_running = running
        self._emit_if_changed(was)

    @staticmethod
    def _is_recent(mtime, now, window):
        """
        Whether a file modified at ``mtime`` counts as "being written now".

        Pure and unit-testable: ``True`` when the log was touched within
        ``window`` seconds of ``now``.

        Parameters
        ----------
        mtime : :obj:`float`
            The file modification time (epoch seconds).
        now : :obj:`float`
            The current time (epoch seconds).
        window : :obj:`float`
            The "recent" window, in seconds.

        Returns
        -------
        :obj:`bool`
            True if ``now - mtime < window``.
        """
        return (now - mtime) < window

    def _state_mtime(self):
        """
        The state file's modification time, or ``None`` if absent/unreadable.

        Returns
        -------
        :obj:`float` or None
            The ``*_state.json`` mtime (epoch seconds).
        """
        if self._state_path is None or not self._state_path.is_file():
            return None
        try:
            return self._state_path.stat().st_mtime
        except OSError:
            return None

    def poll(self):
        """
        Re-check the ``.log`` mtime for an external run (emitting
        :attr:`lockChanged` on a transition) and the ``*_state.json`` mtime
        for a live update (emitting :attr:`stateChanged` while active).
        """
        active = False
        if self._log_path is not None and self._log_path.is_file():
            try:
                mtime = self._log_path.stat().st_mtime
            except OSError:
                active = False
            else:
                active = self._is_recent(mtime, time.time(),
                                         self.ACTIVE_WINDOW_S)
        was = self.is_locked()
        self._external_active = active
        self._emit_if_changed(was)
        # Live monitoring: while a run is active, signal a state-file change
        # so the views refresh.  Track the mtime every tick so a change that
        # happened while idle does not fire a stale update.
        state_mtime = self._state_mtime()
        if state_mtime is not None and self._last_state_mtime is not None \
                and state_mtime > self._last_state_mtime and self.is_locked():
            self.stateChanged.emit()
        self._last_state_mtime = state_mtime

    def _emit_if_changed(self, was_locked):
        """
        Emit :attr:`lockChanged` if the locked state differs from
        ``was_locked``.

        Parameters
        ----------
        was_locked : :obj:`bool`
            The locked state before the change.
        """
        now_locked = self.is_locked()
        if now_locked != was_locked:
            self.lockChanged.emit(now_locked)
