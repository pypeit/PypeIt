"""
The main window for the PypeIt Dashboard.

The window mirrors the ``setup_gui`` ``QMainWindow``-based layout: a shared
header banner at the top and a tab bar with the **Status**, **Calibrations**,
and **Science** tabs below.  It is driven by a
:class:`~pypeit.dashboard.model.DashboardModel`: the header is populated from
``model.header_info`` and the tabs render the model.
"""

from qtpy.QtWidgets import (QMainWindow, QWidget, QVBoxLayout, QTabWidget)

from pypeit.dashboard import model
from pypeit.dashboard.model import HeaderInfo
from pypeit.dashboard.launcher import Launcher
from pypeit.dashboard.runlock import RunLock
from pypeit.dashboard.view.header import HeaderBanner
from pypeit.dashboard.view.status_view import StatusView
from pypeit.dashboard.view.calibrations_view import CalibrationsView
from pypeit.dashboard.view.science_view import ScienceView
from pypeit.dashboard.view.activity import ActivityBar

# Tab labels, in display order.  Kept as a module-level constant so tests can
# assert the exact tab bar (Status | Calibrations | Science).
TAB_LABELS = ('Status', 'Calibrations', 'Science')

# The design's default window size (matches setup_gui and the layout
# sketches).
DEFAULT_SIZE = (1650, 900)


class DashboardMainWindow(QMainWindow):
    """
    The PypeIt Dashboard main window.

    Parameters
    ----------
    model : :class:`~pypeit.dashboard.model.DashboardModel`
        The reduction-state model; supplies the header metadata and drives
        the tab views.
    parent : QWidget, optional
        The parent widget (normally ``None`` for a top-level window).
    """

    def __init__(self, model, parent=None):
        super().__init__(parent=parent)

        self.model = model
        self.header_info = model.header_info \
            if model.header_info is not None \
            else self._fallback_header_info(model)

        # Central widget: header banner stacked above the tab bar.
        central = QWidget(self)
        layout = QVBoxLayout(central)

        self.header = HeaderBanner(self.header_info, parent=central)
        layout.addWidget(self.header)

        # Shared Dashboard status/activity area, the single-run lock, and the
        # subprocess launcher that reports to the activity bar and holds the
        # lock while a run is active.  The lock also drives live monitoring:
        # it watches the .log (active-run detection) and the *_state.json
        # (live-update signal) on one timer.
        self.activity = ActivityBar()
        self.setStatusBar(self.activity)
        self.run_lock = RunLock(log_path=model.log_path,
                                state_path=model.state_path)
        self.launcher = Launcher(self.activity, run_lock=self.run_lock)

        # The tab bar (Status | Calibrations | Science): the state-first
        # Status overview, the per-step Calibrations drill-down with its
        # (Re)Build control, and the per-frame Science view, all driven by
        # the model and sharing the activity bar.
        self.tab_widget = QTabWidget(parent=central)
        self.status_view = StatusView(model, activity=self.activity,
                                      parent=self.tab_widget)
        self.calibrations_view = CalibrationsView(
            model, launcher=self.launcher, run_lock=self.run_lock,
            on_run_finished=self._on_run_finished, parent=self.tab_widget)
        self.science_view = ScienceView(
            model, launcher=self.launcher, run_lock=self.run_lock,
            on_run_finished=self._on_run_finished, parent=self.tab_widget)
        self.tab_widget.addTab(self.status_view, TAB_LABELS[0])
        self.tab_widget.addTab(self.calibrations_view, TAB_LABELS[1])
        self.tab_widget.addTab(self.science_view, TAB_LABELS[2])
        layout.addWidget(self.tab_widget, stretch=1)

        # Clicking a Status-view science-navigator cell switches to the
        # Science tab and selects that frame.  The StatusView instance
        # persists across set_model, so this connection survives refreshes.
        self.status_view.scienceFrameActivated.connect(
            self._on_science_frame_activated)

        self.setCentralWidget(central)

        # Whether the most recently completed (Re)Build run exited with a
        # failure (non-zero code).  Sticky until the next run starts, so the
        # Build channel keeps reporting the failure instead of being reset to
        # "Monitoring…"/"Idle" by the still-active live monitor (the .log
        # mtime can stay "recent" for several seconds after a crash).
        self._last_run_failed = False

        # Engaging/releasing the run lock gates the (Re)Build controls and
        # updates the Build-channel status.
        self.run_lock.lockChanged.connect(self._on_lock_changed)
        # A change to the state file during a run triggers a live refresh of
        # the views.
        self.run_lock.stateChanged.connect(self._on_state_changed)
        # Start the one polling timer behind both signals.
        self.run_lock.start()

        self.setWindowTitle('PypeIt Dashboard')
        self.resize(*DEFAULT_SIZE)
        self.activity.idle()

    @staticmethod
    def _fallback_header_info(model):
        """
        Build a minimal :class:`~pypeit.dashboard.model.HeaderInfo` when the
        model has none (e.g. the ``.pypeit`` file was not found), so the
        banner still renders while the Status view shows the corresponding
        edge-state message.

        Parameters
        ----------
        model : :class:`~pypeit.dashboard.model.DashboardModel`
            The model whose header metadata is missing.

        Returns
        -------
        :class:`~pypeit.dashboard.model.HeaderInfo`
            A placeholder header.
        """
        return HeaderInfo(pypeit_file=model.pypeit_file.name,
                          spectrograph='—', setup='—', path='—',
                          redux_dir=str(model.redux_dir))

    def _refresh_from_state(self):
        """
        Re-acquire the model **from the state file** (the state is loaded,
        never re-derived, while a run may be active) and re-render all views,
        preserving the user's scope and selected step.  Shared by the
        run-completion refresh and the live monitor.
        """
        fresh = model.DashboardModel(
            self.model.pypeit_file, redux_path=str(self.model.redux_dir),
            derive=False)
        # Skip a transient mid-write read (the state file is rewritten on each
        # step transition); keep the last good views rather than flicker.
        if fresh.load_status == model.LOAD_MALFORMED:
            return
        self.model = fresh
        self.status_view.set_model(fresh)
        # refresh() (not set_model) so the selected step/frame is preserved
        # through live-monitoring updates.
        self.calibrations_view.refresh(fresh)
        self.science_view.refresh(fresh)

    def _on_run_finished(self, code):
        """
        Refresh once a (Re)Build run the Dashboard launched completes.

        A non-zero ``code`` means the run failed (e.g. an error crashed the
        reduction): record it and flag the Build channel as failed so the
        user sees it, rather than the run silently returning to "Idle".  The
        on-disk state is re-read regardless, so the failed step now shows red
        (its status is persisted as ``fail`` by the pipeline).

        Parameters
        ----------
        code : int
            The run's exit code (0 = success).
        """
        self._last_run_failed = code != 0
        self._refresh_from_state()
        if self._last_run_failed:
            self.activity.set_build(
                f'(Re)Build failed (exit code {code}) — see the log.',
                busy=False)

    def _on_science_frame_activated(self, frame, det):
        """
        Switch to the Science tab and select a frame, in response to a click
        on the Status-view science navigator.

        Parameters
        ----------
        frame : str
            The exposure basename.
        det : int or tuple
            Detector (int) or mosaic (tuple).
        """
        self.tab_widget.setCurrentWidget(self.science_view)
        self.science_view.select_frame(frame, det)

    def _on_lock_changed(self, locked):
        """
        React to the run lock engaging/releasing: gate the (Re)Build
        controls, show "Monitoring…" on the Build channel while active, and
        do a final refresh + idle when the run ends.

        Parameters
        ----------
        locked : bool
            Whether a run is now active.
        """
        self.calibrations_view.set_locked(locked)
        self.science_view.set_locked(locked)
        if locked:
            # A fresh run is (or may be) starting: clear any prior failure so
            # the Build channel monitors this run rather than the old result.
            self._last_run_failed = False
            self.activity.set_build('Monitoring run — updating live…',
                                    busy=True)
        elif self._last_run_failed:
            # The run that just ended failed: keep the failure message
            # visible instead of resetting to "Idle".
            self._refresh_from_state()
        else:
            # Run finished/stalled: one final refresh, then idle.
            self._refresh_from_state()
            self.activity.set_build('Idle', busy=False)

    def _on_state_changed(self):
        """
        Live-refresh the views when the state file changes during a run.
        """
        self._refresh_from_state()
        # Don't clobber a sticky failure message with "Monitoring…": a
        # crashed run can keep the lock held (recent .log mtime) for a few
        # more polls.
        if not self._last_run_failed:
            self.activity.set_build('Monitoring run — updating live…',
                                    busy=True)
