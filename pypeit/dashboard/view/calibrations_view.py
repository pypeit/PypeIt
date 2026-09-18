"""
The Calibrations view of the PypeIt Dashboard.

The drill-down companion to the Status view: scope selectors, a path-aware
step-button row, and a detail panel for the selected step — metrics, input
files, the output viewer, QA files, and a per-slit/order table.  It launches
the existing inspection tools as subprocesses via the
:class:`~pypeit.dashboard.launcher.Launcher` and stays thin: data via the
model, colors via the palette, commands via :mod:`pypeit.dashboard.inspect`.
"""

from pathlib import Path

from qtpy.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QLabel,
                            QComboBox, QPushButton, QTableWidget,
                            QTableWidgetItem, QHeaderView, QListWidget,
                            QListWidgetItem, QScrollArea, QGroupBox,
                            QAbstractItemView, QMessageBox)
from qtpy.QtCore import Qt
from qtpy.QtGui import QColor

from pypeit import log
from pypeit.dashboard import palette
# Direct function imports: importing the pypeit.dashboard.inspect module by
# its bare name would shadow the standard-library inspect module.
from pypeit.dashboard.inspect import (output_command, output_target,
                                      run_command, view_input_command)
from pypeit.dashboard.view.qa_dialog import QaImageDialog
from pypeit.dashboard.view.util import (text_on, detect_theme, clear_layout,
                                        fmt_value)
from pypeit.state.run_state import same_det


class CalibrationsView(QWidget):
    """
    The Calibrations tab: the per-(group, detector) calibration drill-down.

    Parameters
    ----------
    model : :class:`~pypeit.dashboard.model.DashboardModel`
        The reduction-state model.
    launcher : :class:`~pypeit.dashboard.launcher.Launcher`, optional
        Used to launch inspection tools and (re)build runs; if ``None``, the
        controls are built but do nothing (e.g. in headless tests).
    run_lock : :class:`~pypeit.dashboard.runlock.RunLock`, optional
        The single-run lock; the (Re)Build control is disabled while it is
        locked.
    on_run_finished : callable, optional
        Called as ``on_run_finished(code)`` when a (Re)Build run ends, so the
        window can refresh the state.
    parent : QWidget, optional
        The parent widget.
    """

    #: Steps that carry per-slit/order detail.
    _PER_SLIT_STEPS = ('slits', 'wv_calib', 'tilts', 'flats')

    #: Border color of the selected step button (from the design sketches).
    _MAGENTA = '#D81B60'

    #: Steps whose viewer opens a Ginga window (vs a matplotlib/Qt plot).
    _GINGA_HINT_STEPS = ('bias', 'dark', 'arc', 'tiltimg', 'slits',
                         'scattlight', 'wv_calib')

    def __init__(self, model, launcher=None, run_lock=None,
                 on_run_finished=None, parent=None):
        super().__init__(parent=parent)
        self._model = model
        self._launcher = launcher
        self._run_lock = run_lock
        self._on_run_finished = on_run_finished
        self._theme = 'light'
        self._scope = None              # (group, det)
        self._selected_step = None
        self._group_combo = None
        self._det_combo = None
        self._button_row = None         # QHBoxLayout holding the step buttons
        self._step_buttons = {}         # step -> QPushButton
        self._detail = None             # QVBoxLayout for the detail panel
        self._table = None              # the step-button-row's status lookup
        # The current detail panel's (Re)Build button + whether it has a
        # command, so the lock can enable/disable it without a full rebuild.
        self._rebuild_button = None
        self._rebuild_available = False
        self._rebuild_step = None       # step the (Re)Build button targets

        self._outer = QVBoxLayout(self)
        self._build()

    # -- public API ------------------------------------------------------

    def set_model(self, model):
        """
        Swap in a new model and rebuild.

        Parameters
        ----------
        model : :class:`~pypeit.dashboard.model.DashboardModel`
            The new model.
        """
        self._model = model
        self._scope = None
        self._selected_step = None
        clear_layout(self._outer)
        self._build()

    def refresh(self, model):
        """
        Swap in a new model but **preserve** the current scope and selected
        step (used after a (Re)Build completes, so the rebuilt step stays
        selected rather than resetting to the first step).

        Parameters
        ----------
        model : :class:`~pypeit.dashboard.model.DashboardModel`
            The freshly-acquired model.
        """
        prev_scope = self._scope
        prev_step = self._selected_step
        self._model = model
        self._scope = None
        # Keep the desired selection; _rebuild_button_row honors it if the
        # step is still present.
        self._selected_step = prev_step
        clear_layout(self._outer)
        self._build()
        # _build scopes to the first pair; restore the previous scope if it
        # still exists in the refreshed state.
        if prev_scope is not None and prev_scope in model.calib_det_pairs():
            self._set_scope(prev_scope)

    # -- build -----------------------------------------------------------

    def _build(self):
        """
        Build the view (edge message or full drill-down).
        """
        self._theme = detect_theme(self)
        if not self._model.is_started():
            msg = QLabel('No calibration state to inspect yet.')
            msg.setAlignment(Qt.AlignCenter)
            msg.setStyleSheet('color: grey; font-style: italic; '
                              'font-size: 16px;')
            self._outer.addStretch(1)
            self._outer.addWidget(msg)
            self._outer.addStretch(1)
            return

        self._build_scope_toolbar()

        # The step-button row (rebuilt on scope change).
        self._button_row = QHBoxLayout()
        row_wrap = QHBoxLayout()
        row_wrap.addLayout(self._button_row)
        row_wrap.addStretch(1)
        self._outer.addLayout(row_wrap)

        # The detail panel (rebuilt on step selection), in a scroll area.
        self._detail = QVBoxLayout()
        detail_container = QWidget()
        detail_container.setLayout(self._detail)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(detail_container)
        self._outer.addWidget(scroll, stretch=1)

        pairs = self._model.calib_det_pairs()
        if len(pairs) > 0:
            self._set_scope(pairs[0])

    def _build_scope_toolbar(self):
        """
        Build the calibration-group + detector drop-downs, independent of the
        Status tab's.
        """
        bar = QHBoxLayout()
        bar.addWidget(QLabel('Calibration group:'))
        self._group_combo = QComboBox()
        for g in sorted({g for g, _ in self._model.calib_det_pairs()}):
            self._group_combo.addItem(str(g), g)
        # Changing the group reloads the detector drop-down and rebuilds the
        # step-button row.
        self._group_combo.currentIndexChanged.connect(self._on_group_changed)
        bar.addWidget(self._group_combo)
        bar.addWidget(QLabel('Detector:'))
        self._det_combo = QComboBox()
        # Changing the detector rebuilds the step-button row.
        self._det_combo.currentIndexChanged.connect(self._on_det_changed)
        bar.addWidget(self._det_combo)
        bar.addStretch(1)
        self._outer.addLayout(bar)
        self._reload_det_combo()

    # -- scope handling --------------------------------------------------

    def _reload_det_combo(self, select=None):
        """
        Repopulate the detector drop-down for the current group.

        Parameters
        ----------
        select : int or tuple, optional
            Detector to select after repopulating, if present.
        """
        group = self._group_combo.currentData()
        dets = [d for g, d in self._model.calib_det_pairs() if g == group]
        self._det_combo.blockSignals(True)
        self._det_combo.clear()
        for d in dets:
            self._det_combo.addItem(self._model.det_name(d), d)
        if select is not None:
            idx = self._det_combo.findData(select)
            if idx >= 0:
                self._det_combo.setCurrentIndex(idx)
        self._det_combo.blockSignals(False)

    def _set_scope(self, pair):
        """
        Set the scoped ``(group, det)``, syncing the drop-downs and
        rebuilding the button row.

        Parameters
        ----------
        pair : tuple
            ``(group, det)``.
        """
        self._scope = pair
        group, det = pair
        self._group_combo.blockSignals(True)
        idx = self._group_combo.findData(group)
        if idx >= 0:
            self._group_combo.setCurrentIndex(idx)
        self._group_combo.blockSignals(False)
        self._reload_det_combo(select=det)
        self._rebuild_button_row()

    def _on_group_changed(self, _index):
        """
        Handle a group change: reload detectors, rebuild the button row.

        Parameters
        ----------
        _index : int
            The new combo index (unused).
        """
        self._reload_det_combo()
        self._scope = (self._group_combo.currentData(),
                       self._det_combo.currentData())
        self._rebuild_button_row()

    def _on_det_changed(self, _index):
        """
        Handle a detector change: rebuild the button row.

        Parameters
        ----------
        _index : int
            The new combo index (unused).
        """
        self._scope = (self._group_combo.currentData(),
                       self._det_combo.currentData())
        self._rebuild_button_row()

    # -- step-button row -------------------------------------------------

    def _scope_status(self):
        """
        Return a ``step -> (required, status, in_pipeline)`` lookup for the
        current scope.

        Returns
        -------
        dict
            Per-step status tuples for the selected (group, det).
        """
        group, det = self._scope
        table = self._model.status_table()
        lookup = {}
        for row in table:
            if row['calibration_group'] == group \
                    and same_det(row['detector'], det):
                lookup[row['step']] = (row['required'], row['status'],
                                       row['in_pipeline'])
        return lookup

    def _rebuild_button_row(self):
        """
        Rebuild the path-aware step-button row for the current scope.
        """
        if self._button_row is None:
            return
        clear_layout(self._button_row)
        self._step_buttons = {}
        status = self._scope_status()
        steps = self._model.step_order()
        for i, step in enumerate(steps):
            required, st, in_pipe = status.get(step, (None, 'absent', True))
            style = palette.step_style(required, st, in_pipe,
                                       theme=self._theme)
            button = QPushButton(f'{step}\n{style.glyph}')
            button.setMinimumWidth(86)
            button.setMinimumHeight(46)
            button._dash_style = style       # stash for re-styling on select
            # Clicking a step button selects that step (restyles the row and
            # rebuilds the detail panel).  Qt's clicked signal passes a
            # "checked" bool that must be swallowed, and the ``s=step``
            # default argument binds the current loop variable.
            button.clicked.connect(
                lambda _checked=False, s=step: self._select_step(s))
            self._step_buttons[step] = button
            self._button_row.addWidget(button)
            # Connector conveying precedence (not after the last button).
            if i < len(steps) - 1:
                self._button_row.addWidget(QLabel('→'))
        # Keep or reset the selected step.
        if self._selected_step not in self._step_buttons:
            self._selected_step = steps[0] if len(steps) > 0 else None
        self._restyle_buttons()
        if self._selected_step is not None:
            self._build_detail(self._selected_step)

    def _restyle_buttons(self):
        """
        Apply the palette fill + selected-ring to each step button.
        """
        for step, button in self._step_buttons.items():
            style = button._dash_style
            fg = text_on(style.color)
            if step == self._selected_step:
                # Magenta selection ring (design).
                border = f'3px solid {self._MAGENTA}'
            else:
                border = '1px solid #888'
            button.setStyleSheet(
                f'background-color: {style.color}; color: {fg}; '
                f'border: {border}; border-radius: 4px; padding: 2px;')

    def _select_step(self, step):
        """
        Select a step: re-style the buttons and rebuild the detail panel.

        Parameters
        ----------
        step : str
            The selected step.
        """
        self._selected_step = step
        self._restyle_buttons()
        self._build_detail(step)

    # -- detail panel ----------------------------------------------------

    def _build_detail(self, step):
        """
        Rebuild the detail panel for the selected step.

        Parameters
        ----------
        step : str
            The selected step.
        """
        clear_layout(self._detail)
        # The detail panel is rebuilt, so the old (Re)Build button is gone.
        self._rebuild_button = None
        self._rebuild_available = False
        group, det = self._scope
        entry = self._model.step_entry(step, group, det)

        heading = QLabel(f'{step}')
        heading.setStyleSheet('font-size: 15px; font-weight: bold;')
        self._detail.addWidget(heading)

        # Metrics — entry-level; per-slit metrics are in the table below.
        # Numeric values are shown to a few significant figures.
        metrics = self._model.step_metrics(step, group, det)
        if metrics is not None and len(metrics) > 0:
            text = '   '.join(f'{k}: {fmt_value(v)}'
                              for k, v in metrics.items())
            self._detail.addWidget(QLabel(text))

        # Action row (top): "Inspect output" + "(Re)Build", the step's two
        # primary actions side by side.
        self._add_action_row(step, group, det)

        # Output filename: show it for every step that has one — useful
        # especially for wv_calib, which has no viewer.
        self._add_output_label(entry)

        # QA files.
        self._add_qa_files(entry)

        # Per-slit/order drill-down.
        if step in self._PER_SLIT_STEPS:
            self._add_slit_table(step, group, det)

        # Input files live at the bottom, half-width and short; grouped by
        # role for flats, union otherwise.
        self._detail.addStretch(1)
        self._add_input_files(step, entry)

    def _add_action_row(self, step, group, det):
        """
        Add the action row: "Inspect output" beside "(Re)Build".

        Parameters
        ----------
        step : str
            The step name.
        group : int
            The calibration group.
        det : int or tuple
            The detector.
        """
        row = QHBoxLayout()
        row.addWidget(self._build_output_button(step, group, det))
        row.addWidget(self._build_rebuild_button(step, group, det))
        row.addStretch(1)
        self._detail.addLayout(row)

    def _add_output_label(self, entry):
        """
        Show the step's output filename, so it is visible even when there is
        no viewer (e.g. ``wv_calib``).

        Parameters
        ----------
        entry : object or None
            The step's state entry.
        """
        name = getattr(entry, 'output_file', None) if entry is not None \
            else None
        if name is None or name == '':
            return
        label = QLabel(f'Output: {Path(name).name}')
        label.setStyleSheet('color: grey;')
        # Selectable so the user can copy the filename.
        label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self._detail.addWidget(label)

    def _build_output_button(self, step, group, det):
        """
        Build the "Inspect output" button.

        Parameters
        ----------
        step : str
            The step name.
        group : int
            The calibration group.
        det : int or tuple
            The detector.

        Returns
        -------
        QPushButton
            The inspect-output button.
        """
        argv = output_command(self._model, step, group, det)
        target = output_target(self._model, step, group, det)
        # Enable when a viewer command exists and the file it opens is on
        # disk — not gated on entry.status (which would wrongly disable the
        # viewer for, e.g., skipped flats whose file still exists).
        available = argv is not None and target is not None \
            and target.exists()
        button = QPushButton('Inspect output')
        button.setEnabled(available)
        if available:
            # Bright, distinct (teal) style so it reads as an action like the
            # (Re)Build button; only when enabled (a disabled output greys
            # out normally).
            inspect = palette.inspect_color(self._theme)
            button.setStyleSheet(
                f'background-color: {inspect}; color: {text_on(inspect)}; '
                f'border-radius: 4px; padding: 4px 10px; font-weight: bold;')
            hint = self._viewer_hint(step)
            # Launch the step's output viewer.  Qt's clicked signal passes a
            # "checked" bool that must be swallowed; the default arguments
            # bind this step's command and window hint.
            button.clicked.connect(
                lambda _c=False, a=argv, h=hint:
                self._launch(a, f'view {step}', hint=h))
        return button

    def _build_rebuild_button(self, step, group, det):
        """
        Build the "(Re)Build" control: launch ``pypeit_run_to_calibstep`` for
        the selected step.  Styled with the distinct **action** color (blue),
        and disabled while a run is in progress (the single-run lock).

        Parameters
        ----------
        step : str
            The step name.
        group : int
            The calibration group.
        det : int or tuple
            The detector.

        Returns
        -------
        QPushButton
            The (Re)Build button.
        """
        argv = run_command(self._model, step, group, det)
        button = QPushButton(f'(Re)Build {step}')
        self._rebuild_button = button
        self._rebuild_step = step
        self._rebuild_available = argv is not None
        if self._rebuild_available:
            # Launch the confirmation + (re)build for this step; the lambda
            # swallows Qt's "checked" bool and binds the step/group/det.
            button.clicked.connect(
                lambda _c=False, s=step, g=group, d=det:
                self._on_rebuild(s, g, d))
        # Style for the current lock state (orange + disabled while a run is
        # active, blue + enabled when idle).
        locked = self._run_lock is not None and self._run_lock.is_locked()
        self._style_rebuild_button(locked)
        return button

    def _style_rebuild_button(self, locked):
        """
        Style the (Re)Build button for the lock state: a distinct blue
        **action** when idle, but **orange** and disabled (with a "run in
        progress" label) while *any* PypeIt run is active — the visual clue
        that it is locked.

        Parameters
        ----------
        locked : bool
            Whether a run is in progress.
        """
        button = self._rebuild_button
        if button is None:
            return
        if not self._rebuild_available:
            # No (re)build target for this step: a plain, disabled button.
            button.setEnabled(False)
            button.setStyleSheet('')
            return
        if locked:
            # Orange + disabled: the run is active, so this is not
            # selectable.
            color = palette.LIGHT_COLORS[palette.RUNNING] \
                if self._theme == 'light' \
                else palette.DARK_COLORS[palette.RUNNING]
            glyph = palette.GLYPHS[palette.RUNNING][0]
            button.setText(f'{glyph} Run in progress')
            button.setEnabled(False)
            button.setToolTip('A PypeIt run is in progress; (re)build is '
                              'disabled until it finishes.')
        else:
            color = palette.action_color(self._theme)
            button.setText(f'(Re)Build {self._rebuild_step}')
            button.setEnabled(True)
            button.setToolTip('Regenerate this step (and any preceding '
                              'steps), overwriting its output.')
        button.setStyleSheet(
            f'background-color: {color}; color: {text_on(color)}; '
            f'border-radius: 4px; padding: 4px 10px; font-weight: bold;')

    def set_locked(self, locked):
        """
        Restyle the (Re)Build control when the run lock changes: orange +
        disabled while locked, blue + enabled when idle.

        Parameters
        ----------
        locked : bool
            Whether a run is in progress.
        """
        self._style_rebuild_button(locked)

    def _on_rebuild(self, step, group, det):
        """
        Confirm and launch a (re)build of the selected step.

        Shows a clobber confirmation naming the file(s) that will be
        overwritten; on confirm, moves only the selected step's output(s)
        aside and launches ``pypeit_run_to_calibstep`` through the launcher,
        which holds the lock and refreshes on completion.

        The move-aside is **crash-safe**: each output is renamed to a
        ``.dashboard_bak`` sibling rather than unlinked, and
        :meth:`_after_rebuild` restores it if the run **fails** (so a failed
        (re)build never loses the existing calibration).

        Parameters
        ----------
        step : str
            The step name.
        group : int
            The calibration group.
        det : int or tuple
            The detector.
        """
        # Refuse if a run is already active (defensive; the button is also
        # disabled while locked).
        if self._run_lock is not None and self._run_lock.is_locked():
            return
        argv = run_command(self._model, step, group, det)
        if argv is None:
            return
        existing = self._model.step_output_files(step, group, det)
        if not self._confirm_rebuild(step, existing):
            return
        # Move the selected step's output(s) aside so the reused
        # run_to_calibstep genuinely rebuilds them, keeping a backup to
        # restore if the run fails (never a shell rm).
        backups = []
        for path in existing:
            bak = path.parent / (path.name + '.dashboard_bak')
            try:
                if bak.exists():
                    bak.unlink()
                path.rename(bak)
                backups.append((bak, path))
            except OSError as exc:
                log.warning(f'Could not move {path} aside: {exc}')
        if self._launcher is not None:
            # Launching engages the run lock, which turns the (Re)Build
            # button orange + disabled via set_locked().
            self._launcher.run(
                argv, description=f'(Re)building {step}',
                on_finished=lambda code, s=step, b=backups:
                self._after_rebuild(code, s, b))

    def _after_rebuild(self, code, step, backups):
        """
        Finish a (re)build: drop the backups on success, restore them on
        failure, keep the rebuilt step selected, then trigger the window
        refresh.

        Parameters
        ----------
        code : int
            The run's exit code (0 = success).
        step : str
            The step that was (re)built — kept selected so the refreshed
            view shows it, not the first step.
        backups : list
            ``(backup_path, original_path)`` pairs moved aside before the
            run.
        """
        for bak, original in backups:
            try:
                if code == 0 or original.exists():
                    # Success (or the run wrote a fresh output anyway): the
                    # backup is no longer needed.
                    bak.unlink()
                else:
                    # Failed and nothing was rebuilt: restore the original so
                    # the calibration is not lost.
                    bak.rename(original)
            except OSError as exc:
                log.warning(f'Could not finalize backup {bak}: {exc}')
        # Keep the rebuilt step selected through the refresh.
        self._selected_step = step
        if self._on_run_finished is not None:
            self._on_run_finished(code)

    def _confirm_rebuild(self, step, existing):
        """
        Show the clobber confirmation dialog, naming the file(s) to be
        overwritten.

        Parameters
        ----------
        step : str
            The step name.
        existing : list
            Existing output `Path` files.

        Returns
        -------
        bool
            True if the user confirmed the (re)build.
        """
        if len(existing) > 0:
            files = '\n'.join(f'  • {p.name}' for p in existing)
            text = (f'(Re)build "{step}" — this will overwrite the existing '
                    f'output file(s):\n\n{files}\n\nPreceding calibrations '
                    f'are reused. Continue?')
            icon = QMessageBox.Warning
        else:
            text = (f'(Re)build "{step}"? No existing output for this step '
                    f'was found, so it (and any preceding steps) will be '
                    f'built fresh. Continue?')
            icon = QMessageBox.Question
        box = QMessageBox(self)
        box.setIcon(icon)
        box.setWindowTitle('Confirm (Re)Build')
        box.setText(text)
        box.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
        box.setDefaultButton(QMessageBox.Cancel)
        return box.exec_() == QMessageBox.Ok

    def _viewer_hint(self, step):
        """
        Return where a step's output viewer shows its result, for the
        activity bar's finished message.

        Parameters
        ----------
        step : str
            The step name.

        Returns
        -------
        str
            ``'Ginga window'`` or ``'plot window'``.
        """
        return 'Ginga window' if step in self._GINGA_HINT_STEPS \
            else 'plot window'

    def _add_input_files(self, step, entry):
        """
        Add the input-file list; grouped by role for flats.

        Parameters
        ----------
        step : str
            The step name.
        entry : object or None
            The step's state entry.
        """
        if entry is None:
            return
        listw = QListWidget()
        # Keep the list short; it scrolls if there are many inputs.
        listw.setMaximumHeight(90)
        if step == 'flats':
            # Grouped by role, with the pixel-flat provenance noted.
            groups = [('Pixelflat', getattr(entry, 'pixelflat_files', None)),
                      ('Illum', getattr(entry, 'illumflat_files', None)),
                      ('Lamp-off', getattr(entry, 'lampoff_files', None))]
            source = getattr(entry, 'pixelflat_source', None)
            if source is not None:
                groups[0] = (f'Pixelflat ({source})', groups[0][1])
            for role, files in groups:
                for f in (files if files is not None else []):
                    self._add_file_item(listw, f, prefix=f'[{role}] ')
        else:
            input_files = entry.input_files \
                if entry.input_files is not None else []
            for f in input_files:
                self._add_file_item(listw, f)
        # Open the double-clicked input frame in pypeit_view_fits.
        listw.itemDoubleClicked.connect(self._on_input_double_clicked)
        # Put the list in a titled box so the section reads as
        # self-contained, at half width — the right half is free.
        box = QGroupBox('Input files (double-click to view)')
        box_layout = QVBoxLayout(box)
        box_layout.addWidget(listw)
        wrap = QHBoxLayout()
        wrap.addWidget(box, stretch=1)
        wrap.addStretch(1)
        self._detail.addLayout(wrap)

    def _add_file_item(self, listw, path, prefix=''):
        """
        Add one file row to an input-file list, stashing the path.

        Parameters
        ----------
        listw : QListWidget
            The list to add to.
        path : str
            The file path.
        prefix : str, optional
            Display prefix (e.g. a role tag).
        """
        item = QListWidgetItem(f'{prefix}{Path(path).name}')
        item.setData(Qt.UserRole, str(path))
        listw.addItem(item)

    def _on_input_double_clicked(self, item):
        """
        View the double-clicked input frame via ``pypeit_view_fits``.

        Parameters
        ----------
        item : QListWidgetItem
            The clicked item.
        """
        path = item.data(Qt.UserRole)
        det = self._scope[1] if self._scope is not None else None
        argv = view_input_command(self._model, path, det=det)
        self._launch(argv, 'view input', hint='Ginga window')

    def _add_qa_files(self, entry):
        """
        Add the QA-file list; double-click opens the PNG full-view.

        Parameters
        ----------
        entry : object or None
            The step's state entry.
        """
        qa_files = getattr(entry, 'qa_files', None) if entry is not None \
            else None
        if qa_files is None or len(qa_files) == 0:
            return
        self._detail.addWidget(QLabel('QA files (double-click to open):'))
        listw = QListWidget()
        listw.setMaximumHeight(110)
        for f in qa_files:
            self._add_file_item(listw, f)
        # Open the double-clicked QA PNG in the full-view image dialog.
        listw.itemDoubleClicked.connect(self._on_qa_double_clicked)
        self._detail.addWidget(listw)

    def _on_qa_double_clicked(self, item):
        """
        Open the double-clicked QA PNG in an image dialog.

        Parameters
        ----------
        item : QListWidgetItem
            The clicked item.
        """
        path = Path(item.data(Qt.UserRole))
        if not path.is_absolute() and not path.exists():
            # QA PNGs live under <redux_dir>/QA/PNGs/.
            candidate = self._model.redux_dir / 'QA' / 'PNGs' / path.name
            if candidate.exists():
                path = candidate
        dialog = QaImageDialog(str(path), parent=self)
        dialog.show()

    def _add_slit_table(self, step, group, det):
        """
        Add the per-slit/order table.  For flats the table has per-correction
        mean/RMS columns.

        Parameters
        ----------
        step : str
            The step name.
        group : int
            The calibration group.
        det : int or tuple
            The detector.
        """
        rows = self._model.slit_table(step, group, det)
        if rows is None or len(rows) == 0:
            return
        # For the Echelle pipeline these rows are spectral orders, so label
        # them "Order" rather than "Slit".
        is_echelle = self._model.header_info is not None \
            and self._model.header_info.path == 'Echelle'
        unit = 'order' if is_echelle else 'slit'
        self._detail.addWidget(QLabel(f'Per-{unit} ({len(rows)}):'))

        # Build columns: Slit/Order, Status, then step-specific metric
        # columns.
        if step == 'slits':
            metric_cols = ['center', 'slitord_id']
        elif step in ('wv_calib', 'tilts'):
            metric_cols = ['rms']
        else:  # flats: a mean/rms pair per present correction
            corrections = self._model.step_metrics(
                step, group, det).get('corrections', [])
            metric_cols = []
            for corr in corrections:
                metric_cols += [f'{corr} mean', f'{corr} rms']

        columns = [unit.capitalize(), 'Status'] + metric_cols
        table = QTableWidget()
        table.setColumnCount(len(columns))
        table.setHorizontalHeaderLabels(columns)
        table.setRowCount(len(rows))
        table.verticalHeader().setVisible(False)
        table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        table.setSelectionMode(QAbstractItemView.NoSelection)
        table.setAlternatingRowColors(True)
        # Neutral selection fill, consistent with the science tables.
        table.setStyleSheet(palette.selection_style(self._theme))
        table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)

        for r, row in enumerate(rows):
            slit_item = QTableWidgetItem(str(row['slit']))
            slit_item.setTextAlignment(Qt.AlignCenter)
            table.setItem(r, 0, slit_item)
            style = palette.slit_style(row['status'], theme=self._theme)
            status_item = QTableWidgetItem(f'{style.glyph} {style.label}')
            status_item.setTextAlignment(Qt.AlignCenter)
            if style.category != palette.REQUIRED_UNDONE:
                status_item.setForeground(QColor(style.color))
            table.setItem(r, 1, status_item)
            for c, col in enumerate(metric_cols, start=2):
                cell = QTableWidgetItem(self._slit_cell(row, col))
                cell.setTextAlignment(Qt.AlignCenter)
                table.setItem(r, c, cell)
        # Keep the per-slit table from collapsing for many slits.
        table.setMinimumHeight(160)
        self._detail.addWidget(table)

    def _slit_cell(self, row, col):
        """
        Format one per-slit metric cell.

        Parameters
        ----------
        row : dict
            The slit_table row.
        col : str
            The column key (e.g. ``rms`` or ``pixelflat mean``).

        Returns
        -------
        str
            The cell text (em-dash when missing).
        """
        if col in ('center', 'slitord_id', 'rms'):
            return fmt_value(row.get(col))
        # flats: "<correction> mean" / "<correction> rms"
        corr, _, field = col.rpartition(' ')
        metric = row.get('corrections', {}).get(corr, {})
        return fmt_value(metric.get(field))

    def _launch(self, argv, description, hint='viewer window'):
        """
        Launch a command via the launcher (no-op if no launcher / command).

        Parameters
        ----------
        argv : list or None
            The command argv.
        description : str
            Human description for the activity bar.
        hint : str, optional
            Where the result appears (e.g. ``'Ginga window'``).
        """
        if argv is not None and self._launcher is not None:
            self._launcher.launch(argv, description=description, hint=hint)
