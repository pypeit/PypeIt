"""
The Status (Initialization) view of the PypeIt Dashboard.

This is the landing, *state-first* view.  It renders a
:class:`~pypeit.dashboard.model.DashboardModel` and stays thin: all data comes
from the model, all colors/glyphs from :mod:`pypeit.dashboard.palette`.  Top
to bottom it shows a global summary strip, a scope toolbar
(calibration-group + detector drop-downs, a stale badge, a Refresh button), a
configuration-overview navigator grid, the scoped calibration status table,
and a science-frames section with a clickable per-frame navigator.

.. include:: ../include/links.rst
"""

from qtpy.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
                            QLabel, QComboBox, QPushButton, QTableWidget,
                            QTableWidgetItem, QHeaderView, QAbstractItemView,
                            QFrame, QScrollArea)
from qtpy.QtCore import Qt, Signal
from qtpy.QtGui import QColor

from pypeit.dashboard import palette
from pypeit.dashboard import model
from pypeit.dashboard.view.util import text_on, detect_theme, clear_layout
from pypeit.state.run_state import same_det


class ScienceNavCell(QFrame):
    """
    One clickable science-navigator cell for the Status view: a four-segment
    mini strip (``process`` / ``findobj`` / ``skysub`` / ``extract``, each
    colored+glyphed by the palette) above a compact frame caption.  Clicking
    it emits :attr:`clicked` with the ``(frame, det)`` so the window can
    switch to the Science tab.

    Parameters
    ----------
    frame : str
        The exposure basename.
    det : int or tuple
        Detector (int) or mosaic (tuple/list).
    det_name : str
        Human detector name (e.g. ``DET01``).
    objtype : str
        ``'science'`` or ``'standard'``.
    statuses : dict
        ``{step: status}`` for the four macro-steps.
    theme : str
        ``'light'`` or ``'dark'``.
    parent : QWidget, optional
        The parent widget.
    """

    #: The four science macro-steps shown as segments of the mini strip (and
    #: summarized elsewhere in this view), in reduction order.
    STRIP_STEPS = ('process', 'findobj', 'skysub', 'extract')

    #: Emitted as ``clicked(frame, det)`` when the cell is clicked.
    clicked = Signal(object, object)

    def __init__(self, frame, det, det_name, objtype, statuses, theme='light',
                 parent=None):
        super().__init__(parent=parent)
        self._frame = frame
        self._det = det
        self.setFrameShape(QFrame.StyledPanel)
        self.setCursor(Qt.PointingHandCursor)
        self.setToolTip('\n'.join(
            [f'{frame} ({det_name}) — {objtype}']
            + [f'  {s}: {statuses.get(s, "pending")}'
               for s in self.STRIP_STEPS]))

        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(2)

        # The four-segment status strip.
        strip = QHBoxLayout()
        strip.setSpacing(1)
        for step in self.STRIP_STEPS:
            style = palette.slit_style(statuses.get(step, 'pending'),
                                       theme=theme)
            seg = QLabel(style.glyph)
            seg.setAlignment(Qt.AlignCenter)
            seg.setFixedSize(22, 20)
            seg.setToolTip(f'{step}: {statuses.get(step, "pending")}')
            seg.setStyleSheet(
                f'background-color: {style.color}; '
                f'color: {text_on(style.color)}; border: 1px solid #888;')
            strip.addWidget(seg)
        strip.addStretch(1)
        layout.addLayout(strip)

        # Compact caption: a shortened frame name + a standard tag.
        short = frame if len(frame) <= 18 else frame[:17] + '…'
        tag = ' ⭑' if objtype == 'standard' else ''
        caption = QLabel(f'{short}{tag}')
        caption.setStyleSheet('font-size: 11px;')
        caption.setToolTip(f'{frame} ({det_name})')
        layout.addWidget(caption)

    def mousePressEvent(self, event):
        """
        Emit :attr:`clicked` with ``(frame, det)`` on any mouse press.

        Parameters
        ----------
        event : QMouseEvent
            The mouse event.
        """
        self.clicked.emit(self._frame, self._det)
        super().mousePressEvent(event)


class StatusView(QWidget):
    """
    The Status tab: the state-first overview of a reduction.

    Parameters
    ----------
    model : :class:`~pypeit.dashboard.model.DashboardModel`
        The reduction-state model to render.
    activity : :class:`~pypeit.dashboard.view.activity.ActivityBar`, optional
        The shared activity bar, used to report what a Refresh did.
    parent : QWidget, optional
        The parent widget.
    """

    #: Columns of the scoped calibration status table.
    _TABLE_COLUMNS = ('Step', 'Required', 'Status', 'Output')

    #: The science macro-step columns (for the compact summary; the full
    #: per-frame drill-down is the Science tab).
    _SCI_STEP_COLS = ScienceNavCell.STRIP_STEPS

    #: Message shown when the reduction has not produced a state yet.
    _NOT_STARTED_MSG = ('This reduction has not been started — no '
                        'calibration state found.')

    #: Human-readable messages for each empty/edge load status.
    _EDGE_MESSAGES = {
        model.LOAD_NOT_STARTED: _NOT_STARTED_MSG,
        # A state file that loaded but has no entries is also "not started".
        model.LOAD_STATE_FILE: _NOT_STARTED_MSG,
        model.LOAD_FILE_NOT_FOUND:
            'The PypeIt file was not found.',
        model.LOAD_MALFORMED:
            'The reduction state file could not be read (malformed or '
            'out-of-date).',
        model.LOAD_ERROR:
            'The reduction state could not be derived.',
    }

    #: Emitted as ``scienceFrameActivated(frame, det)`` when the user clicks
    #: a science-navigator cell; the window switches to the Science tab and
    #: selects that frame.
    scienceFrameActivated = Signal(object, object)

    def __init__(self, model, activity=None, parent=None):
        super().__init__(parent=parent)
        self._model = model
        self._activity = activity
        self._theme = 'light'
        # The currently scoped (calibration_group, detector) pair.
        self._scope = None
        # Kept so the navigator/drop-down handlers can update without
        # rebuilding the whole view.
        self._group_combo = None
        self._det_combo = None
        self._table = None

        self._outer = QVBoxLayout(self)
        self._build()

    # -- public API ------------------------------------------------------

    def set_model(self, model):
        """
        Swap in a new model and re-render (used by Refresh).

        Parameters
        ----------
        model : :class:`~pypeit.dashboard.model.DashboardModel`
            The new model.
        """
        self._model = model
        self._scope = None
        clear_layout(self._outer)
        self._build()

    # -- build -----------------------------------------------------------

    def _build(self):
        """
        Build the view for the current model (edge message or full view).
        """
        self._theme = detect_theme(self)
        if not self._model.is_started():
            self._build_edge_message()
            self._outer.addStretch(1)
            return
        self._build_summary()
        self._build_scope_toolbar()
        self._build_navigator()
        self._build_table_section()
        self._build_science_section()
        # Initial scope: the first available pair.
        pairs = self._model.calib_det_pairs()
        if len(pairs) > 0:
            self._set_scope(pairs[0])

    def _build_edge_message(self):
        """
        Render a clear, centered message for an empty/edge state.
        """
        msg = self._EDGE_MESSAGES.get(self._model.load_status,
                                      'No reduction state available.')
        label = QLabel(msg)
        label.setAlignment(Qt.AlignCenter)
        label.setWordWrap(True)
        label.setStyleSheet('color: grey; font-style: italic; '
                            'font-size: 16px;')
        self._outer.addStretch(1)
        self._outer.addWidget(label)

    def _build_summary(self):
        """
        Build the always-visible global summary strip: counts of required
        calibrations succeeded/failed/running/to-do, extended with the
        science-frame counts when any frame has been reduced.
        """
        table = self._model.status_table()
        # Whole-run health: required steps that are in the pipeline.
        req = [row for row in table
               if row['required'] is True and row['in_pipeline']]
        n_req = len(req)
        n_ok = sum(row['status'] == 'success' for row in req)
        n_fail = sum(status == 'fail' for status in table['status'])
        n_run = sum(status == 'running' for status in table['status'])
        n_undone = sum(row['status'] in ('pending', 'absent') for row in req)
        text = (f'Calibrations: {n_ok}/{n_req} required succeeded'
                f'   ·   {n_fail} failed   ·   {n_run} running'
                f'   ·   {n_undone} to-do')
        sci = self._model.science_table()
        if len(sci) > 0:
            n_sci = len(sci)
            n_extracted = sum(row['extract'] == 'success' for row in sci)
            n_sci_fail = sum(any(row[col] == 'fail'
                                 for col in self._SCI_STEP_COLS)
                             for row in sci)
            text += (f'\nScience: {n_extracted}/{n_sci} extracted'
                     f'   ·   {n_sci_fail} failed')
        label = QLabel(text)
        label.setObjectName('summaryStrip')
        label.setStyleSheet('font-weight: bold; padding: 6px;')
        self._outer.addWidget(label)

    def _build_scope_toolbar(self):
        """
        Build the scope toolbar: calibration-group + detector drop-downs, a
        stale badge (shown only when the state looks out of date), and a
        Refresh button.
        """
        bar = QHBoxLayout()
        bar.addWidget(QLabel('Calibration group:'))
        self._group_combo = QComboBox()
        groups = sorted({g for g, _ in self._model.calib_det_pairs()})
        for g in groups:
            self._group_combo.addItem(str(g), g)
        # Changing the group reloads the detector drop-down and re-renders
        # the scoped table.
        self._group_combo.currentIndexChanged.connect(self._on_group_changed)
        bar.addWidget(self._group_combo)

        bar.addWidget(QLabel('Detector:'))
        self._det_combo = QComboBox()
        # Changing the detector re-renders the scoped table.
        self._det_combo.currentIndexChanged.connect(self._on_det_changed)
        bar.addWidget(self._det_combo)

        bar.addStretch(1)

        # Stale badge, shown only when the state looks out of date.
        if self._model.is_stale():
            stale = QLabel('⚠ state may be stale')
            stale.setStyleSheet('color: #EF6C00; font-weight: bold;')
            bar.addWidget(stale)

        refresh = QPushButton('Refresh')
        # The Refresh button re-acquires the state and re-renders the view.
        refresh.clicked.connect(self._on_refresh)
        bar.addWidget(refresh)

        self._outer.addLayout(bar)
        # Populate detectors for the initial group.
        self._reload_det_combo()

    def _build_navigator(self):
        """
        Build the configuration-overview navigator grid: one cell per
        ``(group, detector)``, colored by the worst step status, clickable to
        scope the status table.
        """
        pairs = self._model.calib_det_pairs()
        groups = sorted({g for g, _ in pairs})
        dets = sorted({d for _, d in pairs}, key=str)

        grid = QGridLayout()
        grid.addWidget(QLabel(''), 0, 0)
        for col, det in enumerate(dets, start=1):
            grid.addWidget(QLabel(self._model.det_name(det)), 0, col)
        for rrow, group in enumerate(groups, start=1):
            grid.addWidget(QLabel(f'Group {group}'), rrow, 0)
            for col, det in enumerate(dets, start=1):
                if (group, det) in pairs:
                    grid.addWidget(self._navigator_cell(group, det),
                                   rrow, col)
        # Keep the navigator compact (left-aligned), not stretched wide.
        wrap = QHBoxLayout()
        wrap.addLayout(grid)
        wrap.addStretch(1)
        self._outer.addLayout(wrap)

    def _navigator_cell(self, group, det):
        """
        Build one clickable, status-colored navigator cell.

        Parameters
        ----------
        group : int
            Calibration group ID.
        det : int or tuple
            Detector (int) or mosaic (tuple).

        Returns
        -------
        QPushButton
            The cell button.
        """
        category = self._cell_category(group, det)
        color = palette.LIGHT_COLORS[category] if self._theme == 'light' \
            else palette.DARK_COLORS[category]
        glyph, _ = palette.GLYPHS[category]
        button = QPushButton(glyph)
        button.setFixedSize(40, 40)
        # Color the cell by the worst status; pair it with the glyph so the
        # status also reads without color.
        button.setStyleSheet(
            f'background-color: {color}; border: 1px solid #888;')
        button.setToolTip(f'Group {group}, {self._model.det_name(det)}')
        # Clicking a cell scopes the table to its (group, det).  Qt's
        # clicked signal passes a "checked" bool that must be swallowed, and
        # the default arguments bind the current loop variables.
        button.clicked.connect(
            lambda _checked=False, g=group, d=det: self._set_scope((g, d)))
        return button

    def _cell_category(self, group, det):
        """
        The worst palette category among a navigator cell's steps.

        Parameters
        ----------
        group : int
            Calibration group ID.
        det : int or tuple
            Detector (int) or mosaic (tuple).

        Returns
        -------
        str
            The worst palette category for that ``(group, det)``.
        """
        rows = self._scoped_table(group, det)
        cats = [palette.classify(row['required'], row['status'],
                                 row['in_pipeline'])
                for row in rows]
        return palette.worst_category(cats)

    def _build_table_section(self):
        """
        Build the scoped calibration status table.
        """
        heading = QLabel('Calibrations')
        heading.setStyleSheet('font-size: 15px; font-weight: bold; '
                              'padding-top: 8px;')
        self._outer.addWidget(heading)

        self._table = QTableWidget()
        self._table.setColumnCount(len(self._TABLE_COLUMNS))
        self._table.setHorizontalHeaderLabels(list(self._TABLE_COLUMNS))
        self._table.verticalHeader().setVisible(False)
        self._table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self._table.setSelectionMode(QAbstractItemView.NoSelection)
        self._table.setAlternatingRowColors(True)
        header = self._table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(2, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(3, QHeaderView.Stretch)
        # Neutral selection fill, consistent with the science tables (applied
        # for consistency even though this table is non-selectable).
        self._table.setStyleSheet(palette.selection_style(self._theme))
        # Let the table grow to show its rows (stretch in the column).
        self._outer.addWidget(self._table, stretch=1)

    def _build_science_section(self):
        """
        Build the Science-frames section of the Status view: a one-line count
        plus a science navigator grid mirroring the calibration navigator —
        one clickable four-segment cell per ``(frame, det)``, science frames
        first then standards, scrollable.
        """
        heading = QLabel('Science frames')
        heading.setStyleSheet('font-size: 15px; font-weight: bold; '
                              'padding-top: 8px;')
        self._outer.addWidget(heading)
        sci = self._model.science_table()
        if len(sci) == 0:
            stub = QLabel('No science frames reduced yet.')
            stub.setStyleSheet('color: grey; font-style: italic;')
            self._outer.addWidget(stub)
            return
        n_sci = len(sci)
        n_std = sum(objtype == 'standard' for objtype in sci['objtype'])
        n_extracted = sum(status == 'success' for status in sci['extract'])
        n_fail = sum(any(row[col] == 'fail' for col in self._SCI_STEP_COLS)
                     for row in sci)
        summary = QLabel(
            f'{n_sci} frame(s) ({n_std} standard) · {n_extracted} extracted · '
            f'{n_fail} failed — click a cell for the Science tab.')
        self._outer.addWidget(summary)
        self._build_science_navigator(sci)

    def _build_science_navigator(self, sci):
        """
        Build the science navigator grid: a wrapped grid of
        :class:`ScienceNavCell`, science frames first then standards, in a
        capped scroll area so many frames stay compact.

        Parameters
        ----------
        sci : `astropy.table.Table`_
            The science table
            (:meth:`~pypeit.dashboard.model.DashboardModel.science_table`).
        """
        # Science frames first, then standards (stable within each group).
        rows = list(sci)
        rows.sort(key=lambda r: r['objtype'] == 'standard')

        grid = QGridLayout()
        grid.setSpacing(4)
        n_cols = 4                       # wrap into rows of four cells
        for i, row in enumerate(rows):
            statuses = {s: row[s] for s in ScienceNavCell.STRIP_STEPS}
            cell = ScienceNavCell(
                row['frame'], row['detector'],
                self._model.det_name(row['detector']),
                row['objtype'], statuses, theme=self._theme)
            # Forward a cell click as this view's scienceFrameActivated, so
            # the window can switch to the Science tab and select the frame.
            cell.clicked.connect(self.scienceFrameActivated)
            grid.addWidget(cell, i // n_cols, i % n_cols)

        holder = QWidget()
        holder.setLayout(grid)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(holder)
        # Cap the height so a long list scrolls rather than dominating the
        # view.
        scroll.setMaximumHeight(220)
        self._outer.addWidget(scroll)

    # -- scope handling --------------------------------------------------

    def _scoped_table(self, group, det):
        """
        Return the status table restricted to one ``(group, det)``, ordered
        by the path-aware step order.

        Parameters
        ----------
        group : int
            Calibration group ID.
        det : int or tuple
            Detector (int) or mosaic (tuple).

        Returns
        -------
        list
            The scoped, step-ordered rows (one `astropy.table.Row`_ per
            step).
        """
        table = self._model.status_table()
        order = {s: i for i, s in enumerate(self._model.step_order())}
        # Keep only the rows for this (group, det) whose step is in the
        # button-row order (drops bpm and any step not part of the
        # pipeline's display order), sorted by that order.
        scoped = [row for row in table
                  if row['calibration_group'] == group
                  and same_det(row['detector'], det)
                  and row['step'] in order]
        scoped.sort(key=lambda row: order[row['step']])
        return scoped

    def _set_scope(self, pair):
        """
        Set the scoped ``(group, det)``, syncing the drop-downs and table.

        Parameters
        ----------
        pair : tuple
            The ``(group, det)`` to scope to.
        """
        self._scope = pair
        group, det = pair
        # Sync the drop-downs without retriggering a re-render loop.
        self._group_combo.blockSignals(True)
        idx = self._group_combo.findData(group)
        if idx >= 0:
            self._group_combo.setCurrentIndex(idx)
        self._group_combo.blockSignals(False)
        self._reload_det_combo(select=det)
        self._render_table()

    def _reload_det_combo(self, select=None):
        """
        Repopulate the detector drop-down for the current group.

        Parameters
        ----------
        select : int or tuple, optional
            Detector to select after repopulating, if present.
        """
        if self._det_combo is None or self._group_combo is None:
            return
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

    def _on_group_changed(self, _index):
        """
        Handle a calibration-group change: reload detectors, re-render.

        Parameters
        ----------
        _index : int
            The new combo index (unused).
        """
        self._reload_det_combo()
        self._render_table()

    def _on_det_changed(self, _index):
        """
        Handle a detector change: re-render the table.

        Parameters
        ----------
        _index : int
            The new combo index (unused).
        """
        self._render_table()

    def _on_refresh(self):
        """
        Re-acquire the state via a fresh model and re-render.
        """
        fresh = model.DashboardModel(self._model.pypeit_file,
                                     redux_path=str(self._model.redux_dir))
        self.set_model(fresh)
        # Report what the refresh did to the shared activity area.
        if self._activity is not None:
            if fresh.load_status == model.LOAD_STATE_FILE:
                self._activity.set_build('Reloaded state file.')
            elif fresh.load_status == model.LOAD_DERIVED:
                self._activity.set_build('Re-derived state (no state file).')
            else:
                self._activity.set_build(f'Refreshed: {fresh.load_status}.')

    def _render_table(self):
        """
        Render the scoped status table from the current drop-down selection.
        """
        if self._table is None or self._group_combo is None:
            return
        group = self._group_combo.currentData()
        det = self._det_combo.currentData()
        if group is None or det is None:
            return
        self._scope = (group, det)
        scoped = self._scoped_table(group, det)
        self._table.setRowCount(len(scoped))
        for r, row in enumerate(scoped):
            self._set_table_row(r, row)

    def _set_table_row(self, r, row):
        """
        Populate one table row (Step | Required | Status | Output) with the
        palette color + glyph.

        Parameters
        ----------
        r : int
            The row index.
        row : `astropy.table.Row`_
            A row of the scoped status table.
        """
        required = row['required']
        status = row['status']
        in_pipeline = row['in_pipeline']
        category = palette.classify(required, status, in_pipeline)
        style = palette.step_style(required, status, in_pipeline,
                                   theme=self._theme)

        step_item = QTableWidgetItem(row['step'])
        req_text = 'Yes' if required is True else \
            ('No' if required is False else '—')
        req_item = QTableWidgetItem(req_text)
        status_item = QTableWidgetItem(f'{style.glyph} {style.label}')
        # Color the status text by the palette, but keep "required, not
        # done" (white fill) readable as text via a neutral color; the ○
        # glyph + label still convey it.
        if category != palette.REQUIRED_UNDONE:
            status_item.setForeground(QColor(style.color))
        # output_file is a basename string or missing (None); show an
        # em-dash for anything that is not a non-empty string.
        out = row['output_file']
        out_text = out if isinstance(out, str) and out != '' else '—'
        out_item = QTableWidgetItem(out_text)
        out_item.setToolTip(out_text if out_text != '—' else '')

        # Dim optional (not-required) rows so a pending optional never reads
        # as a problem.
        if required is False:
            for it in (step_item, req_item):
                it.setForeground(QColor('#9E9E9E'))

        for col, it in enumerate((step_item, req_item, status_item,
                                  out_item)):
            # Center the cell contents.
            it.setTextAlignment(Qt.AlignCenter)
            self._table.setItem(r, col, it)
