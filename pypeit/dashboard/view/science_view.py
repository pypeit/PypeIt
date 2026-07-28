"""
The Science view of the PypeIt Dashboard.

The science companion to the Calibrations view: a per-frame **table** (one row
per reduced ``(frame, detector)`` exposure, science and standard) with the
four macro-step statuses (``process`` / ``findobj`` / ``skysub`` /
``extract``) as color+glyph cells, ``nobj``, and product presence; and a
per-frame **detail panel** with per-slit and per-object tables plus the
product viewers and a **(Re)Build** control.  It stays thin: data via the
model, colors via the palette, commands via
:mod:`pypeit.dashboard.inspect`, launches via the shared
:class:`~pypeit.dashboard.launcher.Launcher`.

.. include:: ../include/links.rst
"""

from pathlib import Path

from qtpy.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QLabel,
                            QPushButton, QTableWidget, QTableWidgetItem,
                            QHeaderView, QScrollArea, QAbstractItemView,
                            QMessageBox, QListWidget, QListWidgetItem)
from qtpy.QtCore import Qt
from qtpy.QtGui import QColor

from pypeit.dashboard import palette
# Direct function imports: importing the pypeit.dashboard.inspect module by
# its bare name would shadow the standard-library inspect module.
from pypeit.dashboard.inspect import (spec2d_command, spec1d_command,
                                      science_run_command,
                                      run_pypeit_command,
                                      science_object_name)
from pypeit.dashboard.view.qa_dialog import QaImageDialog
from pypeit.dashboard.view.util import (text_on, detect_theme, clear_layout,
                                        fmt_value)


class ScienceView(QWidget):
    """
    The Science tab: a per-frame table + per-frame detail drill-down.

    Parameters
    ----------
    model : :class:`~pypeit.dashboard.model.DashboardModel`
        The reduction-state model.
    launcher : :class:`~pypeit.dashboard.launcher.Launcher`, optional
        Used to launch product viewers and (Re)Build runs.
    run_lock : :class:`~pypeit.dashboard.runlock.RunLock`, optional
        The single-run lock; the (Re)Build controls are disabled while it is
        locked.
    on_run_finished : callable, optional
        Called as ``on_run_finished(code)`` when a (Re)Build ends.
    parent : QWidget, optional
        The parent widget.
    """

    #: The science macro-steps that get a (Re)Build control.
    _REBUILD_STEPS = ('process', 'findobj', 'extract')

    #: The prerequisite step whose success enables each (Re)Build step
    #: (``process`` has none).
    _PREREQ = {'process': None, 'findobj': 'process', 'extract': 'findobj'}

    #: Columns of the per-frame science table, in display order.
    _COLUMNS = ('Frame', 'Calib', 'Detector', 'Type', 'process', 'findobj',
                'skysub', 'extract', 'nobj', 'spec2d', 'spec1d')

    #: The subset of :attr:`_COLUMNS` that are macro-step status cells.
    _STEP_COLS = ('process', 'findobj', 'skysub', 'extract')

    def __init__(self, model, launcher=None, run_lock=None,
                 on_run_finished=None, parent=None):
        super().__init__(parent=parent)
        self._model = model
        self._launcher = launcher
        self._run_lock = run_lock
        self._on_run_finished = on_run_finished
        self._theme = 'light'
        self._table = None
        self._detail = None
        self._selected = None           # (frame, det)
        self._rows = []                 # row index -> (frame, det)
        # step -> (button, available) for the current detail panel, so the
        # lock can restyle without a full rebuild.
        self._rebuild_buttons = {}
        # Column index of the per-object QA-count cell (set per detail
        # rebuild), and the frame QA-file list it cross-references.
        self._obj_qa_col = None
        self._qa_list = None
        # The view-level "Run PypeIt" button (full reduction).
        self._run_pypeit_button = None
        self._run_pypeit_available = False
        self._outer = QVBoxLayout(self)
        self._build()

    # -- public API ------------------------------------------------------

    def set_model(self, model):
        """
        Swap in a new model and rebuild from scratch.

        Parameters
        ----------
        model : :class:`~pypeit.dashboard.model.DashboardModel`
            The new model.
        """
        self._model = model
        self._selected = None
        clear_layout(self._outer)
        self._build()

    def refresh(self, model):
        """
        Swap in a new model but preserve the selected frame (used by the live
        monitor / completion refresh).

        Parameters
        ----------
        model : :class:`~pypeit.dashboard.model.DashboardModel`
            The new model.
        """
        prev = self._selected
        self._model = model
        self._selected = None
        clear_layout(self._outer)
        self._build()
        if prev is not None:
            self._select_pair(prev)

    def select_frame(self, frame, det):
        """
        Publicly select a ``(frame, det)`` row (used when the Status-view
        science navigator activates a frame).

        Parameters
        ----------
        frame : str
            The exposure basename.
        det : int or tuple
            Detector (int) or mosaic.
        """
        if self._table is not None:
            self._select_pair((frame, det))

    # -- build -----------------------------------------------------------

    def _build(self):
        """
        Build the view (edge message or the table + detail panel).
        """
        self._theme = detect_theme(self)
        self._rebuild_buttons = {}
        self._run_pypeit_button = None
        if not self._model.has_science():
            msg = QLabel('No science frames reduced yet.')
            msg.setAlignment(Qt.AlignCenter)
            msg.setStyleSheet('color: grey; font-style: italic; '
                              'font-size: 16px;')
            self._outer.addStretch(1)
            self._outer.addWidget(msg)
            self._outer.addStretch(1)
            return

        # The full-reduction "Run PypeIt" action, above the table.
        self._add_run_pypeit_button()
        self._build_table()

        self._detail = QVBoxLayout()
        container = QWidget()
        container.setLayout(self._detail)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(container)
        self._outer.addWidget(scroll, stretch=1)

        if len(self._rows) > 0:
            self._select_pair(self._rows[0])

    def _build_table(self):
        """
        Build the per-frame science table: flat ``(frame, det)`` rows with
        the four step statuses as color+glyph cells.
        """
        heading = QLabel('Science frames')
        heading.setStyleSheet('font-size: 15px; font-weight: bold;')
        self._outer.addWidget(heading)

        table = self._model.science_table()
        self._rows = []
        self._table = QTableWidget()
        self._table.setColumnCount(len(self._COLUMNS))
        self._table.setHorizontalHeaderLabels(list(self._COLUMNS))
        self._table.verticalHeader().setVisible(False)
        self._table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self._table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self._table.setSelectionMode(QAbstractItemView.SingleSelection)
        self._table.setAlternatingRowColors(True)
        self._table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeToContents)
        self._table.setRowCount(len(table))
        for r, row in enumerate(table):
            self._rows.append((row['frame'], row['detector']))
            self._set_table_row(r, row)
        self._table.setMinimumHeight(150)
        # Neutral selected-row fill, so a selected frame never inherits the
        # desktop theme's (possibly red) Highlight color.
        self._table.setStyleSheet(palette.selection_style(self._theme))
        # Selecting a row rebuilds the detail panel for that frame.
        self._table.itemSelectionChanged.connect(self._on_selection_changed)
        self._outer.addWidget(self._table, stretch=1)

    def _set_table_row(self, r, row):
        """
        Populate one science table row.

        Parameters
        ----------
        r : int
            Row index.
        row : `astropy.table.Row`_
            A row from
            :meth:`~pypeit.dashboard.model.DashboardModel.science_table`.
        """
        det_name = self._model.det_name(row['detector'])
        cells = {'Frame': str(row['frame']), 'Calib': str(row['calib']),
                 'Detector': det_name, 'Type': str(row['objtype']),
                 'nobj': str(row['nobj']), 'spec2d': str(row['spec2d']),
                 'spec1d': str(row['spec1d'])}
        for c, col in enumerate(self._COLUMNS):
            if col in self._STEP_COLS:
                status = row[col]
                style = palette.slit_style(status, theme=self._theme)
                item = QTableWidgetItem(f'{style.glyph} {style.label}')
                if style.category != palette.REQUIRED_UNDONE:
                    item.setForeground(QColor(style.color))
            else:
                item = QTableWidgetItem(cells[col])
            item.setTextAlignment(Qt.AlignCenter)
            self._table.setItem(r, c, item)

    # -- selection -------------------------------------------------------

    def _on_selection_changed(self):
        """
        Rebuild the detail panel for the selected row.
        """
        rows = self._table.selectionModel().selectedRows()
        if len(rows) == 0:
            return
        idx = rows[0].row()
        if 0 <= idx < len(self._rows):
            self._selected = self._rows[idx]
            self._build_detail(*self._selected)

    def _select_pair(self, pair):
        """
        Select the table row for ``(frame, det)`` (syncs the detail panel).

        Parameters
        ----------
        pair : tuple
            ``(frame, det)``.
        """
        if pair in self._rows:
            self._table.selectRow(self._rows.index(pair))
        elif len(self._rows) > 0:
            self._table.selectRow(0)

    # -- detail panel ----------------------------------------------------

    def _build_detail(self, frame, det):
        """
        Rebuild the detail panel for one science frame: actions + (Re)Build,
        per-slit table, per-object table, and the frame's QA files.

        Parameters
        ----------
        frame : str
            The exposure basename.
        det : int or tuple
            Detector (int) or mosaic.
        """
        if self._detail is None:
            return
        clear_layout(self._detail)
        self._rebuild_buttons = {}
        self._obj_qa_col = None
        self._qa_list = None
        entry = self._model.science_frame_entry(frame, det)

        heading = QLabel(f'{frame}  ({self._model.det_name(det)})')
        heading.setStyleSheet('font-size: 15px; font-weight: bold;')
        self._detail.addWidget(heading)

        self._add_action_row(frame, det, entry)

        # Per-slit table.
        self._add_slit_table(frame, det)
        # Per-object table (double-click a row to view its 1D spectrum; the
        # QA column opens/highlights that object's QA figures instead).
        self._add_object_table(frame, det, entry)
        # All of the frame's QA PNGs (per-object + frame-level), exposed
        # like the Calibrations QA list.
        self._add_qa_files(frame, det)

        self._detail.addStretch(1)

    def _add_action_row(self, frame, det, entry):
        """
        Add the View-2D button and the per-step (Re)Build controls.

        Parameters
        ----------
        frame : str
            The exposure basename.
        det : int or tuple
            Detector (int) or mosaic.
        entry : :class:`~pypeit.state.run_state.ScienceFrameState` or None
            The frame's state entry.
        """
        row = QHBoxLayout()

        # View spec2d (the whole 2D image).
        spec2d = getattr(entry, 'spec2d_file', None) if entry is not None \
            else None
        view2d = QPushButton('View spec2d')
        available2d = spec2d is not None and spec2d != '' \
            and self._product_path(spec2d).exists()
        view2d.setEnabled(available2d)
        if available2d:
            inspect = palette.inspect_color(self._theme)
            view2d.setStyleSheet(
                f'background-color: {inspect}; color: {text_on(inspect)}; '
                f'border-radius: 4px; padding: 4px 10px; font-weight: bold;')
            argv2d = spec2d_command(
                self._product_path(spec2d), det=det)
            # Launch the spec2d viewer; the lambda swallows Qt's "checked"
            # bool and binds the command argv.
            view2d.clicked.connect(
                lambda _c=False, a=argv2d:
                self._launch(a, 'view spec2d', hint='Ginga window'))
        row.addWidget(view2d)

        # (Re)Build controls: one per step, blue + lock-gated, enabled only
        # when the prerequisite step succeeded.
        for step in self._REBUILD_STEPS:
            row.addWidget(self._build_rebuild_button(frame, det, entry, step))
        row.addStretch(1)
        self._detail.addLayout(row)

    def _product_path(self, name):
        """
        Resolve a product filename to its path under ``Science/``.

        Parameters
        ----------
        name : str
            The product basename (or full path).

        Returns
        -------
        `Path`
            The resolved path.
        """
        p = Path(name)
        if p.is_absolute():
            return p
        return self._model.redux_dir / 'Science' / p.name

    def _build_rebuild_button(self, frame, det, entry, step):
        """
        Build a (Re)Build button for a science step: blue action when idle,
        orange + disabled while a run is active; enabled only when its
        prerequisite step succeeded, its calibrations are built, and a raw
        frame is known.

        Parameters
        ----------
        frame : str
            The exposure basename.
        det : int or tuple
            Detector (int) or mosaic.
        entry : :class:`~pypeit.state.run_state.ScienceFrameState` or None
            The frame's state entry.
        step : str
            ``process``, ``findobj``, or ``extract``.

        Returns
        -------
        QPushButton
            The (Re)Build button.
        """
        button = QPushButton(f'(Re)Build {step}')
        argv = science_run_command(self._model, frame, det, step)
        prereq = self._PREREQ[step]
        prereq_ok = prereq is None or (
            entry is not None and getattr(entry, prereq).status == 'success')
        # The science (Re)Build is disabled until the frame's calibrations
        # are built successfully: a science step can only run on finished
        # calibrations.
        calib_id = getattr(entry, 'calib_id', None) if entry is not None \
            else None
        calib_ready = self._model.calibrations_ready(calib_id, det)
        available = argv is not None and prereq_ok and calib_ready
        if available:
            tip = 'Re-run this science step (pypeit_reduce_by_step) for ' \
                  'this frame.'
        elif not calib_ready:
            tip = 'Unavailable: the calibrations for this group are not yet ' \
                  'built successfully.'
        else:
            tip = 'Unavailable: needs the prerequisite step and a known raw ' \
                  'frame.'
        button.setToolTip(tip)
        self._rebuild_buttons[step] = (button, available)
        if available:
            # Launch the confirmation + (re)build for this step; the lambda
            # swallows Qt's "checked" bool and binds the frame/det/step.
            button.clicked.connect(
                lambda _c=False, f=frame, d=det, s=step:
                self._on_rebuild(f, d, s))
        locked = self._run_lock is not None and self._run_lock.is_locked()
        self._style_rebuild_button(step, locked)
        return button

    def _style_rebuild_button(self, step, locked):
        """
        Style one (Re)Build button for the lock state (blue idle; orange +
        disabled while a run is active; plain-disabled if unavailable).

        Parameters
        ----------
        step : str
            The step whose button to style.
        locked : bool
            Whether a run is in progress.
        """
        button, available = self._rebuild_buttons.get(step, (None, False))
        if button is None:
            return
        if not available:
            button.setEnabled(False)
            button.setStyleSheet('')
            return
        if locked:
            color = palette.LIGHT_COLORS[palette.RUNNING] \
                if self._theme == 'light' \
                else palette.DARK_COLORS[palette.RUNNING]
            button.setText('⏳ Run in progress')
            button.setEnabled(False)
        else:
            color = palette.action_color(self._theme)
            button.setText(f'(Re)Build {step}')
            button.setEnabled(True)
        button.setStyleSheet(
            f'background-color: {color}; color: {text_on(color)}; '
            f'border-radius: 4px; padding: 4px 10px; font-weight: bold;')

    def set_locked(self, locked):
        """
        Restyle the (Re)Build controls and the "Run PypeIt" button when the
        run lock changes.

        Parameters
        ----------
        locked : bool
            Whether a run is in progress.
        """
        for step in list(self._rebuild_buttons):
            self._style_rebuild_button(step, locked)
        self._style_run_pypeit_button(locked)

    # -- "Run PypeIt" (full reduction) -----------------------------------

    def _add_run_pypeit_button(self):
        """
        Add the view-level **"Run PypeIt"** button: launch the full reduction
        (``run_pypeit -o``) over all science/standard frames, with an
        overwrite warning, governed by the same single-run lock.
        """
        row = QHBoxLayout()
        button = QPushButton('Run PypeIt')
        self._run_pypeit_button = button
        argv = run_pypeit_command(self._model)
        self._run_pypeit_available = argv is not None
        button.setToolTip(
            'Run the full reduction (run_pypeit -o): process all '
            'science/standard frames, overwriting existing outputs.'
            if self._run_pypeit_available else
            'Unavailable: no .pypeit file.')
        if self._run_pypeit_available:
            # Launch the confirmation + full reduction; the lambda swallows
            # Qt's "checked" bool.
            button.clicked.connect(lambda _c=False: self._on_run_pypeit())
        locked = self._run_lock is not None and self._run_lock.is_locked()
        self._style_run_pypeit_button(locked)
        row.addWidget(button)
        row.addStretch(1)
        self._outer.addLayout(row)

    def _style_run_pypeit_button(self, locked):
        """
        Style the "Run PypeIt" button for the lock state: a distinct
        **indigo** action when idle (set apart from the per-step blue
        (Re)Build), orange + disabled "⏳ Run in progress" while any run is
        active.

        Parameters
        ----------
        locked : bool
            Whether a run is in progress.
        """
        button = self._run_pypeit_button
        if button is None:
            return
        if not self._run_pypeit_available:
            button.setEnabled(False)
            button.setStyleSheet('')
            return
        if locked:
            color = palette.LIGHT_COLORS[palette.RUNNING] \
                if self._theme == 'light' \
                else palette.DARK_COLORS[palette.RUNNING]
            button.setText('⏳ Run in progress')
            button.setEnabled(False)
        else:
            # A distinct indigo so the global run is not mistaken for a
            # per-step (Re)Build (which is blue).
            color = '#4527A0' if self._theme == 'light' else '#7E57C2'
            button.setText('Run PypeIt')
            button.setEnabled(True)
        button.setStyleSheet(
            f'background-color: {color}; color: {text_on(color)}; '
            f'border-radius: 4px; padding: 4px 14px; font-weight: bold;')

    def _on_run_pypeit(self):
        """
        Confirm (overwrite warning) and launch the full reduction
        (``run_pypeit -o``) via the launcher.
        """
        if self._run_lock is not None and self._run_lock.is_locked():
            return
        argv = run_pypeit_command(self._model)
        if argv is None:
            return
        box = QMessageBox(self)
        box.setIcon(QMessageBox.Warning)
        box.setWindowTitle('Confirm Run PypeIt')
        box.setText('Run the full PypeIt reduction (run_pypeit -o)?\n\n'
                    'This processes ALL science and standard frames and '
                    '**overwrites** the existing science outputs (the '
                    'spec2d/spec1d files and reduction products).')
        box.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
        box.setDefaultButton(QMessageBox.Cancel)
        if box.exec_() != QMessageBox.Ok:
            return
        if self._launcher is not None:
            self._launcher.run(argv, description='Running PypeIt (full)',
                               on_finished=self._on_run_finished)

    def _on_rebuild(self, frame, det, step):
        """
        Confirm and launch a science (Re)Build (re-run with a confirmation;
        no move-aside of the existing products).

        Parameters
        ----------
        frame : str
            The exposure basename.
        det : int or tuple
            Detector (int) or mosaic.
        step : str
            The science step to (re)build.
        """
        if self._run_lock is not None and self._run_lock.is_locked():
            return
        argv = science_run_command(self._model, frame, det, step)
        if argv is None:
            return
        box = QMessageBox(self)
        box.setIcon(QMessageBox.Warning)
        box.setWindowTitle('Confirm (Re)Build')
        box.setText(f'(Re)build the "{step}" step for science frame\n'
                    f'{frame} ({self._model.det_name(det)})?\n\n'
                    f'This re-runs pypeit_reduce_by_step and overwrites that '
                    f'step\'s products (the Intermediate/ files, and the '
                    f'spec2d/spec1d for "extract").')
        box.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
        box.setDefaultButton(QMessageBox.Cancel)
        if box.exec_() != QMessageBox.Ok:
            return
        if self._launcher is not None:
            self._launcher.run(argv, description=f'(Re)building {step}',
                               on_finished=self._on_run_finished)

    # -- per-slit / per-object tables ------------------------------------

    def _add_slit_table(self, frame, det):
        """
        Add the per-slit science table (status + nobj).

        Parameters
        ----------
        frame : str
            The exposure basename.
        det : int or tuple
            Detector (int) or mosaic.
        """
        rows = self._model.science_slit_table(frame, det)
        if rows is None or len(rows) == 0:
            return
        self._detail.addWidget(QLabel(f'Per-slit ({len(rows)}):'))
        table = QTableWidget()
        cols = ['Slit', 'Status', 'nobj']
        table.setColumnCount(len(cols))
        table.setHorizontalHeaderLabels(cols)
        table.verticalHeader().setVisible(False)
        table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        table.setSelectionMode(QAbstractItemView.NoSelection)
        table.setAlternatingRowColors(True)
        table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        table.setRowCount(len(rows))
        for r, row in enumerate(rows):
            slit_item = QTableWidgetItem(str(row['slit']))
            slit_item.setTextAlignment(Qt.AlignCenter)
            table.setItem(r, 0, slit_item)
            style = palette.slit_style(row['status'], theme=self._theme)
            st_item = QTableWidgetItem(f'{style.glyph} {style.label}')
            st_item.setTextAlignment(Qt.AlignCenter)
            if style.category != palette.REQUIRED_UNDONE:
                st_item.setForeground(QColor(style.color))
            table.setItem(r, 1, st_item)
            nobj_item = QTableWidgetItem(
                '—' if row['nobj'] is None else str(row['nobj']))
            nobj_item.setTextAlignment(Qt.AlignCenter)
            table.setItem(r, 2, nobj_item)
        table.setMinimumHeight(110)
        self._detail.addWidget(table)

    def _add_object_table(self, frame, det, entry):
        """
        Add the per-object science table; double-click a row to view that
        object's 1D spectrum.

        Parameters
        ----------
        frame : str
            The exposure basename.
        det : int or tuple
            Detector (int) or mosaic.
        entry : :class:`~pypeit.state.run_state.ScienceFrameState` or None
            The frame's state entry.
        """
        rows = self._model.science_object_table(frame, det)
        if rows is None or len(rows) == 0:
            return
        self._detail.addWidget(QLabel(
            f'Objects ({len(rows)}) — double-click a metric cell for the 1D '
            f'spectrum, or the QA cell for the object\'s QA figures:'))
        # The QA column holds the count of this object's QA figures (found
        # by the QA naming convention — no fixed list of QA types, so new
        # figures are picked up automatically); the metric columns open the
        # 1D spectrum.
        cols = ['ObjID', 'Slit', 'Spat', 'FWHM', 'snr_find', 's2n', 'sign',
                'extracted', 'QA']
        keys = ['objid', 'slitid', 'spat_pixpos', 'fwhm', 'snr_find', 's2n',
                'sign', 'extracted']
        self._obj_qa_col = len(keys)
        table = QTableWidget()
        table.setColumnCount(len(cols))
        table.setHorizontalHeaderLabels(cols)
        table.verticalHeader().setVisible(False)
        table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        table.setSelectionBehavior(QAbstractItemView.SelectRows)
        table.setAlternatingRowColors(True)
        table.setStyleSheet(palette.selection_style(self._theme))
        table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        table.setRowCount(len(rows))
        spec1d = getattr(entry, 'spec1d_file', None) if entry is not None \
            else None
        det_name = self._model.det_name(det)
        for r, row in enumerate(rows):
            for c, key in enumerate(keys):
                cell = QTableWidgetItem(fmt_value(row.get(key)))
                cell.setTextAlignment(Qt.AlignCenter)
                # Stash the object name (for --obj) on the first column.
                if c == 0 and spec1d is not None:
                    name = science_object_name(
                        row.get('spat_pixpos'), row.get('slitid'), det_name)
                    cell.setData(Qt.UserRole, name)
                table.setItem(r, c, cell)
            # The QA cell: the count of QA figures discovered for this
            # object; the paths are stashed for the double-click handler.
            qa_paths = self._model.science_object_qa_files(
                frame, det, row.get('slitid'))
            cell = QTableWidgetItem(f'🔍 {len(qa_paths)}' if qa_paths
                                    else '—')
            cell.setTextAlignment(Qt.AlignCenter)
            if qa_paths:
                cell.setData(Qt.UserRole, [str(p) for p in qa_paths])
                cell.setForeground(
                    QColor(palette.inspect_color(self._theme)))
            table.setItem(r, self._obj_qa_col, cell)
        table.setMinimumHeight(120)
        # Route a double-click to the QA PNG or the 1D-spectrum viewer,
        # binding the frame's spec1d product.
        table.itemDoubleClicked.connect(
            lambda item, s=spec1d: self._on_object_double_clicked(item, s))
        self._detail.addWidget(table)

    def _on_object_double_clicked(self, item, spec1d):
        """
        Handle a double-click in the per-object table: the QA cell opens
        that object's QA figure (or, if it has several, highlights them in
        the frame's QA-file list); any other cell opens the object's 1D
        spectrum (``pypeit_show_1dspec --obj``).

        Parameters
        ----------
        item : QTableWidgetItem
            The clicked cell.
        spec1d : str or None
            The frame's ``spec1d`` product.
        """
        # The QA cell stashes the object's QA PNG paths → open/highlight.
        if item.column() == self._obj_qa_col:
            paths = item.data(Qt.UserRole)
            if paths is None or len(paths) == 0:
                return
            if len(paths) == 1:
                QaImageDialog(paths[0], parent=self).show()
            else:
                # Several figures: highlight them in the frame QA-file list
                # below (which opens each on double-click) rather than
                # opening a stack of windows.
                self._highlight_qa_files(paths)
            return
        # Otherwise view the object's 1D spectrum (name stashed on column 0).
        if spec1d is None or spec1d == '':
            return
        row = item.row()
        name_item = item.tableWidget().item(row, 0)
        obj_name = name_item.data(Qt.UserRole) if name_item is not None \
            else None
        argv = spec1d_command(self._product_path(spec1d),
                                           obj_name=obj_name)
        self._launch(argv, 'view 1D spectrum', hint='plot window')

    def _add_qa_files(self, frame, det):
        """
        Add a QA-file list for the frame (all of its
        ``obj_prof``/``obj_trace`` and ``spec_flex_*`` PNGs), mirroring the
        Calibrations QA list; double-click opens the PNG full-view.

        Parameters
        ----------
        frame : str
            The exposure basename.
        det : int or tuple
            Detector (int) or mosaic.
        """
        qa_files = self._model.science_qa_files(frame, det)
        if qa_files is None or len(qa_files) == 0:
            return
        self._detail.addWidget(
            QLabel(f'QA files ({len(qa_files)}) — double-click to open:'))
        listw = QListWidget()
        listw.setMaximumHeight(130)
        # Multi-select so the per-object QA cell can highlight all of one
        # object's figures at once.
        listw.setSelectionMode(QAbstractItemView.ExtendedSelection)
        # Strip the (long) frame prefix so the list reads compactly.
        prefix = f'{frame}_'
        for path in qa_files:
            label = path.name
            if label.startswith(prefix):
                label = label[len(prefix):]
            item = QListWidgetItem(label)
            item.setData(Qt.UserRole, str(path))
            listw.addItem(item)
        # Open the double-clicked QA PNG in the full-view image dialog.
        listw.itemDoubleClicked.connect(
            lambda it: QaImageDialog(it.data(Qt.UserRole), parent=self).show())
        self._detail.addWidget(listw)
        # Kept so the per-object QA cell can highlight an object's figures.
        self._qa_list = listw

    def _highlight_qa_files(self, paths):
        """
        Select (and scroll to) the given QA files in the frame's QA-file
        list, so the user can open each with a double-click.

        Parameters
        ----------
        paths : list
            The QA PNG paths (str) to highlight.
        """
        if self._qa_list is None:
            # No QA list built (should not happen if paths exist); fall
            # back to opening the first figure directly.
            QaImageDialog(paths[0], parent=self).show()
            return
        wanted = set(paths)
        first = None
        for i in range(self._qa_list.count()):
            item = self._qa_list.item(i)
            selected = item.data(Qt.UserRole) in wanted
            item.setSelected(selected)
            if selected and first is None:
                first = item
        if first is not None:
            self._qa_list.scrollToItem(first)
            self._qa_list.setFocus()

    def _launch(self, argv, description, hint='viewer window'):
        """
        Launch a viewer command via the launcher (no-op if unavailable).

        Parameters
        ----------
        argv : list or None
            The command argv.
        description : str
            Human description for the activity bar.
        hint : str, optional
            Where the result appears.
        """
        if argv is not None and self._launcher is not None:
            self._launcher.launch(argv, description=description, hint=hint)
