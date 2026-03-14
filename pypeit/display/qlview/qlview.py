"""
Ginga local plugin that provides an interactive quicklook reduction interface.

Overview
--------
``QLView`` is the top-level Ginga ``LocalPlugin`` that orchestrates the
PypeIt quicklook pipeline from within the Ginga image viewer.  It is
split into several collaborating components

``QLView`` (this module)
    The plugin class itself.  Owns all mutable GUI state, wires callback
    methods, and coordinates the other components.  Inherits from
    ``ginga.GingaPlugin.LocalPlugin`` and therefore follows the Ginga
    plugin lifecycle: ``build_gui`` → ``start`` → [user interaction] →
    ``stop`` / ``close``.

``QLViewState`` (``.state``)
    A lightweight dataclass that groups the current viewer state:
    the active raw filepath, reduced-calibrations filepath, reduction
    output path, loaded ``SlitTraceSet`` objects, slit polygon dict,
    and the active slit key.  Passed between methods to avoid scattering
    mutable state across ``self``.

``QLViewUI`` (``.ui``)
    Builds the entire Ginga widget tree and attaches it to the plugin
    container.  Stores every widget reference on ``self.plugin`` (i.e.,
    on the ``QLView`` instance) so that callback methods can reach them
    without needing to navigate the widget hierarchy.  Keeps UI
    construction completely separate from business logic.

``InstrumentRegistry`` / ``Instrument`` subclasses (``.instruments``)
    A registry of supported instruments.  Each ``Instrument`` knows how to read display-ready
    raw image data (``get_display_image``), extract FITS header metadata
    for the file-browser tree columns (``get_raw_info`` /
    ``get_reduced_info``), and match raw frames to their best calibration
    directory (``recommend_calibrations``).  Swapping instruments at
    runtime rebuilds the tree-view columns via
    ``_rebuild_treeview_columns``.

``FileBrowserController`` (``.file_browser``)
    Translates a directory path and an ``Instrument`` into a Ginga
    tree-view listing dict.  Delegates all filesystem access to the
    injected ``FileBrowserBackend`` so that the same controller works
    against a local disk or a remote server without changes.

``FileBrowserBackend`` / ``ReductionBackend`` (``.backends``)
    Protocol-based backend abstractions with two concrete pairs:

    - ``LocalFileBrowserBackend`` + ``LocalReductionBackend`` — operate
      on the local filesystem and run ``pypeit_ql`` in-process on a
      daemon thread.
    - ``RemoteFileBrowserBackend`` + ``RemoteReductionBackend`` — delegate
      all operations to an HTTP server via ``requests``, enabling use
      cases where the raw data lives on a remote instrument workstation.

    The active backends are selected at startup from ``~/.quicklook.cfg``
    and can be changed at runtime through the Settings dialog.

``SlitOverlay`` (``.slit_overlay``)
    Manages the Ginga ``DrawingCanvas`` layer that renders slit polygons
    and optional label text over the raw image.  Maintains a dict of
    ``slit_key`` → ``Polygon`` and provides ``activate`` / ``deactivate``
    methods to highlight the currently selected slit in blue.

Data flow
---------
1. **File browsing** — ``_browse_and_update`` calls
   ``FileBrowserController.browse``, which uses the active
   ``FileBrowserBackend`` to list a directory and reads per-file header
   metadata via ``Instrument.get_raw_info`` / ``get_reduced_info``.  The
   resulting dict is pushed directly into the Ginga ``TreeView`` widget.

2. **Raw image display** — double-clicking a FITS file calls
   ``open_raw_file``, which delegates mosaic assembly to
   ``Instrument.get_display_image`` and loads the result into the Ginga
   ``AstroImage`` canvas.

3. **Calibration suggestion** — after a raw file is opened,
   ``_suggest_calibrations`` calls ``Instrument.recommend_calibrations``
   on a background thread, then highlights the best-matching calibration
   directory in the reduced tree on the GUI thread via ``fv.gui_do``.

4. **Slit rendering** — ``render_slits_cb`` reads ``Slits_*.fits*`` from
   the calibration ``Calibrations/`` directory, loads them as
   ``SlitTraceSet`` objects, and hands them to ``SlitOverlay.build`` to
   draw polygons on the ``slit_canvas`` ``DrawingCanvas`` layer. Similarly,
   wavelength information is pulled from "Wave*.fits" files when requested,
   and overlaid on the raw frame using the SlitWavelength plugin.

5. **Reduction** — ``reduce_slit_cb`` assembles the ``pypeit_ql``
   argument list and calls ``ReductionBackend.submit`` on a daemon thread.
   A Ginga timer (``_register_reduction_timer``) polls for output files
   every ``reduction_cadence`` seconds on the GUI thread, enabling and
   wiring the "Show" button when the ``spec1d`` file appears.

6. **Trace overlay** — once reduction is complete, ``show_traces_cb``
   loads the ``SpecObjs`` from the ``spec1d`` file on a background thread
   and draws per-object ``TRACE_SPAT`` paths in orange on a per-slit
   ``DrawingCanvas``.  A polling timer watches the paired ``Spec1dView``
   plugin for selection changes and recolors the active trace cyan.

Configuration
-------------
``~/.quicklook.cfg`` (INI format, ``[DEFAULT]`` section) controls startup
paths, file-filter defaults, backend selection, poll cadence, and
reduction timeout.  Path values support ``strftime``-style format codes
(e.g. ``raw_path_template = /data/raw/%Y%m%d``) that are expanded at
startup.  The "Save Default Config" button writes a template for the
user to edit. 
"""

from __future__ import annotations

import configparser
import datetime
import glob
import os
import re
import threading
import time
from pathlib import Path
from typing import Dict, Optional

import numpy as np
from astropy.io import fits
from ginga import GingaPlugin
from ginga.AstroImage import AstroImage
from ginga.canvas.types.layer import DrawingCanvas
from ginga.qtw.QtHelp import QtCore, QtGui

from pypeit.slittrace import SlitTraceSet

from .backends import (
    FileBrowserBackend,
    LocalFileBrowserBackend,
    LocalReductionBackend,
    ReductionBackend,
    RemoteFileBrowserBackend,
    RemoteReductionBackend,
)
from .file_browser import FileBrowserController
from .instruments import InstrumentRegistry
from .slit_overlay import SlitOverlay
from .state import QLViewState
from .ui import QLViewUI


class QLView(GingaPlugin.LocalPlugin):
    def __init__(self, fv, fitsimage):
        """Initialise the QLView plugin and set up all internal state.

        Called by the Ginga plugin machinery when the plugin is first loaded.
        Sets default values for file-filter flags, reduction timers, trace
        overlay state, and calibration-configuration tracking.  Icon assets
        are fetched from the Ginga icon directory and injected into
        :class:`~.file_browser.FileBrowserController`.

        Parameters
        ----------
        fv : ginga.GingaPlugin.GingaPlugin
            The top-level Ginga application object.
        fitsimage : ginga.ImageViewCanvas.ImageViewCanvas
            The FITS image viewer canvas this plugin is attached to.
        """
        super().__init__(fv, fitsimage)

        keywords = [("Object", "OBJECT"), ("Date", "DATE-OBS"), ("Time UT", "UT")]

        prefs = self.fv.get_preferences()
        self.settings = prefs.create_category("plugin_QLView")
        self.settings.add_defaults(
            scan_fits_headers=False,
            scan_limit=100,
            keywords=keywords,
            color_alternate_rows=True,
            max_rows_for_col_resize=5000,
        )
        self.settings.load(onError="silent")

        self.state = QLViewState()
        self.gui_up = False
        self.instrument_registry = InstrumentRegistry(self.logger)
        self.instrument = self.instrument_registry.create("DEIMOS")

        self.file_backend: FileBrowserBackend = LocalFileBrowserBackend()
        self.reduction_backend: ReductionBackend = LocalReductionBackend()

        self.backend_mode = "local"
        self.remote_host = ""
        self.remote_port = ""
        self.remote_api_key = ""

        # File-filter state (shown/edited in the Settings dialog)
        self.raw_filter_fits = True
        self.raw_filter_nonfits = False
        self.raw_filter_dirs = True
        self.reduced_filter_fits = False
        self.reduced_filter_nonfits = False
        self.reduced_filter_dirs = True

        self._raw_name_col_idx = 0
        self._reduced_name_col_idx = 0
        self._compute_name_col_indices()

        self.dc = fitsimage.get_canvas().get_draw_classes()
        self.slit_canvas: Optional[DrawingCanvas] = None
        self.overlay = SlitOverlay(self.dc)
        self.file_browser = FileBrowserController(self.logger, self.settings, self.file_backend)
        self.reduction_timers: Dict[str, object] = {}
        self.reduction_start_times: Dict[str, float] = {}
        self.reduction_timeout: float = 300.0  # 5 minutes
        self.reduction_cadence: float = 5.0
        self.reduction_control_elements: Dict[str, dict] = {}
        self._slit_combo_keys: list = []  # slit keys in combo box order

        # Per-reduction trace overlay state
        self._trace_canvases: Dict[str, DrawingCanvas] = {}
        self._trace_paths: Dict[str, Dict[int, dict]] = {}
        self._trace_timers: Dict[str, object] = {}
        self._trace_last_exten: Dict[str, Optional[int]] = {}

        # Configuration dict read from the currently rendered calibration dir.
        # Keys are PypeIt configuration_key names (e.g. "filter1", "decker");
        # values are strings.  Used to detect setup changes when a new raw
        # image is selected.
        self._rendered_cal_config: Dict[str, str] = {}

        # Suppress tree activation events briefly after settings dialog closes
        self._suppress_activate: bool = False

        # Manual extraction state
        self.manual_extract_mode: bool = False
        self._manual_x: Optional[float] = None
        self._manual_y: Optional[float] = None
        self._manual_det_label: Optional[str] = None
        self._manual_spat_det: Optional[float] = None
        self._manual_slit_key: Optional[str] = None
        self.manual_extract_canvas: Optional[DrawingCanvas] = None

        icondir = self.fv.iconpath
        self.file_browser.folderpb = self.fv.get_icon(icondir, "folder.svg")
        self.file_browser.filepb = self.fv.get_icon(icondir, "file.svg")
        self.file_browser.fitspb = self.fv.get_icon(icondir, "fits.svg")

        self.ui = QLViewUI(self)

    # --- Column management ---

    def _compute_name_col_indices(self) -> None:
        """Recompute the name-column indices for the current instrument's column defs.
        Used to know where to pull a filename from when the tree is rebuilt"""
        self._raw_name_col_idx = 0
        for idx, (_col, attr) in enumerate(self.instrument.columns["raw"]):
            if attr == "name":
                self._raw_name_col_idx = idx
                break

        self._reduced_name_col_idx = 0
        for idx, (_col, attr) in enumerate(self.instrument.columns["reduced"]):
            if attr == "name":
                self._reduced_name_col_idx = idx
                break

    def _rebuild_treeview_columns(self) -> None:
        """Re-setup treeview column headers for the current instrument, then re-browse."""
        if not self.gui_up:
            return
        self._compute_name_col_indices()
        raw_cols = self.instrument.columns["raw"]
        reduced_cols = self.instrument.columns["reduced"]
        self.raw_treeview.setup_table(raw_cols, 1, "name")
        self.reduced_treeview.setup_table(reduced_cols, 1, "name")
        # Re-browse both trees so data rows match the new columns
        raw_base = self._get_tree_base_dir(self.state.raw_filepath)
        if raw_base and os.path.isdir(raw_base):
            self._browse_and_update(raw_base, which_tree="raw")
        reduced_base = self._get_tree_base_dir(self.state.reduced_filepath)
        if reduced_base and os.path.isdir(reduced_base):
            self._browse_and_update(reduced_base, which_tree="reduced")

    def build_gui(self, container):
        """Build the plugin GUI and attach it to *container*.

        Delegates widget construction to :class:`~.ui.QLViewUI`, then marks
        the GUI as ready so that subsequent callbacks that guard on
        ``self.gui_up`` can proceed.

        Parameters
        ----------
        container : ginga.gw.Widgets.Box
            The Ginga container widget provided by the plugin framework.
        """
        self.ui.build(container)
        self.gui_up = True

    # --- Callbacks ---

    def create_config_cb(self, w):
        """Write the current settings to ``~/.quicklook.cfg``.

        Serialises the active raw path, reduced calibrations path, redux
        output path, file-filter flags, and reduction timeout to a
        :mod:`configparser` INI file.  Template keys (``*_path_template``)
        are intentionally omitted so that the saved file always uses plain
        absolute paths; template expansion occurs only when reading the file
        at startup.

        Parameters
        ----------
        w : ginga.gw.Widgets.Button
            The "Save Default Config" button widget (unused).
        """
        config_file = Path.home() / ".quicklook.cfg"
        config = configparser.ConfigParser()

        raw_path = self.raw_text_entry.get_text()
        if os.path.isfile(raw_path):
            raw_path = str(Path(raw_path).parent)

        reduced_path = self.reduced_text_entry.get_text()
        if re.match(r'^[a-z][a-z0-9_]+_[A-Z]$', Path(reduced_path).name):
            reduced_path = str(Path(reduced_path).parent)

        config["DEFAULT"] = {
            "redux_path": self.state.redux_path,
            "raw_path": raw_path,
            "reduced_path": reduced_path,
            "raw_show_fits": str(self.raw_filter_fits),
            "raw_show_nonfits": str(self.raw_filter_nonfits),
            "raw_show_dirs": str(self.raw_filter_dirs),
            "reduced_show_fits": str(self.reduced_filter_fits),
            "reduced_show_nonfits": str(self.reduced_filter_nonfits),
            "reduced_show_dirs": str(self.reduced_filter_dirs),
            "reduction_timeout": str(self.reduction_timeout),
        }
        with open(config_file, "w") as f:
            config.write(f)
        self.logger.info(f"Saved default config to {config_file}")

    def open_settings_dialog(self) -> None:
        """Open the Settings dialog.

        Contains backend configuration, reduction path, poll cadence, and
        per-tree file-filter toggles.  Changes are applied on OK, *not* in real
        time.
        """
        dialog = QtGui.QDialog()
        dialog.setWindowTitle("Settings")
        dialog.setMinimumWidth(420)

        outer = QtGui.QVBoxLayout()

        # ── Backend ──────────────────────────────────────────────────────────
        backend_group = QtGui.QGroupBox("Backend")
        backend_form = QtGui.QFormLayout()
        local_checkbox = QtGui.QCheckBox("Use local backend")
        local_checkbox.setChecked(self.backend_mode == "local")
        backend_form.addRow(local_checkbox)
        host_edit = QtGui.QLineEdit(self.remote_host)
        port_edit = QtGui.QLineEdit(self.remote_port)
        key_edit = QtGui.QLineEdit(self.remote_api_key)
        key_edit.setPlaceholderText("Leave blank if server has no --api-key set")
        backend_form.addRow("Host:", host_edit)
        backend_form.addRow("Port:", port_edit)
        backend_form.addRow("API Key:", key_edit)
        backend_group.setLayout(backend_form)
        outer.addWidget(backend_group)

        def _toggle_backend_fields(checked: bool) -> None:
            """Set all of the backend stuff in one go"""
            host_edit.setEnabled(not checked)
            port_edit.setEnabled(not checked)
            key_edit.setEnabled(not checked)

        _toggle_backend_fields(local_checkbox.isChecked())
        local_checkbox.toggled.connect(_toggle_backend_fields)

        # ── Reduction ─────────────────────────────────────────────────────────
        reduction_group = QtGui.QGroupBox("Reduction")
        reduction_form = QtGui.QFormLayout()
        redux_path_edit = QtGui.QLineEdit(self.state.redux_path)
        redux_path_edit.setPlaceholderText("Path to the reduction output directory")
        reduction_form.addRow("Reduction Path:", redux_path_edit)
        cadence_edit = QtGui.QLineEdit(str(self.reduction_cadence))
        cadence_edit.setPlaceholderText("Seconds between completion checks")
        reduction_form.addRow("Poll Cadence (s):", cadence_edit)
        timeout_edit = QtGui.QLineEdit(str(self.reduction_timeout))
        timeout_edit.setPlaceholderText("Seconds before declaring a reduction timed out")
        reduction_form.addRow("Reduction Timeout (s):", timeout_edit)
        reduction_group.setLayout(reduction_form)
        outer.addWidget(reduction_group)

        # ── File Filters ──────────────────────────────────────────────────────
        filters_group = QtGui.QGroupBox("File Filters")
        filters_layout = QtGui.QVBoxLayout()

        raw_label = QtGui.QLabel("Raw tree:")
        raw_fits_cb = QtGui.QCheckBox("FITS Files")
        raw_fits_cb.setChecked(self.raw_filter_fits)
        raw_nonfits_cb = QtGui.QCheckBox("Non-FITS Files")
        raw_nonfits_cb.setChecked(self.raw_filter_nonfits)
        raw_dirs_cb = QtGui.QCheckBox("Directories")
        raw_dirs_cb.setChecked(self.raw_filter_dirs)

        reduced_label = QtGui.QLabel("Reduced tree:")
        red_fits_cb = QtGui.QCheckBox("FITS Files")
        red_fits_cb.setChecked(self.reduced_filter_fits)
        red_nonfits_cb = QtGui.QCheckBox("Non-FITS Files")
        red_nonfits_cb.setChecked(self.reduced_filter_nonfits)
        red_dirs_cb = QtGui.QCheckBox("Directories")
        red_dirs_cb.setChecked(self.reduced_filter_dirs)

        filters_layout.addWidget(raw_label)
        for cb in (raw_fits_cb, raw_nonfits_cb, raw_dirs_cb):
            filters_layout.addWidget(cb)
        filters_layout.addWidget(reduced_label)
        for cb in (red_fits_cb, red_nonfits_cb, red_dirs_cb):
            filters_layout.addWidget(cb)
        filters_group.setLayout(filters_layout)
        outer.addWidget(filters_group)

        # ── Dialog buttons ───────────────────────────────────────────────────
        buttons = QtGui.QDialogButtonBox(
            QtGui.QDialogButtonBox.Ok | QtGui.QDialogButtonBox.Cancel
        )
        outer.addWidget(buttons)
        dialog.setLayout(outer)
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)

        if dialog.exec_() != QtGui.QDialog.Accepted:
            return

        # ── Apply backend ────────────────────────────────────────────────────
        if local_checkbox.isChecked():
            self.backend_mode = "local"
            self.file_backend = LocalFileBrowserBackend()
            self.reduction_backend = LocalReductionBackend()
            self.remote_host = ""
            self.remote_port = ""
            self.remote_api_key = ""
        else:
            host = host_edit.text().strip()
            port = port_edit.text().strip()
            api_key = key_edit.text().strip() or None
            base_url = f"http://{host}:{port}"

            from pypeit.display.qlview.backends import _requests
            if _requests is None:
                QtGui.QMessageBox.critical(
                    None,
                    "Missing Dependency",
                    "The 'requests' package is required for remote backends.\n\n"
                    "Install it with:  pip install requests",
                )
                return
            try:
                resp = _requests.get(f"{base_url}/api/health", timeout=5)
                resp.raise_for_status()
            except Exception as exc:
                QtGui.QMessageBox.critical(
                    None,
                    "Backend Unreachable",
                    f"Could not connect to {base_url}:\n\n{exc}",
                )
                return

            self.backend_mode = "remote"
            self.remote_host = host
            self.remote_port = port
            self.remote_api_key = api_key or ""
            self.file_backend = RemoteFileBrowserBackend(base_url, api_key=api_key)
            self.reduction_backend = RemoteReductionBackend(base_url, api_key=api_key)

        self.file_browser.backend = self.file_backend

        # ── Apply reduction settings ─────────────────────────────────────────
        self.state.redux_path = redux_path_edit.text().strip()
        try:
            val = float(cadence_edit.text())
            if val > 0:
                self.reduction_cadence = val
        except ValueError:
            pass
        try:
            val = float(timeout_edit.text())
            if val > 0:
                self.reduction_timeout = val
        except ValueError:
            pass

        # ── Apply file filters ────────────────────────────────────────────────
        self.raw_filter_fits = raw_fits_cb.isChecked()
        self.raw_filter_nonfits = raw_nonfits_cb.isChecked()
        self.raw_filter_dirs = raw_dirs_cb.isChecked()
        self.reduced_filter_fits = red_fits_cb.isChecked()
        self.reduced_filter_nonfits = red_nonfits_cb.isChecked()
        self.reduced_filter_dirs = red_dirs_cb.isChecked()

        # Use the parent of whatever state path is set — this is the directory
        # whose contents the tree is currently showing (whether filepath points
        # to the browse root glob, a selected subdir, or a selected file).
        def _filter_base(filepath):
            """Return the parent directory of *filepath* for use as a filter base.

            Uses the *parent* of the stored filepath rather than the filepath
            itself so that the filter base is the directory whose *contents*
            are currently shown in the tree — whether the state path points to
            a browse-root glob (``/dir/*``), a selected subdirectory, or a
            selected file.

            Parameters
            ----------
            filepath : str or None
                Value of ``state.raw_filepath`` or ``state.reduced_filepath``.

            Returns
            -------
            str or None
                Absolute path of the containing directory, or ``None`` when
                *filepath* is empty.
            """
            if not filepath:
                return None
            return os.path.abspath(os.path.join(filepath, os.pardir))

        raw_base = _filter_base(self.state.raw_filepath)
        if raw_base:
            self._apply_file_filter(
                self.raw_treeview, raw_base,
                self.raw_filter_fits, self.raw_filter_nonfits, self.raw_filter_dirs,
                name_col_idx=self._raw_name_col_idx,
            )
        reduced_base = _filter_base(self.state.reduced_filepath)
        if reduced_base:
            self._apply_file_filter(
                self.reduced_treeview, reduced_base,
                self.reduced_filter_fits, self.reduced_filter_nonfits, self.reduced_filter_dirs,
                name_col_idx=self._reduced_name_col_idx,
            )

        # Suppress any spurious tree-activation events (e.g. a stray Enter
        # keypress from the dialog OK button) for a brief window after OK.
        self._suppress_activate = True
        QtCore.QTimer.singleShot(300, lambda: setattr(self, "_suppress_activate", False))

    def hide_reduced_tree_cb(self, w, val):
        """Show or hide the reduced-calibrations tree view.

        Parameters
        ----------
        w : ginga.gw.Widgets.CheckBox
            The "Show Reduced Tree" checkbox widget.
        val : bool
            ``True`` to show the tree, ``False`` to hide it.
        """
        if val:
            self.reduced_treeview.show()
        else:
            self.reduced_treeview.hide()

    def hide_raw_tree_cb(self, w, val):
        """Show or hide the raw-data tree view.

        Parameters
        ----------
        w : ginga.gw.Widgets.CheckBox
            The "Show Raw Tree" checkbox widget.
        val : bool
            ``True`` to show the tree, ``False`` to hide it.
        """
        if val:
            self.raw_treeview.show()
        else:
            self.raw_treeview.hide()

    def instrument_combo_cb(self, *args):
        """Swap the active instrument and refresh both tree views.

        Called when the instrument combo box selection changes.  Instantiates
        a new :class:`~.instruments.Instrument` via
        :class:`~.instruments.InstrumentRegistry` and rebuilds the tree
        column headers and listings for the new instrument vocabulary.

        Parameters
        ----------
        *args
            Positional arguments forwarded by the Ginga ``activated``
            callback; unused.
        """
        selected = self.instrument_combo.get_text()
        self.instrument = self.instrument_registry.create(selected)

        # Rebuild the treeview with the new columns:
        self._rebuild_treeview_columns()

    # --- Manual extraction ---

    def manual_extract_mode_cb(self, w, val):
        """Toggle manual-extraction click mode on/off."""
        self.manual_extract_mode = bool(val)
        if not self.manual_extract_mode:
            # Clear stored position and canvas marker
            self._manual_x = None
            self._manual_y = None
            self._manual_det_label = None
            self._manual_spat_det = None
            self._manual_slit_key = None
            self.state.manual_extract_str = None
            self.manual_extract_params_entry.set_text("")
            self._clear_manual_extract_marker()

    def fwhm_box_changed_cb(self, w):
        """Redraw the aperture marker and rebuild the params string when FWHM changes."""
        if self._manual_x is None:
            return
        self._update_manual_extract(self._manual_x, self._manual_y)

    def manual_extract_params_cb(self, w):
        """User edited the params text entry directly; store the new string and redraw."""
        text = self.manual_extract_params_entry.get_text().strip()
        if not text:
            self.state.manual_extract_str = None
            self._clear_manual_extract_marker()
            return
        # Try to parse it back so we can redraw the marker
        try:
            parts = text.split(":")
            # Format: det:spat:spec:fwhm[:boxcar_rad]
            spat = float(parts[1])
            spec = float(parts[2])
            fwhm = float(parts[3])
            # Recover image-space x from spat_det (best effort — mosaic offset unknown
            # without det_label context, so use stored values if available)
            x_img = self._manual_x if self._manual_x is not None else spat
            y_img = self._manual_y if self._manual_y is not None else spec
            # Re-parse FWHM from string and update the FWHM box to keep them in sync
            self.fwhm_box.set_text(f"{fwhm:.1f}")
            self._draw_manual_extract_marker(x_img, y_img, fwhm)
        except (IndexError, ValueError):
            pass
        self.state.manual_extract_str = text

    def _update_manual_extract(self, x: float, y: float) -> None:
        """Compute the extraction string and redraw the marker for click position (x, y)."""
        if not self.state.slittracesets:
            return

        # Find which detector/slit contains the click
        det_label = None
        spat_det = x
        found_det_idx = None
        found_slit_key = None
        for det_idx, slits in self.state.slittracesets.items():
            if slits is None:
                continue
            offset = (int(det_idx) - 1) * slits.nspat
            spec_row = int(np.clip(np.round(y), 0, slits.nspec - 1))
            left = slits.left_init[spec_row] + offset
            right = slits.right_init[spec_row] + offset
            for i in range(slits.nslits):
                if left[i] < x < right[i]:
                    found_det_idx = det_idx
                    found_slit_key = f"S{slits.spat_id[i]}"
                    det_label = f"{self.instrument.detector_prefix}{det_idx}"
                    spat_det = x - offset
                    break
            if det_label:
                break

        if det_label is None:
            self.logger.warning("Manual extract click did not land inside any slit.")
            return

        self._manual_x = x
        self._manual_y = y
        self._manual_det_label = det_label
        self._manual_spat_det = spat_det
        self._manual_slit_key = found_slit_key

        try:
            fwhm = float(self.fwhm_box.get_text())
        except ValueError:
            fwhm = 3.0
            self.fwhm_box.set_text("3.0")

        extract_str = f"{int(found_det_idx)}:{spat_det:.1f}:{y:.1f}:{fwhm:.1f}"
        self.state.manual_extract_str = extract_str
        self.manual_extract_params_entry.set_text(extract_str)
        self._draw_manual_extract_marker(x, y, fwhm)

    def _draw_manual_extract_marker(self, x: float, y: float, fwhm: float) -> None:
        """Draw a dot at (x, y) and a horizontal line of length fwhm centered on x.
        FWHM is measured in pixels, so this is easy."""
        if self.manual_extract_canvas is None:
            return
        self.manual_extract_canvas.delete_all_objects()

        half = fwhm / 2.0
        # Horizontal FWHM bar
        bar = self.dc.Line(x - half, y, x + half, y, color="yellow", linewidth=2)
        # Vertical tick marks at each end of the bar
        tick_h = 3.0
        tick_left = self.dc.Line(x - half, y - tick_h, x - half, y + tick_h,
                                 color="yellow", linewidth=2)
        tick_right = self.dc.Line(x + half, y - tick_h, x + half, y + tick_h,
                                  color="yellow", linewidth=2)
        # Central dot (small cross)
        dot_r = 3.0
        dot_h = self.dc.Line(x - dot_r, y, x + dot_r, y, color="yellow", linewidth=2)
        dot_v = self.dc.Line(x, y - dot_r, x, y + dot_r, color="yellow", linewidth=2)

        for obj in (bar, tick_left, tick_right, dot_h, dot_v):
            self.manual_extract_canvas.add(obj, redraw=False)
        self.manual_extract_canvas.update_canvas(whence=3)

    def _clear_manual_extract_marker(self) -> None:
        """Remove all manual-extraction marker objects from the canvas.

        No-op when ``manual_extract_canvas`` has not yet been initialised.
        """
        if self.manual_extract_canvas is not None:
            self.manual_extract_canvas.delete_all_objects()
            self.manual_extract_canvas.update_canvas(whence=3)

    def slit_list_box_cb(self, w, res_dict):
        """Combo-box callback: switch the active slit when the user picks a new entry.

        Deactivates (reverts to green) the previously active slit polygon on
        :attr:`slit_canvas`, then activates (highlights blue) the newly
        selected slit.  The slit key is extracted from the first token of the
        combo item text because entries are formatted as ``"S{spat_id} (name)"``
        or simply ``"S{spat_id}"``.

        Parameters
        ----------
        w : ginga widget
            The ``ComboBox`` widget that fired the callback.  ``w.get_text()``
            returns the full display string for the selected item.
        res_dict : dict
            Not used; present because Ginga passes a result dict to all
            ``"activated"`` callbacks.
        """
        if self.state.active_slit_key:
            self.overlay.deactivate(self.state.active_slit_key, self.slit_canvas)
        self.state.active_slit_key = w.get_text().split()[0]
        self.overlay.activate(self.state.active_slit_key, self.slit_canvas)

    def reduce_slit_cb(self, w):
        """Launch a PypeIt QuickLook reduction for the selected slit.

        Validates that a raw file, calibration directory, and redux output
        path are all set, then builds a ``ql.py`` argument list and submits
        it on a **daemon thread** via
        :meth:`~.backends.ReductionBackend.submit` so that the GUI remains
        responsive.  After submission a status row is added to the Reduction
        Control panel via :meth:`_make_reduction_row` and a polling timer is
        started via :meth:`_register_reduction_timer` to watch for output
        files.

        Parameters
        ----------
        w : ginga.gw.Widgets.Button
            The "Reduce Slit" button widget (unused).

        Notes
        -----
        The reduction subprocess is launched inside a :class:`threading.Thread`
        with ``daemon=True``.  All GUI updates that follow must be marshalled
        back to the GUI thread (Ginga handles this for timer callbacks
        automatically).
        """
        if self.state.manual_extract_str and self._manual_slit_key:
            slit_key = self._manual_slit_key
        else:
            slit_key = self.slit_list_box.get_text().split()[0]
        if not slit_key:
            self.logger.error("No slit selected for reduction.")
            return

        # Get the file paths:
        raw_path = self.state.raw_filepath or self.raw_text_entry.get_text()
        reduced_path = self.state.reduced_filepath or self.reduced_text_entry.get_text()
        redux_path = self.state.redux_path

        # Check the paths:
        # TODO: have the raw and reduced paths also open a dialog box?
        if not raw_path or not os.path.isfile(raw_path):
            self.logger.error("Raw file path is invalid or not selected.")
            return
        if not reduced_path or not os.path.isdir(reduced_path):
            self.logger.error("Reduced path is invalid or not selected.")
            return
        if not redux_path or not os.path.isdir(redux_path):
            dialog = QtGui.QMessageBox()
            dialog.setWindowTitle("Reduction Path Not Found")
            dialog.setIcon(QtGui.QMessageBox.Warning)
            dialog.setText(
                f"The reduction output path does not exist:\n\n{redux_path}\n\n"
                "Would you like to create it?"
            )
            create_btn = dialog.addButton("Create Directory", QtGui.QMessageBox.AcceptRole)
            dialog.addButton("Cancel", QtGui.QMessageBox.RejectRole)
            dialog.exec_()
            if dialog.clickedButton() != create_btn:
                return
            try:
                os.makedirs(redux_path)
            except OSError as exc:
                QtGui.QMessageBox.critical(None, "Error", f"Could not create directory:\n{exc}")
                return

        if not self.state.slittracesets:
            self.logger.error("No slit traces loaded. Render slits first.")
            return

        # Figure out which slit we're supposed to reduce and assemble the arguement
        # for it.
        slit_id = slit_key[1:] if slit_key.startswith("S") else slit_key
        det_label: Optional[str] = None
        for det_idx, slittrace in self.state.slittracesets.items():
            if slittrace is None:
                continue
            for spat_id in slittrace.spat_id:
                if slit_key == f"S{spat_id}":
                    det_label = f"{self.instrument.detector_prefix}{det_idx}"
                    break
            if det_label:
                break

        if not det_label:
            self.logger.error(f"Could not resolve detector index for slit {slit_key}.")
            return

        # Generate a unique directory for each reduction, to keep everything clean
        now = datetime.datetime.now()
        run_dir = os.path.join(redux_path, f"{det_label}_{slit_id}_{now.strftime('%H%M%S')}")
        os.makedirs(run_dir, exist_ok=True)

        # Create the logpath. The log is used to monitor the reduction status.
        log_path = os.path.abspath(os.path.join(run_dir, f"{det_label}_{slit_id}.log"))

        # If there's a B frame (i.e. we want to try doing an AB subtraction), 
        # check it.
        b_frame = self.state.ab_partner_filepath
        raw_files = [str(Path(raw_path).name)]
        if b_frame and os.path.isfile(b_frame):
            if Path(b_frame).parent.resolve() == Path(raw_path).parent.resolve():
                raw_files.append(str(Path(b_frame).name))
            else:
                self.logger.warning(
                    "B frame is in a different directory than A frame; ignoring for reduction."
                )

        # Arguments for the ql call. These are parsed directly as if they were
        # command line arguments.
        args = [
            self.instrument.pypeit_name,
            "--raw_files",
            *raw_files,
            "--raw_path",
            str(Path(raw_path).parent.absolute()),
            "--setup_calib_dir",
            f"{Path(reduced_path).absolute()}/Calibrations",
            "--slitspatnum",
            f"{det_label}:{slit_id}",
            "--redux_path",
            run_dir,
            "--skip_display",
            "--log_file",
            log_path,
            "-v", "2",
        ]
        if not self.state.manual_extract_str:
            args += ["--snr_thresh", self.SNR_box.get_text()]
        if self.state.manual_extract_str:
            args += ["--manual_extract", self.state.manual_extract_str]
        if self.coadd2d_box.get_state():
            args += ["--coadd2d"]
        self.logger.info("Launching reduction: {0}".format(" ".join(args)))

        def _run() -> None:
            """Submit the reduction job on the background thread.

            Calls :meth:`~.backends.ReductionBackend.submit` which may block
            for the full duration of the reduction when using
            :class:`~.backends.LocalReductionBackend`.  The GUI polls for
            output files independently via the timer registered in
            :meth:`_register_reduction_timer`.
            """
            self.reduction_backend.submit(args)

        # Launch the reduction in a thread
        threading.Thread(target=_run, daemon=True).start()

        # Create the buttons to view this reduction, and register them with the 
        # timer.
        coadd2d = self.coadd2d_box.get_state()
        self._make_reduction_row(slit_key, raw_path, now.strftime("%H:%M:%S"), coadd2d=coadd2d)
        self._register_reduction_timer(raw_path, run_dir, slit_key, log_path, coadd2d=coadd2d)

    def _register_reduction_timer(
        self, raw_path: str, run_dir: str, slit_key: str, log_path: str,
        coadd2d: bool = False,
    ) -> None:
        """Register a Ginga timer that polls *run_dir* for reduction output.

        Creates (or replaces) a timer entry in ``self.reduction_timers`` keyed
        by ``{raw_stem}_{slit_key}``.  On each expiry the timer fires
        :meth:`_check_reduction_complete`.  The start time is recorded in
        ``self.reduction_start_times`` so that :meth:`_check_reduction_complete`
        can enforce ``self.reduction_timeout``.

        Parameters
        ----------
        raw_path : str
            Absolute path to the raw FITS file used for the reduction.
            Passed through to :meth:`_check_reduction_complete` so that
            trace-overlay callbacks can gate on instrument configuration.
        run_dir : str
            Output directory created for this reduction job (under
            ``self.state.redux_path``).
        slit_key : str
            Slit identifier string (e.g. ``"S1234"``).
        log_path : str
            Absolute path to the reduction log file checked for error strings.
        coadd2d : bool, optional
            When ``True`` the timer uses two-phase polling: phase 1 waits for
            ``*/Science/spec1d*``, phase 2 waits for
            ``*/science_coadd/spec1d*``.
        """
        raw_stem = Path(raw_path).name.split(".fits")[0]
        timer_key = f"{raw_stem}_{slit_key}"

        existing = self.reduction_timers.get(timer_key)
        if existing is not None:
            existing.cancel()

        self.reduction_start_times[timer_key] = time.monotonic()

        timer = self.fitsimage.make_timer()
        timer.add_callback(
            "expired",
            lambda t: self._check_reduction_complete(
                timer_key, run_dir, slit_key, log_path,
                raw_path=raw_path, coadd2d=coadd2d,
            ),
        )
        self.reduction_timers[timer_key] = timer
        timer.set(self.reduction_cadence)

    def _check_reduction_complete(
        self, timer_key: str, run_dir: str, slit_key: str, log_path: str,
        raw_path: str = "", coadd2d: bool = False,
    ) -> None:
        """Timer callback: check for reduction output and update the status row.

        Called by the Ginga timer on the GUI thread every
        ``self.reduction_cadence`` seconds.  Implements two-phase polling
        when *coadd2d* is ``True``:

        - **Phase 1** — glob for ``*/Science/spec1d*.fits*`` inside
          *run_dir*.  Handles both single-frame output
          (``{raw_stem}/Science/``) and AB-pair output
          (``{A_stem}-{B_stem}/Science/``).
        - **Phase 2** (coadd2d only) — once ``spec1d_found`` is set, glob
          for ``*/science_coadd/spec1d*.fits*``.

        On success the relevant "Show" / "Show CoAdd2D" button is enabled.
        On timeout or logged exception the row label is set to an error
        string and the timer is cancelled.

        Also periodically watches the log for relevant strings

        Parameters
        ----------
        timer_key : str
            Key used to look up the timer and start-time entries.
        run_dir : str
            Reduction output directory to search for output files.
        slit_key : str
            Slit identifier (e.g. ``"S1234"``); used to locate the status
            row in ``self.reduction_control_elements``.
        log_path : str
            Absolute path to the reduction log file.
        raw_path : str, optional
            Path to the raw FITS file used for the reduction; forwarded to
            the "Show Traces" button callback lambda so that
            :meth:`show_traces_cb` can compare instrument configurations.
        coadd2d : bool, optional
            Whether to perform two-phase polling for CoAdd2D output.
        """
        timer = self.reduction_timers.get(timer_key)
        control = self.reduction_control_elements.get(slit_key)

        def _stop_timer() -> None:
            if timer is not None:
                timer.cancel()
                self.reduction_timers.pop(timer_key, None)
            self.reduction_start_times.pop(timer_key, None)

        def _fail(msg: str) -> None:
            self.logger.warning(msg)
            if control is not None:
                control["label"].set_text(f"Error {slit_key}")
            _stop_timer()

        # Timeout check
        start = self.reduction_start_times.get(timer_key)
        if start is not None and (time.monotonic() - start) > self.reduction_timeout:
            _fail(f"Reduction timed out for {timer_key} after {self.reduction_timeout:.0f}s.")
            return

        # Exception in log
        if self.file_backend.check_log_for_failure(log_path, "Exception"):
            _fail(f"Exception detected in reduction log for {timer_key}.")
            return

        if self.file_backend.check_log_for_failure(
            log_path, "No science frames found among the files provided."
        ):
            self.logger.warning(f"Reduction failed for {timer_key}: no science frames found.")
            if control is not None:
                control["label"].set_text(f"Failed {slit_key}: not a science frame")
            _stop_timer()
            return

        # Phase 2 (coadd2d only): regular spec1d already found; wait for
        # science_coadd output from the CoAdd2D pipeline.
        if coadd2d and control is not None and control.get("spec1d_found"):
            coadd_files = self.file_backend.glob(run_dir, "*/science_coadd/spec1d*.fits*")
            if coadd_files:
                coadd_path = coadd_files[0]
                self.logger.info(f"CoAdd2D complete for {timer_key}: {coadd_path}")
                control["label"].set_text(f"Reduced {slit_key}")
                control["btn_coadd2d"].set_enabled(True)
                control["btn_coadd2d"].add_callback(
                    "activated",
                    lambda w, p=coadd_path, k=slit_key: self.show_spec1d_cb(
                        w, path=p, slit_key=f"{k}_coadd2d"
                    ),
                )
                _stop_timer()
            else:
                self.logger.info(
                    f"Waiting for science_coadd output for {timer_key}; "
                    f"rechecking in {self.reduction_cadence}s"
                )
                if timer is not None:
                    timer.set(self.reduction_cadence)
            return

        # Phase 1: look for the regular per-frame spec1d in Science/.
        # For single-frame reductions ql.py creates {raw_stem}/Science/; for
        # AB pairs it creates {A_stem}-{B_stem}/Science/, so search one level deep.
        spec1d_files = self.file_backend.glob(run_dir, "*/Science/spec1d*.fits*")
        if spec1d_files:
            spec1d_path = spec1d_files[0]
            self.logger.info(f"Reduction complete for {timer_key}: {spec1d_path}")
            if control is not None:
                control["spec1d_found"] = True
                control["button"].set_enabled(True)
                control["button"].add_callback(
                    "activated",
                    lambda w, p=spec1d_path, k=slit_key: self.show_spec1d_cb(w, path=p, slit_key=k),
                )
                control["btn_traces"].set_enabled(True)
                control["btn_traces"].add_callback(
                    "activated",
                    lambda w, p=spec1d_path, k=slit_key, r=raw_path: self.show_traces_cb(
                        w, slit_key=k, spec1d_path=p, reduction_raw_path=r,
                    ),
                )
                if coadd2d:
                    control["label"].set_text(f"Coadd2D {slit_key}...")
                    if timer is not None:
                        timer.set(self.reduction_cadence)
                    return
                control["label"].set_text(f"Reduced {slit_key}")
            _stop_timer()
            return

        # Reduction finished (pipeline exited) but produced no spec1d — extraction failed.
        if self.file_backend.check_log_for_failure(log_path, "Quicklook execution time"):
            self.logger.warning(
                f"Reduction finished for {timer_key} but no spec1d written — extraction failed."
            )
            if control is not None:
                control["label"].set_text(f"Extraction failed {slit_key}")
            _stop_timer()
            return

        self.logger.info(
            f"Reduction not complete for {timer_key}; "
            f"rechecking in {self.reduction_cadence}s"
        )
        if timer is not None:
            timer.set(self.reduction_cadence)

    def _make_reduction_row(
        self, slit_key: str, raw_path: str, start_time: str, coadd2d: bool = False
    ) -> None:
        """Add a per-slit status row to vbox_redux.

        Displays the source filename, start time, a status label, a Show
        button (initially disabled), an optional Show CoAdd2D button (when
        coadd2d=True), and a Remove button.
        """
        from ginga.gw import Widgets as GWidgets

        filename = Path(raw_path).name

        vbox = GWidgets.VBox()

        # Top line: filename and start time
        hbox_info = GWidgets.HBox()
        hbox_info.add_widget(GWidgets.Label(f"{filename}  |  started {start_time}"), stretch=1)
        vbox.add_widget(hbox_info, stretch=0)

        # Bottom line: status label, Show Traces, Show, [Show CoAdd2D], Remove
        hbox_controls = GWidgets.HBox()
        label = GWidgets.Label(f"Reducing {slit_key}...")
        btn_traces = GWidgets.Button("Show Traces")
        btn_traces.set_enabled(False)
        btn_show = GWidgets.Button("Show")
        btn_show.set_enabled(False)
        btn_remove = GWidgets.Button("Remove")

        hbox_controls.add_widget(label, stretch=1)
        hbox_controls.add_widget(btn_traces, stretch=0)
        hbox_controls.add_widget(btn_show, stretch=0)

        btn_coadd2d = None
        if coadd2d:
            btn_coadd2d = GWidgets.Button("Show CoAdd2D")
            btn_coadd2d.set_enabled(False)
            hbox_controls.add_widget(btn_coadd2d, stretch=0)

        hbox_controls.add_widget(btn_remove, stretch=0)
        vbox.add_widget(hbox_controls, stretch=0)

        self.vbox_redux.add_widget(vbox, stretch=0)

        def _remove(w):
            """Button callback: remove this status row and clean up its resources.

            Removes the VBox row from :attr:`vbox_redux`, deletes the slit key
            from :attr:`reduction_control_elements`, detaches any trace canvas
            from the main Ginga canvas, and cancels any active trace-polling
            timer.

            Parameters
            ----------
            w : ginga widget
                The ``Remove`` button widget that fired the callback.
            """
            self.vbox_redux.remove(vbox)
            self.reduction_control_elements.pop(slit_key, None)
            canvas = self._trace_canvases.pop(slit_key, None)
            if canvas is not None:
                try:
                    self.fitsimage.get_canvas().delete_object(canvas)
                except Exception:
                    pass
            timer = self._trace_timers.pop(slit_key, None)
            if timer is not None:
                timer.cancel()
            self._trace_paths.pop(slit_key, None)
            self._trace_last_exten.pop(slit_key, None)

        btn_remove.add_callback("activated", _remove)

        self.reduction_control_elements[slit_key] = {
            "label": label,
            "button": btn_show,
            "btn_traces": btn_traces,
            "btn_coadd2d": btn_coadd2d,
            "spec1d_found": False,
            "vbox": vbox,
        }

    def show_spec1d_cb(self, w, path: Optional[str] = None, slit_key: Optional[str] = None):
        """Button callback: open a spec1d FITS file in a dedicated Ginga channel.

        Loads *path* into a new channel named ``"Spec1D{slit_key}"`` and
        immediately starts the ``Spec1dView`` Ginga plugin on that channel so
        the user can inspect extracted spectra interactively.

        Called both by the "Show" button in a reduction status row (regular
        spec1d) and by the "Show CoAdd2D" button (coadded spec1d), with the
        *slit_key* suffixed with ``"_coadd2d"`` in the latter case.

        Parameters
        ----------
        w : ginga widget
            The button widget that fired the callback.
        path : str, optional
            Absolute path to the spec1d FITS file to display.  If *None* or
            the file does not exist the method logs an error and returns.
        slit_key : str, optional
            Slit identifier (e.g. ``"S1234"`` or ``"S1234_coadd2d"``); used
            to construct a unique Ginga channel name so multiple slits can be
            viewed simultaneously.
        """
        if not path or not os.path.isfile(path):
            self.logger.error("No spec1d file available to show.")
            return
        self.logger.info(f"Showing reduced spectrum: {path}")
        ch_name = f"Spec1D{slit_key}" if slit_key else "Spec1D"
        self.fv.load_file(path, chname=ch_name)
        self.fv.start_local_plugin(ch_name, "Spec1dView")

    def show_traces_cb(
        self, w, *, slit_key: str, spec1d_path: str, reduction_raw_path: str = ""
    ) -> None:
        """Load object traces from a spec1d file and overlay them on the raw image.

        Refuses to draw if the currently displayed raw image belongs to a
        different calibration group (instrument configuration) than the raw
        frame used for this reduction.

        Runs file I/O on a background thread, then draws on the GUI thread.
        A polling timer watches the paired Spec1dView plugin and highlights
        whichever object is currently selected there.

        This will not be happy if the Spec1DView isn't open!
        """
        current_raw = self.state.raw_filepath
        if reduction_raw_path and current_raw and current_raw != reduction_raw_path:
            current_cfg = self._get_raw_config(current_raw)
            reduction_cfg = self._get_raw_config(reduction_raw_path)
            if current_cfg and reduction_cfg and not self._configs_match(
                current_cfg, reduction_cfg
            ):
                self.logger.warning(
                    f"Cannot draw traces for {slit_key}: the displayed image "
                    f"({Path(current_raw).name}) has a different instrument "
                    f"configuration from the reduction frame "
                    f"({Path(reduction_raw_path).name}). "
                    f"Select the matching raw file to view traces."
                )
                return

        def _load():
            """Background worker: load SpecObjs from *spec1d_path* and schedule drawing.

            Runs on a daemon thread started by :meth:`show_traces_cb`.  On
            success calls :meth:`_draw_trace_objects` via ``fv.gui_do`` so
            that canvas operations happen on the GUI thread.  Any exception is
            logged without propagating.
            """
            try:
                from pypeit.specobjs import SpecObjs
                sobjs = SpecObjs.from_fitsfile(spec1d_path, chk_version=False)
                self.fv.gui_do(self._draw_trace_objects, slit_key, sobjs)
            except Exception as exc:
                self.logger.error(
                    f"Failed to load traces from {spec1d_path}: {exc}", exc_info=True
                )
        threading.Thread(target=_load, daemon=True).start()

    def _draw_trace_objects(self, slit_key: str, sobjs) -> None:
        """Draw per-object extraction traces on the raw image canvas.
        
        Largely looks like parts of show_spec2D"""
        canvas = self._trace_canvases.get(slit_key)
        if canvas is None:
            canvas = DrawingCanvas()
            canvas.enable_draw(False)
            canvas.enable_edit(False)
            canvas.set_surface(self.fitsimage)
            self.fitsimage.get_canvas().add(canvas)
            self._trace_canvases[slit_key] = canvas
        else:
            canvas.delete_all_objects()

        paths: Dict[int, dict] = {}
        step = 10  # subsample spectral rows for canvas performance

        for i in range(sobjs.nobj):
            sobj = sobjs[i]
            trace = sobj.TRACE_SPAT
            if trace is None or len(trace) == 0:
                continue
            nspec = len(trace)

            # Compute spatial offset for multi-detector mosaics (0 for single-det)
            det_str = sobj.DET or "01"
            det_num = int(det_str[-2:]) if len(det_str) >= 2 else 1
            slitset = (
                self.state.slittracesets.get(f"{det_num:02d}")
                if self.state.slittracesets
                else None
            )
            x_offset = (det_num - 1) * (slitset.nspat if slitset is not None else 0)

            spec_rows = np.arange(nspec)[::step]
            trace_x = trace[::step] + x_offset
            pts = list(zip(trace_x.tolist(), spec_rows.tolist()))
            path = self.dc.Path(pts, color="orange", linewidth=1.5)
            canvas.add(path, redraw=False)

            mid = nspec // 2
            name = sobj.NAME or str(i)
            lbl = self.dc.Text(
                float(trace[mid]) + x_offset, float(mid), name,
                color="orange", fontsize=8, rot_deg=90,
            )
            canvas.add(lbl, redraw=False)
            paths[i] = {"path": path, "label": lbl}

        canvas.update_canvas(whence=3)
        self._trace_paths[slit_key] = paths
        self._trace_last_exten[slit_key] = None

        # Start polling the paired Spec1dView channel for selection changes.
        existing = self._trace_timers.get(slit_key)
        if existing is not None:
            existing.cancel()
        timer = self.fitsimage.make_timer()
        timer.add_callback(
            "expired",
            lambda t: self._poll_trace_highlight(slit_key),
        )
        self._trace_timers[slit_key] = timer
        # This is pretty frequent, but doesn't appear to have adverse resource issues
        timer.set(0.5)

    def _poll_trace_highlight(self, slit_key: str) -> None:
        """Poll the paired Spec1dView for its current selection and recolor traces."""
        ch_name = f"Spec1D{slit_key}"
        current: Optional[int] = None
        try:
            plugin = self.fv.get_channel(ch_name).opmon.get_plugin("Spec1dView")
            current = plugin.exten
        except Exception:
            pass

        last = self._trace_last_exten.get(slit_key)
        if current != last:
            canvas = self._trace_canvases.get(slit_key)
            paths = self._trace_paths.get(slit_key, {})
            if canvas is not None:
                if last is not None and last in paths:
                    paths[last]["path"].color = "orange"
                    paths[last]["label"].color = "orange"
                if current is not None and current in paths:
                    paths[current]["path"].color = "cyan"
                    paths[current]["label"].color = "cyan"
                canvas.update_canvas(whence=3)
            self._trace_last_exten[slit_key] = current

        timer = self._trace_timers.get(slit_key)
        if timer is not None:
            timer.set(0.5)

    def show_labels_box_cb(self, w, val):
        """Checkbox callback: toggle slit-label text visibility on the canvas.

        Delegates to :meth:`~pypeit.display.qlview.slit_overlay.SlitOverlay.set_labels_visible`.
        Is a no-op when no slit canvas has been created yet.

        Parameters
        ----------
        w : ginga widget
            The ``CheckBox`` widget that fired the callback.
        val : bool
            ``True`` to show slit labels; ``False`` to hide them.
        """
        if self.slit_canvas is not None:
            self.overlay.set_labels_visible(val, self.slit_canvas)

    def display_slits_box_cb(self, w, val):
        """Checkbox callback: show or hide all slit polygons on the canvas.

        Sets the ``alpha`` (outline) and ``fillalpha`` (interior) of every
        polygon stored in :attr:`state.slit_polys` to their visible values
        (1.0 / 0.05) or to zero, then redraws the canvas.

        Parameters
        ----------
        w : ginga widget
            The ``CheckBox`` widget that fired the callback.
        val : bool
            ``True`` to make slit polygons visible; ``False`` to hide them.
        """
        for poly in self.state.slit_polys.values():
            poly.alpha = 1.0 if val else 0.0
            poly.fillalpha = 0.05 if val else 0.0
        if self.slit_canvas is not None:
            self.slit_canvas.update_canvas(whence=3)

    # ── File selection stuff ──────────────────────────────────────────────────

    def reduced_table_double_click_cb(self, w, res_dict):
        """Tree-view callback: navigate into a directory on double-click.

        Registered on the ``"activated"`` event of :attr:`reduced_treeview`.
        Ignored when :attr:`_suppress_activate` is set (used to prevent
        spurious activations during programmatic tree updates).  Only
        directories trigger navigation; file entries are silently ignored
        because the reduced tree does not open files on double-click.

        Parameters
        ----------
        w : ginga widget
            The ``TreeView`` widget that fired the callback.
        res_dict : dict
            Mapping of row key → row-info object with a ``path`` attribute,
            as provided by the Ginga tree-view selection API.
        """
        if self._suppress_activate:
            return
        paths = [info.path for key, info in res_dict.items()]
        if not paths:
            return
        path = paths[0]
        if os.path.isdir(path):
            self._browse_and_update(path, which_tree="reduced")

    def reduced_table_selected_cb(self, w, res_dict):
        """Tree-view callback: update the reduced-path entry and button state on selection.

        Registered on the ``"selected"`` event of :attr:`reduced_treeview`.
        Writes the selected path to :attr:`reduced_text_entry` and
        :attr:`state.reduced_filepath`, then enables :attr:`reduced_btn` and
        :attr:`show_wavelengths_btn` only when the selected entry is a valid
        calibration directory (i.e., it is a directory with a ``Calibrations/``
        subdirectory directly inside it).

        Parameters
        ----------
        w : ginga widget
            The ``TreeView`` widget that fired the callback.
        res_dict : dict
            Mapping of row key → row-info object with a ``path`` attribute.
        """
        paths = [info.path for key, info in res_dict.items()]
        if not paths:
            return
        path = paths[0]
        self.reduced_text_entry.set_text(path)
        self.state.reduced_filepath = path

        # Enable render/wavelength buttons only when the selected path is a
        # calibration directory (i.e., has a Calibrations/ subdirectory directly
        # inside it).  Use a direct exists check rather than iterating all entries
        # so that permissions errors or non-directory selections are handled safely.
        try:
            enabled = Path(path).is_dir() and Path(path, "Calibrations").is_dir()
        except OSError:
            enabled = False
        self.reduced_btn.set_enabled(enabled)
        self.show_wavelengths_btn.set_enabled(enabled)

    def raw_table_double_click_cb(self, w, res_dict):
        """Tree-view callback: navigate or open a raw file on double-click.

        Registered on the ``"activated"`` event of :attr:`raw_treeview`.
        Ignored when :attr:`_suppress_activate` is set.

        - **Directory** → calls :meth:`_browse_and_update` to navigate into it.
        - **FITS file** → calls :meth:`_check_instrument_match`, then
          :meth:`open_raw_file` and :meth:`_suggest_calibrations` if the
          instrument check passes.

        In both cases :attr:`state.raw_filepath` is updated to the selected
        path.

        Parameters
        ----------
        w : ginga widget
            The ``TreeView`` widget that fired the callback.
        res_dict : dict
            Mapping of row key → row-info object with a ``path`` attribute.
        """
        if self._suppress_activate:
            return
        paths = [info.path for key, info in res_dict.items()]
        if not paths:
            return
        path = paths[0]
        if os.path.isdir(path):
            self._browse_and_update(path, which_tree="raw")
        elif os.path.isfile(path):
            if self._check_instrument_match(path):
                self.open_raw_file(path)
                self._suggest_calibrations(path)
        self.state.raw_filepath = path

    def raw_table_selected_cb(self, w, res_dict):
        """Tree-view callback: update the raw-path entry on single selection.

        Registered on the ``"selected"`` event of :attr:`raw_treeview`.
        Updates :attr:`raw_text_entry` and :attr:`state.raw_filepath` with the
        first selected path; no file is opened and no buttons are toggled.

        Parameters
        ----------
        w : ginga widget
            The ``TreeView`` widget that fired the callback.
        res_dict : dict
            Mapping of row key → row-info object with a ``path`` attribute.
        """
        paths = [info.path for key, info in res_dict.items()]
        if not paths:
            return
        path = paths[0]
        self.raw_text_entry.set_text(path)
        self.state.raw_filepath = path

    # ── Dithering support stuff ───────────────────────────────────────────────

    def b_frame_entry_cb(self, w):
        """User manually typed a B frame path."""
        path = self.b_frame_entry.get_text().strip()
        self.state.ab_partner_filepath = path if path and os.path.isfile(path) else None

    def clear_b_frame_cb(self, w):
        """Clear the B frame field."""
        self.b_frame_entry.set_text("")
        self.state.ab_partner_filepath = None

    def detect_ab_pair_cb(self, w) -> None:
        """Button callback: search for a matching B frame in a background thread."""
        raw_path = self.state.raw_filepath or self.raw_text_entry.get_text().strip()
        if not raw_path or not os.path.isfile(raw_path):
            self.b_frame_entry.set_text("No raw file selected")
            return
        self.b_frame_entry.set_text("Searching…")
        self.state.ab_partner_filepath = None

        def _worker():
            """Background worker: run AB-pair detection and schedule UI update.

            Calls :meth:`_suggest_ab_partner` on a daemon thread started by
            :meth:`detect_ab_pair_cb`, then delivers the result to the GUI
            thread via ``fv.gui_do`` so that widget updates happen safely.
            """
            partner = self._suggest_ab_partner(raw_path)
            self.fv.gui_do(self._apply_b_frame_suggestion, partner)

        threading.Thread(target=_worker, daemon=True).start()

    def _apply_b_frame_suggestion(self, partner: Optional[str]) -> None:
        """GUI-thread callback: apply a detected B-frame path to the UI.

        Called via ``fv.gui_do`` from the ``_worker`` closure inside
        :meth:`detect_ab_pair_cb`.  When *partner* is set, writes the path to
        :attr:`b_frame_entry` and :attr:`state.ab_partner_filepath` and logs
        the discovery.  When *partner* is ``None``, writes a friendly
        "not found" message to the entry and logs accordingly.

        Parameters
        ----------
        partner : str or None
            Absolute path of the best-matching B-position frame, or ``None``
            if no suitable candidate was found.
        """
        if partner:
            self.b_frame_entry.set_text(partner)
            self.state.ab_partner_filepath = partner
            self.logger.info(f"Detected B frame: {partner}")
        else:
            self.b_frame_entry.set_text("No matching B frame found")
            self.logger.info("No matching B frame found in directory.")

    def _suggest_ab_partner(self, raw_path: str) -> Optional[str]:
        """Return the best-matching B frame for *raw_path*, or None.

        Uses PypeIt's own spectrograph metadata system — the same fields that
        ``set_combination_groups`` / ``get_comb_group`` uses — to match on
        instrument configuration keys, target, exptime, and dither pattern.
        Picks the B/B' candidate with the closest MJD.
        """
        try:
            from pypeit.spectrographs.util import load_spectrograph
            spec = load_spectrograph(self.instrument.pypeit_name)
        except Exception as exc:
            self.logger.warning(f"Could not load spectrograph for B frame suggestion: {exc}")
            return None

        def _get(headarr, key):
            try:
                return spec.get_meta_value(headarr, key, ignore_bad_header=True)
            except Exception:
                return None

        try:
            a_headarr = spec.get_headarr(raw_path)
        except Exception:
            return None

        dithpos = _get(a_headarr, "dithpos") or ""
        if dithpos not in ("A", "A'"):
            return None

        a_dithpat = _get(a_headarr, "dithpat") or ""
        a_target  = _get(a_headarr, "target")
        a_obsmode = _get(a_headarr, "obsmode")
        try:
            a_exptime = float(_get(a_headarr, "exptime") or 0)
        except (TypeError, ValueError):
            a_exptime = None
        try:
            a_mjd = float(_get(a_headarr, "mjd") or 0)
        except (TypeError, ValueError):
            a_mjd = 0.0

        config_keys = spec.configuration_keys()
        a_config = {k: _get(a_headarr, k) for k in config_keys}

        raw_dir = Path(raw_path).parent
        best_path: Optional[str] = None
        best_delta = float("inf")

        for candidate in raw_dir.glob("*.fits*"):
            if str(candidate) == raw_path:
                continue
            try:
                c_headarr = spec.get_headarr(str(candidate))
            except Exception:
                continue

            # Must be B/B' dither position
            cpos = _get(c_headarr, "dithpos") or ""
            if cpos not in ("B", "B'"):
                continue

            # Dither pattern must match
            cpat = _get(c_headarr, "dithpat") or ""
            if a_dithpat and cpat and cpat != a_dithpat:
                continue

            # All spectrograph configuration keys must match
            if any(_get(c_headarr, k) != a_config[k] for k in config_keys):
                continue

            # Target must match
            if a_target and _get(c_headarr, "target") != a_target:
                continue

            # Obs mode must match
            if a_obsmode and _get(c_headarr, "obsmode") != a_obsmode:
                continue

            # Exptime must match within 10 %
            if a_exptime:
                try:
                    cexp = float(_get(c_headarr, "exptime") or 0)
                    if abs(cexp - a_exptime) > 0.1 * a_exptime:
                        continue
                except (TypeError, ValueError):
                    pass

            # Pick the candidate closest in time
            try:
                cmjd = float(_get(c_headarr, "mjd") or 0)
            except (TypeError, ValueError):
                cmjd = 0.0
            delta = abs(cmjd - a_mjd)
            if delta < best_delta:
                best_delta = delta
                best_path = str(candidate)

        return best_path

    def _suggest_calibrations(self, raw_path: str) -> None:
        """Auto-populate reduced cals path with the best matching calibration.
        This takes a long time, potentially, so it uns the PypeIt calibration
        search on a background thread and updates the status label and reduced
        tree when finished.
        """
        cal_root = self._get_tree_base_dir(self.state.reduced_filepath)
        if not cal_root:
            cal_root = self.reduced_text_entry.get_text().strip()
        if not cal_root or not os.path.isdir(cal_root):
            if hasattr(self, "cal_status_label"):
                self.cal_status_label.set_text("")
            return

        if hasattr(self, "cal_status_label"):
            self.cal_status_label.set_text("Searching for matching calibrations...")

        instrument = self.instrument  # capture for thread

        def _search() -> None:
            try:
                suggestions = instrument.recommend_calibrations(raw_path, cal_root)
                if suggestions:
                    self.fv.gui_do(self._on_cal_found, suggestions[0])
                else:
                    self.fv.gui_do(self._on_cal_not_found)
            except Exception as exc:
                self.fv.gui_do(self._on_cal_error, str(exc))

        threading.Thread(target=_search, daemon=True).start()

    def _on_cal_found(self, best: str) -> None:
        if hasattr(self, "cal_status_label"):
            self.cal_status_label.set_text("Calibrations matched")

        # ``best`` is the Calibrations/ subdirectory (e.g. keck_mosfire_A/Calibrations).
        # We want to highlight the setup directory (keck_mosfire_A) in its parent listing.
        cal_set_path = Path(best).parent   # keck_mosfire_A/
        parent = str(cal_set_path.parent)  # directory that contains keck_mosfire_A
        dir_name = cal_set_path.name       # "keck_mosfire_A"

        # Navigate to the parent directory if the tree is not already showing it
        current_base = self._get_tree_base_dir(self.state.reduced_filepath)
        if os.path.normpath(current_base or "") != os.path.normpath(parent):
            self._browse_and_update(parent, which_tree="reduced")

        # Highlight the matching calibration set without navigating into it
        self._select_reduced_tree_item(dir_name)
        self.reduced_text_entry.set_text(str(cal_set_path))
        self.state.reduced_filepath = str(cal_set_path)
        self.reduced_btn.set_enabled(True)

    def _on_cal_not_found(self) -> None:
        if hasattr(self, "cal_status_label"):
            self.cal_status_label.set_text("No calibrations found")

    def _on_cal_error(self, msg: str) -> None:
        if hasattr(self, "cal_status_label"):
            self.cal_status_label.set_text(f"Error: {msg[:60]}")
        self.logger.error(f"Calibration suggestion error: {msg}")

    def _select_reduced_tree_item(self, name: str) -> None:
        """Select and scroll to a top-level item in the reduced treeview by name."""
        widget = self.reduced_treeview.widget
        widget.clearSelection()
        for i in range(widget.topLevelItemCount()):
            item = widget.topLevelItem(i)
            if item.text(self._reduced_name_col_idx) == name:
                item.setSelected(True)
                widget.scrollToItem(item)
                break

    def _check_instrument_match(self, path: str) -> bool:
        """Check that the file's INSTRUME keyword matches the selected instrument.

        Returns True when it is safe to proceed with opening the file (either
        the instruments match, the keyword is absent, or the user chose to
        continue / switch instruments).  Returns False when the user cancelled.
        """
        try:
            with fits.open(path) as hdul:
                instrume = hdul[0].header.get("INSTRUME", None)
        except Exception:
            return True  # can't read header; proceed gracefully

        if instrume is None:
            return True

        instrume = instrume.strip().upper()
        expected = self.instrument.instrume_value.upper()

        if not expected or instrume == expected:
            return True

        # Find the registry entry that matches the file's instrument.
        # Use instrume_values() to avoid constructing an instrument instance
        # just to read a class attribute.
        matching_name: Optional[str] = None
        for name, instrume_val in self.instrument_registry.instrume_values():
            if instrume_val.upper() == instrume:
                matching_name = name
                break

        current_name = self.instrument.__class__.__name__
        msg = (
            f"The selected file reports INSTRUME = '{instrume}',\n"
            f"but the active instrument is '{current_name}'."
        )

        dialog = QtGui.QMessageBox()
        dialog.setWindowTitle("Instrument Mismatch")
        dialog.setText(msg)

        switch_btn = None
        if matching_name:
            switch_btn = dialog.addButton(
                f"Switch to {matching_name}", QtGui.QMessageBox.ActionRole
            )
        continue_btn = dialog.addButton("Continue Anyway", QtGui.QMessageBox.ActionRole)
        cancel_btn = dialog.addButton("Cancel", QtGui.QMessageBox.RejectRole)
        dialog.setDefaultButton(cancel_btn)
        dialog.exec_()

        clicked = dialog.clickedButton()
        if switch_btn is not None and clicked == switch_btn:
            idx = self.instrument_registry.names().index(matching_name)
            self.instrument_combo.set_index(idx)
            self.instrument = self.instrument_registry.create(matching_name)
            self._rebuild_treeview_columns()
            return True
        if clicked == continue_btn:
            return True
        return False  # user cancelled

    def render_slits_cb(self, w):
        """Button callback: load slit-trace files and draw slit polygons on the canvas.

        Reads ``Slits_*.fits*`` files from the ``Calibrations/`` subdirectory
        of :attr:`state.reduced_filepath`, builds the slit overlay via
        :attr:`overlay`, and populates :attr:`slit_list_box` with an entry for
        every slit.  Also snapshots the calibration instrument configuration
        into :attr:`_rendered_cal_config` so that :meth:`open_raw_file` can
        clear the overlay when the user subsequently opens a frame from a
        different observing setup.

        Parameters
        ----------
        w : ginga widget
            The ``Button`` widget that fired the callback.
        """
        if not self.state.reduced_filepath:
            self.logger.error("No reduced filepath set.")
            return

        cal_path = self.state.reduced_filepath
        # Strip trailing glob wildcard left by _browse_and_update
        if cal_path.endswith("*"):
            cal_path = os.path.dirname(cal_path)
        cal_path = os.path.join(cal_path, "Calibrations")
        self.logger.info(f"Searching for slit files in: {cal_path}")

        slitsets = self.open_slits_files(cal_path)
        if not slitsets:
            self.logger.error(f"No Slits_*.fits* files found in {cal_path}")
            return
        self.logger.info(f"Loaded {len(slitsets)} slit set(s): {list(slitsets.keys())}")
        self.state.slittracesets = slitsets

        slit_names = self._build_slit_names(slitsets)

        self.slit_canvas.delete_all_objects()
        show_labels = self.show_labels_box.get_state()
        slit_polys = self.overlay.build(slitsets, self.slit_canvas,
                                        slit_names=slit_names, show_labels=show_labels)
        self.state.slit_polys = slit_polys
        self.logger.info(f"Built {len(slit_polys)} slit polygons.")

        self.slit_list_box.clear()
        self._slit_combo_keys = []
        for slit_key in sorted(self.state.slit_polys.keys()):
            display = f"{slit_key} ({slit_names[slit_key]})" if slit_key in slit_names else slit_key
            self.slit_list_box.append_text(display)
            self._slit_combo_keys.append(slit_key)

        # Snapshot the calibration configuration so we can detect setup changes
        # when the user selects a different raw file.
        cal_dir = self.state.reduced_filepath or ""
        if cal_dir.endswith("*"):
            cal_dir = os.path.dirname(cal_dir)
        self._rendered_cal_config = self.instrument._read_pypeit_setup_config(cal_dir)

    def show_wavelengths_cb(self, w):
        """Build a wavelength image from calibration files and display it in a new channel.

        Loads WaveCalib, WaveTilts, and Slits files from the selected calibration
        directory, builds a 2D wavelength image, and displays it using the
        SlitWavelength Ginga plugin so that hovering over the image shows the
        wavelength at each pixel.
        """
        if not self.state.reduced_filepath:
            self.logger.error("No reduced filepath set.")
            return

        cal_path = self.state.reduced_filepath
        if cal_path.endswith("*"):
            cal_path = os.path.dirname(cal_path)
        cal_path = os.path.join(cal_path, "Calibrations")
        self.logger.info(f"Searching for wavelength calibration files in: {cal_path}")

        wv_files = sorted(glob.glob(os.path.join(cal_path, "WaveCalib_*.fits*")))
        tilt_files = sorted(glob.glob(os.path.join(cal_path, "Tilts_*.fits*")))
        slit_files = sorted(glob.glob(os.path.join(cal_path, "Slits_*.fits*")))

        if not wv_files:
            self.logger.error(f"No WaveCalib files found in {cal_path}")
            return
        if not tilt_files:
            self.logger.error(f"No Tilts files found in {cal_path}")
            return
        if not slit_files:
            self.logger.error(f"No Slits files found in {cal_path}")
            return

        def _suffix(filename):
            """Extract the _{setup}_{id}_{det} suffix from a calibration filename."""
            m = re.match(r'(?:WaveCalib|Tilts|Slits)(_[^.]+)', os.path.basename(filename))
            return m.group(1) if m else ""

        # Build a lookup by suffix for tilts and slits.
        tilt_by_suffix = {_suffix(f): f for f in tilt_files}
        slit_by_suffix = {_suffix(f): f for f in slit_files}

        # Find the first WaveCalib that has matching Tilts and Slits.
        matched = None
        for wv_file in wv_files:
            sfx = _suffix(wv_file)
            if sfx in tilt_by_suffix and sfx in slit_by_suffix:
                matched = (wv_file, tilt_by_suffix[sfx], slit_by_suffix[sfx])
                break

        if matched is None:
            self.logger.error(
                "Could not find matching WaveCalib/Tilts/Slits triplet in "
                f"{cal_path}"
            )
            return

        wv_file, tilt_file, slit_file = matched
        self.logger.info(
            f"Loading wavelength calibration:\n"
            f"  WaveCalib: {wv_file}\n"
            f"  Tilts:     {tilt_file}\n"
            f"  Slits:     {slit_file}"
        )

        def _build_and_show():
            try:
                from pypeit.wavecalib import WaveCalib
                from pypeit.wavetilts import WaveTilts
                from pypeit.slittrace import SlitTraceSet

                wvcalib = WaveCalib.from_file(wv_file)
                wavetilts = WaveTilts.from_file(tilt_file)
                slits = SlitTraceSet.from_file(slit_file)

                slitmask = slits.slit_img()
                tilts = wavetilts.fit2tiltimg(slitmask, flexure=wavetilts.spat_flexure)
                waveimg = wvcalib.build_waveimg(tilts, slits).astype(np.float32)

                # Collect per-slit RMS values: list of (spat_id, rms_or_None)
                rms_data = []
                if wvcalib.spat_ids is not None and wvcalib.wv_fits is not None:
                    for spat_id, wvfit in zip(wvcalib.spat_ids, wvcalib.wv_fits):
                        rms = None if (wvfit is None or wvfit.rms is None) else wvfit.rms
                        rms_data.append((int(spat_id), rms))

                self.fv.gui_do(self._display_waveimg, waveimg, rms_data)
            except Exception as exc:
                self.logger.error(f"Failed to build wavelength image: {exc}", exc_info=True)

        threading.Thread(target=_build_and_show, daemon=True).start()

    def _display_waveimg(self, waveimg: np.ndarray, rms_data=None) -> None:
        """Attach a wavelength map to the currently displayed raw image.

        Replaces the current AstroImage in fitsimage with a SlitImage that
        carries the wavelength array.  When hovering over a pixel, Ginga will
        show the wavelength value (in Å) in the info bar instead of RA/Dec.

        Overlays per-slit λ RMS annotations on the slit canvas.
        """
        from pypeit.display.slitwavelength import SlitImage

        current = self.fitsimage.get_image()
        if current is None:
            self.logger.error("No image currently displayed; open a raw file first.")
            return

        raw_data = current.get_data()
        if raw_data.shape != waveimg.shape:
            self.logger.warning(
                f"Raw image shape {raw_data.shape} differs from wavelength map "
                f"shape {waveimg.shape}. Wavelength cursor will only be accurate "
                "within the overlapping region."
            )

        wave_img = SlitImage(wav_np=waveimg, logger=self.logger)
        wave_img.load_data(raw_data)
        self.fitsimage.set_image(wave_img)

        # Draw per-slit λ RMS labels on the slit canvas
        if rms_data and self.slit_canvas is not None:
            nspec = waveimg.shape[0]
            # Place labels in the top-left corner of the image, stacked vertically
            x_pos = 5.0
            line_height = max(12.0, nspec * 0.012)
            y_start = line_height
            for i, (spat_id, rms) in enumerate(rms_data):
                rms_str = f"{rms:.3f}" if rms is not None else "N/A"
                label = f"\u03bb RMS S{spat_id}: {rms_str} \u00c5"
                y_pos = y_start + i * line_height
                txt = self.dc.Text(
                    x_pos, y_pos, label,
                    color="cyan", fontsize=10,
                )
                self.slit_canvas.add(txt, redraw=False)
            self.slit_canvas.update_canvas(whence=3)

    def raw_button_cb(self, w):
        """Button callback: handle the "Go" button next to the raw-data path entry.

        Reads the current text of :attr:`raw_text_entry` and behaves as
        follows:

        - **Directory** → calls :meth:`_browse_and_update` to populate the
          raw tree with the directory's contents.
        - **File** → calls :meth:`_check_instrument_match` and, on success,
          :meth:`open_raw_file` followed by :meth:`_suggest_calibrations`.
        - **Invalid path** → writes ``"Invalid path"`` back into the entry.

        Parameters
        ----------
        w : ginga widget
            The ``Button`` widget that fired the callback.
        """
        path = self.raw_text_entry.get_text()
        if os.path.isdir(path):
            self._browse_and_update(path, which_tree="raw")
        elif os.path.isfile(path):
            if self._check_instrument_match(path):
                self.open_raw_file(path)
                self._suggest_calibrations(path)
        else:
            self.raw_text_entry.set_text("Invalid path")

    def canvas_clicked_cb(self, canvas, pnt, x, y):
        """Cursor-down callback: handle clicks on the main image canvas.

        Registered on the ``"cursor-down"`` event of :attr:`fitsimage` by
        :meth:`start`.  Two modes are supported:

        - **Manual extraction mode** (``self.manual_extract_mode is True``) —
          delegates to :meth:`_update_manual_extract` with the click
          coordinates.
        - **Slit selection mode** — searches :attr:`state.slittracesets` for
          the slit that contains the clicked pixel, updates the slit combo-box
          selection, and highlights the slit polygon on :attr:`slit_canvas`.

        Is a no-op when no slit trace sets have been loaded.

        Parameters
        ----------
        canvas : ginga Canvas
            The image canvas that received the click.
        pnt : object
            Ginga point object (not used directly).
        x : float
            Image-coordinate column of the click.
        y : float
            Image-coordinate row of the click.
        """
        if not self.state.slittracesets:
            return

        if self.manual_extract_mode:
            self._update_manual_extract(float(x), float(y))
            return

        for msc_idx, slits in self.state.slittracesets.items():
            if slits is None:
                continue
            offset = (int(msc_idx) - 1) * slits.nspat
            row = int(np.clip(np.round(y), 0, slits.nspec - 1))
            left_bound_at_y = slits.left_init[row] + offset
            right_bound_at_y = slits.right_init[row] + offset

            for i in range(slits.nslits):
                if left_bound_at_y[i] < x < right_bound_at_y[i]:
                    slit_id = slits.spat_id[i]
                    slit_key = f"S{slit_id}"
                    if slit_key in self._slit_combo_keys:
                        self.slit_list_box.set_index(self._slit_combo_keys.index(slit_key))
                    if self.state.active_slit_key:
                        self.overlay.deactivate(self.state.active_slit_key, self.slit_canvas)
                    self.state.active_slit_key = slit_key
                    self.overlay.activate(slit_key, self.slit_canvas)
                    return

    # --- Helpers ---

    def _browse_and_update(self, path: str, which_tree: str) -> None:
        """Navigate to *path* and refresh the appropriate file-browser tree.

        Delegates to :attr:`file_browser` to obtain the sorted directory
        listing, then pushes it into the raw or reduced ``TreeView`` widget.
        After a navigation into the reduced tree the "Render Slits" and "Show
        Wavelengths" button states are updated based on whether a
        ``Calibrations/`` subdirectory exists in the newly browsed directory.

        Parameters
        ----------
        path : str
            Directory to browse.  Relative paths are resolved by the file
            browser.
        which_tree : str
            ``"reduced"`` to update the reduced calibrations tree and button
            state; any other value updates the raw data tree.

        Raises
        ------
        The method silently catches :class:`ValueError` raised by
        :meth:`~pypeit.display.qlview.file_browser.FileBrowserController.browse`
        when *path* is invalid and sets an error label in the relevant text
        entry instead of propagating.
        """
        mode = "reduced" if which_tree == "reduced" else "raw"
        columns = self.instrument.columns[mode]
        try:
            listing, resize, fullpath = self.file_browser.browse(
                path, self.instrument, columns=columns, mode=mode
            )
        except ValueError:
            if which_tree == "reduced":
                self.reduced_text_entry.set_text("Not a valid path")
            else:
                self.raw_text_entry.set_text("Not a valid path")
            return

        if which_tree == "reduced":
            self.state.reduced_filepath = fullpath
            self.reduced_treeview.set_tree(listing)
            if resize:
                self.reduced_treeview.set_optimal_column_widths()
            self._apply_file_filter(
                self.reduced_treeview,
                self._get_tree_base_dir(self.state.reduced_filepath),
                self.reduced_filter_fits,
                self.reduced_filter_nonfits,
                self.reduced_filter_dirs,
                name_col_idx=self._reduced_name_col_idx,
            )
            # Update render/wavelength button state based on the directory we
            # just navigated into.  fullpath ends with "/*"; dirname is the dir.
            cal_dir = os.path.dirname(fullpath)
            try:
                enabled = bool(cal_dir) and os.path.isdir(
                    os.path.join(cal_dir, "Calibrations")
                )
            except OSError:
                enabled = False
            self.reduced_btn.set_enabled(enabled)
            self.show_wavelengths_btn.set_enabled(enabled)
        else:
            self.state.raw_filepath = fullpath
            self.raw_treeview.set_tree(listing)
            if resize:
                self.raw_treeview.set_optimal_column_widths()
            self._apply_file_filter(
                self.raw_treeview,
                self._get_tree_base_dir(self.state.raw_filepath),
                self.raw_filter_fits,
                self.raw_filter_nonfits,
                self.raw_filter_dirs,
                name_col_idx=self._raw_name_col_idx,
            )

    def _get_raw_config(self, filepath: str) -> Dict[str, str]:
        """Return the instrument configuration_keys() values for a raw FITS file.

        Uses PypeIt's spectrograph class so the result is authoritative and
        consistent with how PypeIt groups frames into calibration sets.  Returns
        an empty dict if the file cannot be read or the spectrograph is unknown.
        """
        if not self.instrument.pypeit_name or not filepath or not os.path.isfile(filepath):
            return {}
        try:
            from pypeit.spectrographs.util import load_spectrograph
            spec = load_spectrograph(self.instrument.pypeit_name)
            config: Dict[str, str] = {}
            for key in spec.configuration_keys():
                try:
                    val = spec.get_meta_value(filepath, key)
                    config[key] = str(val) if val is not None else "N/A"
                except Exception:
                    pass
            return config
        except Exception as exc:
            self.logger.debug(f"Could not read configuration from {filepath}: {exc}")
            return {}

    def _configs_match(self, cfg_a: Dict[str, str], cfg_b: Dict[str, str]) -> bool:
        """Return True if two configuration dicts are compatible.

        Compares only keys that appear in both dicts.  Numeric values are
        compared with a generous tolerance to absorb floating-point formatting
        differences between pypeit-file strings and live header reads.
        """
        if not cfg_a or not cfg_b:
            return True
        common = set(cfg_a) & set(cfg_b)
        if not common:
            return True
        for k in common:
            a, b = cfg_a[k].strip(), cfg_b[k].strip()
            if a == b:
                continue
            try:
                if abs(float(a) - float(b)) < 0.5:
                    continue
            except ValueError:
                pass
            return False
        return True

    def _get_tree_base_dir(self, path: Optional[str]) -> Optional[str]:
        """Return the directory that the tree is currently browsing.

        The internal path stored in :attr:`state.raw_filepath` or
        :attr:`state.reduced_filepath` after a :meth:`_browse_and_update` call
        ends with ``"/*"`` (the glob wildcard appended by the file browser).
        This helper strips the wildcard and resolves the parent.

        Parameters
        ----------
        path : str or None
            The path stored in plugin state (may end with ``"/*"``), a plain
            directory, a file path, or ``None``.

        Returns
        -------
        str or None
            Absolute path of the directory currently shown in the tree, or
            ``None`` when *path* is falsy.

        Examples
        --------
        >>> self._get_tree_base_dir("/data/raw/keck_deimos_A/*")
        '/data/raw/keck_deimos_A'
        >>> self._get_tree_base_dir("/data/raw/frame.fits")
        '/data/raw'
        >>> self._get_tree_base_dir("/data/raw")
        '/data/raw'
        """
        if not path:
            return None
        if path.endswith("*"):
            return os.path.abspath(os.path.join(path, os.pardir))
        if os.path.isdir(path):
            return path
        return os.path.abspath(os.path.join(path, os.pardir))

    def _apply_file_filter(
        self,
        treeview,
        base_dir: Optional[str],
        show_fits: bool,
        show_nonfits: bool,
        show_dirs: bool,
        name_col_idx: int = 0,
    ) -> None:
        """Hide or show tree-view rows according to the current filter settings.

        Iterates every top-level item in *treeview*'s underlying
        ``QTreeWidget`` and sets its hidden state based on whether the
        corresponding filesystem entry is a directory, a ``.fits`` /
        ``.fits.gz`` file, or some other file type.

        Called by :meth:`_browse_and_update` after each navigation and also
        when the user toggles filter checkboxes in the Settings dialog.

        Parameters
        ----------
        treeview : ginga TreeView widget
            The tree-view to filter.  The underlying Qt widget is accessed via
            ``treeview.widget``.
        base_dir : str or None
            Directory whose contents are currently shown in *treeview*.  Used
            to resolve each row's display name to an absolute path for the
            ``os.path.isdir`` test.
        show_fits : bool
            Whether to show ``.fits`` and ``.fits.gz`` files.
        show_nonfits : bool
            Whether to show all other non-directory entries.
        show_dirs : bool
            Whether to show subdirectory entries.
        name_col_idx : int, optional
            Column index that holds the filename in the ``QTreeWidgetItem``.
            Defaults to ``0``.
        """
        if not base_dir or treeview is None:
            return

        tree_iterator = QtGui.QTreeWidgetItemIterator(treeview.widget)
        while tree_iterator.value():
            item = tree_iterator.value()
            name = item.text(name_col_idx)
            child_path = os.path.join(base_dir, name)
            if os.path.isdir(child_path):
                item.setHidden(not show_dirs)
            else:
                is_fits = name.lower().endswith(".fits") or name.lower().endswith(".fits.gz")
                if is_fits:
                    item.setHidden(not show_fits)
                else:
                    item.setHidden(not show_nonfits)
            tree_iterator += 1

    def _build_slit_names(self, slitsets: Dict[str, SlitTraceSet]) -> Dict[str, str]:
        """Return a dict mapping slit key → OBJNAME from mask design, if available."""
        names: Dict[str, str] = {}
        for slittrace in slitsets.values():
            maskdef_id_arr = getattr(slittrace, "maskdef_id", None)
            maskdef_designtab = getattr(slittrace, "maskdef_designtab", None)
            if maskdef_id_arr is None or maskdef_designtab is None:
                continue
            design = maskdef_designtab
            for idx, spat_id in enumerate(slittrace.spat_id):
                maskdef_id = maskdef_id_arr[idx]
                rows = design[design['MASKDEF_ID'] == maskdef_id]
                if len(rows) > 0:
                    objname = str(rows['OBJNAME'][0]).strip()
                    if objname:
                        names[f"S{spat_id}"] = objname
        return names

    def open_slits_files(self, path: str) -> Dict[str, SlitTraceSet]:
        """Read all slit-trace FITS files under *path* and return them keyed by mosaic index.

        Searches *path* for files matching ``Slits_*.fits*`` and loads each
        one as a :class:`~pypeit.slittrace.SlitTraceSet`.  The mosaic detector
        index is extracted from the filename stem: for a file named
        ``Slits_A_0_DET01.fits``, the last two characters of the last
        ``_``-separated token before the extension give the index (``"01"``).

        Parameters
        ----------
        path : str
            Directory to search; usually the ``Calibrations/`` subdirectory
            of a PypeIt reduction output directory.

        Returns
        -------
        dict
            Mapping of mosaic index string (e.g. ``"01"``) →
            :class:`~pypeit.slittrace.SlitTraceSet`.  An empty dict is
            returned when no matching files exist.
        """
        p = Path(path)
        slit_files = p.glob("Slits_*.fits*")
        slit_dict: Dict[str, SlitTraceSet] = {}
        for slit_file in slit_files:
            msc_idx = slit_file.stem.split(".")[0].split("_")[-1][-2:]
            slit_dict[msc_idx] = SlitTraceSet.from_file(slit_file)
        return slit_dict

    def open_raw_file(self, path: str) -> None:
        """Load a raw FITS file and display it in the main Ginga viewer.

        Delegates image assembly to
        :meth:`~pypeit.display.qlview.instruments.Instrument.get_display_image`,
        wraps the result in a Ginga :class:`~ginga.AstroImage.AstroImage`, and
        pushes it into :attr:`fitsimage`.

        After loading, all existing trace overlays are cleared (they are
        spatially meaningless on a different frame) and their polling timers
        are cancelled.  If slits have previously been rendered and the new
        file's instrument configuration differs from the calibration that
        produced them, the slit overlay is also cleared.

        Is a no-op when *path* does not exist or does not contain ``".fits"``
        in its name.

        Parameters
        ----------
        path : str
            Absolute path to the raw FITS file to open.
        """
        if not os.path.isfile(path):
            return
        p = Path(path)
        if ".fits" not in p.name:
            return
        img_data = self.instrument.get_display_image(path)
        img = AstroImage(logger=self.logger)
        img.load_data(img_data)
        self.fitsimage.set_image(img)
        self.state.raw_filepath = path

        # --- Clear trace overlays ---
        # Traces drawn for a previous raw image are spatially meaningless on a
        # different frame; remove them from the canvas and stop polling timers.
        for canvas in self._trace_canvases.values():
            try:
                canvas.delete_all_objects()
            except Exception:
                pass
        for timer in self._trace_timers.values():
            timer.cancel()
        self._trace_timers.clear()
        self._trace_paths.clear()
        self._trace_last_exten.clear()

        # --- Optionally clear the slit overlay ---
        # If slits have been rendered, compare the new file's instrument
        # configuration against the cal setup that produced those slits.  If
        # the setups differ (different filter, grating, mask, …), the slit
        # polygons no longer correspond to this image.
        if self._rendered_cal_config:
            new_config = self._get_raw_config(path)
            if new_config and not self._configs_match(new_config, self._rendered_cal_config):
                self.logger.info(
                    f"New raw file has a different instrument configuration — "
                    f"clearing slit overlay."
                )
                self.slit_canvas.delete_all_objects()
                self.state.slittracesets = None
                self.state.slit_polys = {}
                self.slit_list_box.clear()
                self._slit_combo_keys = []
                self._rendered_cal_config = {}

    # --- Ginga plugin lifecycle ---

    def close(self):
        """Stop this Ginga local plugin.

        Calls the Ginga framework method
        ``fv.stop_local_plugin(chname, plugin_name)`` to deactivate the plugin
        and remove its UI from the panel.  Ginga subsequently calls
        :meth:`stop` to release all timers and canvases.
        """
        self.fv.stop_local_plugin(self.chname, str(self))

    @staticmethod
    def _expand_path(raw: str) -> str:
        """Expand strftime-style format codes in a path using the current datetime.

        A path with no ``%`` characters is returned unchanged, so absolute
        paths that happen not to contain any format codes are unaffected.

        Examples
        --------
        ``/data/raw/%Y%m%d``  →  ``/data/raw/20260313``
        ``/hqdrpdata/outputs`` →  ``/hqdrpdata/outputs``  (unchanged)
        """
        return datetime.datetime.now().strftime(raw)

    def start(self):
        """Ginga plugin lifecycle hook: initialize state and populate the UI.

        Called by the Ginga framework when the plugin is started.  Reads
        ``~/.quicklook.cfg`` (if it exists) to seed the raw-data path,
        reduced-calibrations path, reduction output path, file-filter flags,
        and reduction timeout.  Template keys (``raw_path_template``,
        ``reduced_path_template``, ``redux_path_template``) take precedence
        over their plain counterparts and support ``strftime``-style format
        codes expanded at startup via :meth:`_expand_path`.

        After applying configuration values, populates both file-browser
        trees, creates the slit-overlay and manual-extraction ``DrawingCanvas``
        objects, and registers :meth:`canvas_clicked_cb` on the main image
        canvas.
        """
        config_file = Path.home() / ".quicklook.cfg"
        raw_path = os.getcwd()
        reduced_path = os.getcwd()
        redux_path = self.state.redux_path

        if config_file.exists():
            config = configparser.ConfigParser()
            config.read(config_file)
            defaults = config["DEFAULT"]

            # Template keys (raw_path_template, reduced_path_template,
            # redux_path_template) take precedence over the static keys when
            # present.  They are expanded with strftime so that format codes
            # such as %Y, %m, %d, %H, %M are replaced with the current date/
            # time at startup.  Absolute paths without any % characters work
            # identically with both the template and static keys.
            raw_tmpl = defaults.get("raw_path_template", "")
            raw_path = self._expand_path(raw_tmpl) if raw_tmpl else defaults.get("raw_path", raw_path)

            red_tmpl = defaults.get("reduced_path_template", "")
            reduced_path = self._expand_path(red_tmpl) if red_tmpl else defaults.get("reduced_path", reduced_path)

            rdx_tmpl = defaults.get("redux_path_template", "")
            redux_path = self._expand_path(rdx_tmpl) if rdx_tmpl else defaults.get("redux_path", redux_path)

            self.raw_filter_fits = defaults.getboolean("raw_show_fits", True)
            self.raw_filter_nonfits = defaults.getboolean("raw_show_nonfits", False)
            self.raw_filter_dirs = defaults.getboolean("raw_show_dirs", True)
            self.reduced_filter_fits = defaults.getboolean("reduced_show_fits", False)
            self.reduced_filter_nonfits = defaults.getboolean("reduced_show_nonfits", False)
            self.reduced_filter_dirs = defaults.getboolean("reduced_show_dirs", True)
            try:
                self.reduction_timeout = float(defaults.get("reduction_timeout", str(self.reduction_timeout)))
            except ValueError:
                pass
            self.logger.info(f"Loaded config from {config_file}")

        self.state.redux_path = redux_path
        self.raw_text_entry.set_text(raw_path)
        self.reduced_text_entry.set_text(reduced_path)

        self._browse_and_update(raw_path, which_tree="raw")
        self._browse_and_update(reduced_path, which_tree="reduced")
        self.instrument_combo.make_callback("activated")

        self.slit_canvas = DrawingCanvas()
        self.slit_canvas.enable_draw(False)
        self.slit_canvas.enable_edit(False)
        self.slit_canvas.set_surface(self.fitsimage)
        self.fitsimage.get_canvas().add(self.slit_canvas)

        self.manual_extract_canvas = DrawingCanvas()
        self.manual_extract_canvas.enable_draw(False)
        self.manual_extract_canvas.enable_edit(False)
        self.manual_extract_canvas.set_surface(self.fitsimage)
        self.fitsimage.get_canvas().add(self.manual_extract_canvas)

        self.fitsimage.add_callback("cursor-down", self.canvas_clicked_cb)

    def stop(self):
        """Ginga plugin lifecycle hook: cancel all timers and clean up canvases.

        Called by the Ginga framework when the plugin is stopped or the window
        is closed.  Cancels every active reduction-polling timer and every
        trace-highlight polling timer, then removes the slit-overlay canvas,
        manual-extraction canvas, and all per-slit trace canvases from the
        main Ginga image canvas.  Sets :attr:`gui_up` to ``False`` to signal
        that widget access is no longer safe.
        """
        for timer in self.reduction_timers.values():
            timer.cancel()
        self.reduction_timers.clear()
        for timer in self._trace_timers.values():
            timer.cancel()
        self._trace_timers.clear()
        if self.slit_canvas is not None:
            try:
                self.fitsimage.get_canvas().delete_object(self.slit_canvas)
            except Exception:
                pass
        if self.manual_extract_canvas is not None:
            try:
                self.fitsimage.get_canvas().delete_object(self.manual_extract_canvas)
            except Exception:
                pass
        for canvas in self._trace_canvases.values():
            try:
                self.fitsimage.get_canvas().delete_object(canvas)
            except Exception:
                pass
        self._trace_canvases.clear()
        self.gui_up = False

    def __str__(self):
        return "qlview"
