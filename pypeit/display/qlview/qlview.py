from __future__ import annotations

import configparser
import datetime
import os
import re
import threading
from pathlib import Path
from typing import Dict, Optional

import numpy as np
from astropy.io import fits
from ginga import GingaPlugin
from ginga.AstroImage import AstroImage
from ginga.canvas.types.layer import DrawingCanvas
from ginga.qtw.QtHelp import QtGui

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
        """Recompute the name-column indices for the current instrument's column defs."""
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
        self.ui.build(container)
        self.gui_up = True

    # --- Callbacks ---

    def create_config_cb(self, w):
        config_file = Path.home() / ".quicklook.cfg"
        config = configparser.ConfigParser()

        raw_path = self.raw_text_entry.get_text()
        if os.path.isfile(raw_path):
            raw_path = str(Path(raw_path).parent)

        reduced_path = self.reduced_text_entry.get_text()
        if re.match(r'^.+_[A-Za-z]$', Path(reduced_path).name):
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
        per-tree file-filter toggles.  Changes are applied on OK.
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

            try:
                from pypeit.display.qlview.backends import _requests
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

        raw_base = self._get_tree_base_dir(self.state.raw_filepath)
        if raw_base:
            self._apply_file_filter(
                self.raw_treeview, raw_base,
                self.raw_filter_fits, self.raw_filter_nonfits, self.raw_filter_dirs,
                name_col_idx=self._raw_name_col_idx,
            )
        reduced_base = self._get_tree_base_dir(self.state.reduced_filepath)
        if reduced_base:
            self._apply_file_filter(
                self.reduced_treeview, reduced_base,
                self.reduced_filter_fits, self.reduced_filter_nonfits, self.reduced_filter_dirs,
                name_col_idx=self._reduced_name_col_idx,
            )

    def hide_reduced_tree_cb(self, w, val):
        if val:
            self.reduced_treeview.show()
        else:
            self.reduced_treeview.hide()

    def hide_raw_tree_cb(self, w, val):
        if val:
            self.raw_treeview.show()
        else:
            self.raw_treeview.hide()

    def instrument_combo_cb(self, *args):
        selected = self.instrument_combo.get_text()
        self.instrument = self.instrument_registry.create(selected)
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
        for det_idx, slits in self.state.slittracesets.items():
            if slits is None:
                continue
            offset = (int(det_idx) - 1) * slits.nspat
            spec_row = int(np.clip(np.round(y), 0, slits.nspec - 1))
            left = slits.left_init[spec_row] + offset
            right = slits.right_init[spec_row] + offset
            for i in range(slits.nslits):
                if left[i] < x < right[i]:
                    det_label = f"{self.instrument.detector_prefix}{det_idx}"
                    spat_det = x - offset
                    slit_key = f"S{slits.spat_id[i]}"
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
        self._manual_slit_key = slit_key

        try:
            fwhm = float(self.fwhm_box.get_text())
        except ValueError:
            fwhm = 3.0
            self.fwhm_box.set_text("3.0")

        extract_str = f"{int(det_idx)}:{spat_det:.1f}:{y:.1f}:{fwhm:.1f}"
        self.state.manual_extract_str = extract_str
        self.manual_extract_params_entry.set_text(extract_str)
        self._draw_manual_extract_marker(x, y, fwhm)

    def _draw_manual_extract_marker(self, x: float, y: float, fwhm: float) -> None:
        """Draw a dot at (x, y) and a horizontal line of length fwhm centered on x."""
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
        if self.manual_extract_canvas is not None:
            self.manual_extract_canvas.delete_all_objects()
            self.manual_extract_canvas.update_canvas(whence=3)

    def reduce_slit_cb(self, w):
        if self.state.manual_extract_str and self._manual_slit_key:
            slit_key = self._manual_slit_key
        else:
            slit_key = self.slit_list_box.get_text().split()[0]
        if not slit_key:
            self.logger.error("No slit selected for reduction.")
            return

        raw_path = self.state.raw_filepath or self.raw_text_entry.get_text()
        reduced_path = self.state.reduced_filepath or self.reduced_text_entry.get_text()
        redux_path = self.state.redux_path

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

        now = datetime.datetime.now()
        run_dir = os.path.join(redux_path, f"{det_label}_{slit_id}_{now.strftime('%H%M%S')}")
        os.makedirs(run_dir, exist_ok=True)

        log_path = os.path.abspath(os.path.join(run_dir, f"{det_label}_{slit_id}.log"))

        args = [
            self.instrument.pypeit_name,
            "--raw_files",
            str(Path(raw_path).name),
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
            "-v",
            str(2)
        ]
        if not self.state.manual_extract_str:
            args += ["--snr_thresh", self.SNR_box.get_text()]
        if self.state.manual_extract_str:
            args += ["--manual_extract", self.state.manual_extract_str]
        self.logger.info("Launching reduction: {0}".format(" ".join(args)))

        def _run() -> None:
            self.reduction_backend.submit(args)

        threading.Thread(target=_run, daemon=True).start()

        self._make_reduction_row(slit_key, raw_path, now.strftime("%H:%M:%S"))
        self._register_reduction_timer(raw_path, run_dir, slit_key, log_path)

    def _register_reduction_timer(
        self, raw_path: str, redux_path: str, slit_key: str, log_path: str
    ) -> None:
        raw_stem = Path(raw_path).name.split(".fits")[0]
        timer_key = f"{raw_stem}_{slit_key}"
        science_dir = os.path.join(redux_path, raw_stem, "Science")

        existing = self.reduction_timers.get(timer_key)
        if existing is not None:
            existing.cancel()

        import time as _time
        self.reduction_start_times[timer_key] = _time.monotonic()

        timer = self.fitsimage.make_timer()
        timer.add_callback(
            "expired",
            lambda t: self._check_reduction_complete(timer_key, science_dir, slit_key, log_path),
        )
        self.reduction_timers[timer_key] = timer
        timer.set(self.reduction_cadence)

    def _check_reduction_complete(
        self, timer_key: str, science_dir: str, slit_key: str, log_path: str
    ) -> None:
        import time as _time
        timer = self.reduction_timers.get(timer_key)

        def _fail(msg: str) -> None:
            self.logger.warning(msg)
            control = self.reduction_control_elements.get(slit_key)
            if control is not None:
                control["label"].set_text(f"Error {slit_key}")
            if timer is not None:
                timer.cancel()
                self.reduction_timers.pop(timer_key, None)
            self.reduction_start_times.pop(timer_key, None)

        # Timeout check
        start = self.reduction_start_times.get(timer_key)
        if start is not None and (_time.monotonic() - start) > self.reduction_timeout:
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
            control = self.reduction_control_elements.get(slit_key)
            if control is not None:
                control["label"].set_text(f"Failed {slit_key}: not a science frame")
            if timer is not None:
                timer.cancel()
                self.reduction_timers.pop(timer_key, None)
            self.reduction_start_times.pop(timer_key, None)
            return

        spec1d_files = self.file_backend.glob(science_dir, "spec1d*.fits*")
        if spec1d_files:
            spec1d_path = spec1d_files[0]
            self.logger.info(f"Reduction complete for {timer_key}: {spec1d_path}")
            control = self.reduction_control_elements.get(slit_key)
            if control is not None:
                control["label"].set_text(f"Reduced {slit_key}")
                control["button"].set_enabled(True)
                control["button"].add_callback(
                    "activated",
                    lambda w, p=spec1d_path, k=slit_key: self.show_spec1d_cb(w, path=p, slit_key=k),
                )
                control["btn_traces"].set_enabled(True)
                control["btn_traces"].add_callback(
                    "activated",
                    lambda w, p=spec1d_path, k=slit_key: self.show_traces_cb(w, slit_key=k, spec1d_path=p),
                )
            if timer is not None:
                timer.cancel()
                self.reduction_timers.pop(timer_key, None)
            self.reduction_start_times.pop(timer_key, None)
            return

        # Reduction finished (pipeline exited) but produced no spec1d — extraction failed.
        if self.file_backend.check_log_for_failure(log_path, "Quicklook execution time"):
            self.logger.warning(f"Reduction finished for {timer_key} but no spec1d written — extraction failed.")
            control = self.reduction_control_elements.get(slit_key)
            if control is not None:
                control["label"].set_text(f"Extraction failed {slit_key}")
                control["button"].set_enabled(True)
            if timer is not None:
                timer.cancel()
                self.reduction_timers.pop(timer_key, None)
            self.reduction_start_times.pop(timer_key, None)
            return

        self.logger.info(
            f"Reduction not complete for {timer_key}; "
            f"rechecking in {self.reduction_cadence}s"
        )
        if timer is not None:
            timer.set(self.reduction_cadence)

    def slit_list_box_cb(self, w, res_dict):
        if self.state.active_slit_key:
            self.overlay.deactivate(self.state.active_slit_key, self.slit_canvas)
        self.state.active_slit_key = w.get_text().split()[0]
        self.overlay.activate(self.state.active_slit_key, self.slit_canvas)

    def _make_reduction_row(self, slit_key: str, raw_path: str, start_time: str) -> None:
        """Add a per-slit status row to vbox_redux.

        Displays the source filename, start time, a status label, a Show
        button (initially disabled), and a Remove button.
        """
        from ginga.gw import Widgets as GWidgets

        filename = Path(raw_path).name

        vbox = GWidgets.VBox()

        # Top line: filename and start time
        hbox_info = GWidgets.HBox()
        hbox_info.add_widget(GWidgets.Label(f"{filename}  |  started {start_time}"), stretch=1)
        vbox.add_widget(hbox_info, stretch=0)

        # Bottom line: status label, Show Traces button, Show button, Remove button
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
        hbox_controls.add_widget(btn_remove, stretch=0)
        vbox.add_widget(hbox_controls, stretch=0)

        self.vbox_redux.add_widget(vbox, stretch=0)

        def _remove(w):
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
            "vbox": vbox,
        }

    def show_spec1d_cb(self, w, path: Optional[str] = None, slit_key: Optional[str] = None):
        if not path or not os.path.isfile(path):
            self.logger.error("No spec1d file available to show.")
            return
        self.logger.info(f"Showing reduced spectrum: {path}")
        ch_name = f"Spec1D{slit_key}" if slit_key else "Spec1D"
        self.fv.load_file(path, chname=ch_name)
        self.fv.start_local_plugin(ch_name, "Spec1dView")

    def show_traces_cb(self, w, *, slit_key: str, spec1d_path: str) -> None:
        """Load object traces from a spec1d file and overlay them on the raw image.

        Runs file I/O on a background thread, then draws on the GUI thread.
        A polling timer watches the paired Spec1dView plugin and highlights
        whichever object is currently selected there.
        """
        def _load():
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
        """Draw per-object extraction traces on the raw image canvas."""
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
        if self.slit_canvas is not None:
            self.overlay.set_labels_visible(val, self.slit_canvas)

    def display_slits_box_cb(self, w, val):
        for poly in self.state.slit_polys.values():
            poly.alpha = 1.0 if val else 0.0
            poly.fillalpha = 0.05 if val else 0.0
        if self.slit_canvas is not None:
            self.slit_canvas.update_canvas(whence=3)

    def reduced_table_double_click_cb(self, w, res_dict):
        paths = [info.path for key, info in res_dict.items()]
        if not paths:
            return
        path = paths[0]
        if os.path.isdir(path):
            self._browse_and_update(path, which_tree="reduced")

    def reduced_table_selected_cb(self, w, res_dict):
        paths = [info.path for key, info in res_dict.items()]
        if not paths:
            return
        path = paths[0]
        self.reduced_text_entry.set_text(path)
        self.state.reduced_filepath = path

        if "keck_" in path:
            p = Path(path)
            subdirs = [x for x in p.iterdir() if x.is_dir()]
            if "Calibrations" in [x.name for x in subdirs]:
                self.reduced_btn.set_enabled(True)
                self.show_wavelengths_btn.set_enabled(True)
                return
        self.reduced_btn.set_enabled(False)
        self.show_wavelengths_btn.set_enabled(False)

    def raw_table_double_click_cb(self, w, res_dict):
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
        paths = [info.path for key, info in res_dict.items()]
        if not paths:
            return
        path = paths[0]
        self.raw_text_entry.set_text(path)
        self.state.raw_filepath = path

    def _suggest_calibrations(self, raw_path: str) -> None:
        """Auto-populate reduced cals path with the best matching calibration.

        Runs the PypeIt calibration search on a background thread and updates
        the status label and reduced tree when finished.
        """
        cal_root = self._get_tree_base_dir(self.state.reduced_filepath)
        if not cal_root:
            cal_root = self.reduced_text_entry.get_text().strip()
        if not cal_root or not os.path.isdir(cal_root):
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

        # Find the registry entry that matches the file's instrument
        matching_name: Optional[str] = None
        for name in self.instrument_registry.names():
            inst = self.instrument_registry.create(name)
            if inst.instrume_value.upper() == instrume:
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

        import glob as _glob
        import re

        wv_files = sorted(_glob.glob(os.path.join(cal_path, "WaveCalib_*.fits*")))
        tilt_files = sorted(_glob.glob(os.path.join(cal_path, "Tilts_*.fits*")))
        slit_files = sorted(_glob.glob(os.path.join(cal_path, "Slits_*.fits*")))

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
        if not self.state.slittracesets:
            return

        if self.manual_extract_mode:
            self._update_manual_extract(float(x), float(y))
            return

        for msc_idx, slits in self.state.slittracesets.items():
            if slits is None:
                continue
            offset = (int(msc_idx) - 1) * slits.nspat
            left_bound_at_y = slits.left_init[np.round(y).astype(int)] + offset
            right_bound_at_y = slits.right_init[np.round(y).astype(int)] + offset

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

    def _get_tree_base_dir(self, path: Optional[str]) -> Optional[str]:
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
            if slittrace.maskdef_id is None or slittrace.maskdef_designtab is None:
                continue
            design = slittrace.maskdef_designtab
            for idx, spat_id in enumerate(slittrace.spat_id):
                maskdef_id = slittrace.maskdef_id[idx]
                rows = design[design['MASKDEF_ID'] == maskdef_id]
                if len(rows) > 0:
                    objname = str(rows['OBJNAME'][0]).strip()
                    if objname:
                        names[f"S{spat_id}"] = objname
        return names

    def open_slits_files(self, path: str) -> Dict[str, SlitTraceSet]:
        p = Path(path)
        slit_files = p.glob("Slits_*.fits*")
        slit_dict: Dict[str, SlitTraceSet] = {}
        for slit_file in slit_files:
            msc_idx = slit_file.stem.split(".")[0].split("_")[-1][-2:]
            slit_dict[msc_idx] = SlitTraceSet.from_file(slit_file)
        return slit_dict

    def open_raw_file(self, path: str) -> None:
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

    # --- Ginga plugin lifecycle ---

    def close(self):
        self.fv.stop_local_plugin(self.chname, str(self))

    def start(self):
        config_file = Path.home() / ".quicklook.cfg"
        raw_path = os.getcwd()
        reduced_path = os.getcwd()
        redux_path = self.state.redux_path

        if config_file.exists():
            config = configparser.ConfigParser()
            config.read(config_file)
            defaults = config["DEFAULT"]
            raw_path = defaults.get("raw_path", raw_path)
            reduced_path = defaults.get("reduced_path", reduced_path)
            redux_path = defaults.get("redux_path", redux_path)
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
