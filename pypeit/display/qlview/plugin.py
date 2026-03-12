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


class QLTEST(GingaPlugin.LocalPlugin):
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

        self._raw_name_col_idx = 0
        self._reduced_name_col_idx = 0
        self._compute_name_col_indices()

        self.dc = fitsimage.get_canvas().get_draw_classes()
        self.slit_canvas: Optional[DrawingCanvas] = None
        self.overlay = SlitOverlay(self.dc)
        self.file_browser = FileBrowserController(self.logger, self.settings, self.file_backend)
        self.reduction_timers: Dict[str, object] = {}
        self.reduction_cadence: float = 5.0
        self.reduction_control_elements: Dict[str, dict] = {}
        self._slit_combo_keys: list = []  # slit keys in combo box order

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
            "redux_path": self.redux_path_entry.get_text(),
            "raw_path": raw_path,
            "reduced_path": reduced_path,
            "raw_show_fits": str(self.raw_show_fits.get_state()),
            "raw_show_nonfits": str(self.raw_show_nonfits.get_state()),
            "raw_show_dirs": str(self.raw_show_dirs.get_state()),
            "reduced_show_fits": str(self.reduced_show_fits.get_state()),
            "reduced_show_nonfits": str(self.reduced_show_nonfits.get_state()),
            "reduced_show_dirs": str(self.reduced_show_dirs.get_state()),
        }
        with open(config_file, "w") as f:
            config.write(f)
        self.logger.info(f"Saved default config to {config_file}")

    def open_backend_dialog(self) -> None:
        dialog = QtGui.QDialog()
        dialog.setWindowTitle("Backend Configuration")

        layout = QtGui.QVBoxLayout()

        local_checkbox = QtGui.QCheckBox("Use local backend")
        local_checkbox.setChecked(self.backend_mode == "local")
        layout.addWidget(local_checkbox)

        form_layout = QtGui.QFormLayout()
        host_edit = QtGui.QLineEdit(self.remote_host)
        port_edit = QtGui.QLineEdit(self.remote_port)
        key_edit = QtGui.QLineEdit(self.remote_api_key)
        key_edit.setPlaceholderText("Leave blank if server has no --api-key set")
        form_layout.addRow("Host:", host_edit)
        form_layout.addRow("Port:", port_edit)
        form_layout.addRow("API Key:", key_edit)
        layout.addLayout(form_layout)

        def _toggle_fields(checked: bool) -> None:
            host_edit.setEnabled(not checked)
            port_edit.setEnabled(not checked)
            key_edit.setEnabled(not checked)

        _toggle_fields(local_checkbox.isChecked())
        local_checkbox.toggled.connect(_toggle_fields)

        buttons = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Ok | QtGui.QDialogButtonBox.Cancel)
        layout.addWidget(buttons)

        dialog.setLayout(layout)
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)

        if dialog.exec_() != QtGui.QDialog.Accepted:
            return

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

            # Verify the server is reachable before switching backends.
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

    def raw_filter_cb(self, w, val):
        base_dir = self._get_tree_base_dir(self.state.raw_filepath)
        if base_dir:
            self._apply_file_filter(
                self.raw_treeview, base_dir,
                self.raw_show_fits.get_state(),
                self.raw_show_nonfits.get_state(),
                self.raw_show_dirs.get_state(),
                name_col_idx=self._raw_name_col_idx,
            )

    def reduced_filter_cb(self, w, val):
        base_dir = self._get_tree_base_dir(self.state.reduced_filepath)
        if base_dir:
            self._apply_file_filter(
                self.reduced_treeview, base_dir,
                self.reduced_show_fits.get_state(),
                self.reduced_show_nonfits.get_state(),
                self.reduced_show_dirs.get_state(),
                name_col_idx=self._reduced_name_col_idx,
            )

    def instrument_combo_cb(self, *args):
        selected = self.instrument_combo.get_text()
        self.instrument = self.instrument_registry.create(selected)
        self._rebuild_treeview_columns()

    def reduce_slit_cb(self, w):
        slit_key = self.slit_list_box.get_text().split()[0]
        if not slit_key:
            self.logger.error("No slit selected for reduction.")
            return

        raw_path = self.state.raw_filepath or self.raw_text_entry.get_text()
        reduced_path = self.state.reduced_filepath or self.reduced_text_entry.get_text()
        redux_path = self.redux_path_entry.get_text() or self.state.redux_path

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
            "--snr_thresh",
            self.SNR_box.get_text(),
            "--log_file",
            log_path,
            "-v",
            str(2)
        ]
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
        timer = self.reduction_timers.get(timer_key)

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
            if timer is not None:
                timer.cancel()
                self.reduction_timers.pop(timer_key, None)
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

        # Bottom line: status label, Show button, Remove button
        hbox_controls = GWidgets.HBox()
        label = GWidgets.Label(f"Reducing {slit_key}...")
        btn_show = GWidgets.Button("Show")
        btn_show.set_enabled(False)
        btn_remove = GWidgets.Button("Remove")

        hbox_controls.add_widget(label, stretch=1)
        hbox_controls.add_widget(btn_show, stretch=0)
        hbox_controls.add_widget(btn_remove, stretch=0)
        vbox.add_widget(hbox_controls, stretch=0)

        self.vbox_redux.add_widget(vbox, stretch=0)

        def _remove(w):
            self.vbox_redux.remove(vbox)
            self.reduction_control_elements.pop(slit_key, None)

        btn_remove.add_callback("activated", _remove)

        self.reduction_control_elements[slit_key] = {
            "label": label,
            "button": btn_show,
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

    def cadence_entry_cb(self, w):
        try:
            val = float(w.get_text())
            if val <= 0:
                raise ValueError
            self.reduction_cadence = val
        except ValueError:
            self.logger.error("Invalid cadence value; must be a positive number.")
            w.set_text(str(self.reduction_cadence))

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
                return
        self.reduced_btn.set_enabled(False)

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
            if hasattr(self, "reduced_show_fits"):
                self._apply_file_filter(
                    self.reduced_treeview,
                    self._get_tree_base_dir(self.state.reduced_filepath),
                    self.reduced_show_fits.get_state(),
                    self.reduced_show_nonfits.get_state(),
                    self.reduced_show_dirs.get_state(),
                    name_col_idx=self._reduced_name_col_idx,
                )
        else:
            self.state.raw_filepath = fullpath
            self.raw_treeview.set_tree(listing)
            if resize:
                self.raw_treeview.set_optimal_column_widths()
            if hasattr(self, "raw_show_fits"):
                self._apply_file_filter(
                    self.raw_treeview,
                    self._get_tree_base_dir(self.state.raw_filepath),
                    self.raw_show_fits.get_state(),
                    self.raw_show_nonfits.get_state(),
                    self.raw_show_dirs.get_state(),
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
        with fits.open(path) as hdul:
            img_data = self.instrument.get_display_image(hdul)
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
            self.raw_show_fits.set_state(defaults.getboolean("raw_show_fits", True))
            self.raw_show_nonfits.set_state(defaults.getboolean("raw_show_nonfits", False))
            self.raw_show_dirs.set_state(defaults.getboolean("raw_show_dirs", True))
            self.reduced_show_fits.set_state(defaults.getboolean("reduced_show_fits", False))
            self.reduced_show_nonfits.set_state(defaults.getboolean("reduced_show_nonfits", False))
            self.reduced_show_dirs.set_state(defaults.getboolean("reduced_show_dirs", True))
            self.logger.info(f"Loaded config from {config_file}")

        self.raw_text_entry.set_text(raw_path)
        self.reduced_text_entry.set_text(reduced_path)
        self.redux_path_entry.set_text(redux_path)

        self._browse_and_update(raw_path, which_tree="raw")
        self._browse_and_update(reduced_path, which_tree="reduced")
        self.instrument_combo.make_callback("activated")

        self.slit_canvas = DrawingCanvas()
        self.slit_canvas.enable_draw(False)
        self.slit_canvas.enable_edit(False)
        self.slit_canvas.set_surface(self.fitsimage)
        self.fitsimage.get_canvas().add(self.slit_canvas)

        self.fitsimage.add_callback("cursor-down", self.canvas_clicked_cb)

    def stop(self):
        for timer in self.reduction_timers.values():
            timer.cancel()
        self.reduction_timers.clear()
        if self.slit_canvas is not None:
            try:
                self.fitsimage.get_canvas().delete_object(self.slit_canvas)
            except Exception:
                pass
        self.gui_up = False

    def __str__(self):
        return "qltest"
