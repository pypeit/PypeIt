from __future__ import annotations

from typing import TYPE_CHECKING

from ginga.gw import Widgets

if TYPE_CHECKING:
    from .qlview import QLView


class QLViewUI:
    def __init__(self, plugin: "QLView") -> None:
        """
        Parameters
        ----------
        plugin : QLView
            The parent ``QLView`` Ginga plugin instance.  A reference is
            stored as ``self.plugin`` and used throughout :meth:`build` to
            attach widgets and callbacks.
        """
        self.plugin = plugin

    def build(self, container) -> None:
        """Construct the full Ginga widget tree and attach it to *container*.

        The method builds four logical sections and wires all callbacks:

        - **Raw Data frame** — tree view of raw files, path entry, "Go"
          button, and optional B-frame entry for A-B dithered sky subtraction.
        - **Reduced Calibrations frame** — tree view of reduced calibration
          products, path entry, "Render Slits", and "Show Wavelengths" buttons.
        - **Reduction Control frame** — slit selector, "Reduce Slit" button,
          SNR threshold, manual extraction toggle with FWHM entry, CoAdd2D
          checkbox, and slit/label display toggles.
        - **Button bar** — Close, Help, Save Default Config, and Settings
          buttons.

        All widgets are stored as attributes on ``self.plugin`` so that the
        plugin's callback methods can reference them directly.  Callbacks are
        registered here via ``add_callback``.

        Parameters
        ----------
        container : ginga widget
            The top-level container widget supplied by the Ginga plugin
            framework into which the entire UI is inserted.

        Returns
        -------
        None
        """
        top = Widgets.VBox()
        top.set_border_width(4)

        vbox = Widgets.VBox()
        vbox.set_border_width(4)
        vbox.set_spacing(2)

        config_hbox = Widgets.HBox()
        config_hbox.add_widget(Widgets.Label("Instrument:"), stretch=0)
        self.plugin.instrument_combo = Widgets.ComboBox()
        self.plugin.instrument_combo.set_tooltip(
            "Select the active instrument — controls FITS header keywords "
            "and file-browser column layout"
        )
        for name in self.plugin.instrument_registry.names():
            self.plugin.instrument_combo.append_text(name)
        self.plugin.instrument_combo.add_callback("activated", self.plugin.instrument_combo_cb)
        config_hbox.add_widget(self.plugin.instrument_combo, stretch=0)

        vbox_show = Widgets.VBox()
        self.plugin.hide_reduced_tree = Widgets.CheckBox("Show Reduced Tree")
        self.plugin.hide_reduced_tree.set_state(True)
        self.plugin.hide_reduced_tree.set_tooltip("Show or collapse the reduced calibrations file browser")
        self.plugin.hide_reduced_tree.add_callback("activated", self.plugin.hide_reduced_tree_cb)
        self.plugin.hide_raw_tree = Widgets.CheckBox("Show Raw Tree")
        self.plugin.hide_raw_tree.set_state(True)
        self.plugin.hide_raw_tree.set_tooltip("Show or collapse the raw data file browser")
        self.plugin.hide_raw_tree.add_callback("activated", self.plugin.hide_raw_tree_cb)
        vbox_show.add_widget(self.plugin.hide_raw_tree, stretch=0)
        vbox_show.add_widget(self.plugin.hide_reduced_tree, stretch=0)

        config_hbox.add_widget(vbox_show, stretch=0)
        vbox.add_widget(config_hbox, stretch=0)

        color_alternate = self.plugin.settings.get("color_alternate_rows", True)

        self.plugin.reduced_treeview = Widgets.TreeView(
            sortable=True,
            selection="multiple",
            use_alt_row_color=color_alternate,
            dragable=True,
        )
        self.plugin.raw_treeview = Widgets.TreeView(
            sortable=True,
            selection="multiple",
            use_alt_row_color=color_alternate,
            dragable=True,
        )

        # --- Raw Data frame ---
        fr = Widgets.Frame("Raw Data")
        fr_vbox = Widgets.VBox()
        raw_cols = self.plugin.instrument.columns["raw"]
        self.plugin.raw_treeview.setup_table(raw_cols, 1, "name")
        self.plugin.raw_treeview.add_callback("selected", self.plugin.raw_table_selected_cb)
        self.plugin.raw_treeview.add_callback("activated", self.plugin.raw_table_double_click_cb)
        fr_vbox.add_widget(self.plugin.raw_treeview, stretch=1)

        hbox_raw = Widgets.HBox()
        hbox_raw.add_widget(Widgets.Label("Raw Data Path:"), stretch=0)
        self.plugin.raw_text_entry = Widgets.TextEntry()
        self.plugin.raw_text_entry.set_tooltip(
            "Directory or file to browse. Supports strftime format codes "
            "(e.g. /data/raw/%Y%m%d is expanded to today's date at startup)."
        )
        hbox_raw.add_widget(self.plugin.raw_text_entry, stretch=0)
        self.plugin.raw_btn = Widgets.Button("Go")
        self.plugin.raw_btn.set_tooltip(
            "Browse to the entered directory, or open the file if a FITS path is given"
        )
        self.plugin.raw_btn.add_callback("activated", self.plugin.raw_button_cb)
        hbox_raw.add_widget(self.plugin.raw_btn, stretch=0)
        fr_vbox.add_widget(hbox_raw, stretch=0)

        hbox_bframe = Widgets.HBox()
        hbox_bframe.add_widget(Widgets.Label("B frame (optional):"), stretch=0)
        self.plugin.b_frame_entry = Widgets.TextEntry()
        self.plugin.b_frame_entry.set_tooltip(
            "Path to the B-position frame for A-B dithered sky subtraction."
        )
        self.plugin.b_frame_entry.add_callback("activated", self.plugin.b_frame_entry_cb)
        hbox_bframe.add_widget(self.plugin.b_frame_entry, stretch=1)
        btn_detect_ab = Widgets.Button("Detect AB Pair")
        btn_detect_ab.set_tooltip(
            "Scan the raw directory for the best-matching B-position frame"
        )
        btn_detect_ab.add_callback("activated", self.plugin.detect_ab_pair_cb)
        hbox_bframe.add_widget(btn_detect_ab, stretch=0)
        btn_clear_b = Widgets.Button("Clear")
        btn_clear_b.set_tooltip("Clear the B frame path")
        btn_clear_b.add_callback("activated", lambda w: self.plugin.clear_b_frame_cb(w))
        hbox_bframe.add_widget(btn_clear_b, stretch=0)
        fr_vbox.add_widget(hbox_bframe, stretch=0)

        fr.set_widget(fr_vbox)
        vbox.add_widget(fr, stretch=0)

        # --- Reduced Calibrations frame ---
        fr = Widgets.Frame("Reduced Calibrations")
        fr_vbox = Widgets.VBox()

        # Calibration-suggestion status label
        self.plugin.cal_status_label = Widgets.Label("")
        self.plugin.cal_status_label.set_tooltip(
            "Status of the automatic calibration recommendation"
        )
        fr_vbox.add_widget(self.plugin.cal_status_label, stretch=0)

        reduced_cols = self.plugin.instrument.columns["reduced"]
        self.plugin.reduced_treeview.setup_table(reduced_cols, 1, "name")
        self.plugin.reduced_treeview.add_callback("selected", self.plugin.reduced_table_selected_cb)
        self.plugin.reduced_treeview.add_callback("activated", self.plugin.reduced_table_double_click_cb)
        fr_vbox.add_widget(self.plugin.reduced_treeview, stretch=1)

        hbox_reduced = Widgets.HBox()
        hbox_reduced.add_widget(Widgets.Label("Reduced Cals Path:"), stretch=0)
        self.plugin.reduced_text_entry = Widgets.TextEntry()
        self.plugin.reduced_text_entry.set_tooltip(
            "Path to a PypeIt reduction output directory that contains a Calibrations/ subdirectory"
        )
        hbox_reduced.add_widget(self.plugin.reduced_text_entry, stretch=0)
        self.plugin.reduced_btn = Widgets.Button("Render Slits")
        self.plugin.reduced_btn.set_enabled(False)
        self.plugin.reduced_btn.set_tooltip(
            "Load Slits_*.fits* from the selected Calibrations/ directory "
            "and draw slit polygons over the raw image"
        )
        self.plugin.reduced_btn.add_callback("activated", self.plugin.render_slits_cb)
        hbox_reduced.add_widget(self.plugin.reduced_btn, stretch=0)
        self.plugin.show_wavelengths_btn = Widgets.Button("Show Wavelengths")
        self.plugin.show_wavelengths_btn.set_enabled(False)
        self.plugin.show_wavelengths_btn.set_tooltip(
            "Build a 2D wavelength map from WaveCalib/Tilts/Slits files — "
            "hover over the image to see the wavelength at each pixel"
        )
        self.plugin.show_wavelengths_btn.add_callback("activated", self.plugin.show_wavelengths_cb)
        hbox_reduced.add_widget(self.plugin.show_wavelengths_btn, stretch=0)
        fr_vbox.add_widget(hbox_reduced, stretch=0)
        fr.set_widget(fr_vbox)
        vbox.add_widget(fr, stretch=0)

        # --- Reduction Control frame ---
        fr = Widgets.Frame("Reduction Control")
        self.plugin.vbox_redux = Widgets.VBox()

        hbox = Widgets.HBox()
        self.plugin.slit_list_box = Widgets.ComboBox()
        self.plugin.slit_list_box.set_tooltip(
            "Select the slit to reduce or highlight — you can also click directly on a slit in the image"
        )
        self.plugin.slit_list_box.add_callback("activated", self.plugin.slit_list_box_cb)
        hbox.add_widget(self.plugin.slit_list_box)
        self.plugin.btn_reduce = Widgets.Button("Reduce Slit")
        self.plugin.btn_reduce.set_tooltip("Reduce the selected slit")
        self.plugin.btn_reduce.add_callback("activated", self.plugin.reduce_slit_cb)
        hbox.add_widget(self.plugin.btn_reduce, stretch=0)

        self.plugin.SNR_box = Widgets.TextEntry()
        self.plugin.SNR_box.set_text("25")
        self.plugin.SNR_box.set_tooltip("SNR Threshold for object extraction")
        label = Widgets.Label("SNR Threshold:")
        label.set_tooltip("SNR Threshold for object extraction")
        hbox.add_widget(label, stretch=0)
        hbox.add_widget(self.plugin.SNR_box, stretch=0)
        self.plugin.vbox_redux.add_widget(hbox, stretch=0)

        # Manual extraction row: toggle checkbox + FWHM entry
        hbox_manual = Widgets.HBox()
        self.plugin.manual_extract_checkbox = Widgets.CheckBox("Set Manual Extraction")
        self.plugin.manual_extract_checkbox.set_state(False)
        self.plugin.manual_extract_checkbox.set_tooltip(
            "Enable click-to-extract mode: click on the raw image to place a manual extraction aperture"
        )
        self.plugin.manual_extract_checkbox.add_callback(
            "activated", self.plugin.manual_extract_mode_cb
        )
        hbox_manual.add_widget(self.plugin.manual_extract_checkbox, stretch=0)

        fwhm_label = Widgets.Label("FWHM (px):")
        fwhm_label.set_tooltip("Extraction FWHM in pixels for the manual aperture")
        hbox_manual.add_widget(fwhm_label, stretch=0)
        self.plugin.fwhm_box = Widgets.TextEntry()
        self.plugin.fwhm_box.set_text("3.0")
        self.plugin.fwhm_box.set_tooltip("Extraction FWHM in pixels")
        self.plugin.fwhm_box.add_callback("activated", self.plugin.fwhm_box_changed_cb)
        hbox_manual.add_widget(self.plugin.fwhm_box, stretch=0)
        self.plugin.vbox_redux.add_widget(hbox_manual, stretch=0)

        # Manual extraction params display/edit row
        hbox_params = Widgets.HBox()
        params_label = Widgets.Label("Extract params:")
        params_label.set_tooltip("Manual extraction string passed to --manual_extract")
        hbox_params.add_widget(params_label, stretch=0)
        self.plugin.manual_extract_params_entry = Widgets.TextEntry()
        self.plugin.manual_extract_params_entry.set_text("")
        self.plugin.manual_extract_params_entry.set_tooltip(
            "det:spat:spec:fwhm — updates on each click or FWHM change; editable directly"
        )
        self.plugin.manual_extract_params_entry.add_callback(
            "activated", self.plugin.manual_extract_params_cb
        )
        hbox_params.add_widget(self.plugin.manual_extract_params_entry, stretch=1)
        self.plugin.vbox_redux.add_widget(hbox_params, stretch=0)

        self.plugin.coadd2d_box = Widgets.CheckBox("CoAdd2D")
        self.plugin.coadd2d_box.set_state(False)
        self.plugin.coadd2d_box.set_tooltip(
            "Pass --coadd2d to ql.py: co-add the 2D spectra after reduction "
            "(required for dithered A-B pairs with multiple frames)"
        )
        self.plugin.vbox_redux.add_widget(self.plugin.coadd2d_box, stretch=0)

        self.plugin.display_slits_box = Widgets.CheckBox("Display Slits")
        self.plugin.display_slits_box.set_state(True)
        self.plugin.display_slits_box.set_tooltip("Show or hide slit polygon outlines on the raw image")
        self.plugin.display_slits_box.add_callback("activated", self.plugin.display_slits_box_cb)
        self.plugin.vbox_redux.add_widget(self.plugin.display_slits_box, stretch=0)

        self.plugin.show_labels_box = Widgets.CheckBox("Show Labels")
        self.plugin.show_labels_box.set_state(False)
        self.plugin.show_labels_box.set_tooltip(
            "Overlay slit ID and object name labels on each slit polygon"
        )
        self.plugin.show_labels_box.add_callback("activated", self.plugin.show_labels_box_cb)
        self.plugin.vbox_redux.add_widget(self.plugin.show_labels_box, stretch=0)

        fr.set_widget(self.plugin.vbox_redux)
        vbox.add_widget(fr, stretch=0)

        top.add_widget(vbox, stretch=0)
        top.add_widget(Widgets.Label(""), stretch=1)

        btns = Widgets.HBox()
        btns.set_spacing(3)
        btn = Widgets.Button("Close")
        btn.add_callback("activated", lambda w: self.plugin.close())
        btns.add_widget(btn, stretch=0)
        btn = Widgets.Button("Help")
        btn.add_callback("activated", lambda w: self.plugin.help())
        btns.add_widget(btn, stretch=0)
        btn = Widgets.Button("Save Default Config")
        btn.set_tooltip("Create a config file with the current settings")
        btn.add_callback("activated", lambda w: self.plugin.create_config_cb(w))
        btns.add_widget(btn, stretch=0)
        btn = Widgets.Button("Settings")
        btn.set_tooltip("Configure backend, reduction path, poll cadence, and file filters")
        btn.add_callback("activated", lambda w: self.plugin.open_settings_dialog())
        btns.add_widget(btn, stretch=0)
        btns.add_widget(Widgets.Label(""), stretch=1)
        top.add_widget(btns, stretch=0)

        container.add_widget(top, stretch=1)
