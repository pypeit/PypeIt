from __future__ import annotations

from typing import TYPE_CHECKING

from ginga.gw import Widgets

if TYPE_CHECKING:
    from .plugin import QLTEST


class QLViewUI:
    def __init__(self, plugin: "QLTEST") -> None:
        self.plugin = plugin

    def build(self, container) -> None:
        top = Widgets.VBox()
        top.set_border_width(4)

        vbox = Widgets.VBox()
        vbox.set_border_width(4)
        vbox.set_spacing(2)

        config_hbox = Widgets.HBox()
        config_hbox.add_widget(Widgets.Label("Instrument:"), stretch=0)
        self.plugin.instrument_combo = Widgets.ComboBox()
        for name in self.plugin.instrument_registry.names():
            self.plugin.instrument_combo.append_text(name)
        self.plugin.instrument_combo.add_callback("activated", self.plugin.instrument_combo_cb)
        config_hbox.add_widget(self.plugin.instrument_combo, stretch=0)

        vbox_show = Widgets.VBox()
        self.plugin.hide_reduced_tree = Widgets.CheckBox("Show Reduced Tree")
        self.plugin.hide_reduced_tree.set_state(True)
        self.plugin.hide_reduced_tree.add_callback("activated", self.plugin.hide_reduced_tree_cb)
        self.plugin.hide_raw_tree = Widgets.CheckBox("Show Raw Tree")
        self.plugin.hide_raw_tree.set_state(True)
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
        hbox_raw.add_widget(self.plugin.raw_text_entry, stretch=0)
        self.plugin.raw_btn = Widgets.Button("Go")
        self.plugin.raw_btn.add_callback("activated", self.plugin.raw_button_cb)
        hbox_raw.add_widget(self.plugin.raw_btn, stretch=0)
        fr_vbox.add_widget(hbox_raw, stretch=0)
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
        hbox_reduced.add_widget(self.plugin.reduced_text_entry, stretch=0)
        self.plugin.reduced_btn = Widgets.Button("Render Slits")
        self.plugin.reduced_btn.set_enabled(False)
        self.plugin.reduced_btn.add_callback("activated", self.plugin.render_slits_cb)
        hbox_reduced.add_widget(self.plugin.reduced_btn, stretch=0)
        self.plugin.show_wavelengths_btn = Widgets.Button("Show Wavelengths")
        self.plugin.show_wavelengths_btn.set_enabled(False)
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

        self.plugin.display_slits_box = Widgets.CheckBox("Display Slits")
        self.plugin.display_slits_box.set_state(True)
        self.plugin.display_slits_box.add_callback("activated", self.plugin.display_slits_box_cb)
        self.plugin.vbox_redux.add_widget(self.plugin.display_slits_box, stretch=0)

        self.plugin.show_labels_box = Widgets.CheckBox("Show Labels")
        self.plugin.show_labels_box.set_state(False)
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
