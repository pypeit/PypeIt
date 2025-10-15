"""
Spec1dView is a plugin for the Ginga viewer that provides functionality for
visualizing and analyzing 1D spectra from FITS files. The plugin allows users
to plot spectra, identify spectral lines from various line lists, and
customize the display according to different parameters.

**Plugin Type: Local**

Spec1dView is a local plugin, which means it is associated with a specific
channel in the Ginga viewer. An instance of the plugin can be opened for
each channel, allowing for multiple spectra to be analyzed simultaneously.

**Usage**
- Load and visualize 1D spectra from FITS files.
- Customize the display by selecting different line lists, extraction types,
  and flux/mask settings.
- Update the redshift to shift the spectral lines accordingly.

**Editing**
Users can modify the visualization by:
- Choosing from a variety of line lists to identify spectral features.
- Selecting different types of extraction methods (OPT, BOX).
- Applying or removing flux calibration and masking options.
- Updating the redshift value to reflect the observed wavelengths.

**UI**
The user interface provides controls for:
- Selecting the line list from a combobox.
- Entering a redshift value to shift the spectrum.
- Choosing the extraction type, flux calibration, and masking options via comboboxes.
- Buttons to load a FITS file and clear the current selection.

**Buttons**
- Update z: Updates the redshift value and refreshes the spectrum plot.
- Enter: Loads the specified FITS file for analysis.
- Clear: Clears the current inputs and resets the UI settings.

**Tips**
- Use the comboboxes to switch between different line lists and adjust the spectrum display settings.
- Ensure that the correct FITS file path is entered before attempting to load the data.
"""
import time
import glob
import os
import numpy as np
from pathlib import Path
import subprocess
import multiprocessing
import threading
import asyncio
from ginga import GingaPlugin
from ginga.misc import Bunch
from ginga.gw import Widgets
from ginga.table.AstroTable import AstroTable
from ginga.AstroImage import AstroImage
from ginga.plot.Plotable import Plotable
from ginga.canvas.CanvasObject import get_canvas_types
from ginga.canvas.types.basic import Polygon

from ginga.qtw.QtHelp import QtGui, QtCore


from pypeit import specobjs
from pypeit import utils
from pypeit.slittrace import SlitTraceSet
from pypeit.scripts import ql

from astropy.io import fits
# TODO: need to make PypeIt LineList class to deprecate linetools
from linetools.lists.linelist import LineList

__all__ = ['QLView']


class QLView(GingaPlugin.LocalPlugin):

    def __init__(self, fv, fitsimage):

        """Constructor for the plugin."""
        # superclass defines some variables for us, like logger
        super().__init__(fv, fitsimage)

        # get QLView preferences
        keywords = [('Object', 'OBJECT'),
                    ('Date', 'DATE-OBS'),
                    ('Time UT', 'UT')]
        # columns = [('Type', 'icon'),
        #            ('Name', 'name'),
        #            ('Size', 'st_size_str'),
        #            ('Mode', 'st_mode_oct'),
        #            ('Last Changed', 'st_mtime_str')]
        columns = [
            ('Type', 'icon'),
            ('Frame No', 'FRAMENO'),
            ('Name', 'name'),
            ('Object', 'OBJECT'),
            ('Type', 'IMTYPE'),
            ('Mask Name', 'MASKNAME'),
            ('Obs Mode', 'OBSMODE'),
            ('Exposure Time', 'EXPTIME'),
            ('Last Changed', 'st_mtime_str')
        ]
        
        prefs = self.fv.get_preferences()
        self.settings = prefs.create_category('plugin_QLView')
        self.settings.add_defaults( scan_fits_headers=False,
                                    scan_limit=100,
                                    keywords=keywords,
                                    columns=columns,
                                    color_alternate_rows=True,
                                    max_rows_for_col_resize=5000)
        self.settings.load(onError='silent')

        self.scan_limit = self.settings.get('scan_limit', 100)
        self.keywords = self.settings.get('keywords', keywords)
        self.columns = self.settings.get('columns', columns)

        self.raw_filepath = None
        self.reduced_filepath = None
        self.redux_path = "/hqdrpdata/DRP_TESTING/outputs"
        self.reduction_control_elements = {}
        # self.data_source = LocalDataSource(self.logger, "DEIMOS")


        # Make icons
        icondir = self.fv.iconpath
        self.folderpb = self.fv.get_icon(icondir, 'folder.svg')
        self.filepb = self.fv.get_icon(icondir, 'file.svg')
        self.fitspb = self.fv.get_icon(icondir, 'fits.svg')

        
        # dictionary of plotable types
        self.dc = get_canvas_types()
        # self.plot = None
        self.gui_up = False
        self.slit_canvas = None
        self.slittracesets = None
        self.active_slit = None
        self.slit_polys = {}

        # Multiprocessing pool for running reductions
        self.pool = multiprocessing.Pool(processes=2)

        self.instrument = DEIMOS(self.logger) # Gotta start somewhere

    def build_gui(self, container):
        """Construct the UI in the plugin container.

        This get's called just prior to the ``start`` method when the plugin is
        activated.
        """
        top = Widgets.VBox()
        top.set_border_width(4)

        vbox = Widgets.VBox()
        vbox.set_border_width(4)
        vbox.set_spacing(2)


        config_hbox = Widgets.HBox()
        config_hbox.add_widget(Widgets.Label("Instrument:"), stretch=0)
        self.instrument_combo = Widgets.ComboBox()
        self.instrument_combo.append_text("DEIMOS")
        self.instrument_combo.append_text("MOSFIRE")
        self.instrument_combo.add_callback('activated', self.instrument_combo_cb)
        config_hbox.add_widget(self.instrument_combo, stretch=0)

        vbox_show = Widgets.VBox()
        self.hide_reduced_tree = Widgets.CheckBox("Show Reduced Tree")
        self.hide_reduced_tree.set_state(True)
        self.hide_reduced_tree.add_callback('activated', self.hide_reduced_tree_cb)
        self.hide_raw_tree = Widgets.CheckBox("Show Raw Tree")
        self.hide_raw_tree.set_state(True)
        self.hide_raw_tree.add_callback('activated', self.hide_raw_tree_cb)
        vbox_show.add_widget(self.hide_raw_tree, stretch=0)
        vbox_show.add_widget(self.hide_reduced_tree, stretch=0)

        config_hbox.add_widget(vbox_show, stretch=0)
        vbox.add_widget(config_hbox, stretch=0)

        # # # # # # # # # #
        # File selection  #
        # # # # # # # # # #

        # configure the table:
        color_alternate = self.settings.get('color_alternate_rows', True)

        # Create the file tree for the directory with the reduced calibrations
        self.reduced_treeview = Widgets.TreeView(sortable=True, selection='multiple',
                                 use_alt_row_color=color_alternate,
                                 dragable=True)
        
        # Similarly, create the file tree for the directory with the raw data
        self.raw_treeview = Widgets.TreeView(sortable=True, selection='multiple',
                            use_alt_row_color=color_alternate,
                            dragable=True)

        # Set the column headers for both file tree tables
        # col = 0
        # self._name_idx = 0
        # for hdr, attrname in self.settings.get('columns'):
        #     if attrname == 'name':
        #         self._name_idx = col
        #     col += 1

        # Finalize the raw file tree and add it to the UI
        fr = Widgets.Frame("Raw Data")
        fr_vbox = Widgets.VBox()
        self.raw_treeview.setup_table(self.settings.get('columns'), 1, 'name')
        self.raw_treeview.add_callback('selected', self.raw_table_selected_cb)
        self.raw_treeview.add_callback('activated', self.raw_table_double_click_cb)
        fr_vbox.add_widget(self.raw_treeview, stretch=1)
        
        hbox_raw = Widgets.HBox()
        hbox_raw.add_widget(Widgets.Label("Raw Data Path:"), stretch=0)
        self.raw_text_entry = Widgets.TextEntry()
        hbox_raw.add_widget(self.raw_text_entry, stretch=0)
        self.raw_btn = Widgets.Button("Go")
        self.raw_btn.add_callback('activated', self.raw_button_cb)
        hbox_raw.add_widget(self.raw_btn, stretch=0)
        fr_vbox.add_widget(hbox_raw, stretch=0)
        fr.set_widget(fr_vbox)
        vbox.add_widget(fr, stretch=0)

        # Finalize the reduced file tree and add it to the UI
        fr = Widgets.Frame("Reduced Calibrations")
        fr_vbox = Widgets.VBox()       
        self.reduced_treeview.setup_table(self.settings.get('columns'), 1, 'name')
        self.reduced_treeview.add_callback('selected', self.reduced_table_selected_cb)
        self.reduced_treeview.add_callback('activated', self.reduced_table_double_click_cb)
        fr_vbox.add_widget(self.reduced_treeview, stretch=1)
        
        hbox_reduced = Widgets.HBox()
        hbox_reduced.add_widget(Widgets.Label("Reduced Cals Path:"), stretch=0)
        self.reduced_text_entry = Widgets.TextEntry()
        hbox_reduced.add_widget(self.reduced_text_entry, stretch=0)
        self.reduced_btn = Widgets.Button("Render Slits")
        self.reduced_btn.set_enabled(False)
        self.reduced_btn.add_callback('activated', self.render_slits_cb)
        hbox_reduced.add_widget(self.reduced_btn, stretch=0)
        fr_vbox.add_widget(hbox_reduced, stretch=0)
        fr.set_widget(fr_vbox)
        vbox.add_widget(fr, stretch=0)


        # # # # # # # # # # #
        # Reduction Control #
        # # # # # # # # # # #
        fr = Widgets.Frame("Reduction Control")
        self.vbox_redux = Widgets.VBox()

        hbox = Widgets.HBox()
        self.redux_path_entry = Widgets.TextEntry()
        self.redux_path_entry.set_text(self.redux_path)
        self.redux_path_entry.set_tooltip("Path to the reduction directory")
        hbox.add_widget(Widgets.Label("Reduction Path:"), stretch=0)
        hbox.add_widget(self.redux_path_entry, stretch=1)
        self.vbox_redux.add_widget(hbox, stretch=0)

        hbox = Widgets.HBox()
        self.slit_list_box = Widgets.ComboBox()
        self.slit_list_box.add_callback('activated', self.slit_list_box_cb)
        hbox.add_widget(self.slit_list_box)
        self.btn_reduce = Widgets.Button("Reduce Slit")
        self.btn_reduce.set_tooltip("Reduce the selected slit")
        self.btn_reduce.add_callback('activated', self.reduce_slit_cb)
        hbox.add_widget(self.btn_reduce, stretch=0)

        self.SNR_box = Widgets.TextEntry()
        self.SNR_box.set_text("25")
        self.SNR_box.set_tooltip("SNR Threshold for object extraction")
        label = Widgets.Label("SNR Threshold:")
        label.set_tooltip("SNR Threshold for object extraction")
        hbox.add_widget(label, stretch=0)
        hbox.add_widget(self.SNR_box, stretch=0)

        self.vbox_redux.add_widget(hbox, stretch=0)

        self.display_slits_box = Widgets.CheckBox("Display Slits")
        self.display_slits_box.set_state(True)
        self.display_slits_box.add_callback('activated', self.display_slits_box_cb)
        self.vbox_redux.add_widget(self.display_slits_box, stretch=0)



        fr.set_widget(self.vbox_redux)
        vbox.add_widget(fr, stretch=0)

        top.add_widget(vbox, stretch=0)

        spacer = Widgets.Label('')
        top.add_widget(spacer, stretch=1)

        btns = Widgets.HBox()
        btns.set_spacing(3)

        btn = Widgets.Button("Close")
        btn.add_callback('activated', lambda w: self.close())
        btns.add_widget(btn, stretch=0)
        btn = Widgets.Button("Help")
        btn.add_callback('activated', lambda w: self.help())
        btns.add_widget(btn, stretch=0)
        btn = Widgets.Button("Save Default Config")
        btn.set_tooltip("Create a config file with the current settings")
        btn.add_callback('activated', lambda w: self.create_config_cb())
        btns.add_widget(btn, stretch=0)
        btns.add_widget(Widgets.Label(''), stretch=1)
        top.add_widget(btns, stretch=0)

        container.add_widget(top, stretch=1)


        self.gui_up = True

    # # # # # # #
    # Callbacks #
    # # # # # # #

    def create_config_cb(self):
        """Callback for the "Create Config" button.
        
        This method is called when the user clicks the "Create Config" button.
        It creates a configuration file with the current settings and saves it
        to the user's home directory.
        """
        config_file = Path.home() / ".quicklook.cfg"
        self.logger.info(f"Creating config file: {config_file}")
        import configparser
        config = configparser.ConfigParser()
        config["DEFAULT"] = {
            "redux_path": self.redux_path_entry.get_text(),
            "raw_path": self.raw_text_entry.get_text(),
            "reduced_path": self.reduced_text_entry.get_text()
        }
        with open(config_file, 'w') as f:
            config.write(f)

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

    def reduce_slit_cb(self, w):
        # Launch the reduction using popen:

        msc = self.slittracesets
        for msc_idx in self.slittracesets.keys():
            
            slittrace = self.slittracesets[msc_idx]
            spatial_ids = slittrace.spat_id

            for idx, spat_id in enumerate(spatial_ids):
                if self.slit_list_box.get_text() == f"S{spat_id}":
                    msc = f"MSC{msc_idx}"

        command = ["pypeit_ql"]
        command.append(self.instrument.pypeit_name)
        command.append("--raw_files")
        command.append(str(Path(self.raw_filepath).name))
        command.append("--raw_path")
        command.append(str(Path(self.raw_filepath).parent.absolute()))
        command.append("--setup_calib_dir")
        command.append(f"{Path(self.reduced_filepath).absolute()}/Calibrations")
        command.append("--slitspatnum")
        slitspatnum = f"{msc}:{self.slit_list_box.get_text()[1:]}"
        command.append(slitspatnum)
        command.append("--redux_path")
        command.append(self.redux_path_entry.get_text())
        command.append("--skip_display")
        command.append("--snr_thresh")
        command.append(self.SNR_box.get_text())

        self.logger.info("Launching command: {0}".format(" ".join(command)))

        old_method = False
        if old_method:
            # popen_w_cb(self.show_reduced_spec, command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            logfile = open(os.path.join(self.redux_path_entry.get_text(), f"{slitspatnum}.log"), "w")
            p = subprocess.Popen(command, stdout=logfile, stderr=logfile)

            print("Launched reduction")

            hbox, elements = self.make_reduced_slit_hbox()
            elements['name'].set_text(f"Reducing {self.slit_list_box.get_text()}...")

            self.vbox_redux.add_widget(hbox, stretch=0)
            self.reduction_control_elements[self.slit_list_box.get_text()] = {"elements": elements, "hbox": hbox}
            
            # need to sleep for just long enough for the directory to be created
            # redux_path = self.redux_path_from_filename(self.raw_filepath)
            # while not os.path.exists(redux_path):
                # self.logger.info(f"Waiting for {redux_path} to be created...")
                # time.sleep(.5)
            # self.add_dirs_to_watcher(self.redux_path_from_filename(self.raw_filepath))

            # out, err = p.communicate()

            # while p.poll() is None:
            #     time.sleep(0.1)
        else:
            print("Using new method to launch reduction")
            # args = ql.QL.get_parser().parse_args(command[1:])
            # print(f"Args to pass to QL.main: {args}")
            # def a_cool_callback(out):
            #     self.logger.info("Reduction complete as seen in callback~!")
            #     self.raw_btn.set_enabled(False)
            # # self.pool.apply_async(ql.QL.main, args, callback=a_cool_callback)
            # # ql.QL.main(args)
            print(self.logger.name)
            self.pool.apply(sleepy, args=())
            self.pool.join()
            hbox, elements = self.make_reduced_slit_hbox()
            elements['name'].set_text(f"Reducing {self.slit_list_box.get_text()}...")

            self.vbox_redux.add_widget(hbox, stretch=0)
            self.reduction_control_elements[self.slit_list_box.get_text()] = {"elements": elements, "hbox": hbox}


    def slit_list_box_cb(self, w, res_dict):
        self.deactivate_slit()
        slit = self.slit_polys[w.get_text()]
        self.active_slit = slit
        self.activate_slit()


    def display_slits_box_cb(self, w, val):
        """Callback for the "Display Slits" checkbox in the Reduction Control frame.
        
        This method is called when the user toggles the "Display Slits" checkbox.
        """
        if self.slit_canvas:
            self.slit_canvas.viewer.viewable(val)
            tag = self.fitsimage.get_canvas().lookup_object_tag(self.slit_canvas)
            if not val:
                self.fitsimage.get_canvas().delete_object_by_tag(tag)
            else:
                self.fitsimage.get_canvas().add(self.slit_canvas)

    def reduced_table_double_click_cb(self, w, res_dict):
        """Callback for double-clicking on a row in the REDUCED data file tree.
        
        Double-clicking on a row in the table will attempt to load the FITS file,
        process it, and display it in the viewer, along with attempting to load
        the slit edges for overlay.
        """
        self.logger.debug(f"table_double_click: {res_dict}")
        
        paths = [info.path for key, info in res_dict.items()]
        path = paths[0]
        if os.path.isdir(path):
            self.browse(path, "reduced")
        elif os.path.isfile(path):
            self.logger.info("File path entered, but we don't open those here")
            pass
        else:
            self.logger.error("Invalid path entered.")
            self.reduced_text_entry.set_text("Invalid path")
        # Do whatever you need with the LocalDataSource object

        # TODO: Handle the double-click event to load slit edges

    def reduced_table_selected_cb(self, w, res_dict):
        """Callback for selecting a row in the REDUCED data file tree.
        
        Selecting a row in the table will update the displayed information in the
        plugin UI.
        """
        self.logger.debug(f"table_selected: {res_dict}")
        paths = [info.path for key, info in res_dict.items()]
        if not paths:
            return
        path = paths[0]

        # If the selected item is a PypeIt directory, enable the "Load Slits" button
        self.reduced_text_entry.set_text(path)

        if "keck_" in self.reduced_text_entry.get_text():
            self.logger.debug("We have (presumably) entered a PypeIt reduced directory.")
            self.reduced_filepath = path
            p = Path(path)
            subdirs = [x for x in p.iterdir() if x.is_dir()]
            if "Calibrations" in [x.name for x in subdirs]:
                self.logger.info("Calibrations directory found.")
                self.reduced_btn.set_enabled(True)
        else:
            self.logger.debug("Not a PypeIt reduced directory, disabling render_slits button.")
            self.reduced_btn.set_enabled(False)

        # Don't load the file, but do display some metadata!

        # TODO: Display the metadata for the selected file in the info box
    

    def raw_table_double_click_cb(self, w, res_dict):
        """Callback for double-clicking on a row in RAW data file tree.
        
        Double-clicking on a row in the table will attempt to load the FITS file,
        process it, and display it in the viewer, along with attempting to load
        the slit edges for overlay.
        """
        self.logger.debug(f"table_double_click: {res_dict}")
        paths = [info.path for key, info in res_dict.items()]
        path = paths[0]
        if os.path.isdir(path):
            self.browse(path, "raw")
        elif os.path.isfile(path):
            self.open_raw_file(path)

        else:
            self.logger.error("Invalid path entered.")
            self.raw_text_entry.set_text("Invalid path")
        self.raw_filepath = path


        # TODO: Handle the double-click event to load the FITS file

    def raw_table_selected_cb(self, w, res_dict):
        """Callback for selecting a row in the RAW data file tree.
        
        Selecting a row in the table will update the displayed information in the
        plugin UI.
        """
        self.logger.debug(f"table_selected: {res_dict}")
        paths = [info.path for key, info in res_dict.items()]
        if not paths:
            return
        path = paths[0]
        self.raw_text_entry.set_text(path)
        self.raw_filepath = path

        # Don't load the file, but do display some metadata!

        # TODO: Display the metadata for the selected file in the info box

    def render_slits_cb(self, w):
        """Callback for reduced data directory text entry widget. This will render
        slits on the viewer, if the user has selected a valid directory.
        
        This method is called when the user enters text into a text entry widget.
        """
        self.logger.debug(f"Render Slits called with reduced path {self.reduced_filepath}")

        self.slittracesets = self.open_slits_files(os.path.join(self.reduced_filepath, "Calibrations"))
        
        if self.fitsimage.get_canvas() is not None:
            if self.display_slits_box.get_state():
                try:
                    self.logger.info("Trying to remove old slits")
                    tag = self.fitsimage.get_canvas().lookup_object_tag(self.slit_canvas)
                    self.logger.info(f"Tag: {tag}")
                    self.fitsimage.get_canvas().delete_object_by_tag(tag)
                    self.logger.info("Removed old slits (I think...)")
                except KeyError:
                    self.logger.info("No slits found to remove")
            self.slit_canvas = self.construct_slits(self.slittracesets)
            self.fitsimage.get_canvas().add(self.slit_canvas)
        else:
            self.logger.error("No image has been rendered yet. Cannot add slits.")

    def raw_button_cb(self, w):
        """Callback for raw data directory text entry widget.
        
        This method is called when the user enters text into a text entry widget.
        """
        self.logger.debug(f"text_entry: {w}")
        if os.path.isdir(self.raw_text_entry.get_text()):
            self.browse(self.raw_text_entry.get_text(), "raw")
        elif os.path.isfile(self.raw_text_entry.get_text()):
            self.logger.info("File path entered. Loading file.")
            self.open_raw_file(self.raw_text_entry.get_text())
        else:
            self.logger.error("Invalid path entered.")
            self.raw_text_entry.set_text("Invalid path")
    

    def instrument_combo_cb(self, *args):
        selected_inst_str = self.instrument_combo.get_text()
        match selected_inst_str:
            case "DEIMOS":
                    self.instrument = DEIMOS(self.logger)
            case "MOSFIRE":
                    self.instrument = MOSFIRE(self.logger)
            case _:
                    self.logger.error(f"Instrument not recognized: {selected_inst_str}")

    def canvas_clicked_cb(self, canvas, pnt, x, y):
        self.logger.info("Canvas clicked! X: {0}, Y: {1}".format(x, y))

        if self.slittracesets is None:
            self.logger.error("No slits have been loaded yet.")
            return
        for msc_idx in self.slittracesets.keys():
            slits = self.slittracesets[msc_idx]
            if slits is None:
                continue
            offset = (int(msc_idx) - 1) * slits.nspat
            left_bound_at_y = slits.left_init[np.round(-1 * y + (2048 * 4)).astype(int)] + offset
            right_bound_at_y = slits.right_init[np.round(-1 * y + (2048 * 4)).astype(int)] + offset

            for i in range(slits.nslits):
                if left_bound_at_y[i] < x < right_bound_at_y[i]:
                    slit_id = slits.spat_id[i]
                    self.logger.info("Found slit {0}".format(slit_id))
                    self.slit_list_box.show_text(f"S{slit_id}")
                    try:
                        self.deactivate_slit()
                        self.active_slit = self.slit_polys[f"S{slit_id}"]
                        self.activate_slit()
                    except KeyError:
                        self.logger.error(f"Slit {slit_id} not found in rendered slits!")
                    break

    # # # # # # #
    # File Tree #
    # # # # # # #

    def update_tree_and_text(self):
        """Update the treeview and text entry widgets in the plugin.

        This method is called when the user loads a new FITS file or changes the
        selection in the treeview.
        """
        self.reduced_treeview.set_tree(self.tree_dict)
        self.reduced_text_entry.set_text(self.filepath)

    def browse(self, path, which_tree):
        """Browse the directory at the given path and update the file tree."""   
        if which_tree not in ["reduced", "raw"]:
            raise ValueError(f"Invalid tree type: {which_tree}")
        
        self.logger.info(f"Browse {which_tree} path: {path}")
        
        if os.path.isdir(path):
            dirname = path
            globname = None
        else:
            # dirname, globname = os.path.split(path)
            self.logger.error("Attempting to browse a file path")
            if which_tree == "reduced":
                self.reduced_text_entry.set_text("Not a valid path")
            else:
                self.raw_text_entry.set_text("Not a valid path")
            return
        dirname = os.path.abspath(dirname)

        if not globname:
            globname = '*'
        path = os.path.join(dirname, globname)

        # Make a directory listing
        self.logger.debug("globbing path: %s" % (path))
        filelist = list(glob.glob(path))
        filelist.sort(key=lambda s: s.lower())
        filelist.insert(0, os.path.join(dirname, '..'))

        jumpinfo = list(map(self.get_info, filelist))
        if which_tree == "reduced":
            self.reduced_filepath = path
        else:
            self.raw_filepath = path
        # if self.settings.get('scan_fits_headers', False):
        #     num_files = len(jumpinfo)
        #     if num_files <= self.scan_limit:
        #         self.scan_fits()
        #     else:
        #         self.logger.warning(
        #             "Number of files (%d) is greater than scan limit (%d)"
        #             "--skipping header scan" % (num_files, self.scan_limit))

        listing, resize = self.makelisting(jumpinfo)

        if which_tree == "reduced":
            self.reduced_treeview.set_tree(listing)
            # Check to see if there are appropriate PypeIt Slits files in this directory
            # If there are, enable only those

            if "keck_" in self.reduced_filepath:
                self.logger.info("We have (presumably) entered a PypeIt reduced directory.")

                tree_iterator = QtGui.QTreeWidgetItemIterator(self.reduced_treeview.widget)
                while tree_iterator.value():
                    item = tree_iterator.value()

                    # Takes the raw filepath and removes the last character (which is a wildcard), then appends the item text
                    child_path = os.path.join(self.reduced_filepath[:-1], item.text(1))
                    if not os.path.isdir(child_path):
                        if "Slits" not in item.text(1):
                            item.setDisabled(True)
                    tree_iterator += 1


            if resize:
                self.reduced_treeview.set_optimal_column_widths()
                self.logger.debug("Resized columns on reduced treeview")
        else:
            self.raw_treeview.set_tree(listing)

            # After we set the table, iterate over each item and disable it if it's not a FITS file or directory
            tree_iterator = QtGui.QTreeWidgetItemIterator(self.raw_treeview.widget)
            while tree_iterator.value():
                item = tree_iterator.value()

                # Takes the raw filepath and removes the last character (which is a wildcard), then appends the item text
                child_path = os.path.join(self.raw_filepath[:-1], item.text(1))

                if not os.path.isdir(child_path):
                    if "fits" not in item.text(2):
                        print(f"hiding {item.text(2)}")
                        item.setHidden(True)
                tree_iterator += 1
            if resize:
                self.raw_treeview.set_optimal_column_widths()
                self.logger.debug("Resized columns on raw treeview")


    def get_info(self, path):
        dirname, filename = os.path.split(path)
        name, ext = os.path.splitext(filename)
        ftype = 'file'
        if os.path.isdir(path):
            ftype = 'dir'
        elif os.path.islink(path):
            ftype = 'link'
        elif ext.lower() == '.fits':
            ftype = 'fits'

        header_dict = {}
        if ftype == 'fits':
            self.logger.info(f"Reading FITS header for {path}")
            header_dict = self.instrument.get_raw_info(path)
        #     with fits.open(path) as hdul:
        #         header = hdul[0].header
        #         header_dict['OBJECT'] = header.get('OBJECT', 'N/A')
        #         header_dict['FRAMENO'] = header.get('FRAMENO', 'N/A')
        #         header_dict['IMTYPE'] = header.get('KOAIMTYP', 'N/A') # Just use KOA, since that's the only consistent one
        #         header_dict['MASKNAME'] = header.get('MASKNAME', 'N/A')
        #         header_dict['OBSMODE'] = header.get('OBSMODE', 'N/A')
        #         header_dict['EXPTIME'] = header.get('EXPTIME', None)
        #         if header_dict['EXPTIME'] is None:
        #             header_dict['EXPTIME'] = header.get('TTIME', None)
        #         if header_dict['EXPTIME'] is None:
        #             header_dict['EXPTIME'] = header.get('ITIME', None)
        #         if header_dict['EXPTIME'] is None:
        #             header_dict['EXPTIME'] = header.get('ETIME', None)
        #         if header_dict['EXPTIME'] is None:
        #             header_dict['EXPTIME'] = header.get('ELAPTIME', "N/A")



        na_dict = {attrname: 'N/A' for colname, attrname in self.settings.get('columns')}
        bnch = Bunch.Bunch(na_dict)
        try:
            filestat = os.stat(path)
            if ftype == 'fits':
                self.logger.info(f"Parsing header dict for {path}")
                bnch.update(dict(path=path, name=filename, type=ftype,
                            st_mode=filestat.st_mode,
                            st_mode_oct=oct(filestat.st_mode),
                            st_size=filestat.st_size,
                            st_size_str=str(filestat.st_size),
                            st_mtime=filestat.st_mtime,
                            st_mtime_str=time.ctime(filestat.st_mtime),
                            OBJECT=header_dict['OBJECT'],
                            FRAMENO=header_dict['FRAMENO'],
                            IMTYPE=header_dict['IMTYPE'],
                            EXPTIME=header_dict['EXPTIME'],
                            MASKNAME=header_dict['MASKNAME'],
                            OBSMODE=header_dict['OBSMODE'],
                            ))
            else:
                bnch.update(dict(path=path, name=filename, type=ftype,
                            st_mode=filestat.st_mode,
                            st_mode_oct=oct(filestat.st_mode),
                            st_size=filestat.st_size,
                            st_size_str=str(filestat.st_size),
                            st_mtime=filestat.st_mtime,
                            st_mtime_str=time.ctime(filestat.st_mtime),
                            OBJECT="",
                            FRAMENO="",
                            IMTYPE="",
                            EXPTIME="",
                            MASKNAME="",
                            OBSMODE="",))
        except OSError as e:
            # TODO: identify some kind of error with this path
            self.logger.error("Error getting file info for %s: %s" % (path, str(e)))
            bnch.update(dict(path=path, name=filename, type=ftype,
                            st_mode=0, st_size=0,
                            st_mtime=0))

        return bnch

    def makelisting(self, jumpinfo):
        def file_icon(bnch):
            if bnch.type == 'dir':
                pb = self.folderpb
            elif bnch.type == 'fits':
                pb = self.fitspb
            else:
                pb = self.filepb
            return pb

        tree_dict = {}
        for bnch in jumpinfo:
            icon = file_icon(bnch)
            bnch.setvals(icon=icon)
            entry_key = bnch.name

            if entry_key is None:
                raise Exception("No key for tuple")

            tree_dict[entry_key] = bnch

        # Do we need to resize column widths?
        n_rows = len(tree_dict)
        if n_rows < self.settings.get('max_rows_for_col_resize', 5000):
            resize_table_columns = True
        else:
            resize_table_columns = False

        return tree_dict, resize_table_columns

    def open_slits_files(self, path):
        """Takes a path to a Calibration directory and loads the Slits files within

        Parameters
        ----------
        path : str
            Path to the directory containing the Slits files (not the Slits files themselves)
        """

        # Find all the Slits files
        # For each one in the mosaic, open it using SlitTraceSet.from_file. Note that Mosaic index is one-indexed
        # Save each SlitTraceSet to a dictionary with the key being the mosaic index
        # return the dictionary

        p = Path(path)
        slit_files = p.glob('Slits_*.fits*')
        slit_dict = {}
        for slit_file in slit_files:
            msc_idx = slit_file.stem.split('.')[0].split('_')[-1][-2:]
            slit_dict[msc_idx] = SlitTraceSet.from_file(slit_file)
        
        return slit_dict
    

    def construct_slits(self, slittraceset_dict):
        """Creates a compound object will polygons representing each slit

        Parameters
        ----------
        slittraceset_list : dict
            Dict of SlitTraceSet objects to be rendered, where the keys are MSC keys (as strings)
        """

        slit_polys = lines = self.dc.Canvas()

        for msc_idx in slittraceset_dict.keys():
            
            slittrace = slittraceset_dict[msc_idx]
            spatial_ids = slittrace.spat_id
            left_init = slittrace.left_init.T
            sampling = 200
            y_values_left = np.arange(slittrace.nspec)[::sampling]
            right_init = slittrace.right_init.T
            y_values_right = np.arange(slittrace.nspec)[::-sampling]

            for idx, spat_id in enumerate(spatial_ids):
                self.slit_list_box.append_text(f"S{spat_id}")

                x_vals = np.concatenate((left_init[idx][::sampling], right_init[idx][::-sampling]), axis=0) + ((int(msc_idx) - 1) * slittrace.nspat)
                y_vals = (-1 * np.concatenate((y_values_left, y_values_right), axis=0)) + (2048 * 4)
                slit_boundard_coords = (x_vals, y_vals)
                poly = Polygon(list(zip(slit_boundard_coords[0], slit_boundard_coords[1])), color='green', linewidth=1, fill=True, fillalpha=.05)
                # Ideally we would directly use a polygon callback to handle clicks, but that doesn't seem to be working
                # poly.add_callback('cursor-down', self.test_cb, slit_id=spat_id)
                slit_polys.add(poly)
                self.slit_polys[f"S{spat_id}"] = poly

        return slit_polys

    def activate_slit(self):
        if self.active_slit is None:
            self.logger.error("No active slit selected.")
            return
        self.logger.info(f"Activating slit {self.active_slit}")
        active_slit = self.slit_canvas.get_object_by_tag(self.active_slit.tag)
        active_slit.color = 'blue'
        self.slit_canvas.delete_object_by_tag(self.active_slit.tag)
        self.slit_canvas.add(active_slit)
    
    def deactivate_slit(self):
        if self.active_slit is None:
            self.logger.error("No active slit selected.")
            return
        self.logger.info(f"Deactivating slit {self.active_slit}")
        active_slit = self.slit_canvas.get_object_by_tag(self.active_slit.tag)
        active_slit.color = 'green'
        self.slit_canvas.delete_object_by_tag(self.active_slit.tag)
        self.slit_canvas.add(active_slit)
        
        

    def open_raw_file(self, path):
        """Open a FITS file at the given path, and render into the viewer.
        """
        self.logger.info(f"Call to opn file: {path}")
        if not os.path.isfile(path):
            self.logger.error(f"File is not a directory: {path}")
            return
        
        p = Path(path)

        if '.fits' in p.name:
            hdul = fits.open(path)
            img_data = None
            img_data = self.instrument.get_display_image(hdul)
            
            img = AstroImage(logger=self.logger)
            img.load_data(img_data)
            self.fitsimage.set_image(img)

            # Set the raw filepath, now that we know its a FITS file
            self.raw_filepath = path
            self.logger.info(f"Raw file loaded: {path}")
        else:
            self.logger.info("Not a FITS file. Ignoring.")

            

    # # # # # # # # # #
    # Ginga Functions #
    # # # # # # # # # #


    def plot_lines(self):
        """Plot the line list.

        Lines are made into a single compound object so that it is easier
        to remove them as a group if the line list is changed.
        """
        pass
        # self.plot.make_callback('modified')

    def replot(self):
        """Replot the plot and line list.
        """
        pass

    def display_image(self):
        pass
        # self.channel.add_image()

    def recalc(self):
        """Reprocess the chosen extension, based on current choices for extraction
        method, fluxing and masking.

        Replot everything as a result.
        """
        pass

        self.replot()

    def close(self):
        """Method called to close the plugin when the Close button is pressed."""
        self.pool.close()
        self.fv.stop_local_plugin(self.chname, str(self))

    def start(self):
        """Method called right after `build_gui` when the plugin is activated.

        Simply attempt to process the latest FITS file loaded in the channel.
        """
        self.browse(os.getcwd(), "reduced")
        self.browse(os.getcwd(), "raw")

        # Set the instrument based on the combobox selection
        self.instrument_combo.make_callback('activated')

        self.fitsimage.add_callback('cursor-down', self.canvas_clicked_cb)

        # Load config, if there is one
        config_file = Path.home() / ".quicklook.cfg"
        if config_file.exists():
            self.logger.info(f"Loading config file: {config_file}")
            import configparser
            config = configparser.ConfigParser()
            try:
                config.read(Path.home() / ".quicklook.cfg")

                self.redux_path_entry.set_text(config["DEFAULT"]["redux_path"])
                self.redux_path = config["DEFAULT"]["redux_path"]
                self.logger.info(f"Redux path set to {self.redux_path}")

                self.raw_text_entry.set_text(config["DEFAULT"]["raw_path"])
                self.raw_filepath = config["DEFAULT"]["raw_path"]
                self.logger.info(f"Raw path set to {self.raw_filepath}")

                self.reduced_text_entry.set_text(config["DEFAULT"]["reduced_path"])
                self.reduced_filepath = config["DEFAULT"]["reduced_path"]
                self.logger.info(f"Reduced path set to {self.reduced_filepath}")

                if os.path.isdir(self.raw_filepath):
                    self.logger.info(f"Raw path exists: {self.raw_filepath}")
                    self.browse(self.raw_filepath, "raw")
                else:
                    self.logger.error(f"Raw path does not exist: {self.raw_filepath}")
                if os.path.isdir(self.reduced_filepath):
                    self.logger.info(f"Reduced path exists: {self.raw_filepath}")
                    self.browse(self.reduced_filepath, "reduced")
                else:
                    self.logger.error(f"Reduced path does not exist: {self.reduced_filepath}")
                
            except KeyError as e:
                self.logger.error(f"Error reading config file: {e}")
        else:
            self.logger.info("No config file found. Using defaults.")

        self.watcher = QtCore.QFileSystemWatcher()
        self.watcher.addPath(self.redux_path)
        self.watcher.directoryChanged.connect(self.dir_change)
        self.watcher.fileChanged.connect(self.dir_change)

        self.redo()
    
    def dir_change(self, path):
        print("Directory changed: ", path)
        p = Path(path)
        if path == self.redux_path_entry.get_text():
            self.logger.info("New directory detected in redux path, adding to watcher")
            for child in p.iterdir():
                self.add_dirs_to_watcher([str(child.absolute())])
        
        p = p / "Science"
        if not p.exists():
            print("Not a science directory")
            return
        
        science_files = list(p.glob("*.fits"))
        if len(science_files) > 0:
            self.logger.info("Reduced spectrum found.")
            science_file = None
            for file in science_files:
                if file.name.startswith("spec1d"):
                    science_file = file
            if not science_file:
                self.logger.error("No spec1D file found in Science directory.")
                return
            reduction_control = self.reduction_control_elements.get(self.slit_list_box.get_text(), None)
            if reduction_control:
                elements = reduction_control["elements"]
                elements['button'].set_enabled(True)
                elements['button'].add_callback('activated', lambda x: self.show_reduced_spec(science_file))
                old_text = elements['name'].get_text()
                elements['name'].set_text(old_text.replace("Reducing", "Reduced").replace("...", ""))
            else:
                self.logger.error("No Reduction slits to enable!")
        else:
            # monitor the science directory for changes
            self.logger.info("No reduced spectrum found. Adding Science directory to watcher")
            self.add_dirs_to_watcher([str(p.absolute())])
            # For good measure, re-check the science directory
            self.dir_change(path)

        # for directory in p.iterdir():
        #     if directory.is_dir():
        #         print(f"\tFound directory: {directory.name}")
        #         science_name = Path(self.raw_filepath).name.split(".fits")[0]
        #         print("\tIf there's a science directory, it should look like this: ", science_name)
        #         print(f"\tLooing for {science_name} in {directory.name}")
        #         if science_name in directory.name:
        #             print("\t\tOK, this is a science directory...")
        #             for child in directory.iterdir():
        #                 if "Science" in child.name:
        #                     print("\t\t\tOK, this is a science directory with data inside")
        #                     self.logger.info("Science directory found.")
        #                     for file in child.iterdir():
        #                         if file.suffix == ".fits" and "spec1d" in file.name:
        #                             print("\t\t\t\tOK, there is a spec1d file")
        #                             self.logger.info("Reduced spectrum found.")
        #                             reduction_control = self.reduction_control_elements.get(self.slit_list_box.get_text(), None)
        #                             if reduction_control:
        #                                 print("\t\t\t\t\tOK, I'm ready to show this slit")
        #                                 elements = reduction_control["elements"]
        #                                 elements['button'].set_enabled(True)
        #                                 elements['button'].add_callback('activated', self.show_reduced_spec)
        #                             else:
        #                                 print("\t\t\t\t\t\tNo reduction control found")


    def stop(self):
        """Method called when the plugin is deactivated.

        Clean up instance variables so we don't hang on to any large data
        structures.
        """
        self.gui_up = False

    def redo(self):
        """Method called when a new FITS image or extension is loaded into the
        channel.
        """
        pass

    def show_reduced_spec(self, args):
        """Show the reduced spectrum in a new channel"""
        self.logger.info(f"Showing reduced spectrum: {args}")
        new_ch_name = 'Spec1D' + self.slit_list_box.get_text()
        self.fv.load_file(str(args), chname=new_ch_name)
        self.fv.start_local_plugin(new_ch_name, 'Spec1dView')

        
    
    def make_reduced_slit_hbox(self):

        # def timer_cb():
        #     elements['timer_left'].set_text(f"Time Remaining: {elements['timer'].remainingTime() / 1000} s")

        elements = {}
        hbox = Widgets.HBox()
        elements['name'] = Widgets.Label("Slit:")
        elements['button'] = Widgets.Button("Show")
        elements['button'].set_enabled(False)

        # elements['timer_left'] = Widgets.Label("Time Remainig:")
        
        # elements['timer'] = QtCore.QTimer()
        # elements['timer'].timeout.connect(timer_cb)

        hbox.add_widget(elements['name'], stretch=1)
        hbox.add_widget(elements['button'], stretch=0)
        # hbox.add_widget(elements['timer_left'], stretch=0)
        return hbox, elements
    
    def add_dirs_to_watcher(self, paths):
        for path in paths:
            if len(self.watcher.directories()) > 255:
                self.logger.error("Too many directories being watched. Cannot add more.")
                return
            self.logger.info(f"Adding {path} to watcher")
            if Path(path).exists():
                watching = self.watcher.addPath(path)
                if not watching:
                    self.logger.error(f"Failed to add {path} to watcher!")
            else:
                self.logger.error(f"Path {path} does not exist. Cannot add to watcher.")

    def redux_path_from_filename(self, filename):
        p = Path(filename)
        return os.path.join(self.redux_path_entry.get_text(),  p.name.split(".fits")[0])

    def __str__(self):
        # necessary to identify the plugin and provide correct operation in Ginga
        return 'qlview'


# # # # # # # # # # # # # # # #
# Instrument-specific classes #
# # # # # # # # # # # # # # # #

class Instrument():
    """
    Class that generalizes instrument-specific information for the QL viewer.
    This includes how to get a single 2D array representation of the entire image
    (i.e. mosaicing), what columns we want to display in both the raw and reduced
    filetrees, and how to get that information from raw FITS files.

    This could conceivably be part of the Spectrograph class in PypeIt, but
    for now we are keeping it separate.
    """
    def __init__(self, logger) -> None:
        self.logger = logger
        self.raw_columns = []

    def get_display_image(self, hdul) -> np.ndarray:
        """Abstract method to return a mosaiced image from the HDUList.

        Parameters
        ----------
        HDUList : astropy.io.fits.HDUList
            The HDUList object containing the raw data.

        Returns
        -------
        np.ndarray
            mosaiced image, in coordinates that are appropriate for the viewer.
        """
        raise NotImplementedError("get_mosaic method not implemented for this instrument")

    def get_raw_info(self, path):
        """Get the raw information from the file at the given path.

        Parameters
        ----------
        path : str
            Path to the file.

        Returns
        -------
        ginga.misc.Bunch
            A Bunch object containing the raw column information.
        """
        raise NotImplementedError("get_raw_info method not implemented for this instrument")

class DEIMOS(Instrument):

    def __init__(self, logger) -> None:
        super().__init__(logger)
        self.pypeit_name = "keck_deimos"

    def get_display_image(self, hdul) -> np.ndarray:
        """Return a mosaiced image from the DEIMOS HDUList.

        Parameters
        ----------
        HDUList : astropy.io.fits.HDUList
            The HDUList object containing the raw data.

        Returns
        -------
        np.ndarray
            mosaiced image, in coordinates that are appropriate for the viewer.
        """
        
        hdr0 = hdul[0].header
        
        ext = np.arange(1,9)
        
        binning = hdr0['BINNING'].split(',')
        
        precol =   int(hdr0['PRECOL'])   // int(binning[0])
        postpix =  int(hdr0['POSTPIX'])  // int(binning[0])
        # preline =  int(hdr0['PRELINE'])  // int(binning[1])
        # postline = int(hdr0['POSTLINE']) // int(binning[1])
        
        alldata = []
        for i in ext:
            data = hdul[i].data
            
            height, width = hdul[i].shape
            # get bias from overscan region
            x1 = 0
            x2 = height
            y1 = width - postpix
            y2 = width
            bias = np.median(data[x1:x2, y1:y2], axis = 1)
            bias = np.array(bias, dtype = np.int64)            
            
            # bias subtraction
            data = data - bias[:, None]       
            
            # remove the overscan region
            data = data[:, precol: width - postpix]
            
            # append all 8 biased arrays into a list
            alldata.append(data)
        
        # creating CCD mosaic rows
        r0 = np.concatenate(alldata[:4], axis=1)
        r0 = np.flipud(r0)
        r1 = []
        for arr in alldata[4:]:
            arr = np.fliplr(arr)
            r1.append(arr)
        r1 = np.concatenate(r1, axis=1)

        fulldata = np.concatenate((r1, r0), axis=0)
        # fulldata = np.rot90(fulldata)
        
        return fulldata
    
    def get_raw_info(self, path):

        # dirname, filename = os.path.split(path)
        # name, ext = os.path.splitext(filename)
        # ftype = 'file'
        # if os.path.isdir(path):
        #     ftype = 'dir'
        # elif os.path.islink(path):
        #     ftype = 'link'
        # elif ext.lower() == '.fits':
        #     ftype = 'fits'

        header_dict = {}
        with fits.open(path) as hdul:
            header = hdul[0].header
            header_dict['OBJECT'] = header.get('OBJECT', 'N/A')
            header_dict['FRAMENO'] = header.get('FRAMENO', 'N/A')
            header_dict['IMTYPE'] = header.get('KOAIMTYP', 'N/A') # Just use KOA, since that's the only consistent one
            header_dict['MASKNAME'] = header.get('MASKNAME', 'N/A')
            header_dict['OBSMODE'] = header.get('OBSMODE', 'N/A')
            header_dict['EXPTIME'] = header.get('EXPTIME', None)
            if header_dict['EXPTIME'] is None:
                header_dict['EXPTIME'] = header.get('TTIME', None)
            if header_dict['EXPTIME'] is None:
                header_dict['EXPTIME'] = header.get('ITIME', None)
            if header_dict['EXPTIME'] is None:
                header_dict['EXPTIME'] = header.get('ETIME', None)
            if header_dict['EXPTIME'] is None:
                header_dict['EXPTIME'] = header.get('ELAPTIME', "N/A")

        return header_dict

        # # na_dict = {attrname: 'N/A' for colname, attrname in self.settings.get('columns')}
        # # bnch = Bunch.Bunch(na_dict)
        # try:
        #     filestat = os.stat(path)
        #     try:
        #         bnch.update(dict(path=path,
        #                     name=filename,
        #                     type=ftype,
        #                     st_mode=filestat.st_mode,
        #                     st_mode_oct=oct(filestat.st_mode),
        #                     st_size=filestat.st_size,
        #                     st_size_str=str(filestat.st_size),
        #                     st_mtime=filestat.st_mtime,
        #                     st_mtime_str=time.ctime(filestat.st_mtime),
        #                     OBJECT=header_dict['OBJECT'],
        #                     FRAMENO=header_dict['FRAMENO'],
        #                     IMTYPE=header_dict['IMTYPE'],
        #                     EXPTIME=header_dict['EXPTIME'],
        #                     MASKNAME=header_dict['MASKNAME'],
        #                     OBSMODE=header_dict['OBSMODE'],
        #                     ))
        #     except KeyError as e:
        #         self.logger.error(f"KeyError in DEIMOS get_raw_info: {e}")
        #         # bnch.update(dict(path=path, name=filename, type=ftype,
        #         #             st_mode=filestat.st_mode,
        #         #             st_mode_oct=oct(filestat.st_mode),
        #         #             st_size=filestat.st_size,
        #         #             st_size_str=str(filestat.st_size),
        #         #             st_mtime=filestat.st_mtime,
        #         #             st_mtime_str=time.ctime(filestat.st_mtime),))
        # except OSError as e:
        #     # TODO: identify some kind of error with this path
        #     bnch.update(dict(path=path, name=filename, type=ftype,
        #                     st_mode=0, st_size=0,
        #                     st_mtime=0))
        
        # return bnch

class MOSFIRE(Instrument):

    def __init__(self, logger) -> None:
        super().__init__(logger)
        self.pypeit_name = "keck_mosfire"

        self.raw_columns = [
            ('Type', 'icon'),
            ('Frame No', 'FRAMENO'),
            ('Name', 'name'),
            ('Object', 'OBJECT'),
            ('Type', 'IMTYPE'),
            ('Mask Name', 'MASKNAME'),
            ('Obs Mode', 'OBSMODE'),
            ('Exposure Time', 'EXPTIME'),
            ('Last Changed', 'st_mtime_str')
        ]
    
    def get_raw_info(self, path):

        # dirname, filename = os.path.split(path)
        # name, ext = os.path.splitext(filename)
        # ftype = 'file'
        # if os.path.isdir(path):
        #     ftype = 'dir'
        # elif os.path.islink(path):
        #     ftype = 'link'
        # elif ext.lower() == '.fits':
        #     ftype = 'fits'

        header_dict = {}
        with fits.open(path) as hdul:
            header = hdul[0].header
            header_dict['OBJECT'] = header.get('OBJECT', 'N/A')
            header_dict['FRAMENO'] = header.get('FRAMENO', 'N/A')
            header_dict['IMTYPE'] = header.get('KOAIMTYP', 'N/A') # Just use KOA, since that's the only consistent one
            header_dict['MASKNAME'] = header.get('MASKNAME', 'N/A')
            header_dict['OBSMODE'] = header.get('OBSMODE', 'N/A')
            header_dict['EXPTIME'] = header.get('EXPTIME', None)
            if header_dict['EXPTIME'] is None:
                header_dict['EXPTIME'] = header.get('TTIME', None)
            if header_dict['EXPTIME'] is None:
                header_dict['EXPTIME'] = header.get('ITIME', None)
            if header_dict['EXPTIME'] is None:
                header_dict['EXPTIME'] = header.get('ETIME', None)
            if header_dict['EXPTIME'] is None:
                header_dict['EXPTIME'] = header.get('ELAPTIME', "N/A")

        return header_dict

        # na_dict = {attrname: 'N/A' for colname, attrname in self.settings.get('columns')}
        # bnch = Bunch.Bunch(na_dict)
        # try:
        #     filestat = os.stat(path)
        #     try:
        #         bnch.update(dict(path=path, name=filename, type=ftype,
        #                     st_mode=filestat.st_mode,
        #                     st_mode_oct=oct(filestat.st_mode),
        #                     st_size=filestat.st_size,
        #                     st_size_str=str(filestat.st_size),
        #                     st_mtime=filestat.st_mtime,
        #                     st_mtime_str=time.ctime(filestat.st_mtime),
        #                     OBJECT=header_dict['OBJECT'],
        #                     FRAMENO=header_dict['FRAMENO'],
        #                     IMTYPE=header_dict['IMTYPE'],
        #                     EXPTIME=header_dict['EXPTIME'],
        #                     MASKNAME=header_dict['MASKNAME'],
        #                     OBSMODE=header_dict['OBSMODE'],
        #                     ))
        #     except KeyError as e:
        #         bnch.update(dict(path=path, name=filename, type=ftype,
        #                     st_mode=filestat.st_mode,
        #                     st_mode_oct=oct(filestat.st_mode),
        #                     st_size=filestat.st_size,
        #                     st_size_str=str(filestat.st_size),
        #                     st_mtime=filestat.st_mtime,
        #                     st_mtime_str=time.ctime(filestat.st_mtime),))
        # except OSError as e:
        #     # TODO: identify some kind of error with this path
        #     bnch.update(dict(path=path, name=filename, type=ftype,
        #                     st_mode=0, st_size=0,
        #                     st_mtime=0))
        
        # return bnch

    def get_display_image(self, hdul) -> np.ndarray:
        """Return a mosaiced image from the MOSFIRE HDUList.

        Parameters
        ----------
        HDUList : astropy.io.fits.HDUList
            The HDUList object containing the raw data.

        Returns
        -------
        np.ndarray
            image, in coordinates that are appropriate for the viewer.
        """
        
        fulldata = hdul[0].data
        
        return np.rot90(fulldata)


def sleepy():
    print("Starting sleepy...")
    time.sleep(5)
    print("Done sleeping, now running")

# class popen_w_cb():
#     """Mirrors Popen, but with a callback for when the process finishes.
#     This doesn't actually work since the callback is executed in the thread, though
#     """

#     def __init__(self, onExit, *popenArgs, **popenKWArgs):
#         thread = threading.Thread(target=self.runInThread, args=(onExit, popenArgs, popenKWArgs))
#         thread.start()

#     def runInThread(self, onExit, popenArgs, popenKWArgs):
#         self.proc = subprocess.Popen(*popenArgs, **popenKWArgs)
#         self.proc.wait()
#         self.proc = None
#         onExit()

#     def stop(self):
#         if self.proc:
#             self.proc.kill()