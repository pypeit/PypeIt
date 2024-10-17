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
import numpy as np

from ginga import GingaPlugin
from ginga.misc import Bunch
from ginga.gw import Widgets
from ginga.table.AstroTable import AstroTable
from ginga.plot.Plotable import Plotable
from ginga.canvas.CanvasObject import get_canvas_types

from pypeit import specobjs
from pypeit import utils
# TODO: need to make PypeIt LineList class to deprecate linetools
from linetools.lists.linelist import LineList

__all__ = ['Spec1dView']


class Spec1dView(GingaPlugin.LocalPlugin):

    def __init__(self, fv, fitsimage):
        """Constructor for the plugin."""
        # superclass defines some variables for us, like logger
        super().__init__(fv, fitsimage)

        # get Spec1dView preferences
        prefs = self.fv.get_preferences()
        self.settings = prefs.create_category('plugin_Spec1dView')
        self.settings.add_defaults(lines="error",
                                   extraction='OPT', fluxed=False, masked=False)
        self.settings.load(onError='silent')

        # will be set if we are invoked
        self.data = Bunch.Bunch()
        self.sobjs = []
        self.num_exten = 0
        self.exten = 0

        self.w = None
        # selected line list
        self.line_list = 'None'
        # allowed line lists
        self.line_lists = ['None', 'ISM', 'Strong', 'HI', 'H2', 'CO',
                           'EUV', 'Galaxy', 'AGN']
        self.llist = None   # the actual line list object

        # redshift
        self.z = 0.0

        self.extraction_types = ('OPT', 'BOX')
        self.extraction = self.settings.get('extraction', 'OPT')

        self.fluxed_options = (True, False)
        self.fluxed = self.settings.get('fluxed', False)

        self.masked_options = (True, False)
        self.masked = self.settings.get('masked', False)

        # dictionary of plotable types
        self.dc = get_canvas_types()
        self.plot = None
        self.gui_up = False

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

        fr = Widgets.Frame("File")

        captions = (("1D file:", 'label', 'filepath', 'entryset'),
                    ("Extensions:", 'label', 'extensions', 'llabel'),
                    ("Extension:", 'label', 'exten', 'combobox'),
                    ("Extraction:", 'label', 'extraction', 'combobox'),
                    ("Fluxed:", 'label', 'fluxed', 'combobox'),
                    ("Masked:", 'label', 'masked', 'combobox'),
                    )
        w, b = Widgets.build_info(captions, orientation='vertical')
        self.w = b
        b.filepath.add_callback('activated', self.set_filepath_cb)
        b.exten.set_tooltip("Choose extension of file to view")
        b.exten.add_callback('activated', self.set_exten_cb)
        b.exten.set_enabled(False)

        combobox = b.extraction
        for name in self.extraction_types:
            combobox.append_text(name)
        index = self.extraction_types.index(self.extraction)
        combobox.set_index(index)
        combobox.add_callback('activated', self.set_extraction_cb)

        combobox = b.fluxed
        for name in self.fluxed_options:
            combobox.append_text(str(name))
        index = self.fluxed_options.index(self.fluxed)
        combobox.set_index(index)
        combobox.add_callback('activated', self.set_fluxed_cb)

        combobox = b.masked
        for name in self.masked_options:
            combobox.append_text(str(name))
        index = self.masked_options.index(self.masked)
        combobox.set_index(index)
        combobox.add_callback('activated', self.set_masked_cb)

        fr.set_widget(w)
        vbox.add_widget(fr, stretch=0)

        fr = Widgets.Frame("Lines")

        captions = (("Z:", 'label', 'redshift', 'entryset'),
                    ("Line lists:", 'label', 'lists', 'combobox'),
                    )

        w, b = Widgets.build_info(captions, orientation="vertical")
        self.w.update(b)

        b.redshift.set_text(str(self.z))
        b.redshift.add_callback('activated', self.set_z_cb)
        b.redshift.set_tooltip("Set the redshift (Z) value")

        combobox = b.lists
        for name in self.line_lists:
            combobox.append_text(name)
        index = self.line_lists.index(self.line_list)
        combobox.set_index(index)
        combobox.set_tooltip("Line list to plot on this spectrum")
        combobox.add_callback('activated', self.set_line_list_cb)

        fr.set_widget(w)
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
        btns.add_widget(Widgets.Label(''), stretch=1)
        top.add_widget(btns, stretch=0)

        container.add_widget(top, stretch=1)
        self.gui_up = True

    def set_z_cb(self, w):
        """Callback for setting the zed value in the plugin.

        Replot the lines as a result.
        """
        z = float(w.get_text())
        self.logger.info(f"z = {z}")
        self.z = z

        self.fv.gui_do(self.plot_lines)

    def set_line_list_cb(self, w, idx):
        """Callback for setting the line list in the plugin.

        Construct a new line list (or None) and replot the lines as a result.
        """
        self.line_list = self.line_lists[idx]
        if self.line_list == 'None':
            self.llist = None
        else:
            self.llist = LineList(self.line_list)
        self.logger.info(f"Loaded line list: '{self.line_list}'")

        self.fv.gui_do(self.plot_lines)

    def set_extraction_cb(self, w, idx):
        """Callback for changing the `extraction` option in the plugin.

        Redo the data extraction and replot everything as a result.
        """
        self.extraction = self.extraction_types[idx]
        self.logger.debug(f"Selected extraction type: {self.extraction}")
        self.recalc()

    def set_fluxed_cb(self, w, idx):
        """Callback for changing the `fluxed` option in the plugin.

        Redo the data extraction and replot everything as a result.
        """
        self.fluxed = self.fluxed_options[idx]
        self.logger.debug(f"Selected fluxed option: {self.fluxed}")
        self.recalc()

    def set_masked_cb(self, w, idx):
        """Callback for changing the `masked` option in the plugin.

        Redo the data extraction and replot everything as a result.
        """
        self.masked = self.masked_options[idx]
        self.logger.debug(f"Selected masked option: {self.masked}")
        self.recalc()

    def plot_lines(self):
        """Plot the line list.

        Lines are made into a single compound object so that it is easier
        to remove them as a group if the line list is changed.
        """
        canvas = self.plot.get_canvas()
        canvas.delete_object_by_tag('lines', redraw=False)

        lines = self.dc.Canvas()
        if self.llist is not None:
            # NOTE: adapted straight from linetools
            x_min, x_max = self.data.x_min, self.data.x_max
            y_min, y_max = self.data.y_min, self.data.y_max

            z = self.z
            wvobs = np.array((1 + z) * self.llist.wrest)
            ylbl_pos = y_max - 0.2 * (y_max - y_min)
            gdwv = np.where((wvobs > x_min) & (wvobs < x_max))[0]

            for kk in range(len(gdwv)):
                jj = gdwv[kk]
                wrest = self.llist.wrest[jj].value
                lbl = self.llist.name[jj]
                # Plot
                x_data = wrest * np.array([z + 1, z + 1])
                y_data = (y_min, y_max)
                lines.add(self.dc.Line(x_data[0], y_data[0],
                                       x_data[1], y_data[1],
                                       linewidth=1,
                                       #linestyle='dotted',
                                       linestyle='solid',
                                       color='blue'), redraw=False)
                # Label
                x, y = wrest * (z + 1), ylbl_pos
                lines.add(self.dc.Text(x, y, text=lbl, rot_deg=90,
                                       color='blue', fontsize=10),
                          redraw=False)

        # will be empty if there were no lines
        canvas.add(lines, tag='lines', redraw=False)

        # this causes the plot viewer to redraw itself
        self.plot.make_callback('modified')

    def replot(self):
        """Replot the plot and line list.
        """
        if self.plot is None:
            return
        # plot flux vs. wavelength
        self.plot.clear()
        canvas = self.plot.get_canvas()

        # plot flux vs. wavelen
        points = np.array((self.data.wave, self.data.flux)).T
        p1 = self.dc.Path(points, linewidth=1, linestyle='solid',
                          alpha=0.7, color='black')
        canvas.add(p1, tag='spectrum', redraw=False)

        # additionally, plot error vs. wavelength
        points = np.array((self.data.wave, self.data.sig)).T
        p2 = self.dc.Path(points, linewidth=2, linestyle='solid',
                          alpha=1.0, color='red')
        canvas.add(p2, tag='error', redraw=False)

        self.plot.set_titles(title="Spec1dView",
                             x_axis="Wavelength (Ang)", y_axis="Flux")
        self.plot.set_grid(True)

        self.plot_lines()

    def process_file(self, filepath):
        """Process `filepath` creating `SpecObjs` (a series of extensions),
        which can then have data extracted and plotted.
        """
        self.w.filepath.set_text(filepath)

        self.sobjs = specobjs.SpecObjs.from_fitsfile(filepath,
                                                     chk_version=False)
        self.num_exten = len(self.sobjs)
        self.w.extensions.set_text(str(self.num_exten))

        self.w.exten.clear()
        for name in self.sobjs.NAME:
            self.w.exten.append_text(name)
        self.w.exten.set_enabled(True)
        self.exten = 0
        self.w.exten.set_index(self.exten)

        # create a new plotable to contain the plot
        self.plot = Plotable(logger=self.logger)
        self.plot.set(name="plot-{}".format(str(time.time())),
                      path=None, nothumb=True)

        self.recalc()

        # add the plot to this channel
        self.channel.add_image(self.plot)

    def recalc(self):
        """Reprocess the chosen extension, based on current choices for extraction
        method, fluxing and masking.

        Replot everything as a result.
        """
        specobj = self.sobjs[self.exten]
        if specobj['OPT_WAVE'] is None:
            self.fv.show_error("Spectrum not extracted with OPT.  Try --extract BOX")
            return

        wave, flux, ivar, gpm = specobj.to_arrays(extraction=self.extraction,
                                                  fluxed=self.fluxed)
        sig = np.sqrt(utils.inverse(ivar))
        wave_gpm = wave > 1.0
        wave, flux, sig, gpm = (wave[wave_gpm], flux[wave_gpm],
                                sig[wave_gpm], gpm[wave_gpm])
        if self.masked:
            flux = flux*gpm
            sig = sig*gpm

        self.data.x_min, self.data.x_max = np.nanmin(wave), np.nanmax(wave)
        self.data.y_min, self.data.y_max = np.nanmin(flux), np.nanmax(flux)
        self.data.wave = wave
        self.data.flux = flux
        self.data.sig = sig

        self.replot()

    def close(self):
        """Method called to close the plugin when the Close button is pressed."""
        self.fv.stop_local_plugin(self.chname, str(self))

    def start(self):
        """Method called right after `build_gui` when the plugin is activated.

        Simply attempt to process the latest FITS file loaded in the channel.
        """
        self.redo()

    def stop(self):
        """Method called when the plugin is deactivated.

        Clean up instance variables so we don't hang on to any large data
        structures.
        """
        self.sobjs = 0
        self.exten = 0
        self.num_exten = 0
        self.data = Bunch.Bunch()
        self.plot = None
        self.gui_up = False

    def redo(self):
        """Method called when a new FITS image or extension is loaded into the
        channel.

        We do some minimal checks to make sure that it is a table, then call the
        routine to process the file that is behind this table.
        """
        if not self.gui_up:
            return

        dataobj = self.channel.get_current_image()
        if not isinstance(dataobj, AstroTable):
            # NOTE: do we need a better test here for a 1D data item?
            return

        path = dataobj.get('path', None)
        if path is None:
            self.fv.show_error(
                "Cannot open dataobj: no value for metadata key 'path'")
            return

        self.process_file(path)

    def set_filepath_cb(self, w):
        """Callback for changing the path in the UI.

        Try to process the new file.
        """
        filepath = w.get_text().strip()
        self.process_file(filepath)

    def set_exten_cb(self, w, val):
        """Callback for changing the extension in the UI.

        Try to process the new extension.
        """
        self.exten = val
        self.recalc()

    def __str__(self):
        # necessary to identify the plugin and provide correct operation in Ginga
        return 'spec1dview'
