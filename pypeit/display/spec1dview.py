"""
Spec1dView is a plugin for the Ginga viewer that provides functionality for
visualizing and analyzing 1D spectra from FITS files. The plugin allows users
to plot spectra, identify spectral lines from various line lists, and
customize the display according to different parameters.

Plugin Type: Local
==================

Spec1dView is a local plugin, which means it is associated with a specific
channel in the Ginga viewer. An instance of the plugin can be opened for
each channel, allowing for multiple spectra to be analyzed simultaneously.

Usage
-----

- Load and visualize 1D spectra from FITS files.

- Customize the display by selecting different line lists, extraction types,
  and flux/mask settings.

- Update the redshift to shift the spectral lines accordingly.

Editing
-------

Users can modify the visualization by:

- Choosing from a variety of line lists to identify spectral features.

- Selecting different types of extraction methods (OPT, BOX).

- Applying or removing flux calibration and masking options.

- Updating the redshift value to reflect the observed wavelengths.

UI
--

The user interface provides controls for:

- Selecting the line list from a combobox.

- Entering a redshift value to shift the spectrum.

- Choosing the extraction type, flux calibration, and masking options via
  comboboxes.

- Buttons to load a FITS file and clear the current selection.

Buttons
-------

- Update z: Updates the redshift value and refreshes the spectrum plot.

- Enter: Loads the specified FITS file for analysis.

- Clear: Clears the current inputs and resets the UI settings.

Tips
----

- Use the comboboxes to switch between different line lists and adjust the
  spectrum display settings.

- Ensure that the correct FITS file path is entered before attempting to load
  the data.

"""
import time

import numpy as np
# for gaussian line fitting
from astropy.modeling import models, fitting

from ginga import GingaPlugin
from ginga.misc import Bunch
from ginga.gw import Widgets
from ginga.table.AstroTable import AstroTable
from ginga.plot.Plotable import Plotable
from ginga.canvas.CanvasObject import get_canvas_types
from ginga.util.syncops import Shelf

from pypeit import specobjs
from pypeit import utils

__all__ = ['Spec1dView']


class Spec1dView(GingaPlugin.LocalPlugin):

    def __init__(self, fv, fitsimage):
        """Constructor for the plugin."""
        # superclass defines some variables for us, like logger
        super().__init__(fv, fitsimage)

        # get Spec1dView preferences
        prefs = self.fv.get_preferences()
        self.settings = prefs.create_category('plugin_Spec1dView')
        self.settings.add_defaults(lines="error", start_ext=0,
                                   extraction='OPT', fluxed=False, masked=False,
                                   plot_error=True, autozoom=True)
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
        self.line_lists = ['None'] + utils.get_line_list_names()
        self.llist = None   # the actual line list object
        self.ext_name = ''

        # redshift
        self.z = 0.0
        self.start_ext = self.settings.get('start_ext', 0)

        self.extraction_types = ('OPT', 'BOX')
        self.extraction = self.settings.get('extraction', 'OPT')

        self.fluxed_options = (True, False)
        self.fluxed = self.settings.get('fluxed', False)

        self.masked_options = (True, False)
        self.masked = self.settings.get('masked', False)

        # smoothing factor
        self.nsmooth = 0.0

        # dictionary of plotable types
        self.dc = get_canvas_types()
        self.plot = None
        self.plot_shelf = Shelf()
        self.plot_stocker = self.plot_shelf.get_stocker()

        # holds regions of interest in the X axis
        self.region_dct = dict(mark=[None, None])
        # holds markers in the X axis
        self.marker_dct = dict()

        viewer = self.channel.get_viewer('Ginga Plot')
        viewer.add_callback('range-set', self.range_changed_cb)
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
                    ("Plot Error", 'checkbox'),
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

        b.plot_error.set_state(self.settings.get('plot_error', True))
        b.plot_error.set_tooltip("Include the error plot")
        b.plot_error.add_callback('activated', self.plot_error_cb)

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

        tbar = Widgets.Toolbar(orientation='horizontal')
        vbox.add_widget(tbar, stretch=0)

        btn = tbar.add_action("No smooth")
        btn.add_callback('activated', lambda w: self.no_smooth())
        btn.set_enabled(self.nsmooth > 0.0)
        self.w.no_smooth = btn
        btn.set_tooltip("Turn off smoothing")
        btn = tbar.add_action("Smooth-")
        btn.add_callback('activated', lambda w: self.smooth_less())
        btn.set_enabled(self.nsmooth > 0.0)
        btn.set_tooltip("Smooth a bit less")
        self.w.smooth_less = btn
        btn = tbar.add_action("Smooth+")
        btn.add_callback('activated', lambda w: self.smooth_more())
        btn.set_tooltip("Smooth a bit more")
        self.w.smooth_more = btn

        fr = Widgets.Frame("Result")
        self.w.text = Widgets.TextArea(wrap=True)
        fr.set_widget(self.w.text)
        vbox.add_widget(fr, stretch=1)

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
            self.llist = utils.get_line_list(self.line_list)
        self.logger.info(f"Loaded line list: '{self.line_list}'")

        self.fv.gui_do(self.plot_lines)

    def fit_y(self):
        """Zoom to fit Y axis"""
        viewer = self.channel.get_viewer('Ginga Plot')
        viewer.zoom_fit(axis='y')

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
        self.fit_y()

    def set_masked_cb(self, w, idx):
        """Callback for changing the `masked` option in the plugin.

        Redo the data extraction and replot everything as a result.
        """
        self.masked = self.masked_options[idx]
        self.logger.debug(f"Selected masked option: {self.masked}")
        self.recalc()

    def no_smooth(self):
        """Completely remove all smoothing from the plot."""
        self.nsmooth = 0.0
        self.w.no_smooth.set_enabled(False)
        self.w.smooth_less.set_enabled(False)
        self.recalc()

    def smooth_more(self):
        """Increase the level of smoothing in the plot."""
        self.nsmooth += (0.5 if self.nsmooth > 0.0 else 1.0)
        self.w.no_smooth.set_enabled(True)
        self.w.smooth_less.set_enabled(True)
        self.recalc()

    def smooth_less(self):
        """Decrease the level of smoothing in the plot."""
        self.nsmooth -= (0.5 if self.nsmooth > 1.0 else 1.0)
        self.nsmooth = max(self.nsmooth, 0.0)
        self.w.no_smooth.set_enabled(self.nsmooth > 0.0)
        self.w.smooth_less.set_enabled(self.nsmooth > 0.0)
        self.recalc()

    def set_region_side(self, name, side, value):
        """Define one side (left, right, top, bottom) of a named region.

        Parameters
        ----------
        name : str
            The name of a region (e.g. 'mark')

        side : str ('left', 'right', 'top', 'bottom')
            The side of the region to set

        value : float
            The value to set for this side of the region
        """
        region = self.region_dct.setdefault(name, [None, None])
        if side in ('left', 'bottom'):
            region[0] = value
            if region[1] is not None and region[1] <= value:
                region[1] = None
        elif side in ('right', 'top'):
            region[1] = value
            if region[0] is not None and region[0] >= value:
                region[0] = None

        self.plot_lines()

    def clear_region(self, name):
        """Clear the region named `name`."""
        self.region_dct[name] = [None, None]

        self.plot_lines()

    def clear_markers(self):
        """Clear all the markers."""
        self.marker_dct = dict()

        self.plot_lines()

    def cut_x_range(self, name):
        """Zoom the X range of the plot to the region named `name`."""
        x_lo, x_hi = self.region_dct[name]
        if None in (x_lo, x_hi):
            self.fv.show_error("Please set left and right markers to delineate region to be cut!")
            return

        viewer = self.channel.get_viewer('Ginga Plot')
        viewer.set_ranges(x_range=(x_lo, x_hi))
        self.clear_region(name)

    def fit_range(self, name):
        """Perform a fitting on the flux data within the region named `name`."""
        x_lo, x_hi = self.region_dct[name]
        if None in (x_lo, x_hi):
            self.fv.show_error("Please set left and right markers to delineate region to be fitted!")
            return

        try:
            # find indices of low and high X values
            idx_lo = np.searchsorted(self.data.wave, x_lo, side="left")
            idx_hi = np.searchsorted(self.data.wave, x_hi, side="right")

            # fit a gaussian
            model, fitter = fit_gaussian_line(self.data.wave, self.data.flux,
                                              idx_lo, idx_hi)
            continuum = model[0]
            line = model[1]

        except Exception as e:
            errtxt = f"error doing fitting: {e}; traceback in viewer log"
            self.w.text.set_text(errtxt)
            self.logger.error(f"error doing fitting: {e}", exc_info=True)
            return

        # collect some results and display
        line_center = line.mean.value
        line_sigma = line.stddev.value
        line_fwhm = 2.35482 * line_sigma
        line_flux = line.amplitude.value * np.sqrt(2 * np.pi) * line_sigma

        # Plot the fitted gaussian
        # Dense grid for smooth model curve
        x_fitted = np.linspace(self.data.wave[idx_lo],
                              self.data.wave[idx_hi], 100)  # 1000?
        y_fitted = model(x_fitted)

        result = dict(x=line_center, sigma=line_sigma, fwhm=line_fwhm,
                      flux=line_flux, x_fitted=x_fitted, y_fitted=y_fitted)

        # TODO: probably want to plot multiple unique lines and label them
        # name = str(line_center)  # ?? what is a better unique name
        # for now, just plot the last result
        # self.add_marker(name, result)
        self.marker_dct = dict(fitting=result)

        text = """
        Center: {x:.4f}
        Sigma: {sigma:.6f}
        FWHM: {fwhm:.6f}
        Flux: {flux:.4f}
        """.format(**result)
        self.w.text.set_text(text)

        self.plot_lines()

    def add_marker(self, name, dct):
        self.marker_dct[name] = dct

        self.plot_lines()

    def clear_marker(self, name):
        if name in self.marker_dct:
            del self.marker_dct[name]

            self.plot_lines()

    def clear_markers(self):
        self.marker_dct = dict()
        self.plot_lines()

    def plot_markers(self):
        """Plot the regions and markers.
        """
        if self.plot is None:
            return

        canvas = self.plot.get_canvas()
        objs = canvas.get_objects_by_tag_pfx('mark_')
        canvas.delete_objects(objs)

        x_lo, x_hi = self.region_dct['mark']
        if None not in (x_lo, x_hi):
            canvas.add(self.dc.XRange(x_lo, x_hi, fillcolor='aquamarine',
                                      fillalpha=0.25),
                       tag='mark_%region', redraw=False)

        viewer = self.channel.get_viewer('Ginga Plot')
        (x_lo, x_hi), (y_lo, y_hi) = viewer.get_ranges()

        for name, dct in self.marker_dct.items():
            x = dct['x']
            line = self.dc.Line(x, y_lo, x, y_hi, linewidth=2,
                                linestyle='dotted', arrow='start',
                                color='green')
            canvas.add(line, tag=f'mark_{name}', redraw=False)

            if 'x_fitted' in dct:
                # <-- there is a fit to be plotted
                points = np.array((dct['x_fitted'], dct['y_fitted'])).T
                path = self.dc.Path(points, linewidth=2, linestyle='solid',
                                    alpha=1.0, color='green')
                canvas.add(path, tag='mark_{name}_fit', redraw=False)


    def plot_lines(self):
        """Plot the line list + regions + markers.

        Lines are made into a single compound object so that it is easier
        to remove them as a group if the line list is changed.
        """
        if self.plot is None:
            return
        canvas = self.plot.get_canvas()
        canvas.delete_object_by_tag('lines', redraw=False)

        lines = self.dc.Canvas()
        if self.llist is not None:
            # NOTE: adapted straight from linetools
            x_min, x_max = self.data.x_min, self.data.x_max
            y_min, y_max = self.data.y_min, self.data.y_max

            z = self.z
            wvobs = np.array((1 + z) * self.llist['wrest'])
            ylbl_pos = y_max - 0.2 * (y_max - y_min)
            gdwv = np.where((wvobs > x_min) & (wvobs < x_max))[0]

            viewer = self.channel.get_viewer('Ginga Plot')
            (x_lo, x_hi), (y_lo, y_hi) = viewer.get_ranges()
            # make label always sit about 3/4 up the Y range
            y_lbl = y_lo + (y_hi - y_lo) * 0.75
            # line should reach to the label at least
            y_max = max(y_lbl, y_max)

            for kk in range(len(gdwv)):
                jj = gdwv[kk]
                wrest = self.llist['wrest'][jj]
                x_lbl = wrest * (z + 1)
                if not (x_lo < x_lbl < x_hi):
                    # skip plotting lines that are not visible
                    continue

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
                lbl = self.llist['name'][jj]
                lines.add(self.dc.Text(x_lbl, y_lbl, text=lbl, rot_deg=90,
                                       bgcolor='white', bgalpha=1.0,
                                       color='blue', fontsize=10),
                          redraw=False)

        # will be empty if there were no lines
        canvas.add(lines, tag='lines', redraw=False)

        self.plot_markers()

        # this causes the plot viewer to redraw itself
        with self.plot_stocker:
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

        # optionally, plot error vs. wavelength
        if self.settings.get('plot_error', False):
            points = np.array((self.data.wave, self.data.sig)).T
            p2 = self.dc.Path(points, linewidth=2, linestyle='solid',
                              alpha=1.0, color='red')
            canvas.add(p2, tag='error', redraw=False)

        self.plot.set_titles(title=f"Spec1dView: {self.ext_name}",
                             x_axis="Wavelength (Ang)", y_axis="Flux")
        self.plot.set_grid(True)

        self.plot_lines()

    def range_changed_cb(self, viewer, ranges):
        """Called when the plot viewers displayed range has changed.
        We use this to replot lines of interest.
        """
        if not self.plot_shelf.is_blocked():
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
        self.exten = min(self.start_ext, self.num_exten - 1)
        self.w.exten.set_index(self.exten)

        # create a new plotable to contain the plot
        self.plot = Plotable(logger=self.logger)
        self.plot.set(name="plot-{}".format(str(time.time())),
                      path=None, nothumb=True)
        self.recalc()

        # add the plot to this channel
        self.channel.add_image(self.plot)

        self.autozoom_plot()

    def recalc(self):
        """Reprocess the chosen extension, based on current choices for
        extraction method, fluxing and masking.

        Replot everything as a result.
        """
        specobj = self.sobjs[self.exten]

        # check if we have BOX or OPT extractions and adjust UI
        # accordingly
        if specobj['OPT_WAVE'] is None:
            self.w.extraction.set_text('BOX')
            self.extraction = 'BOX'
            self.w.extraction.set_enabled(False)
        elif specobj['BOX_WAVE'] is None:
            self.w.extraction.set_text('OPT')
            self.extraction = 'OPT'
            self.w.extraction.set_enabled(False)
        else:
            self.w.extraction.set_enabled(True)

        # look for OPT_FLAM_IVAR or BOX_FLAM_IVAR
        # if don't have, then fluxed==True cannot be used
        if specobj[f'{self.extraction}_FLAM_IVAR'] is None:
            self.w.fluxed.set_text('False')
            self.fluxed = False
            self.w.fluxed.set_enabled(False)
        else:
            self.w.fluxed.set_enabled(True)

        wave, flux, ivar, gpm = specobj.to_arrays(extraction=self.extraction,
                                                  fluxed=self.fluxed)
        if self.fluxed and ivar is None:
            # <-- fluxed=True cannot be used
            self.fv.show_error("fluxed=True cannot be used with this data")
            return

        sig = np.sqrt(utils.inverse(ivar))
        wave_gpm = wave > 1.0
        wave, flux, sig, gpm = (wave[wave_gpm], flux[wave_gpm],
                                sig[wave_gpm], gpm[wave_gpm])
        if self.masked:
            flux = flux*gpm
            sig = sig*gpm

        if self.nsmooth > 0.0:
            # do smoothing if requested
            flux = utils.convolve_psf(flux, self.nsmooth)
            # TODO: I don't think this is right.  I think it should be:
            #   sig = np.sqrt(utils.convolve_psf(sig**2, self.nsmooth))
            sig = utils.convolve_psf(sig, self.nsmooth)

        self.data.x_min, self.data.x_max = np.nanmin(wave), np.nanmax(wave)
        self.data.y_min, self.data.y_max = np.nanmin(flux), np.nanmax(flux)
        self.data.wave = wave
        self.data.flux = flux
        self.data.sig = sig

        self.ext_name = self.sobjs.NAME[self.exten]
        self.replot()

    def close(self):
        """Method called to close the plugin when the Close button is pressed."""
        self.fv.stop_local_plugin(self.chname, str(self))

    def start(self):
        """Method called right after `build_gui` when the plugin is activated.

        Simply attempt to process the latest FITS file loaded in the channel.
        """
        viewer = self.channel.get_viewer('Ginga Plot')
        bd = viewer.get_bindings()
        bd.set_mode(viewer, 'spec1d', mode_type='locked')

        self.redo()

    def stop(self):
        """Method called when the plugin is deactivated.

        Clean up instance variables so we don't hang on to any large data
        structures.
        """
        viewer = self.channel.get_viewer('Ginga Plot')
        viewer.clear()
        self.sobjs = 0
        self.exten = 0
        self.ext_name = ''
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

        self.autozoom_plot()

    def plot_error_cb(self, w, val):
        """Callback for toggling the "Plot Error" checkbox in the UI.
        """
        self.settings.set(plot_error=val)
        self.recalc()

    def autozoom_plot(self):
        if self.settings.get('autozoom', False):
            viewer = self.channel.get_viewer('Ginga Plot')
            viewer.zoom_fit()

    def set_params(self, ext=None, extraction=None, masked=None, fluxed=None):
        """Used to set up defaults from command line args to pypeit_show_1dspec script."""
        self.logger.info(f"ext={ext} extraction={extraction} masked={masked} fluxed={fluxed}")
        if ext is not None:
            self.start_ext = ext
        if extraction is not None:
            self.extraction = extraction
        if masked is not None:
            self.masked = masked
        if fluxed is not None:
            self.fluxed = fluxed

        if self.gui_up:
            index = self.extraction_types.index(self.extraction)
            self.w.extraction.set_index(index)
            index = self.fluxed_options.index(self.fluxed)
            self.w.fluxed.set_index(index)
            index = self.masked_options.index(self.masked)
            self.w.masked.set_index(index)

    def __str__(self):
        # necessary to identify the plugin and provide correct operation in Ginga
        return 'spec1dview'


# TODO: This should be moved somewhere more general.
def fit_gaussian_line(x, flux, i_low, i_high):
    """
    Fit a Gaussian spectral line with a local linear continuum.

    Parameters
    ----------
    x : array_like
        Monotonically increasing spectral coordinate (wavelength, frequency, etc.).
    flux : array_like
        Flux values.
    i_low : int
        Lower index (inclusive) of fitting window.
    i_high : int
        Upper index (exclusive) of fitting window.

    Returns
    -------
    model : astropy.modeling.CompoundModel
        Best-fit (Gaussian1D + Linear1D) model.
    fitter : astropy.modeling.fitting.LevMarLSQFitter
        Fitter instance containing fit diagnostics.
    """

    x = np.asarray(x)
    flux = np.asarray(flux)

    if x.ndim != 1 or flux.ndim != 1:
        raise ValueError("x and flux must be 1D arrays")
    if len(x) != len(flux):
        raise ValueError("x and flux must have the same length")
    if not (0 <= i_low < i_high <= len(x)):
        raise ValueError("Invalid index range")

    # Extract window
    x_fit = x[i_low:i_high]
    f_fit = flux[i_low:i_high]

    # --- Continuum initialization ---
    # Estimate continuum from window edges
    n_edge = max(1, int(0.1 * len(x_fit)))
    cont_level = np.median(
        np.concatenate([f_fit[:n_edge], f_fit[-n_edge:]])
    )

    cont_init = models.Linear1D(
        slope=0.0,
        intercept=cont_level
    )

    # --- Gaussian initialization ---
    peak_index = np.argmax(np.abs(f_fit - cont_level))
    mean_init = x_fit[peak_index]
    amplitude_init = f_fit[peak_index] - cont_level

    stddev_init = 0.25 * (x_fit[-1] - x_fit[0])

    gauss_init = models.Gaussian1D(
        amplitude=amplitude_init,
        mean=mean_init,
        stddev=stddev_init
    )

    # Compound model: continuum + line
    model_init = cont_init + gauss_init

    fitter = fitting.LevMarLSQFitter()
    model_fit = fitter(model_init, x_fit, f_fit)

    return model_fit, fitter
