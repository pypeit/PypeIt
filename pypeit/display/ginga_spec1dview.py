# This is open-source software licensed under a BSD license.
# Please see the file LICENSE.txt for details.
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

import numpy as np
from ginga import GingaPlugin, colors
from ginga.gw import Widgets
from ginga.util import viewer as gviewer
from pypeit import specobjs
from pypeit import utils
from linetools.lists import linelist


__all__ = ['Spec1dView']


class Spec1dView(GingaPlugin.LocalPlugin):

    def __init__(self, fv, fitsimage):
        # superclass defines some variables for us, like logger
        super(Spec1dView, self).__init__(fv, fitsimage)

        # get Spec1dView preferences
        prefs = self.fv.get_preferences()
        self.settings = prefs.create_category('plugin_Spec1dView')
        self.settings.add_defaults(lines="error")
        
        self.settings.load(onError='silent')

        self.w = None
        self.plot_viewer = None
        self.line_list = None
        self.linetypes = ('None', 'ISM', 'Strong', 'HI', 'H2', 'CO', 'EUV', 'Galaxy', 'AGN')

        self.lines = 'None'
        self.wave = None
        self.flux = None

        self.z = 0.0

        self.extraction_types = ('OPT', 'BOX')
        self.extraction = 'OPT'

        self.fluxed_options = (True, False)
        self.fluxed = False

        self.masked_options = (True, False)
        self.masked = True

        self.gui_up = False


    def build_gui(self, container):
        top = Widgets.VBox()
        top.set_border_width(4)

        vbox = Widgets.VBox()
        vbox.set_border_width(4)
        vbox.set_spacing(2)

        fr = Widgets.Frame("Spec1dView")

        captions = (("Lines:", 'label', 'Lines', 'combobox'),
                    ("Redshift z:", 'label', 'redshift', 'entry', "Update z", 'button'),
                    ("Extraction:", 'label', 'extraction', 'combobox'),
                    ("Fluxed:", 'label', 'fluxed', 'combobox'),
                    ("Masked:", 'label', 'masked', 'combobox'),
                    )
        
        w, b = Widgets.build_info(captions, orientation="vertical")
        self.w = b

        combobox = b.lines
        for name in self.linetypes:
            combobox.append_text(name)
        index = self.linetypes.index(self.lines)
        combobox.set_index(index)

        combobox.add_callback('activated', self.set_line_list_cb)
        b.update_z.add_callback('activated', self.update_z_cb)

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

        fr = Widgets.Frame("1D File")

        captions = (("FITS file:", 'label', 'fits_file', 'entry'),
                    ("Object:", 'label', 'Objects', 'spinbutton'),
                    ("Enter", 'button', "Clear", 'button'),
                    )

        w, b = Widgets.build_info(captions, orientation="vertical")
        self.w.update(b)

        b.enter.add_callback('activated', self.enter_cb)

        btn = b.clear
        btn.add_callback('activated', lambda w: self.clear())
        btn.set_tooltip("Clear the filename")

        spinbox = b.objects
        spinbox.set_limits(0, 0, incr_value=1)
        spinbox.add_callback('value-changed', self.set_objects_cb)

        fr.set_widget(w)
        vbox.add_widget(fr, stretch=0)

        spacer = Widgets.Label('')
        vbox.add_widget(spacer, stretch=1)

        top.add_widget(vbox, stretch=1)

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

    def enter_cb(self, w):
    
        fits_file_value = self.w.fits_file.get_text()
        print(f"Entered FITS file: {fits_file_value}")
 
        self.fv.load_file(fits_file_value, chname=self.channel.name)

    def set_objects_cb(self, w, val):
        self.set_objects(val)

    def update_z_cb(self, w):
        try:
            z_text = self.w.redshift.get_text().strip()
            if z_text == "":
                self.z = 0.0
                self.w.redshift.set_text("0.0")
            else:
                self.z = float(z_text)
            print(f"Updated redshift z: {self.z}")
            self.redo()
        except ValueError as e:
            self.logger.error(f"Invalid redshift value: {str(e)}")
            self.fv.show_error("Invalid redshift value")

    def set_objects(self, objectslimit):
        if hasattr(self, 'objectsmax') and objectslimit > self.objectsmax:
            raise Exception("Limit exceeds maximum value of %d" % (self.objectsmax))
        self.objectslimit = objectslimit
        print(f"Number of Objects: {objectslimit}")

    def set_line_list_cb(self, w, idx):
        self.lines = self.linetypes[idx]
        if self.lines == 'None':
            self.line_list = 'None'
            self.w.redshift.set_text("")
        else:
            self.line_list = linelist.LineList(self.lines)
        print(f"Loaded line list: {self.line_list}")
        
    def set_extraction_cb(self, w, idx):
        self.extraction = self.extraction_types[idx]        
        print(f"Selected extraction type: {self.extraction}")
        self.redo()

    def set_fluxed_cb(self, w, idx):
        self.fluxed = self.fluxed_options[idx]
           
        print(f"Selected fluxed option: {self.fluxed}")
        self.redo()

    def set_masked_cb(self, w, idx):
        self.masked = self.masked_options[idx]
        
        print(f"Selected masked option: {self.masked}")
        self.redo()
    
    def process_file(self, fits_file_value):   
        try:  
            print(f"Processing FITS file: {fits_file_value}")
            sobjs = specobjs.SpecObjs.from_fitsfile(fits_file_value, chk_version=False)
            
            extraction = self.extraction
            fluxed = self.fluxed     
            masked = self.masked 

            num_objects = len(sobjs)
            self.objectsmax = num_objects
            self.w.objects.set_limits(0, num_objects, incr_value=1)
            
            spinbutton_value = self.w.objects.get_value()

            if spinbutton_value > 0:
                exten = spinbutton_value - 1
            else:
                exten = 0

            wave, flux, ivar, gpm = sobjs[exten].to_arrays(extraction, fluxed)
                
            sig = np.sqrt(utils.inverse(ivar))
            wave_gpm = wave > 1.0
            wave, flux, sig, gpm = wave[wave_gpm], flux[wave_gpm], sig[wave_gpm], gpm[wave_gpm]
            if masked:
                flux = flux*gpm 
                sig = sig*gpm


            print("Processed data:")
            print("Wavelengths:", wave)
            print("Flux:", flux)
            print("Sigma:", sig)
            print("Good Pixel Mask:", gpm)

            self.fv.gui_do(self.plot_file, wave, flux, sig, self.z)
        
        except Exception as e:
            print(f"Error processing file: {str(e)}")
            self.fv.show_error("Failed to process file")

    
    def plot_file(self, wave, flux, sig, z):

        self.wave = wave
        self.flux = flux

        self.plot_viewer.clear_plot()

        if self.line_list != 'None':
            self.plot_spectral_lines(wave, flux, sig, z)

        self.plot_viewer.set_titles(x_axis="Wavelength (Ang)", y_axis="Flux")
        self.plot_viewer.set_grid(True)
        
        self.plot_viewer.plot_line(wave, flux, linewidth=1,
                                linestyle='-', color='black')

        self.plot_viewer.plot_line(wave, sig, linewidth=2,
                                linestyle='--', color='red')
        
        self.plot_viewer.redraw()

    def plot_spectral_lines(self, wave, flux, sig, z):
        if self.line_list is None:
            print("No line list selected")
            return

        llist = self.line_list
        y_min = min(np.min(flux), np.min(sig))
        y_max = max(np.max(flux), np.max(sig))

        x_min, x_max = np.min(wave), np.max(wave)

        print(f"Plotting lines for list: {self.lines}")
        print(f"Line list details: {llist}")
        print(f"Redshift: {z}")

        ylbl = y_max - 0.2 * (y_max - y_min)
        wvobs = np.array((1 + z) * llist.wrest)
        print(f"Observed wavelengths: {wvobs}")
        print(f"Plotting range x: [{x_min}, {x_max}]")

        gdwv = np.where((wvobs > x_min) & (wvobs < x_max))[0]
        print(f"Indices within range: {gdwv}, len: {len(gdwv)}")

        for kk in range(len(gdwv)):
            jj = gdwv[kk]
            wrest = llist.wrest[jj].value
            lbl = llist.name[jj]
            print(f"Plotting linelist: {lbl} at {wrest * (z + 1)}")
            self.plot_viewer.plot_line([wrest * (z + 1), wrest * (z + 1)], [y_min, y_max], color='blue', linestyle='--')
            self.plot_viewer.plot_text(wrest * (z + 1), ylbl, lbl, color='blue', rot_deg=90, fontsize=10)
    
    def clear(self):
        if self.gui_up:
            self.w.fits_file.set_text("")
            self.w.objects.set_value(0)
            self.w.lines.set_index(0)
            self.w.redshift.set_text("")
            self.lines = self.linetypes[0]
            self.line_list = self.linetypes[0]
            self.z = 0.0
            self.w.extraction.set_index(0)
            self.w.fluxed.set_index(1)
            self.w.masked.set_index(0)
            self.extraction = self.extraction_types[0]
            self.fluxed = self.fluxed_options[1]
            self.masked = self.masked_options[0]
    
    def close(self):
        self.fv.stop_local_plugin(self.chname, str(self))


    def start(self):
        self.plot_viewer = self.channel.get_viewer("Ginga Plot")
        vinfo = gviewer.get_vinfo("Ginga Plot")
        vinfo.priority = -1
        self.redo()

    def stop(self):
        self.gui_up = False

    def redo(self):
        if not self.gui_up:
            return
        dataobj = self.channel.get_current_image()
        if not self.plot_viewer.viewable(dataobj):
            self.fv.show_error("Failed to process file")
            return
        path = dataobj.get('path', None)
        if path is None:
            self.fv.show_error(
                "Cannot open dataobj: no value for metadata key 'path'")
            return
        self.w.fits_file.set_text(path)
        self.plot_viewer.clear_plot()
        self.process_file(path)
        

    
    def __str__(self):
        return 'spec1dview'

# END
