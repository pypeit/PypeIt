""" Module for finding patterns in arc line spectra
"""
import os
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.transforms as mtransforms
from matplotlib.widgets import TextBox

matplotlib.use('Qt5Agg')

#from pypeit import utils
#from pypeit.core.wavecal import autoid
from pypeit import msgs


class Identify(object):
    """
    Manually identify arc lines.

    """

    def __init__(self, canvas, ax, spec, detns, lines):
        """
        Description to go here
        """
        self.ax = ax
        # Initialise the spectrum properties
        self.spec = spec
        self.specdata = spec.get_ydata()
        self.specx = np.arange(self.specdata.size)
        # Detections, linelist, and line IDs
        self._detns = detns
        self._detnsy = self.get_ann_ypos()  # Get the y locations of the annotations
        self._lines = lines
        self._lineids = np.zeros(self._detns.size, dtype=np.float)
        self._lineflg = np.zeros(self._detns.size, dtype=np.int)  # Flags: 0=no ID, 1=user ID, 2=auto ID
        # Fitting properties
        self._fitdict = dict(polyorder=1,
                             scale=self.specdata.size-1,
                             coeff=None)
        # Initialise the annotations
        self.annlines = []
        self.anntexts = []

        # Initialise the main canvas tools
        canvas.mpl_connect('draw_event', self.draw_callback)
        canvas.mpl_connect('button_press_event', self.button_press_callback)
        canvas.mpl_connect('key_press_event', self.key_press_callback)
        canvas.mpl_connect('button_release_event', self.button_release_callback)
        self.canvas = canvas

        axbox = plt.axes([0.4, 0.05, 0.2, 0.05])
        self._id_entry = TextBox(axbox, 'Wavelength', initial="Select a line")
        self._id_entry.on_submit(self.update_line_id)

        # Interaction variables
        self._detns_idx = -1
        self._fitr = None  # Matplotlib shaded fit region (for refitting lines)
        self._fitregions = np.zeros(self.specdata.size, dtype=np.int)  # Mask of the pixels to be included in a fit
        self._addsub = 0   # Adding a region (1) or removing (0)
        # TODO :: This should be deleted
        self._changes = False

        # Draw the spectrum
        self.canvas.draw()

    def draw_lines(self):
        """ Draw the lines and annotate with their IDs
        """
        # annotations = [child for child in self.ax.get_children() if isinstance(child, matplotlib.text.Annotation)]
        for i in self.annlines: i.remove()
        for i in self.anntexts: i.remove()
        self.annlines = []
        self.anntexts = []
        # Plot the lines
        xmn, xmx = self.ax.get_xlim()
        ymn, ymx = self.ax.get_ylim()
        w = np.where((self._detns > xmn) & (self._detns < xmx))[0]
        for i in range(w.size):
            if self._lineids[w[i]] == 0.0:
                if i == self._detns_idx:
                    self.annlines.append(self.ax.axvline(self._detns[w[i]], color='r'))
                else:
                    self.annlines.append(self.ax.axvline(self._detns[w[i]], color='grey', alpha=0.5))
                continue
            else:
                if i == self._detns_idx:
                    self.annlines.append(self.ax.axvline(self._detns[w[i]], color='r'))
                else:
                    self.annlines.append(self.ax.axvline(self._detns[w[i]], color='b'))
            txt = "{0:.2f}".format(self._lineids[w[i]])
            self.anntexts.append(
                self.ax.annotate(txt, (self._detns[w[i]], self._detnsy[w[i]]), rotation=90.0,
                                 color='b', ha='center', va='bottom'))
        return

    def draw_callback(self, event):
        """ Draw the lines and annotate with their IDs

        Parameters
        ----------
        event : mpl event
          A matplotlib event instance
        """
        # Get the background
        self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        # Set the axis transform
        trans = mtransforms.blended_transform_factory(self.ax.transData, self.ax.transAxes)
        self.draw_fitregions(trans)
        self.ax.draw_artist(self.spec)
        self.draw_lines()

    def draw_fitregions(self, trans):
        """ Refresh the fit regions

        Parameters
        ----------
        trans : axis transform
          A matplotlib axis transform from data to axes coordinates
        """
        if self._fitr is not None:
            self._fitr.remove()
        # Find all regions
        regwhr = np.copy(self._fitregions == 1)
        # Fudge to get the leftmost pixel shaded in too
        regwhr[np.where((self._fitregions[:-1] == 0) & (self._fitregions[1:] == 1))] = True
        self._fitr = self.ax.fill_between(self.specx, 0, 1, where=regwhr, facecolor='green',
                                          alpha=0.5, transform=trans)

    def get_ann_ypos(self, scale=1.02):
        """ Calculate the y locations of the annotated IDs

        Parameters
        ----------
        scale : float
          scale the location relative to the maximum value of the spectrum

        Returns
        -------
        ypos : ndarray
          y locations of the annotations
        """
        ypos = np.zeros(self._detns.size)
        for xx in range(self._detns.size):
            wmin = np.argmin(np.abs(self.specx-self._detns[xx]))
            ypos[xx] = scale * np.max(self.specdata[wmin-1:wmin+2])
        return ypos

    def get_detns(self):
        """ Get the index of the detection closest to the cursor
        """
        return np.argmin(np.abs(self._detns-self.specx[self._end]))

    def get_ind_under_point(self, event):
        """
        Get the index of the line closest to the cursor
        """
        ind = np.argmin(np.abs(self.specx - event.xdata))
        return ind

    def button_press_callback(self, event):
        """
        whenever a mouse button is pressed
        """
        if event.inaxes is None:
            return
        if self.canvas.toolbar.mode != "":
            return
        if event.button == 1:
            self._addsub = 1
        elif event.button == 3:
            self._addsub = 0
        self._start = self.get_ind_under_point(event)

    def button_release_callback(self, event):
        """
        whenever a mouse button is released
        """
        if event.inaxes is None:
            return
        if self.canvas.toolbar.mode != "":
            return
        self._end = self.get_ind_under_point(event)
        if self._end != self._start:
            if self._start > self._end:
                tmp = self._start
                self._start = self._end
                self._end = tmp
            self.update_regions()
        elif self._end == self._start:
            self._detns_idx = self.get_detns()
            if self._fitdict['coeff'] is not None:
                txtval = "{0:.5f}".format(self.fitsol_value())
            else:
                txtval = ""
            self._id_entry.set_val(txtval)
        trans = mtransforms.blended_transform_factory(self.ax.transData, self.ax.transAxes)
        self.canvas.restore_region(self.background)
        self.draw_fitregions(trans)
        self.draw_lines()
        self.canvas.draw()

    def key_press_callback(self, event):
        """
        whenever a key is pressed
        """
        if not event.inaxes:
            return
        if event.key == '?':
            print("===============================================================")
            print("       MAIN OPERATIONS")
            print("cursor  : Select lines (LMB click)")
            print("          Select regions (LMB drag = add, RMB drag = remove)")
            print("          Navigate (LMB drag = pan, RMB drag = zoom)")
            print("p       : toggle pan/zoom with the cursor")
            print("w       : write the line IDs and spectrum")
            print("q       : continue PypeIt reduction")
            print("---------------------------------------------------------------")
            print("       ARC LINE OPERATIONS")
            print("a       : Automatically identify lines using current solution")
            print("c       : Clear automatically identified lines")
            print("d       : Delete all line identifications (start from scratch)")
            print("f       : Fit the wavelength solution")
            print("r       : Refit a line")
            print("z       : Delete a single line identification")
            print("---------------------------------------------------------------")
        elif event.key == 'a':
            msgs.work("Feature not yet implemented")
        elif event.key == 'c':
            msgs.work("Feature not yet implemented")
        elif event.key == 'd':
            msgs.work("Feature not yet implemented")
            self._fitdict['coeff'] = None
        elif event.key == 'f':
            msgs.work("Feature not yet implemented")
        elif event.key == 'r':
            if self._detns_idx == -1:
                msgs.info("You must select a line first")
            elif self._fitr is None:
                msgs.info("You must select a fitting region first")
            else:
                msgs.work("Feature not yet implemented")
        elif event.key == 'w':
            self.write()
        elif event.key == 'q':
            if self._changes:
                print("WARNING: There are unsaved changes!!")
                print("Press q again to exit")
                self._changes = False
            else:
                sys.exit()
        elif event.key == 'z':
            msgs.work("Feature not yet implemented")
        self.canvas.draw()

    def fitsol_value(self):
        return np.polyval(self._fitdict["coeff"], self._detns[self._detns_idx]/self._fitdict["scale"])

    def fitsol_fit(self):
        ord = self._fitdict["polyorder"]
        wfit = np.where(self._lineflg == 1)  # User IDs only!
        xpix = self._detns[wfit]/self._fitdict["scale"]
        ylam = self._lineids[wfit]
        self._fitdict["coeff"] = np.polyfit(xpix, ylam, ord)

    def update_line_id(self, text):
        try:
            wdata = float(text)
            # Find the nearest wavelength in the linelist
            idx = np.argmin(np.abs(self._lines - wdata))
            self._lineids[idx] = self._lines[idx]
            self._lineflg[idx] = 1
        except ValueError:
            msgs.info("Invalid entry in Line ID box")

    def update_regions(self):
        self._fitregions[self._start:self._end] = self._addsub

    def write(self):
        # TODO :: Write this to the master calibrations
        return


def initialise(spec, detns, lines):
    spec = Line2D(np.arange(spec.size), spec, linewidth=1, linestyle='solid', color='k', drawstyle='steps', animated=True)

    fig, ax = plt.subplots(figsize=(16, 9), facecolor="white")
    ax.add_line(spec)
    reg = Identify(fig.canvas, ax, spec, detns, lines)

    ax.set_title("Press '?' to list the available options in a terminal")
    ax.set_ylim((0.0, 1.1*spec.get_ydata().max()))
    plt.show()


if __name__ == '__main__':
    # Start with a linelist and detections
    dir = "/Users/rcooke/Work/Research/vmp_DLAs/observing/WHT_ISIS_2019B/N1/wht_isis_blue_U/"
    spec = np.loadtxt(dir+"spec.dat")
    detns = np.loadtxt(dir+"detns.dat")
    lines = np.loadtxt(dir+"lines.dat")
    initialise(spec, detns, lines)
