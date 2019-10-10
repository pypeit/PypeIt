""" Module for finding patterns in arc line spectra
"""
import os
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
import matplotlib.transforms as mtransforms
from matplotlib.widgets import Button, Slider

matplotlib.use('Qt5Agg')

#from pypeit import utils
#from pypeit.core.wavecal import autoid
from pypeit import msgs


class Identify(object):
    """
    Manually identify arc lines.

    """

    def __init__(self, canvas, ax, axr, axi, spec, specres, detns, lines):
        """
        Description to go here
        """
        self.ax = ax
        self.axr = axr
        self.axi = axi
        # Initialise the spectrum properties
        self.spec = spec
        self.specres = specres   # Residual information
        self.specdata = spec.get_ydata()
        self.specx = np.arange(self.specdata.size)
        # Detections, linelist, and line IDs
        self._detns = detns
        self._detnsy = self.get_ann_ypos()  # Get the y locations of the annotations
        self._lines = lines
        self._lineids = np.zeros(self._detns.size, dtype=np.float)
        self._lineflg = np.zeros(self._detns.size, dtype=np.int)  # Flags: 0=no ID, 1=user ID, 2=auto ID, 3=flag reject
        # Fitting properties
        self._fitdict = dict(polyorder=1,
                             scale=self.specdata.size-1,
                             coeff=None)
        # Initialise the residuals colormap
        residcmap = LinearSegmentedColormap.from_list("my_list", ['grey', 'blue', 'orange', 'red'], N=4)
        self.residmap = ScalarMappable(norm=Normalize(vmin=0, vmax=3), cmap=residcmap)
        # Initialise the annotations
        self.annlines = []
        self.anntexts = []

        # Unset some of the matplotlib keymaps
        matplotlib.pyplot.rcParams['keymap.fullscreen'] = ''        # toggling fullscreen (Default: f, ctrl+f)
        matplotlib.pyplot.rcParams['keymap.home'] = ''              # home or reset mnemonic (Default: h, r, home)
        matplotlib.pyplot.rcParams['keymap.back'] = ''              # forward / backward keys to enable (Default: left, c, backspace)
        matplotlib.pyplot.rcParams['keymap.forward'] = ''           # left handed quick navigation (Default: right, v)
        #matplotlib.pyplot.rcParams['keymap.pan'] = ''              # pan mnemonic (Default: p)
        matplotlib.pyplot.rcParams['keymap.zoom'] = ''              # zoom mnemonic (Default: o)
        matplotlib.pyplot.rcParams['keymap.save'] = ''              # saving current figure (Default: s)
        matplotlib.pyplot.rcParams['keymap.quit'] = ''              # close the current figure (Default: ctrl+w, cmd+w)
        matplotlib.pyplot.rcParams['keymap.grid'] = ''              # switching on/off a grid in current axes (Default: g)
        matplotlib.pyplot.rcParams['keymap.yscale'] = ''            # toggle scaling of y-axes ('log'/'linear') (Default: l)
        matplotlib.pyplot.rcParams['keymap.xscale'] = ''            # toggle scaling of x-axes ('log'/'linear') (Default: L, k)
        matplotlib.pyplot.rcParams['keymap.all_axes'] = ''          # enable all axes (Default: a)

        # Initialise the main canvas tools
        canvas.mpl_connect('draw_event', self.draw_callback)
        canvas.mpl_connect('button_press_event', self.button_press_callback)
        canvas.mpl_connect('key_press_event', self.key_press_callback)
        canvas.mpl_connect('button_release_event', self.button_release_callback)
        self.canvas = canvas

        # Setup slider for the linelist
        self._slideval = 0  # Default starting point for the linelist slider
        self.linelist_init()

        # Interaction variables
        self._detns_idx = -1
        self._fitr = None  # Matplotlib shaded fit region (for refitting lines)
        self._fitregions = np.zeros(self.specdata.size, dtype=np.int)  # Mask of the pixels to be included in a fit
        self._addsub = 0   # Adding a region (1) or removing (0)
        self._respreq = [False, None]  # Does the user need to provide a response before any other operation will be permitted? Once the user responds, the second element of this array provides the action to be performed.
        self._qconf = False  # Confirm quit message
        # TODO :: This should be deleted
        self._changes = False

        # Draw the spectrum
        self.canvas.draw()

    def replot(self):
        self.canvas.restore_region(self.background)
        self.draw_lines()
        self.draw_residuals()
        self.canvas.draw()

    def linelist_update(self, val):
        val = int(round(val))
        self._slidell.label.set_text("{0:.4f}".format(self._lines[val]))
        self._slideval = val

    def linelist_select(self, event):
        if event.button == 1:
            self.update_line_id()
            self._detns_idx = -1
            # Try to perform a fit
            self.fitsol_fit()
            # Now replot everything
            self.replot()

    def linelist_init(self):
        axcolor = 'lightgoldenrodyellow'
        # Slider
        self.axl = plt.axes([0.15, 0.87, 0.7, 0.04], facecolor=axcolor)
        self._slidell = Slider(self.axl, "{0:.4f}".format(self._lines[self._slideval]), self._slideval,
                               self._lines.size-1, valinit=0, valstep=1)
        self._slidell.valtext.set_visible(False)
        self._slidell.on_changed(self.linelist_update)
        # Select button
        selax = plt.axes([0.86, 0.87, 0.1, 0.04])
        self._select = Button(selax, 'Select Line', color=axcolor, hovercolor='y')
        self._select.on_clicked(self.linelist_select)

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
                if w[i] == self._detns_idx:
                    self.annlines.append(self.ax.axvline(self._detns[w[i]], color='r'))
                else:
                    self.annlines.append(self.ax.axvline(self._detns[w[i]], color='grey', alpha=0.5))
                continue
            else:
                if w[i] == self._detns_idx:
                    self.annlines.append(self.ax.axvline(self._detns[w[i]], color='r'))
                else:
                    self.annlines.append(self.ax.axvline(self._detns[w[i]], color='b'))
            txt = "{0:.2f}".format(self._lineids[w[i]])
            self.anntexts.append(
                self.ax.annotate(txt, (self._detns[w[i]], self._detnsy[w[i]]), rotation=90.0,
                                 color='b', ha='right', va='bottom'))
        return

    def draw_residuals(self):
        if self._fitdict["coeff"] is None:
            nid = np.where(self._lineflg==1)[0].size
            msg = "Cannot plot residuals until more lines have been identified\n" +\
                  "Polynomial order = {0:d}, Number of line IDs = {1:d}".format(self._fitdict["polyorder"], nid)
            self.update_infobox(message=msg, yesno=False)
        else:
            # NOTE: self.specres = [respts, resfit, resres]

            # Perform model fits
            wavefit = self.fitsol_value(xfit=self.specx)
            wavevals = self.fitsol_value()
            resvals = (self._lineids - wavevals.copy()) / self.fitsol_deriv()

            # Set identified line wavelengths
            wid = self._lineids != 0.0
            wavevals[wid] = self._lineids[wid]

            # Pixel vs wavelength
            self.specres[0].set_offsets(np.vstack((self._detns, self._lineids)).T)
            self.specres[1].set_ydata(wavefit)
            self.axr[1].set_ylim((np.min(wavevals), np.max(wavevals)))
            self.specres[0].set_color(self.residmap.to_rgba(self._lineflg))

            # Pixel residuals
            self.specres[2].set_offsets(np.vstack((self._detns, resvals)).T)
            self.axr[0].set_ylim((np.min(resvals), np.max(resvals)))
            self.specres[2].set_color(self.residmap.to_rgba(self._lineflg))

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

    def get_axisID(self, event):
        if event.inaxes == self.ax:
            return 0
        elif event.inaxes == self.axr[0]:
            return 1
        elif event.inaxes == self.axr[1]:
            return 2
        elif event.inaxes == self.axi:
            return 3
        return None

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
        axisID = self.get_axisID(event)
        self._start = self.get_ind_under_point(event)

    def button_release_callback(self, event):
        """
        whenever a mouse button is released
        """
        if event.inaxes is None:
            return
        if event.inaxes == self.axi:
            if (event.xdata > 0.8) and (event.xdata < 0.9):
                answer = "y"
            elif event.xdata >= 0.9:
                answer = "n"
            else:
                return
            self.operations(answer, -1, -1)
            self.update_infobox(default=True)
            return
        elif self._respreq[0]:
            # The user is trying to do something before they have responded to a question
            return
        if self.canvas.toolbar.mode != "":
            return
        # Draw an actor
        axisID = self.get_axisID(event)
        if axisID is not None:
            if axisID <= 2:
                self._end = self.get_ind_under_point(event)
                if self._end == self._start:
                    # The mouse button was pressed (not dragged)
                    self._detns_idx = self.get_detns()
                    if self._fitdict['coeff'] is not None:
                        # Find closest line
                        waveest = self.fitsol_value(idx=self._detns_idx)
                        widx = np.argmin(np.abs(waveest-self._lines))
                        self.linelist_update(widx)
                elif self._end != self._start:
                    # The mouse button was dragged
                    if axisID == 0:
                        if self._start > self._end:
                            tmp = self._start
                            self._start = self._end
                            self._end = tmp
                        self.update_regions()
        # Now plot
        trans = mtransforms.blended_transform_factory(self.ax.transData, self.ax.transAxes)
        self.canvas.restore_region(self.background)
        self.draw_fitregions(trans)
        self.draw_lines()
        self.canvas.draw()

    def key_press_callback(self, event):
        """
        whenever a key is pressed
        """
        # Check that the event is in an axis...
        if not event.inaxes:
            return
        # ... but not the information box!
        if event.inaxes == self.axi:
            return
        axisID = self.get_axisID(event)
        self.operations(event.key, axisID)

    def operations(self, key, axisID):
        # Check if the user really wants to quit
        if key == 'q' and self._qconf:
            if self._changes:
                self.update_infobox(message="WARNING: There are unsaved changes!!\nPress q again to exit", yesno=False)
                self._qconf = True
            else:
                msgs.bug("Need to change this to kill and return the results to PypeIt")
                plt.close()
        elif self._qconf:
            self.update_infobox(default=True)
            self._qconf = False

        # Manage responses from questions posed to the user.
        if self._respreq[0]:
            if key != "y" and key != "n":
                return
            else:
                # Switch off the required response
                self._respreq[0] = False
                # Deal with the response
                if self._respreq[1] == "write":
                    # First remove the old file, and save the new one
                    msgs.work("Not implemented yet!")
                    self.write()
                else:
                    return
            # Reset the info box
            self.update_infobox(default=True)
            return

        if key == '?':
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
            print("+/-     : Raise/Lower the order of the fitting polynomial")
            print("---------------------------------------------------------------")
        elif key == 'a':
            msgs.work("Feature not yet implemented")
        elif key == 'c':
            msgs.work("Feature not yet implemented")
        elif key == 'd':
            msgs.work("Feature not yet implemented")
            self._fitdict['coeff'] = None
        elif key == 'f':
            self.fitsol_fit()
        elif key == 'q':
            if self._changes:
                self.update_infobox(message="WARNING: There are unsaved changes!!\nPress q again to exit", yesno=False)
                self._qconf = True
            else:
                plt.close()
        elif key == 'r':
            if self._detns_idx == -1:
                msgs.info("You must select a line first")
            elif self._fitr is None:
                msgs.info("You must select a fitting region first")
            else:
                msgs.work("Feature not yet implemented")
        elif key == 'w':
            self.write()
        elif key == 'z':
            msgs.work("Feature not yet implemented")
        elif key == '+':
            if self._fitdict["polyorder"] < 10:
                self._fitdict["polyorder"] += 1
                self.update_infobox(message="Polynomial order = {0:d}".format(self._fitdict["polyorder"]), yesno=False)
                self.fitsol_fit()
                self.replot()
            else:
                self.update_infobox(message="Polynomial order must be <= 10", yesno=False)
        elif key == '-':
            if self._fitdict["polyorder"] > 1:
                self._fitdict["polyorder"] -= 1
                self.update_infobox(message="Polynomial order = {0:d}".format(self._fitdict["polyorder"]), yesno=False)
                self.fitsol_fit()
                self.replot()
            else:
                self.update_infobox(message="Polynomial order must be >= 1", yesno=False)
        self.canvas.draw()

    def fitsol_value(self, xfit=None, idx=None):
        if xfit is None:
            xfit = self._detns
        if self._fitdict['coeff'] is not None:
            if idx is None:
                return np.polyval(self._fitdict["coeff"], xfit / self._fitdict["scale"])
            else:
                return np.polyval(self._fitdict["coeff"], xfit[idx] / self._fitdict["scale"])
        else:
            msgs.bug("Cannot predict wavelength value - no fit has been performed")
            return None

    def fitsol_deriv(self, xfit=None, idx=None):
        if xfit is None:
            xfit = self._detns
        if self._fitdict['coeff'] is not None:
            cder = np.polyder(self._fitdict["coeff"])
            if idx is None:
                return np.polyval(cder, xfit / self._fitdict["scale"]) / self._fitdict["scale"]
            else:
                return np.polyval(cder, xfit[idx] / self._fitdict["scale"]) / self._fitdict["scale"]
        else:
            msgs.bug("Cannot predict wavelength value - no fit has been performed")
            return None

    def fitsol_fit(self):
        ord = self._fitdict["polyorder"]
        wfit = np.where(self._lineflg == 1)  # Use the user IDs only!
        if ord+1 <= wfit[0].size:
            xpix = self._detns[wfit]/self._fitdict["scale"]
            ylam = self._lineids[wfit]
            self._fitdict["coeff"] = np.polyfit(xpix, ylam, ord)
        else:
            msg = "Polynomial order must be >= number of line IDs\n" +\
                  "Polynomial order = {0:d}, Number of line IDs = {1:d}".format(ord, wfit[0].size)
            self.update_infobox(message=msg, yesno=False)

    def update_infobox(self, message="Press '?' to list the available options",
                       yesno=True, default=False):
        self.axi.clear()
        if default:
            self.axi.text(0.5, 0.5, "Press '?' to list the available options", transform=self.axi.transAxes,
                          horizontalalignment='center', verticalalignment='center')
            self.canvas.draw()
            return
        # Display the message
        self.axi.text(0.5, 0.5, message, transform=self.axi.transAxes,
                      horizontalalignment='center', verticalalignment='center')
        if yesno:
            self.axi.fill_between([0.8, 0.9], 0, 1, facecolor='green', alpha=0.5, transform=self.axi.transAxes)
            self.axi.fill_between([0.9, 1.0], 0, 1, facecolor='red', alpha=0.5, transform=self.axi.transAxes)
            self.axi.text(0.85, 0.5, "YES", transform=self.axi.transAxes,
                          horizontalalignment='center', verticalalignment='center')
            self.axi.text(0.95, 0.5, "NO", transform=self.axi.transAxes,
                          horizontalalignment='center', verticalalignment='center')
        self.axi.set_xlim((0, 1))
        self.axi.set_ylim((0, 1))
        self.canvas.draw()

    def update_line_id(self):
        # Find the nearest wavelength in the linelist
        if self._detns_idx != -1:
            self._lineids[self._detns_idx] = self._lines[self._slideval]
            self._lineflg[self._detns_idx] = 1

    def update_regions(self):
        self._fitregions[self._start:self._end] = self._addsub

    def write(self):
        # TODO :: Write this to the master calibrations
        return


def initialise(arcspec, detns, lines):
    spec = Line2D(np.arange(arcspec.size), arcspec, linewidth=1, linestyle='solid', color='k', drawstyle='steps', animated=True)

    # Add the main figure axis
    fig, ax = plt.subplots(figsize=(16, 9), facecolor="white")
    plt.subplots_adjust(bottom=0.05, top=0.85, left=0.05, right=0.65)
    ax.add_line(spec)
    ax.set_ylim((0.0, 1.1*spec.get_ydata().max()))

    # Add two residual fitting axes
    axr = [fig.add_axes([0.7, .1, .28, 0.35]),
           fig.add_axes([0.7, .5, .28, 0.35])]
    # Residuals
    residcmap = LinearSegmentedColormap.from_list("my_list", ['grey', 'blue', 'orange', 'red'], N=4)
    resres = axr[0].scatter(detns, np.zeros(detns.size), marker='x',
                            c=np.zeros(detns.size), cmap=residcmap, norm=Normalize(vmin=0.0, vmax=3.0))
    axr[0].axhspan(-0.1, 0.1, alpha=0.5, color='grey')  # Residuals of 0.1 pixels
    axr[0].axhline(0.0, color='r', linestyle='-')  # Zero level
    #axr[0].add_line(resres)
    axr[0].set_xlim((0, arcspec.size-1))
    axr[0].set_ylim((-0.3, 0.3))

    # pixel vs wavelength
    respts = axr[1].scatter(detns, np.zeros(detns.size), marker='x',
                            c=np.zeros(detns.size), cmap=residcmap, norm=Normalize(vmin=0.0, vmax=3.0))
    resfit = Line2D(np.arange(arcspec.size), np.zeros(arcspec.size), linewidth=1, linestyle='-', color='r')
    axr[1].add_line(resfit)
    axr[1].set_xlim((0, arcspec.size-1))
    axr[1].set_ylim((-0.3, 0.3))  # This will get updated as lines are identified

    # Add an information GUI axis
    axi = fig.add_axes([0.15, .92, .7, 0.07])
    axi.get_xaxis().set_visible(False)
    axi.get_yaxis().set_visible(False)
    axi.text(0.5, 0.5, "Press '?' to list the available options", transform=axi.transAxes,
             horizontalalignment='center', verticalalignment='center')
    axi.set_xlim((0, 1))
    axi.set_ylim((0, 1))
    specres = [respts, resfit, resres]

    # Initialise the identify window and display to screen
    reg = Identify(fig.canvas, ax, axr, axi, spec, specres, detns, lines)
    plt.show()

    # Now return the results
    if reg is None:
        return
    return reg.get_results()


if __name__ == '__main__':
    # Start with a linelist and detections
    dir = "/Users/rcooke/Work/Research/vmp_DLAs/observing/WHT_ISIS_2019B/N1/wht_isis_blue_U/"
    spec = np.loadtxt(dir+"spec.dat")
    detns = np.loadtxt(dir+"detns.dat")
    lines = np.loadtxt(dir+"lines.dat")
    initialise(spec, detns, lines)
