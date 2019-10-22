import os
import copy
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
import matplotlib.transforms as mtransforms
from matplotlib.widgets import Button, Slider

matplotlib.use('Qt5Agg')

from pypeit.par import pypeitpar
from pypeit.core.wavecal import fitting, waveio, wvutils
from pypeit import msgs

operations = dict({'cursor': "Select lines (LMB click)\n" +
                    "         Select regions (LMB drag = add, RMB drag = remove)\n" +
                    "         Navigate (LMB drag = pan, RMB drag = zoom)",
                   'p' : "Toggle pan/zoom with the cursor",
                   '?' : "Display the available options",
                   'q' : "Close Object ID window and continue PypeIt reduction",
                   })


class ObjFindGUI(object):
    """
    GUI to interactively identify object traces. The GUI can be run within
    PypeIt during data reduction, or as a standalone script outside of
    PypeIt. To initialise the GUI, call the initialise() function in this
    file.
    """

    def __init__(self, canvas, image, axes):
        """Controls for the interactive Object ID tasks in PypeIt.

        The main goal of this routine is to interactively add/delete
        object apertures.

        Args:
            canvas (Matploltib figure canvas): The canvas on which all axes are contained
            axes (dict): Dictionary of four Matplotlib axes instances (Main spectrum panel, two for residuals, one for information)
        """
        # Store the axes
        self.image = image
        self.axes = axes
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

        # Interaction variables
        self._respreq = [False, None]  # Does the user need to provide a response before any other operation will be permitted? Once the user responds, the second element of this array provides the action to be performed.
        self._qconf = False  # Confirm quit message
        self._changes = False

        # Draw the spectrum
        self.canvas.draw()

    def print_help(self):
        """Print the keys and descriptions that can be used for Identification
        """
        keys = operations.keys()
        print("===============================================================")
        print("       OBJECT ID OPERATIONS")
        for key in keys:
            print("{0:6s} : {1:s}".format(key, operations[key]))
        print("---------------------------------------------------------------")

    def replot(self):
        """Redraw the entire canvas
        """
        self.canvas.restore_region(self.background)
        self.draw_objtraces()
        self.canvas.draw()

    def draw_objtraces(self):
        """Draw the lines and annotate with their IDs
        """
        for i in self.annlines: i.remove()
        for i in self.anntexts: i.remove()
        self.annlines = []
        self.anntexts = []
        # Plot the lines
        xmn, xmx = self.axes['main'].get_xlim()
        ymn, ymx = self.axes['main'].get_ylim()
        w = np.where((self._detns > xmn) & (self._detns < xmx))[0]
        for i in range(w.size):
            if self._lineflg[w[i]] != 1:
                if w[i] == self._detns_idx:
                    self.annlines.append(self.axes['main'].axvline(self._detns[w[i]], color='r'))
                else:
                    self.annlines.append(self.axes['main'].axvline(self._detns[w[i]], color='grey', alpha=0.5))
                continue
            else:
                if w[i] == self._detns_idx:
                    self.annlines.append(self.axes['main'].axvline(self._detns[w[i]], color='r'))
                else:
                    self.annlines.append(self.axes['main'].axvline(self._detns[w[i]], color='b'))
                txt = "{0:.2f}".format(self._lineids[w[i]])
                self.anntexts.append(
                    self.axes['main'].annotate(txt, (self._detns[w[i]], self._detnsy[w[i]]), rotation=90.0,
                                     color='b', ha='right', va='bottom'))

    def draw_callback(self, event):
        """Draw the lines and annotate with their IDs

        Args:
            event (Event): A matplotlib event instance
        """
        # Get the background
        self.background = self.canvas.copy_from_bbox(self.axes['main'].bbox)
        self.draw_objtraces()

    def get_ind_under_point(self, event):
        """Get the index of the line closest to the cursor

        Args:
            event (Event): Matplotlib event instance containing information about the event

        Returns:
            ind (int): Index of the spectrum where the event occurred
        """
        ind = np.argmin(np.abs(self.specx - event.xdata))
        return ind

    def get_axisID(self, event):
        """Get the ID of the axis where an event has occurred

        Args:
            event (Event): Matplotlib event instance containing information about the event

        Returns:
            axisID (int, None): Axis where the event has occurred
        """
        if event.inaxes == self.axes['main']:
            return 0
        elif event.inaxes == self.axes['info']:
            return 1
        return None

    def button_press_callback(self, event):
        """What to do when the mouse button is pressed

        Args:
            event (Event): Matplotlib event instance containing information about the event
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
        """What to do when the mouse button is released

        Args:
            event (Event): Matplotlib event instance containing information about the event

        Returns:
            None
        """
        if event.inaxes is None:
            return
        if event.inaxes == self.axes['info']:
            if (event.xdata > 0.8) and (event.xdata < 0.9):
                answer = "y"
            elif event.xdata >= 0.9:
                answer = "n"
            else:
                return
            self.operations(answer, -1)
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
                    pass
                elif self._end != self._start:
                    # The mouse button was dragged
                    if axisID == 0:
                        if self._start > self._end:
                            tmp = self._start
                            self._start = self._end
                            self._end = tmp
                        # Now do something
                        pass
        # Now plot
        self.canvas.restore_region(self.background)
        self.draw_objtraces()
        self.canvas.draw()

#self.image.set_clim(vmin=10, vmax=100)

    def key_press_callback(self, event):
        """What to do when a key is pressed

        Args:
            event (Event): Matplotlib event instance containing information about the event

        Returns:
            None
        """
        # Check that the event is in an axis...
        if not event.inaxes:
            return
        # ... but not the information box!
        if event.inaxes == self.axes['info']:
            return
        axisID = self.get_axisID(event)
        self.operations(event.key, axisID)

    def operations(self, key, axisID):
        """Canvas operations

        Args:
            key (str): Which key has been pressed
            axisID (int): The index of the axis where the key has been pressed (see get_axisID)
        """
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
                else:
                    return
            # Reset the info box
            self.update_infobox(default=True)
            return

        if key == '?':
            self.print_help()
        elif key == 'q':
            if self._changes:
                self.update_infobox(message="WARNING: There are unsaved changes!!\nPress q again to exit", yesno=False)
                self._qconf = True
            else:
                plt.close()
        self.canvas.draw()

    def update_infobox(self, message="Press '?' to list the available options",
                       yesno=True, default=False):
        """Send a new message to the information window at the top of the canvas

        Args:
            message (str): Message to be displayed
        """
        self.axes['info'].clear()
        if default:
            self.axes['info'].text(0.5, 0.5, "Press '?' to list the available options", transform=self.axes['info'].transAxes,
                          horizontalalignment='center', verticalalignment='center')
            self.canvas.draw()
            return
        # Display the message
        self.axes['info'].text(0.5, 0.5, message, transform=self.axes['info'].transAxes,
                      horizontalalignment='center', verticalalignment='center')
        if yesno:
            self.axes['info'].fill_between([0.8, 0.9], 0, 1, facecolor='green', alpha=0.5, transform=self.axes['info'].transAxes)
            self.axes['info'].fill_between([0.9, 1.0], 0, 1, facecolor='red', alpha=0.5, transform=self.axes['info'].transAxes)
            self.axes['info'].text(0.85, 0.5, "YES", transform=self.axes['info'].transAxes,
                          horizontalalignment='center', verticalalignment='center')
            self.axes['info'].text(0.95, 0.5, "NO", transform=self.axes['info'].transAxes,
                          horizontalalignment='center', verticalalignment='center')
        self.axes['info'].set_xlim((0, 1))
        self.axes['info'].set_ylim((0, 1))
        self.canvas.draw()

def initialise(frame, slit_left, slit_righ, sobjs, slit_ids=None):
    """Initialise the 'ObjFindGUI' window for interactive object tracing

        Args:
            frame (ndarray): Sky subtracted science image
            sobjs (SpecObjs): SpecObjs Class

        Returns:
            ofgui (ObjFindGUI): Returns an instance of the ObjFindGUI class
    """
    # This allows the input lord and rord to either be (nspec, nslit) arrays or a single
    # vectors of size (nspec)
    if slit_left.ndim == 2:
        nslit = slit_left.shape[1]
        lordloc = slit_left
        rordloc = slit_righ
    else:
        nslit = 1
        lordloc = slit_left.reshape(slit_left.size,1)
        rordloc = slit_righ.reshape(slit_righ.size,1)

    if slit_ids is None:
        slit_ids = [str(slit) for slit in np.arange(nslit)]

    # Determine the scale of the image
    med = np.median(frame)
    mad = np.median(np.abs(frame-med))
    vmin = med-3*mad
    vmax = med+3*mad

    # Add the main figure axis
    fig, ax = plt.subplots(figsize=(16, 9), facecolor="white")
    plt.subplots_adjust(bottom=0.05, top=0.85, left=0.05, right=0.85)
    image = ax.imshow(frame, aspect=frame.shape[1]/frame.shape[0], cmap = 'Greys', vmin=vmin, vmax=vmax)

    # Add an information GUI axis
    axinfo = fig.add_axes([0.15, .92, .7, 0.07])
    axinfo.get_xaxis().set_visible(False)
    axinfo.get_yaxis().set_visible(False)
    axinfo.text(0.5, 0.5, "Press '?' to list the available options", transform=axinfo.transAxes,
                horizontalalignment='center', verticalalignment='center')
    axinfo.set_xlim((0, 1))
    axinfo.set_ylim((0, 1))

    axes = dict(main=ax, info=axinfo)
    # Initialise the object finding window and display to screen
    fig.canvas.set_window_title('PypeIt - Object Tracing')
    ofgui = ObjFindGUI(fig.canvas, image, axes)
    plt.show()

    return ofgui


if __name__ == '__main__':
    dirname = "/Users/rcooke/Work/Research/vmp_DLAs/observing/WHT_ISIS_2019B/N1/wht_isis_blue_U/"
    frame = np.load(dirname+"frame.npy")
    slit_left = np.load(dirname + "slitleft.npy")
    slit_righ = np.load(dirname + "slitrigh.npy")
    initialise(frame, slit_left, slit_righ, None)
