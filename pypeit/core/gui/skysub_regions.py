"""
This script allows the user to manually select the sky background regions
"""

import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import matplotlib.transforms as mtransforms

from pypeit import msgs
from pypeit.core import skysub
from pypeit.images import buildimage

operations = dict({'cursor': "Add sky region (LMB drag)\n" +
                   "         Remove sky region (RMB drag)\n" +
                   "         Note: If you would like to pan or zoom, you need to activate\n" +
                   "         the pan/zoom tool with the 'p' key, or by selecting the\n" +
                   "         pan/zoom tool on the Matplotlib navigation tool menu. You\n" +
                   "         can also zoom using the magnifying glass (select this option\n" +
                   "         from the Matplotlib navigation tool menu). While the pan/zoom\n" +
                   "         feature is enabled, you will not be able to update sky regions.\n",
                   'c': "Center the window at the location of the mouse",
                   'd': "Delete all sky regions and start again",
                   'h/r': "Return zoom to the original plotting limits",
                   'p': "Toggle pan/zoom with the cursor",
                   '?': "Display the available options",
                   })


class SkySubGUI:
    """
    GUI to interactively define the sky regions. The GUI can be run within
    PypeIt during data reduction, or as a standalone script outside of
    PypeIt. To initialize the GUI, call the initialize() function in this
    file.
    """

    def __init__(self, canvas, image, frame, outname, det, slits, axes, pypeline, spectrograph, printout=False,
                 runtime=False, resolution=None, initial=False, flexure=None, overwrite=False):
        """Controls for the interactive sky regions definition tasks in PypeIt.

        The main goal of this routine is to interactively select sky background
        regions.

        Parameters
        ----------
        canvas : Matploltib figure canvas
            The canvas on which all axes are contained
        image : AxesImage
            The image plotted to screen
        frame : ndarray
            The image data
        outname : str
            The output filename to save the sky regions mask
        det : int
            Detector to add a slit on
        slits : :class:`~pypeit.slittrace.SlitTraceSet`
            Object with the image coordinates of the slit edges
        axes : dict
            Dictionary of four Matplotlib axes instances (Main
            spectrum panel, two for residuals, one for information)
        pypeline : str
            Name of the instrument pipeline
        spectrograph : str
            Name of the spectrograph
        printout : bool
            Should the results be printed to screen
        initial : bool, optional
            To use the initial edges regardless of the presence of
            the tweaked edges, set this to True.
        flexure : float, optional
            If provided, offset each slit by this amount
        runtime : bool
            Is the GUI being launched during data reduction?
        resolution : int
            The resolution of the skysub definitions. It is the
            number of pixels to divide the slit width by (i.e. 1000
            pixels means a resolution of 0.1% of the slit width).
        """
        # Store the axes
        self._det = det
        self.image = image
        self.frame = frame
        self.pypeline = pypeline
        self.spectrograph = spectrograph
        self._outname = outname
        self._overwrite = overwrite
        self.nspec, self.nspat = frame.shape[0], frame.shape[1]
        self._spectrace = np.arange(self.nspec)
        self._printout = printout
        self._runtime = runtime
        self.axes = axes
        self._currslit = -1
        self.slits = slits
        self._nslits = slits.nslits
        self._maxslitlength = np.max(self.slits.get_slitlengths(initial=initial))
        self._resolution = int(10.0 * self._maxslitlength) if resolution is None else int(resolution)
        self._allreg = np.zeros(int(self._resolution), dtype=bool)
        self._specx = np.arange(int(self._resolution))
        self._start = [0, 0]
        self._end = [0, 0]

        # Unset some of the matplotlib keymaps
        for key in plt.rcParams.keys():
            if 'keymap' in key:
                plt.rcParams[key] = []
        # Enable some useful ones, though
        matplotlib.pyplot.rcParams['keymap.home'] = ['h', 'r', 'home']
        matplotlib.pyplot.rcParams['keymap.pan'] = ['p']

        # Initialise the main canvas tools
        canvas.mpl_connect('draw_event', self.draw_callback)
        canvas.mpl_connect('button_press_event', self.button_press_callback)
        canvas.mpl_connect('key_press_event', self.key_press_callback)
        canvas.mpl_connect('button_release_event', self.button_release_callback)
        canvas.mpl_connect('motion_notify_event',  self.mouse_move_callback)
        self.canvas = canvas

        # Interaction variables

        # Does the user need to provide a response before any other
        # operation will be permitted? Once the user responds, the
        # second element of this array provides the action to be
        # performed.
        self._respreq = [False, None]

        self._qconf = False  # Confirm quit message
        self._changes = False
        self._use_updates = True
        self._inslit = -1  # Which slit is the mouse in
        self.mmx, self.mmy = 0, 0
        self._fitr = []  # Matplotlib shaded fit region
        self._fita = None

        self.slits_left, self.slits_right, _ = slits.select_edges(initial=initial, flexure=flexure)
        self.initialize_menu()
        self.reset_regions()

        # Draw the spectrum
        self.canvas.draw()

#        self.reset_regions()

    @classmethod
    def initialize(cls, det, frame, slits, pypeline, spectrograph, outname="skyregions.fits",
                   overwrite=False, initial=False,
                   flexure=None, runtime=False, printout=False):
        """
        Initialize the 'ObjFindGUI' window for interactive object tracing

        Parameters
        ----------
        det : int
            Detector index
        frame : ndarray
            Sky subtracted science image
        slits : :class:`~pypeit.slittrace.SlitTraceSet`
            Object with the image coordinates of the slit edges
        pypeline : str
            Name of the reduction pipeline
        spectrograph : str
            Name of the spectrograph
        printout : bool
            Should the results be printed to screen
        runtime : bool
            Is this GUI being launched during a data reduction?

        Returns
        -------
        srgui : :class:`SkySubGUI`
            Returns an instance of the :class:`SkySubGUI` class
        """
        # NOTE: SlitTraceSet objects always store the left and right
        # traces as 2D arrays, even if there's only one slit.
        nslit = slits.nslits
        lordloc, rordloc, _ = slits.select_edges(initial=initial, flexure=flexure)

        # Determine the scale of the image
        med = np.median(frame)
        mad = np.median(np.abs(frame - med))
        vmin = med - 3 * mad
        vmax = med + 3 * mad

        # Add the main figure axis
        fig, ax = plt.subplots(figsize=(16, 9), facecolor="white")
        plt.subplots_adjust(bottom=0.05, top=0.85, left=0.05, right=0.8)
        image = ax.imshow(frame, aspect='auto', cmap='Greys', vmin=vmin, vmax=vmax)

        # Overplot the slit traces
        specarr = np.arange(lordloc.shape[0])
        for sl in range(nslit):
            ax.plot(lordloc[:, sl], specarr, 'g-')
            ax.plot(rordloc[:, sl], specarr, 'b-')

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
        fig.canvas.manager.set_window_title('PypeIt - Sky regions')
        srgui = SkySubGUI(fig.canvas, image, frame, outname, det, slits, axes, pypeline, spectrograph,
                          printout=printout, runtime=runtime, initial=initial, flexure=flexure, overwrite=overwrite)
        plt.show()

        return srgui

    def finalize(self):
        plt.rcdefaults()
        plt.close()

    def region_help(self):
        print("You can enter the regions in the text box, as a comma separated")
        print("list of percentages. For example, typing  :10,35:65,80:  in the")
        print("text box and pressing enter will add sky regions to the left 10%,")
        print("the inner 30%, and the right 20% of each slit.")
        print("")

    def print_help(self):
        """Print the keys and descriptions that can be used for Identification
        """
        keys = operations.keys()
        print("===============================================================")
        print("Define the sky background regions in each slit by using the left")
        print("mouse button to click and drag over the sky background region.")
        print("Use the right mouse button (click and drag) to delete a region.")
        print("If you click 'Continue (and save changes)' the sky background")
        print("regions file will be saved to the Calibrations directory.")
        print("")
        print("To assign regions to all slits simultaneously, click and drag")
        print("over the gray regions on the right toolbar. Alternatively,")
        self.region_help()
        print("thin green/blue lines  = slit edges")
        print("thin green/blue lines  = slit edges")
        print("")
        print("thin green/blue lines  = slit edges")
        print("shaded red regions     = selected sky regions")
        print("===============================================================")
        print("       OTHER OPERATIONS")
        for key in keys:
            print("{0:6s} : {1:s}".format(key, operations[key]))
        print("---------------------------------------------------------------")

    def initialize_menu(self):
        """Initialize the menu buttons
        """
        axcolor = 'lightgoldenrodyellow'
        # Continue with reduction (using updated specobjs)
        ax_cont = plt.axes([0.82, 0.85, 0.15, 0.05])
        self._ax_cont = Button(ax_cont, "Continue (and save changes)", color=axcolor,
                               hovercolor='y')
        self._ax_cont.on_clicked(self.button_cont)
        # Continue with reduction (using original specobjs)
        ax_exit = plt.axes([0.82, 0.79, 0.15, 0.05])
        self._ax_exit = Button(ax_exit, "Continue (don't save changes)", color=axcolor,
                               hovercolor='y')
        self._ax_exit.on_clicked(self.button_exit)
        # Frame for the axis
        self.axes['allslitreg'] = plt.axes([0.82, 0.68, 0.15, 0.04], facecolor='black',
                                           title="Assign sky regions to all slits")
        self.axes['allslitreg'].get_xaxis().set_ticks([])
        self.axes['allslitreg'].get_yaxis().set_ticks([])
        self.axes['allslitreg'].axvspan(0, self._resolution-1, color='lightgrey')
        self.axes['allslitreg'].set_xlim(-self._resolution/10,
                                         self._resolution + self._resolution/10)
        self.axes['allslitreg'].set_ylim(0, 1)
        # Text box
        ax_regb = plt.axes([0.82, 0.62, 0.15, 0.05])
        self._ax_regb = Button(ax_regb, "Enter regions", color=axcolor, hovercolor='y')
        self._ax_regb.on_clicked(self.button_regb)

    def button_regb(self, event):
        """Allow for text to be entered in the terminal. If the
        text is valid, and then apply to all slits. The text should
        be a comma separated list of percentages to apply to all slits
        Example: The following string   :10, 35:65, 80:
        would select the first 10%, the inner 30%, and the final 20% of all slits

        Args:
            event : Event
                A matplotlib event instance
        """
        self.update_infobox(message='Enter regions in the terminal (see terminal for help)', yesno=False)
        print("")
        self.region_help()
        print("To exit this tool, enter no text, and press enter.")
        while True:
            print("")
            text = input("Enter the regions: ")
            status, reg = skysub.read_userregions(text, self._nslits, self._maxslitlength)
            if status == 0:
                print("Regions parsed successfully")
                for rr in range(len(reg)):
                    # Apply to all slits
                    for sl in range(self._nslits):
                        self._skyreg[sl] = reg[rr].copy()
                self.replot()
                break
            elif status == 1:
                print('Region definition should be a comma-separated list of percentages (see help)')
            elif status == 2:
                print('Region definition should contain a semi-colon (see help)')
            # Break out of the loop if no text is entered
            if len(text.strip()) == 0: break
        return

    def button_cont(self, event):
        """What to do when the 'exit and save' button is clicked

        Args:
            event : Event
                A matplotlib event instance
        """
        self._respreq = [True, "exit_update"]
        self.update_infobox(message='Are you sure you want to exit and save the newly defined '
                                    'sky regions?', yesno=True)

    def button_exit(self, event):
        """What to do when the 'exit and do not save changes' button is clicked

        Args:
            event : Event
                A matplotlib event instance
        """
        self._respreq = [True, "exit_restore"]
        self.update_infobox(message='Are you sure you want to exit without saving the  sky '
                                    'regions?', yesno=True)

    def replot(self):
        """Redraw the entire canvas
        """
        self.canvas.restore_region(self.background)
        self.draw_regions()
        self.canvas.draw()

    def draw_regions(self):
        """Refresh the fit regions
        """
        # Remove the regions and reset the patches
        for rr in range(len(self._fitr)):
            self._fitr[rr].remove()
        if self._fita is not None:
            self._fita.remove()
        self._fitr = []
        # Loop through all slits:
        for sl in range(self._nslits):
            # Fill fraction of the slit
            diff = self.slits_right[:, sl] - self.slits_left[:,sl]
            tmp = np.zeros(self._resolution+2)
            tmp[1:-1] = self._skyreg[sl]
            wl = np.where(tmp[1:] > tmp[:-1])[0]
            wr = np.where(tmp[1:] < tmp[:-1])[0]
            for rr in range(wl.size):
                left = self.slits_left[:, sl] + wl[rr]*diff/(self._resolution-1.0)
                righ = self.slits_left[:, sl] + wr[rr]*diff/(self._resolution-1.0)
                self._fitr.append(self.axes['main'].fill_betweenx(self._spectrace, left, righ,
                                                                  facecolor='red', alpha=0.5))
        # Plot the region on top of the "all slits" panel
        trans = mtransforms.blended_transform_factory(self.axes['allslitreg'].transData,
                                                      self.axes['allslitreg'].transAxes)
        self._fita = self.axes['allslitreg'].fill_between(self._specx, 0, 1, transform=trans,
                                                          where=self._allreg, facecolor='red',
                                                          alpha=0.5, zorder=10)

    def draw_callback(self, event):
        """Draw callback (i.e. everytime the canvas is being drawn/updated)

        Args:
            event : Event
                A matplotlib event instance
        """
        # Get the background
        self.background = self.canvas.copy_from_bbox(self.axes['main'].bbox)
        self.draw_regions()

    def get_current_slit(self, event):
        """Get the index of the slit closest to the cursor

        Args:
            event : Event
                Matplotlib event instance containing information about the event
        """
        # Find the current slit
        self._currslit = -1
        yv = np.argmin(np.abs(event.ydata-self._spectrace))
        wsl = np.where((event.xdata > self.slits_left[yv, :]) &
                       (event.xdata < self.slits_right[yv, :]))[0]
        # Double check there's only one solution
        if wsl.size == 1:
            self._currslit = int(wsl[0])
        return

    def get_axisID(self, event):
        """Get the ID of the axis where an event has occurred

        Args:
            event : Event
                Matplotlib event instance containing information about the event

        Returns:
            int, None: Axis where the event has occurred
        """
        if event.inaxes == self.axes['main']:
            return 0
        elif event.inaxes == self.axes['info']:
            return 1
        elif event.inaxes == self.axes['allslitreg']:
            return 2
        return None

    def mouse_move_callback(self, event):
        """Store the locations of mouse as it moves across the canvas
        """
        if event.inaxes is None:
            return
        if event.inaxes == self.axes['main']:
            self.mmx, self.mmy = event.xdata, event.ydata

    def button_press_callback(self, event):
        """What to do when the mouse button is pressed

        Args:
            event : Event
                Matplotlib event instance containing information about the event
        """
        if event.inaxes is None:
            return
        if self.canvas.toolbar.mode != "":
            return
        if event.button == 1:
            self._addsub = 1
            self.get_current_slit(event)
        elif event.button == 3:
            self._addsub = 0
            self.get_current_slit(event)
        if event.inaxes == self.axes['main']:
            self._start = [event.xdata, event.ydata]
        elif event.inaxes == self.axes['allslitreg']:
            self._start = [int(round(event.xdata)), event.ydata]

    def button_release_callback(self, event):
        """What to do when the mouse button is released

        Args:
            event : Event
                Matplotlib event instance containing information about the event
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
            if axisID == 2:
                # Set the slide l value
                self._end = [int(round(event.xdata)), event.ydata]
                self.add_region_all()
            else:
                self._end = [event.xdata, event.ydata]
                if (self._end[0] == self._start[0]) and (self._end[1] == self._start[1]):
                    # The mouse button was pressed (not dragged)
                    pass
                elif self._end != self._start:
                    # The mouse button was dragged
                    if axisID == 0:
                        if self._start[0] > self._end[0]:
                            self._start[0], self._end[0] = self._end[0], self._start[0]
                        # Now do something
                        self.add_region()
        self.replot()

    def key_press_callback(self, event):
        """What to do when a key is pressed

        Args:
            event : Event
                Matplotlib event instance containing information about the event
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
            key : str
                Which key has been pressed
            axisID : int
                The index of the axis where the key has been pressed (see get_axisID)
        """
        # Check if the user really wants to quit
        if key == 'q' and self._qconf:
            if self._changes:
                self.update_infobox(message='WARNING: There are unsaved changes!!\nPress q '
                                            'again to exit', yesno=False)
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
                if self._respreq[1] == "exit_update" and key == "y":
                    self._use_updates = True
                    self.operations("qu", None)
                elif self._respreq[1] == "exit_restore" and key == "y":
                    self._use_updates = False
                    self.operations("qr", None)
                else:
                    return
            # Reset the info box
            self.update_infobox(default=True)
            return

        if key == '?':
            self.print_help()
        elif key == 'd':
            if axisID == 0:
                # If this is pressed on the main window
                self.reset_regions()
        elif key == 'c':
            if axisID == 0:
                # If this is pressed on the main window
                self.recenter()
        elif key == 'qu' or key == 'qr':
            if self._changes:
                self.update_infobox(message='WARNING: There are unsaved changes!!\nPress q '
                                            'again to exit', yesno=False)
                self._qconf = True
            else:
                plt.close()
        self.replot()

    def get_result(self):
        """Save a mask containing the skysub regions, and print information
        for what the user should include in their .pypeit file
        """
        # Only do this if the user wishes to save the result
        if self._use_updates:
            # Generate the mask
            inmask = skysub.generate_mask(self.pypeline, self._skyreg, self.slits, self.slits_left, self.slits_right)
            # Save the mask
            outfil = self._outname
            if os.path.exists(self._outname) and not self._overwrite:
                outfil = 'temp.fits'
                msgs.warn("File exists:\n{0:s}\nSaving regions to 'temp.fits'")
                self._overwrite = True
            msskyreg = buildimage.SkyRegions(image=inmask.astype(float), PYP_SPEC=self.spectrograph)
            msskyreg.to_file(file_path=outfil)

    def recenter(self):
        xlim = self.axes['main'].get_xlim()
        ylim = self.axes['main'].get_ylim()
        xmin = self.mmx - 0.5*(xlim[1]-xlim[0])
        xmax = self.mmx + 0.5*(xlim[1]-xlim[0])
        ymin = self.mmy - 0.5*(ylim[1]-ylim[0])
        ymax = self.mmy + 0.5*(ylim[1]-ylim[0])
        self.axes['main'].set_xlim([xmin, xmax])
        self.axes['main'].set_ylim([ymin, ymax])

    def update_infobox(self, message="Press '?' to list the available options",
                       yesno=True, default=False):
        """Send a new message to the information window at the top of the canvas

        Args:
            message : str
                Message to be displayed
            yesno : bool
                Is a yes/no option desired?
            default : bool
                Would you like to refresh the info box and just display the default message
        """
        self.axes['info'].clear()
        if default:
            self.axes['info'].text(0.5, 0.5, "Press '?' to list the available options",
                                   transform=self.axes['info'].transAxes,
                                   horizontalalignment='center', verticalalignment='center')
            self.canvas.draw()
            return
        # Display the message
        self.axes['info'].text(0.5, 0.5, message, transform=self.axes['info'].transAxes,
                      horizontalalignment='center', verticalalignment='center')
        if yesno:
            self.axes['info'].fill_between([0.8, 0.9], 0, 1, facecolor='green', alpha=0.5,
                                           transform=self.axes['info'].transAxes)
            self.axes['info'].fill_between([0.9, 1.0], 0, 1, facecolor='red', alpha=0.5,
                                           transform=self.axes['info'].transAxes)
            self.axes['info'].text(0.85, 0.5, "YES", transform=self.axes['info'].transAxes,
                          horizontalalignment='center', verticalalignment='center')
            self.axes['info'].text(0.95, 0.5, "NO", transform=self.axes['info'].transAxes,
                          horizontalalignment='center', verticalalignment='center')
        self.axes['info'].set_xlim((0, 1))
        self.axes['info'].set_ylim((0, 1))
        self.canvas.draw()

    def add_region(self):
        """ Add/subtract a defined region
        """
        # Figure out the locations of the start values
        ys = np.argmin(np.abs(self._start[1]-self._spectrace))
        difs = self.slits_right[ys, self._currslit] - self.slits_left[ys, self._currslit]
        sval = (self._start[0]-self.slits_left[ys, self._currslit]) / difs
        sidx = int(round(self._resolution*sval))
        # Figure out the locations of the start values
        yf = np.argmin(np.abs(self._end[1]-self._spectrace))
        diff = self.slits_right[yf, self._currslit] - self.slits_left[yf, self._currslit]
        fval = (self._end[0]-self.slits_left[yf, self._currslit]) / diff
        fidx = int(round(self._resolution*fval))
        # Switch the indices if needed
        if sidx > fidx:
            sidx, fidx = fidx, sidx
        # Check that we are within bounds
        if sidx < 0:
            sidx = 0
        if fidx > self._resolution:
            fidx = self._resolution
        # Assign the sky regions
        self._skyreg[self._currslit][sidx:fidx] = self._addsub
        # If some regions are removed, remove this from the "all slits" regions, as well
        if self._addsub == 0:
            self._allreg[sidx:fidx] = 0

    def add_region_all(self):
        """ Set the sky regions for all slits simultaneously
        """
        # Do some checks
        xmin, xmax = self._start[0], self._end[0]
        if xmax < xmin:
            xmin, xmax = xmax, xmin
        if xmin < 0:
            xmin = 0
        if xmax > self._resolution:
            xmax = self._resolution
        # Apply to all slits
        for sl in range(self._nslits):
            self._skyreg[sl][xmin:xmax] = self._addsub
        # Set the all regions parameter
        self._allreg[xmin:xmax] = self._addsub

    def reset_regions(self):
        """ Reset the sky regions for all slits simultaneously
        """
        self._skyreg = [np.zeros(self._resolution, dtype=bool) for all in range(self._nslits)]
        self._allreg[:] = False

