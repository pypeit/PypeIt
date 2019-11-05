"""
This script allows the user to add/delete/modify object traces

.. todo::

    Implement color scaling with RMB click+drag

"""

import os, sys
import copy
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, RadioButtons
from scipy.interpolate import RectBivariateSpline

matplotlib.use('Qt5Agg')

from pypeit import specobjs
from pypeit import msgs

operations = dict({'cursor': "Select object trace (LMB click)\n" +
                    "         Navigate (LMB drag = pan, RMB drag = zoom)",
                   'a': "Add a new object trace using the selected method",
                   'd': "Delete selected object trace",
                   'h/r': "Return zoom to the original plotting limits",
                   'p': "Toggle pan/zoom with the cursor",
                   '?': "Display the available options",
                   })
#                   Commands that could be used for manual object location
#                   'c': "Clear the anchor points and start again",
#                   'm': "Insert a fitting anchor for manual object tracing",
#                   'n': "Delete the fitting anchor nearest to the cursor",
#                   '+/-': "Raise/Lower the order of the fitting polynomial"


class ObjFindGUI(object):
    """
    GUI to interactively identify object traces. The GUI can be run within
    PypeIt during data reduction, or as a standalone script outside of
    PypeIt. To initialise the GUI, call the initialise() function in this
    file.
    """

    def __init__(self, canvas, image, frame, sobjs, trace_dict, axes, profdict, slit_ids=None, printout=False):
        """Controls for the interactive Object ID tasks in PypeIt.

        The main goal of this routine is to interactively add/delete/modify
        object apertures.

        Args:
            canvas (Matploltib figure canvas): The canvas on which all axes are contained
            image (AxesImage): The image plotted to screen
            frame (ndarray): The image data
            sobjs (SpecObjs, None): An instance of the SpecObjs class
            trace_dict (dict): A dictionary containing information about the object traces
            axes (dict): Dictionary of four Matplotlib axes instances (Main spectrum panel, two for residuals, one for information)
            profdict (dict): Dictionary containing profile information (profile data, and the left/right lines displayinf the FWHM)
            slit_ids (list, None): List of slit ID numbers
            printout (bool): Should the results be printed to screen
        """
        # Store the axes
        self.image = image
        self.frame = frame
        self.nspec, self.nspat = frame.shape[0], frame.shape[1]
        self.profile = profdict
        self._printout = printout
        self._slit_ids = slit_ids
        self._trcdict = trace_dict
        self.axes = axes
        self.specobjs = sobjs
        self.objtraces = []
        self.anchors = []
        self._obj_idx = -1
        self._spatpos = np.arange(frame.shape[1])[np.newaxis, :].repeat(frame.shape[0], axis=0)  # Spatial coordinate (as the frame shape)
        self.empty_mantrace()

        # Unset some of the matplotlib keymaps
        matplotlib.pyplot.rcParams['keymap.fullscreen'] = ''        # toggling fullscreen (Default: f, ctrl+f)
        #matplotlib.pyplot.rcParams['keymap.home'] = ''              # home or reset mnemonic (Default: h, r, home)
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
        canvas.mpl_connect('motion_notify_event',  self.mouse_move_callback)
        self.canvas = canvas

        # Interaction variables
        self._respreq = [False, None]  # Does the user need to provide a response before any other operation will be permitted? Once the user responds, the second element of this array provides the action to be performed.
        self._qconf = False  # Confirm quit message
        self._changes = False
        self._use_updates = True
        self._trcmthd = 'object'
        self.mmx, self.mmy = 0, 0

        # Draw the spectrum
        self.canvas.draw()

        # Initialise buttons and menu options
        self._ax_meth_default = 'Object'
        self._methdict = dict({'Object': [0, 'object'],
                               'Standard Star': [1, 'std'],
                               'Slit Edges': [2, 'slit']})
#                               'Manual': [3, 'manual']
        self.initialise_menu()

    def print_help(self):
        """Print the keys and descriptions that can be used for Identification
        """
        keys = operations.keys()
        print("===============================================================")
        print("       OBJECT ID OPERATIONS")
        for key in keys:
            print("{0:6s} : {1:s}".format(key, operations[key]))
        print("---------------------------------------------------------------")

    def initialise_menu(self):
        """Initialise the menu buttons
        """
        axcolor = 'lightgoldenrodyellow'
        # Continue with reduction (using updated specobjs)
        ax_cont = plt.axes([0.82, 0.85, 0.15, 0.05])
        self._ax_cont = Button(ax_cont, "Continue (and save changes)", color=axcolor, hovercolor='y')
        self._ax_cont.on_clicked(self.button_cont)
        # Continue with reduction (using original specobjs)
        ax_exit = plt.axes([0.82, 0.79, 0.15, 0.05])
        self._ax_exit = Button(ax_exit, "Continue (don't save changes)", color=axcolor, hovercolor='y')
        self._ax_exit.on_clicked(self.button_exit)
        # Button to select trace method
        rax = plt.axes([0.82, 0.59, 0.15, 0.15], facecolor=axcolor)
        rax.set_title("Select trace method:")
        self._ax_meth = RadioButtons(rax, ('Object', 'Standard Star', 'Slit Edges'))#, 'Manual'))
        self._ax_meth.on_clicked(self.radio_meth)
        # Determine the best default to use:
        if self._trcdict["trace_model"]["object"]["trace_model"] is not None:
            self._ax_meth_default = 'Object'
        elif self._trcdict["trace_model"]["std"]["trace_model"] is not None:
            self._ax_meth_default = 'Standard Star'
        elif self._trcdict["trace_model"]["slit"]["trace_model"] is not None:
            self._ax_meth_default = 'Slit Edges'
#        elif self._trcdict["trace_model"]["manual"]["trace_model"] is not None:
#            self._ax_meth_default = 'Manual'
        # Set the active method
        self._ax_meth.set_active(self._methdict[self._ax_meth_default][0])

    def radio_meth(self, label):
        """Tell the code what to do when a different trace method is selected

        Args:
            label (str): The label of the radio button that was clicked
        """
        # Update the radio button
        if self._methdict[label][1]:
            self._trcmthd = self._methdict[label][1]
        # Check if the method is available, if not, change to the default (manual is always allowed)
        if self._trcmthd != "manual":
            if self._trcdict["trace_model"][self._trcmthd]["trace_model"] is None:
                self.update_infobox(message="That option is not available - changing to default", yesno=False)
                self._ax_meth.set_active(self._methdict[self._ax_meth_default][0])
                self._trcmthd = self._methdict[self._ax_meth_default][1]

    def button_cont(self, event):
        """What to do when the 'exit and save' button is clicked
        """
        self._respreq = [True, "exit_update"]
        self.update_infobox(message="Are you sure you want to exit and use the updated object traces?", yesno=True)

    def button_exit(self, event):
        """What to do when the 'exit and do not save changes' button is clicked
        """
        self._respreq = [True, "exit_restore"]
        self.update_infobox(message="Are you sure you want to exit and use the original object traces?", yesno=True)

    def replot(self):
        """Redraw the entire canvas
        """
        self.canvas.restore_region(self.background)
        self.draw_objtraces()
        self.draw_anchors()
        self.canvas.draw()

    def draw_objtraces(self):
        """Draw the object traces
        """
        for i in self.objtraces: i.pop(0).remove()
        self.objtraces = []
        # Plot the object traces
        for iobj in range(self.specobjs.nobj):
            if iobj == self._obj_idx:
                self.objtraces.append(self.axes['main'].plot(self.specobjs[iobj].trace_spat,
                                                             self.specobjs[iobj].trace_spec,
                                                             'r-', linewidth=4, alpha=0.5))
            else:
                self.objtraces.append(self.axes['main'].plot(self.specobjs[iobj].trace_spat,
                                                             self.specobjs[iobj].trace_spec,
                                                             'r--', linewidth=2, alpha=0.5))

    def draw_anchors(self):
        """Draw the anchors for manual tracing
        """
        for i in self.anchors: i.pop(0).remove()
        self.anchors = []
        # Plot the best fitting trace, if it exists
        if self._mantrace["spat_trc"] is not None:
            self.anchors.append(self.axes['main'].plot(self._mantrace["spat_trc"], self._mantrace["spec_trc"],
                                                       'g-', linewidth=3, alpha=0.5))
        # Plot the anchor points on top
        self.anchors.append(self.axes['main'].plot(self._mantrace["spat_a"], self._mantrace["spec_a"], 'ro', alpha=0.5))

    def draw_profile(self):
        """Draw the object profile
        """
        if self._obj_idx == -1:
            sz = self.profile['profile'].get_xdata.size
            self.profile['profile'].set_ydata(np.zeros(sz))
        else:
            # Plot the extent of the FWHM
            self.profile['fwhm'][0].set_xdata(-self.specobjs[self._obj_idx].fwhm/2.0)
            self.profile['fwhm'][1].set_xdata(+self.specobjs[self._obj_idx].fwhm/2.0)
            # Update the data shown
            objprof = self.make_objprofile()
            self.profile['profile'].set_ydata(objprof)
            self.axes['profile'].set_xlim([-self.specobjs[self._obj_idx].fwhm, +self.specobjs[self._obj_idx].fwhm])
            omin, omax = objprof.min(), objprof.max()
            self.axes['profile'].set_ylim([omin-0.1*(omax-omin), omax+0.1*(omax-omin)])

    def draw_callback(self, event):
        """Draw callback (i.e. everytime the canvas is being drawn/updated)

        Args:
            event (Event): A matplotlib event instance
        """
        # Get the background
        self.background = self.canvas.copy_from_bbox(self.axes['main'].bbox)
        self.draw_objtraces()

    def get_ind_under_point(self, event):
        """Get the index of the object trace closest to the cursor

        Args:
            event (Event): Matplotlib event instance containing information about the event
        """
        mindist = self._spatpos.shape[0]**2
        self._obj_idx = -1
        for iobj in range(self.specobjs.nobj):
            dist = (event.xdata-self.specobjs[iobj].trace_spat)**2 + (event.ydata-self.specobjs[iobj].trace_spec)**2
            if np.min(dist) < mindist:
                mindist = np.min(dist)
                self._obj_idx = iobj
        if self._obj_idx != -1:
            self.draw_profile()
        return

    def get_axisID(self, event):
        """Get the ID of the axis where an event has occurred

        Args:
            event (Event): Matplotlib event instance containing information about the event

        Returns:
            int, None: Axis where the event has occurred
        """
        if event.inaxes == self.axes['main']:
            return 0
        elif event.inaxes == self.axes['info']:
            return 1
        elif event.inaxes == self.axes['profile']:
            return 2
        return None

    def mouse_move_callback(self, event):
        """Store the locations of mouse as it moves across the canvas
        """
        if event.inaxes is None:
            return
        axisID = self.get_axisID(event)
        if event.inaxes == self.axes['main']:
            self.mmx, self.mmy = event.xdata, event.ydata

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
        if event.inaxes == self.axes['main']:
            self._start = [event.x, event.y]
        elif event.inaxes == self.axes['profile']:
            self._start = [event.x, event.y]

    def button_release_callback(self, event):
        """What to do when the mouse button is released

        Args:
            event (Event): Matplotlib event instance containing information about the event
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
        elif event.inaxes == self.axes['profile']:
            if (event.x == self._start[0]) and (event.y == self._start[1]):
                self.set_fwhm(event.xdata)
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
                self._end = [event.x, event.y]
                if (self._end[0] == self._start[0]) and (self._end[1] == self._start[1]):
                    # The mouse button was pressed (not dragged)
                    self.get_ind_under_point(event)
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

    def key_press_callback(self, event):
        """What to do when a key is pressed

        Args:
            event (Event): Matplotlib event instance containing information about the event
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
                if self._respreq[1] == "delete_object" and key == "y":
                    self.delete_object()
                elif self._respreq[1] == "clear_anchors" and key == "y":
                    self.empty_mantrace()
                    self.replot()
                elif self._respreq[1] == "exit_update" and key == "y":
                    self._use_updates = True
                    self.print_pypeit_info()
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
        elif key == 'a':
            self.add_object()
#        elif key == 'c':
#            self._respreq = [True, "clear_anchors"]
#            self.update_infobox(message="Are you sure you want to clear the anchors", yesno=True)
        elif key == 'd':
            if self._obj_idx != -1:
                self._respreq = [True, "delete_object"]
                self.update_infobox(message="Are you sure you want to delete this object trace", yesno=True)
#        elif key == 'm':
#            if self._trcmthd != 'manual':
#                self.update_infobox(message="To add an anchor point, set the 'manual' trace method", yesno=False)
#            else:
#                self.add_anchor()
#        elif key == 'n':
#            self.remove_anchor()
        elif key == 'qu' or key == 'qr':
            if self._changes:
                self.update_infobox(message="WARNING: There are unsaved changes!!\nPress q again to exit", yesno=False)
                self._qconf = True
            else:
                plt.close()
#        elif key == '+':
#            if self._mantrace["polyorder"] < 10:
#                self._mantrace["polyorder"] += 1
#                self.update_infobox(message="Polynomial order = {0:d}".format(self._mantrace["polyorder"]), yesno=False)
#                self.fit_anchors()
#            else:
#                self.update_infobox(message="Polynomial order must be <= 10", yesno=False)
#        elif key == '-':
#            if self._mantrace["polyorder"] > 1:
#                self._mantrace["polyorder"] -= 1
#                self.update_infobox(message="Polynomial order = {0:d}".format(self._mantrace["polyorder"]), yesno=False)
#                self.fit_anchors()
#            else:
#                self.update_infobox(message="Polynomial order must be >= 1", yesno=False)
        self.replot()

    def add_anchor(self):
        """Add a manual anchor point
        """
        self._mantrace['spat_a'].append(self.mmx)
        self._mantrace['spec_a'].append(self.mmy)
        self.fit_anchors()

    def remove_anchor(self):
        """Remove a manual anchor point
        """
        # Find the anchor closest to the mouse position
        if len(self._mantrace['spat_a']) != 0:
            mindist = (self._mantrace['spat_a'][0]-self.mmx)**2 + (self._mantrace['spec_a'][0]-self.mmy)**2
            minidx = 0
            for ianc in range(1, len(self._mantrace['spat_a'])):
                dist = (self._mantrace['spat_a'][ianc] - self.mmx) ** 2 + (self._mantrace['spec_a'][ianc] - self.mmy) ** 2
                if dist < mindist:
                    mindist = dist
                    minidx = ianc
            del self._mantrace['spat_a'][minidx]
            del self._mantrace['spec_a'][minidx]
        self.fit_anchors()

    def fit_anchors(self):
        """Fit the manual anchor points
        """
        if len(self._mantrace['spat_a']) <= self._mantrace['polyorder']:
            self.update_infobox(message="You need to select more trace points before manually adding\n" +
                                        "a manual object trace. To do this, use the 'm' key", yesno=False)
        else:
            # Fit a polynomial to the anchor points
            coeff = np.polyfit(self._mantrace['spec_a'], self._mantrace['spat_a'], self._mantrace['polyorder'])
            self._mantrace['spat_trc'] = np.polyval(coeff, self._mantrace['spec_trc'])
        # Replot, regardless of whether a fit is done (a point might have been added/removed)
        self.replot()

    def add_object(self):
        """Add an object to specobjs
        """
        if self._trcmthd == 'manual' and len(self._mantrace['spat_a']) <= self._mantrace['polyorder']:
            self.update_infobox(message="You need to select more trace points before manually adding\n" +
                                        "a manual object trace. To do this, use the 'm' key", yesno=False)
            return
        # Add an object trace
        spec_vec = self._mantrace['spec_trc']
        if self._trcmthd == 'manual':
            trace_model = self._mantrace['spat_trc'].copy()
            # Now empty the manual tracing
            self.empty_mantrace()
        else:
            trace_model = self._trcdict["trace_model"][self._trcmthd]["trace_model"].copy()
            spat_0 = np.interp(self.mmy, spec_vec, trace_model)
            shift = self.mmx - spat_0
            trace_model += shift
        xsize = self._trcdict["edges_r"] - self._trcdict["edges_l"]
        nsamp = np.ceil(xsize.max())
        # Extract the SpecObj parameters
        par = self._trcdict['sobj_par']
        # Create a SpecObj
        thisobj = specobjs.SpecObj(par['frameshape'], par['slit_spat_pos'], par['slit_spec_pos'],
                                   det=par['det'], setup=par['setup'], slitid=par['slitid'],
                                   orderindx=par['orderindx'], objtype=par['objtype'])
        thisobj.hand_extract_spat = self.mmx
        thisobj.hand_extract_spec = self.mmy
        thisobj.hand_extract_det = par['det']
        thisobj.hand_extract_fwhm = None
        thisobj.hand_extract_flag = True
        f_ximg = RectBivariateSpline(spec_vec, np.arange(self.nspat), par["ximg"])
        thisobj.spat_fracpos = f_ximg(thisobj.hand_extract_spec, thisobj.hand_extract_spat,
                                      grid=False)  # interpolate from ximg
        thisobj.smash_peakflux = np.interp(thisobj.spat_fracpos * nsamp, np.arange(nsamp),
                                           self._trcdict['profile'])  # interpolate from fluxconv
        # assign the trace
        thisobj.trace_spat = trace_model
        thisobj.trace_spec = spec_vec
        thisobj.spat_pixpos = thisobj.trace_spat[self.nspec//2]
        thisobj.set_idx()
        if len(self.specobjs) != 0:
            thisobj.fwhm = self.specobjs[0].fwhm
        else:  # Otherwise just use the fwhm parameter input to the code (or the default value)
            thisobj.fwhm = 2
        # Finally, add new object
        self.specobjs.add_sobj(thisobj)

    def delete_object(self):
        """Delete a specobj
        """
        self.specobjs.remove_sobj(self._obj_idx)
        self._obj_idx = -1
        self.replot()

    def print_pypeit_info(self):
        """print text that the user should insert into their .pypeit file
        """
        print("TO BE IMPLEMENTED!")

    def make_objprofile(self):
        """Generate an object profile from the traces
        """
        coords = self._spatpos - self.specobjs[self._obj_idx].trace_spat[:, np.newaxis]
        ww = np.where(np.abs(coords) < 4*self.specobjs[self._obj_idx].fwhm)
        bincent = self.profile['profile'].get_xdata()
        offs = 0.5*(bincent[1]-bincent[0])
        edges = np.append(bincent[0]-offs, bincent+offs)
        prof, _ = np.histogram(coords[ww], bins=edges, weights=self.frame[ww])
        return prof/ww[0].size

    def set_fwhm(self, xdata):
        """Set the FWHM using the available panel

        Args:
            xdata (float): The x coordinate selected by the user
        """
        self.specobjs[self._obj_idx].fwhm = 2.0*np.abs(xdata)
        self.draw_profile()
        self.replot()
        return

    def get_specobjs(self):
        """Get the updated version of SpecObjs

        Returns:
            SpecObjs: SpecObjs Class
        """
        if self._use_updates:
            return self.specobjs
        else:
            return None

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

    def empty_mantrace(self):
        """Generate an empty dictionary for the manual tracing
        """
        self._mantrace = dict(spat_a=[], spec_a=[], spat_trc=None, spec_trc=np.arange(self.nspec), polyorder=0)
        return


def initialise(frame, trace_dict, sobjs, slit_ids=None, printout=False):
    """Initialise the 'ObjFindGUI' window for interactive object tracing

        Args:
            frame (ndarray): Sky subtracted science image
            trace_dict (dict, None): Dictionary containing slit and object trace information
            sobjs (SpecObjs, None): SpecObjs Class
            slit_ids (list, None): List of slit ID numbers
            printout (bool): Should the results be printed to screen

        Returns:
            ObjFindGUI: Returns an instance of the ObjFindGUI class
    """
    if sobjs is None:
        print("Running the interactive object tracing tool outside of a reduction")
        print("NOT YET IMPLEMENTED")
        sys.exit()

    # This allows the input lord and rord to either be (nspec, nslit) arrays or a single
    # vectors of size (nspec)
    if trace_dict['slit_left'].ndim == 2:
        nslit = trace_dict['slit_left'].shape[1]
        lordloc = trace_dict['slit_left']
        rordloc = trace_dict['slit_righ']
    else:
        nslit = 1
        lordloc = trace_dict['slit_left'].reshape((trace_dict['slit_left'].size, 1))
        rordloc = trace_dict['slit_righ'].reshape((trace_dict['slit_righ'].size, 1))

    # Assign the slit IDs if none were provided
    if slit_ids is None:
        slit_ids = [str(slit) for slit in np.arange(nslit)]

    # Determine the scale of the image
    med = np.median(frame)
    mad = np.median(np.abs(frame-med))
    vmin = med-3*mad
    vmax = med+3*mad

    # Add the main figure axis
    fig, ax = plt.subplots(figsize=(16, 9), facecolor="white")
    plt.subplots_adjust(bottom=0.05, top=0.85, left=0.05, right=0.8)
    image = ax.imshow(frame, aspect='auto', cmap = 'Greys', vmin=vmin, vmax=vmax)

    # Overplot the slit traces
    specarr = np.arange(lordloc.shape[0])
    for sl in range(nslit):
        ax.plot(lordloc[:, sl], specarr, 'g-')
        ax.plot(rordloc[:, sl], specarr, 'b-')

    # Add an object profile axis
    axprof = fig.add_axes([0.82, 0.05, .15, 0.30])
    profx = np.arange(-10, 10.1, 0.1)
    profy = np.zeros(profx.size)
    profile = axprof.plot(profx, profy, 'k-')
    vlinel = axprof.axvline(-1, color='r')
    vliner = axprof.axvline(+1, color='r')
    axprof.set_title("Object profile")
    axprof.set_xlim((-3, 3))
    axprof.set_ylim((0, 1))
    axprof.set_yticks([])

    # Add an information GUI axis
    axinfo = fig.add_axes([0.15, .92, .7, 0.07])
    axinfo.get_xaxis().set_visible(False)
    axinfo.get_yaxis().set_visible(False)
    axinfo.text(0.5, 0.5, "Press '?' to list the available options", transform=axinfo.transAxes,
                horizontalalignment='center', verticalalignment='center')
    axinfo.set_xlim((0, 1))
    axinfo.set_ylim((0, 1))

    axes = dict(main=ax, profile=axprof, info=axinfo)
    profdict = dict(profile=profile[0], fwhm=[vlinel, vliner])
    # Initialise the object finding window and display to screen
    fig.canvas.set_window_title('PypeIt - Object Tracing')
    ofgui = ObjFindGUI(fig.canvas, image, frame, sobjs, trace_dict, axes, profdict,
                       slit_ids=slit_ids, printout=printout)
    plt.show()

    return ofgui
