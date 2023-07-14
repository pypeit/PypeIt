"""
This script allows the user to add/delete/modify object traces

.. todo::

    Implement color scaling with RMB click+drag

"""

import os
import sys
import copy
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, RadioButtons
from scipy.interpolate import RectBivariateSpline

# TODO: Commented out the import of specobjs because it was only used by
# the (presumably) defunct method.
#from pypeit import specobjs
from pypeit import msgs

# TODO: No globals please
operations = dict({'cursor': "Select object trace (LMB click)\n" +
                   "         Navigate (LMB drag = pan, RMB drag = zoom)\n" +
                   "         Note: In order to pan/zoom, you need to first activate\n" +
                   "         the pan/zoom tool with the 'p' key, or by selecting the\n" +
                   "         pan/zoom tool on the Matplotlib navigation tool menu. You\n" +
                   "         can also zoom using the magnifying glass (select this option\n" +
                   "         from the Matplotlib navigation tool menu)\n",
                   'a': "Add a new object trace using the selected method",
                   'c': "Center the window at the location of the mouse",
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


class ObjectTraces:
    """
    Simple class to store the object traces
    """

    def __init__(self):
        self._det = []
        self._add_rm = []  # Flag to say if an object was been added (1), removed (-1), or was there originally (0)
        self._pos_spat = []
        self._pos_spec = []
        self._trace_spat = []
        self._trace_spec = []
        self._fwhm = []

    @property
    def nobj(self):
        return len(self._det)

    def from_specobj(self, sobjs, det):
        """Fill the object traces from a SpecObjs class

        Args:
            sobjs (SpecObjs): An instance of the SpecObjs class
            det (int): Detector to which this applies
        """
        nobj = sobjs.nobj
        for ii in range(nobj):
            pos_spat = sobjs[ii].trace_spat[sobjs[ii].trace_spat.size//2]
            pos_spec = sobjs[ii].trace_spec[sobjs[ii].trace_spec.size//2]
            self.add_object(det, pos_spat, pos_spec, sobjs[ii].trace_spat, sobjs[ii].trace_spec, sobjs[ii].fwhm, addrm=0)

    def from_dict(self, obj_dict, det):
        """Fill the object traces from a SpecObjs class

        Args:
            obj_dict (dict): A dictionary containing the spatial object traces and FWHM
            det (int): Detector to which this applies
        """
        nobj = len(obj_dict['traces'])
        for ii in range(nobj):
            spec_trace = np.arange(obj_dict['traces'][ii].size)
            pos_spat = obj_dict['traces'][ii][obj_dict['traces'][ii].size//2]
            pos_spec = spec_trace[spec_trace.size//2]
            self.add_object(det, pos_spat, pos_spec, obj_dict['traces'][ii], spec_trace, obj_dict['fwhm'][ii], addrm=0)

    def add_object(self, det, pos_spat, pos_spec, trc_spat, trc_spec, fwhm, addrm=1):
        """Add an object trace

        Args:
            det (int): Detector to add a slit on
            pos_spat (float): Spatial pixel position
            pos_spec (float): Spectral pixel position
            trc_spat (ndarray): Spatial trace of object
            trc_spec (ndarray): Spectral trace of object
            fwhm (float): FWHM of the object
            addrm (int): Flag to say if an object was been added (1), removed (-1), or was an auto found slit (0)
        """
        self._det.append(det)
        self._add_rm.append(addrm)
        self._pos_spat.append(pos_spat)
        self._pos_spec.append(pos_spec)
        self._trace_spat.append(trc_spat)
        self._trace_spec.append(trc_spec)
        self._fwhm.append(fwhm)

    def delete_object(self, ind):
        """Delete an object trace

        Args:
            ind (int): Index of object trace to remove
        """
        if self._add_rm[ind] == 1:
            # If an object was added by hand, the user can delete it from the list
            del self._det[ind]
            del self._add_rm[ind]
            del self._pos_spat[ind]
            del self._pos_spec[ind]
            del self._trace_spat[ind]
            del self._trace_spec[ind]
            del self._fwhm[ind]
        else:
            # otherwise, just flag that it needs to be removed
            self._add_rm[ind] = -1

    def get_pypeit_string(self):
        """Construct a string that can be placed in the .pypeit file to define new object traces
        """
        strall = ""
        for ii in range(self.nobj):
            if self._add_rm[ii] == 1:
                strall += ",{0:d}:{1:.1f}:{2:.1f}:{3:.1f}".format(self._det[ii], self._pos_spat[ii], self._pos_spec[ii], self._fwhm[ii])
        strall += "\n"
        return strall[1:]


class ObjFindGUI:
    """
    GUI to interactively identify object traces. The GUI can be run within
    PypeIt during data reduction, or as a standalone script outside of
    PypeIt. To initialise the GUI, call the initialise() function in this
    file.
    """
    def __init__(self, canvas, image, frame, det, sobjs, left, right, obj_trace, trace_models,
                 axes, profdict, slit_ids=None, printout=False, runtime=False):
        """Controls for the interactive Object ID tasks in PypeIt.

        The main goal of this routine is to interactively add/delete/modify
        object apertures.

        Args:
            canvas (Matploltib figure canvas):
                The canvas on which all axes are contained
            image (AxesImage):
                The image plotted to screen
            frame (ndarray):
                The image data
            det (int):
                Detector to add a slit on
            sobjs (SpecObjs, None):
                An instance of the SpecObjs class
            left (numpy.ndarray):
                Slit left edges
            right (numpy.ndarray):
                Slit right edges
            obj_trace (dict):
                Result of
                :func:`pypeit.scripts.object_finding.parse_traces`.
            trace_models (dict):
                Dictionary with the object, standard star, and slit
                trace models
            axes (dict):
                Dictionary of four Matplotlib axes instances (Main
                spectrum panel, two for residuals, one for information)
            profdict (dict):
                Dictionary containing profile information (profile data,
                and the left/right lines displayinf the FWHM)
            slit_ids (list, None):
                List of slit ID numbers
            printout (bool):
                Should the results be printed to screen
            runtime (bool):
                Is the GUI being launched during data reduction?
        """
        # Store the axes
        self._det = det
        self.image = image
        self.frame = frame
        self.nspec, self.nspat = frame.shape[0], frame.shape[1]
        self._spectrace = np.arange(self.nspec)
        self.profile = profdict
        self._printout = printout
        self._runtime = runtime
        self._slit_ids = slit_ids

        self.left = left
        self.right = right
        self.obj_trace = obj_trace
        self.trace_models = trace_models

        self.axes = axes
        self.specobjs = sobjs
        self.objtraces = []
        self._object_traces = ObjectTraces()
        self.anchors = []
        self._obj_idx = -1
        self._spatpos = np.arange(frame.shape[1])[np.newaxis, :].repeat(frame.shape[0], axis=0)  # Spatial coordinate (as the frame shape)
        self.empty_mantrace()
        if sobjs is None:
            self._object_traces.from_dict(self.obj_trace, det)
        else:
            self._object_traces.from_specobj(sobjs, det)

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
        self._respreq = [False, None]  # Does the user need to provide a response before any other operation will be permitted? Once the user responds, the second element of this array provides the action to be performed.
        self._qconf = False  # Confirm quit message
        self._changes = False
        self._use_updates = True
        self._trcmthd = 'object'
        self.mmx, self.mmy = 0, 0
        self._inslit = 0  # Which slit is the mouse in

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
        print("Add/remove object traces until you are happy with the resulting")
        print("traces. When you've finished, click one of the exit buttons on")
        print("the right side of the page. If you click 'Continue (and save changes)'")
        print("the object traces will be printed to the terminal, where you can")
        print("copy them into your .pypeit file.")
        print("")
        print("thick coloured dashed lines = object traces")
        print("thick coloured solid line   = currently selected object trace")
        print("thin green/blue lines       = slit edges")
        print("")
        print("Meanings of the different coloured dashed lines:")
        print(" green = user-defined object trace")
        print(" blue  = trace automatically generated with PypeIt")
        print(" red   = trace automatically generated with PypeIt (deleted)")
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
        if self.trace_models['object']['trace_model'] is not None:
            self._ax_meth_default = 'Object'
        elif self.trace_models['std']['trace_model'] is not None:
            self._ax_meth_default = 'Standard Star'
        elif self.trace_models['slit']['trace_model'] is not None:
            self._ax_meth_default = 'Slit Edges'
#        elif self._trcdict["trace_model"]["manual"]["trace_model"] is not None:
#            self._ax_meth_default = 'Manual'
        # Set the active method
        self._ax_meth.set_active(self._methdict[self._ax_meth_default][0])

    def radio_meth(self, label, infobox=True):
        """Tell the code what to do when a different trace method is selected

        Args:
            label (str): The label of the radio button that was clicked
        """
        # Update the radio button
        if self._methdict[label][1]:
            self._trcmthd = self._methdict[label][1]
        # Check if the method is available, if not, change to the default (manual is always allowed)
        if self._trcmthd != "manual":
            if self.trace_models[self._trcmthd]['trace_model'] is None:
                self.update_infobox(message="That option is not available - changing to default",
                                    yesno=False)
                self._ax_meth.set_active(self._methdict[self._ax_meth_default][0])
                self._trcmthd = self._methdict[self._ax_meth_default][1]
            else:
                if infobox:
                    self.update_infobox(message="Trace method set to: {0:s}".format(label), yesno=False)
        else:
            if infobox:
                self.update_infobox(message="Trace method set to: {0:s}".format(label), yesno=False)

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
        allcols = ['DodgerBlue', 'LimeGreen', 'r']  # colors mean: [original, added, deleted]
        for iobj in range(self._object_traces.nobj):
            color = allcols[self._object_traces._add_rm[iobj]]
            if iobj == self._obj_idx:
                self.objtraces.append(self.axes['main'].plot(self._object_traces._trace_spat[iobj],
                                                             self._object_traces._trace_spec[iobj],
                                                             color=color,
                                                             linestyle='-', linewidth=4, alpha=0.5))
            else:
                self.objtraces.append(self.axes['main'].plot(self._object_traces._trace_spat[iobj],
                                                             self._object_traces._trace_spec[iobj],
                                                             color=color,
                                                             linestyle='--', linewidth=3, alpha=0.5))

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
            self.profile['fwhm'][0].set_xdata(-self._object_traces._fwhm[self._obj_idx]/2.0)
            self.profile['fwhm'][1].set_xdata(+self._object_traces._fwhm[self._obj_idx]/2.0)
            # Update the data shown
            objprof = self.make_objprofile()
            self.profile['profile'].set_ydata(objprof)
            self.axes['profile'].set_xlim([-self._object_traces._fwhm[self._obj_idx],
                                           +self._object_traces._fwhm[self._obj_idx]])
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
        for iobj in range(self._object_traces.nobj):
            dist = (event.xdata-self._object_traces._trace_spat[iobj])**2 +\
                   (event.ydata-self._object_traces._trace_spec[iobj])**2
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
                self.update_infobox(message="FWHM updated for the selected object", yesno=False)
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
                    self.update_infobox(message="Object selected", yesno=False)
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
        elif key == 'c':
            if axisID == 0:
                # If this is pressed on the main window
                self.recenter()
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

    def get_slit(self):
        """Find the slit that the mouse is currently in
        """
        ypos = int(self.mmy)
        for sl in range(self.left.shape[1]):
            if (self.mmx > self.left[ypos, sl]) and (self.mmx < self.right[ypos, sl]):
                self._inslit = sl
                return

    def add_object(self):
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
            if self._trcmthd == 'slit':
                self.get_slit()
                trace_model = self.trace_models[self._trcmthd]['trace_model'][:,self._inslit].copy()
            else:
                trace_model = self.trace_models[self._trcmthd]['trace_model'].copy()
            spat_0 = np.interp(self.mmy, spec_vec, trace_model)
            shift = self.mmx - spat_0
            trace_model += shift
        # Determine the FWHM
        if self._object_traces.nobj != 0:
            fwhm = self._object_traces._fwhm[0]
        else:  # Otherwise just use the fwhm parameter input to the code (or the default value)
            fwhm = 2

        # Finally, add this object to the list
        self._object_traces.add_object(self._det, self.mmx, self.mmy, trace_model, self._spectrace.copy(), fwhm)

# TODO: This method must have been defunct because it uses elements of
# _trcdict that don't exist (as far as I can tell) and instantiates a
# SpecObj from the specobjs module, which hasn't existed for while
#
#    def add_object_sobj(self):
#        """Add an object to specobjs
#        """
#        if self._trcmthd == 'manual' and len(self._mantrace['spat_a']) <= self._mantrace['polyorder']:
#            self.update_infobox(message="You need to select more trace points before manually adding\n" +
#                                        "a manual object trace. To do this, use the 'm' key", yesno=False)
#            return
#        # Add an object trace
#        spec_vec = self._mantrace['spec_trc']
#        if self._trcmthd == 'manual':
#            trace_model = self._mantrace['spat_trc'].copy()
#            # Now empty the manual tracing
#            self.empty_mantrace()
#        else:
#            trace_model = self._trcdict["trace_model"][self._trcmthd]["trace_model"].copy()
#            spat_0 = np.interp(self.mmy, spec_vec, trace_model)
#            shift = self.mmx - spat_0
#            trace_model += shift
#        xsize = self._trcdict["slit_righ"] - self._trcdict["slit_left"]
#        nsamp = np.ceil(xsize.max())
#        # Extract the SpecObj parameters
#        par = self._trcdict['sobj_par']
#        # Create a SpecObj
#        thisobj = specobjs.SpecObj(par['frameshape'], par['slit_spat_pos'], par['slit_spec_pos'],
#                                   det=par['det'], setup=par['setup'], slitid=par['slitid'],
#                                   orderindx=par['orderindx'], objtype=par['objtype'])
#        thisobj.hand_extract_spat = self.mmx
#        thisobj.hand_extract_spec = self.mmy
#        thisobj.hand_extract_det = par['det']
#        thisobj.hand_extract_fwhm = None
#        thisobj.hand_extract_flag = True
#        f_ximg = RectBivariateSpline(spec_vec, np.arange(self.nspat), par["ximg"])
#        thisobj.spat_fracpos = f_ximg(thisobj.hand_extract_spec, thisobj.hand_extract_spat,
#                                      grid=False)  # interpolate from ximg
#        thisobj.smash_peakflux = np.interp(thisobj.spat_fracpos * nsamp, np.arange(nsamp),
#                                           self._trcdict['profile'])  # interpolate from fluxconv
#        # assign the trace
#        thisobj.trace_spat = trace_model
#        thisobj.trace_spec = spec_vec
#        thisobj.spat_pixpos = thisobj.trace_spat[self.nspec//2]
#        thisobj.set_idx()
#        if self._object_traces.nobj != 0:
#            thisobj.fwhm = self._object_traces._fwhm[0]
#        else:  # Otherwise just use the fwhm parameter input to the code (or the default value)
#            thisobj.fwhm = 2
#        # Finally, add new object
#        self.specobjs.add_sobj(thisobj)

    def delete_object(self):
        """Delete an object trace
        """
        self._object_traces.delete_object(self._obj_idx)
        self._obj_idx = -1
        self.replot()

    def delete_object_sobj(self):
        """Delete a specobj
        """
        self.specobjs.remove_sobj(self._obj_idx)
        self._obj_idx = -1

    def print_pypeit_info(self):
        """print text that the user should insert into their .pypeit file
        """
        if 1 in self._object_traces._add_rm:
            msgs.info("Include the following info in the manual_extract column in your .pypeit file:\n")
            print(self._object_traces.get_pypeit_string())

    def recenter(self):
        xlim = self.axes['main'].get_xlim()
        ylim = self.axes['main'].get_ylim()
        xmin = self.mmx - 0.5*(xlim[1]-xlim[0])
        xmax = self.mmx + 0.5*(xlim[1]-xlim[0])
        ymin = self.mmy - 0.5*(ylim[1]-ylim[0])
        ymax = self.mmy + 0.5*(ylim[1]-ylim[0])
        self.axes['main'].set_xlim([xmin, xmax])
        self.axes['main'].set_ylim([ymin, ymax])

    def make_objprofile(self):
        """Generate an object profile from the traces
        """
        coords = self._spatpos - self._object_traces._trace_spat[self._obj_idx][:, np.newaxis]
        ww = np.where(np.abs(coords) < 4*self._object_traces._fwhm[self._obj_idx])
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
        self._object_traces._fwhm[self._obj_idx] = 2.0*np.abs(xdata)
        self.draw_profile()
        self.replot()
        return

    def get_specobjs(self):
        """Get the updated version of SpecObjs

        Returns:
            SpecObjs: SpecObjs Class
        """
        if self._use_updates:
            msgs.work("Have not updated SpecObjs yet")
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


def initialise(det, frame, left, right, obj_trace, trace_models, sobjs, slit_ids=None,
               runtime=False, printout=False):
    """
    Initialise the 'ObjFindGUI' window for interactive object tracing

    Args:
        det (int):
            1-indexed detector number
        frame (numpy.ndarray):
            Sky subtracted science image
        left (numpy.ndarray):
            Slit left edges
        right (numpy.ndarray):
            Slit right edges
        obj_trace (dict):
            Result of
            :func:`pypeit.scripts.object_finding.parse_traces`.
        trace_models (dict):
            Dictionary with the object, standard star, and slit trace
            models
        sobjs (SpecObjs, None):
            SpecObjs Class
        det (int):
            Detector index
        slit_ids (list, None):
            List of slit ID numbers
        runtime (bool):
            Is this GUI being launched during a data reduction?
        printout (bool):
            Should the results be printed to screen

    Returns:
        :class:`ObjFindGUI`: Returns an instance of the ObjFindGUI class

    """

    # This allows the input lord and rord to either be (nspec, nslit) arrays or a single
    # vectors of size (nspec)
    if left.ndim == 2:
        nslit = left.shape[1]
        _left = left.copy()
        _right = right.copy()
    else:
        nslit = 1
        _left = left.reshape(-1,1)
        _right = right.reshape(-1,1)

    # TODO: This edits the input array!
    if trace_models['slit']['trace_model'].ndim == 1:
        trace_models['slit']['trace_model'] = trace_models['slit']['trace_model'].reshape(-1,1)

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
    specarr = np.arange(_left.shape[0])
    for sl in range(nslit):
        ax.plot(_left[:, sl], specarr, 'g-')
        ax.plot(_right[:, sl], specarr, 'b-')

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
    fig.canvas.manager.set_window_title('PypeIt - Object Tracing')
    ofgui = ObjFindGUI(fig.canvas, image, frame, det, sobjs, _left, _right, obj_trace,
                       trace_models, axes, profdict, slit_ids=slit_ids, printout=printout,
                       runtime=runtime)
    plt.show()

    return ofgui
