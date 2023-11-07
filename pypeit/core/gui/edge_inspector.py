"""
Implements a matplotlib GUI to inspect and interact with the slit edge tracing.

.. include:: ../include/links.rst
"""

from IPython import embed

import numpy as np
from matplotlib import pyplot, widgets

from pypeit import msgs
from pypeit.core.gui import gui_util


class EdgeInspectorGUI:
    """
    matplotlib-based GUI to analyze/manipulate edge traces.

    Args:
        edges (:class:`~pypeit.edgetrace.EdgeTraceSet`):
            Edge tracing to inspect/edit.  Note this object is edited directly!
    """
    def __init__(self, edges):

        self.edges = edges

        # Remove all matplotlib plotting window key bindings to make sure they
        # don't overlap with any defined by the function.
        gui_util.clean_pyplot_keymap()

        # Layout 
        sx = 0.07
        sy = 0.15
        dx = 0.78
        dy = 0.78

        # Figure
        w,h = pyplot.figaspect(1)
        self.fig = pyplot.figure(figsize=(2*w,2*h))

        # Axes
        image_ax = self.fig.add_axes([sx, sy, dx, dy])
        image_ax.tick_params(which='both', top=True, right=True)
        image_ax.minorticks_on()
        image_cax = self.fig.add_axes([sx + dx + 0.02, sy, 0.01, dy])
        image_slider_ax = self.fig.add_axes([sx, sy-0.08, dx/1.2, 0.02])
        image_change_button_ax = self.fig.add_axes([sx + dx - 0.1, sy - 0.1, 0.1, 0.05])
        save_button_ax = self.fig.add_axes([sx + dx + 0.01, sy - 0.1, 0.1, 0.05])

        # Extent
        ny, nx = edges.traceimg.image.shape
        extent = [-0.5, nx-0.5, -0.5, ny-0.5]

        # Images that can be toggled
        self.images = [edges.traceimg.image, edges.sobelsig]

        # Instantiate by plotting the trace image
        lim = [np.amin(self.images[0]), np.amax(self.images[0])]
        self.image_plot = image_ax.imshow(self.images[0], origin='lower', interpolation='nearest',
                                          vmin=lim[0], vmax=lim[1], cmap='gray', extent=extent,
                                          aspect='auto')
        # Get the reference row, used to set where to predict new traces.
        self.has_pca = edges.pcatype is not None
        if self.has_pca:
            self.reference_row = self.edges.left_pca.reference_row \
                                    if self.edges.par['left_right_pca'] \
                                    else self.edges.pca.reference_row
        else:
            msgs.warn('Edges object does not include a PCA decomposition of the traces.')
            self.reference_row = self.edges.nspec // 2
        # NOTE: line properties match what is used for the Pointer
        self.ref_row_line = image_ax.axhline(self.reference_row, color='C1', lw=0.5)

        # Add the colorbar
        self.image_cb = self.fig.colorbar(self.image_plot, cax=image_cax)
        # And the colorbar slider
        self.image_slider = gui_util.UpdateableRangeSlider(image_slider_ax, 'Data Range',
                                                           lim[0], lim[1], valinit=lim)
        # And the ability to update the image
        self.image_updater = gui_util.UpdateableImage(self.images, self.image_plot,
                                                      self.image_slider)
        # And the button that will toggle the images
        self.image_change_button = widgets.Button(image_change_button_ax, 'Image', color='0.7',
                                                  hovercolor='0.9')
        self.image_change_button.on_clicked(self.image_updater.next_image)
        # And the button that will update the Edges object (rm, add, save
        # traces) and update the line plot
        self.save_button = widgets.Button(save_button_ax, 'Save', color='0.7', hovercolor='0.9')
        self.save_button.on_clicked(self.update_traces)

        # Initialize and plot the traces
        self._reset_traces()
        self.plot_traces()

        # Create the pointer
        self.img_pointer = gui_util.Pointer(self.image_plot.axes, name='image', horizOn=False,
                                            color='C1', lw=0.5)
        # Add its functions
        self.img_pointer.register('m', self.move, descr='Move the nearest trace to the cursor')
        self.img_pointer.register('d', self.delete, descr='Delete the nearest trace')
        self.img_pointer.register('l', self.add_left, descr='Add a left trace')
        self.img_pointer.register('r', self.add_right, descr='Add a right trace')
        self.img_pointer.register('U', self.undo, descr='Undo all changes since last save')
        # Add the help dialog
        self.img_pointer.build_help()

    def close(self):
        """
        Close the GUI.  Restores matplotlib RC defaults and closes the figure.
        """
        pyplot.rcdefaults()
        self.fig.clear()
        pyplot.close(self.fig)

    def _trace_color(self, side):
        """
        Return the color to use for each edge sided.
        
        Args:
            side (:obj:`int`):
                Edge side.  Either -1 for left or 1 for right.

        Returns:
            :obj:`str`: The color to use for the trace.
        """
        return 'C2' if side < 0 else 'C4'

    def plot_traces(self):
        """
        Plot the traces from scratch.
        """
        self.trace_plot = []
        for i in range(self.ntrace):
            if not self.good[i]:
                self.trace_plot += [None]
                continue
            color = self._trace_color(self.side[i])
            # Trace data should *not* be masked so that calling plot only
            # produces a single Line2D object
            self.trace_plot += [self.image_plot.axes.plot(self.trace_cen[:,i], self.spec_pix,
                                                          color=color, lw=2)[0]]

    def update_traces(self, *args):
        """
        Update the underlying trace data and edges object.

        Changes to the traces are kept by the internals until we're ready to
        "update" the edges object.  This performs the update and replots the
        trace data.

        All arguments to this function are accepted but ignored.  The reason is
        because this is the function passed to
        `matplotlib.widgets.Button.on_clicked`_, but this function does not need
        any of the input from the ``Button`` event.
        """
        # Offsets are "applied" by first removing the initial traces and then
        # adding new ones with the locations offset.
        _offset = (np.absolute(self.offset) > 0)
        # Remove traces
        _remove = self.remove | _offset
        if np.any(_remove):
            self.edges.remove_traces(_remove[:self.edges.ntrace], rebuild_pca=self.has_pca)
        # Add traces
        _add = self.add | _offset
        if np.any(_add):
            new_traces = self.trace_cen[:,_add] + self.offset[None,_add]
            self.edges.insert_traces(self.side[_add], new_traces, nudge=False)
        # Re-sync traces (make this optional?)
        if np.any(_remove) or np.any(_add):
            success = self.edges.sync()
            if not success:
                msgs.warn('Unable to synchronize left-right traces!')

        # Remove the trace lines from the plot
        # TODO: There may be an easier way to do this, but I couldn't find it.
        while len(self.trace_plot) > 0:
            for i, l in enumerate(self.image_plot.axes.lines):
                if l is self.trace_plot[0]:
                    break
            self.image_plot.axes.lines[i].remove()
            del self.trace_plot[0]

        # Reset the traces using the updated edges object
        self._reset_traces()
        # Plot the traces
        self.plot_traces()
        # Update the plot
        pyplot.draw()
            
    def _reset_traces(self):
        """
        Re-initialize the internals that perform the trace bookkeeping.
        """
        self.spec_pix = np.arange(self.edges.nspec)
        self.ntrace = self.edges.ntrace
        self.good = self.edges.good_traces(include_box=True)
        self.side = np.clip(self.edges.traceid, -1, 1) 
        self.add = np.zeros(self.ntrace, dtype=bool)
        self.remove = np.zeros(self.ntrace, dtype=bool)
        self.offset = np.zeros(self.ntrace, dtype=float)
        self.trace_cen = self.edges.edge_cen if self.edges.edge_fit is None \
                            else self.edges.edge_fit
        
    def undo(self, *args):
        """
        Undo all operations since the last time the ``edges`` object was
        updated.

        The function call includes ``*args`` because it is used as an event
        call-back function that must support arguments passed by the event.  But
        all of these arguments are ignored.
        """
        self._reset_traces()
        self.update_traces()
        
    def nearest_trace(self, spat):
        """
        Find the trace nearest to the spatial position of the cursor.

        Args:
            spat (:obj:`float`):
                The spatial (x) position of the cursor
        """
        # Only select from the traces that are good and not being removed
        indx = self.good & np.logical_not(self.remove)
        diff = np.absolute(self.trace_cen[self.reference_row][indx] + self.offset[indx] - spat)
        return np.arange(self.ntrace)[indx][np.argmin(diff)]
        
    def move(self, pos):
        """
        Move the trace nearest the cursor spatial (x) position to the cursor's
        location.

        This is an event call-back function that must accept an array-like
        object giving the pointer coordinates.

        Args:
            pos (array-like):
                List with x and y position of the cursor at the time of the
                window event.
        """
        # Find the nearest trace
        i = self.nearest_trace(pos[0])
        # Get the spatial offset
        self.offset[i] = pos[0] - self.trace_cen[self.reference_row,i]
        # Report
        msgs.info(f'Offsetting trace {i} by {self.offset[i]} pixels')
        # Reset the line data
        self.trace_plot[i].set_data((self.trace_cen[:,i] + self.offset[i], self.spec_pix))
        # Update the plot
        pyplot.draw()

    def delete(self, pos):
        """
        Delete the trace nearest the cursor spatial (x) position.

        This is an event call-back function that must accept an array-like
        object giving the pointer coordinates.

        Args:
            pos (array-like):
                List with x and y position of the cursor at the time of the
                window event.
        """
        # Find the nearest trace
        i = self.nearest_trace(pos[0])
        # Set to remove it
        self.remove[i] = True
        # For now, simply set the line to be invisible
        self.trace_plot[i].set_visible(False)
        # Update the plot
        pyplot.draw()

    def add_trace(self, spat, side):
        """
        Add a new trace nearest the cursor spatial (x) position.

        Args:
            spat (:obj:`float`):
                The spatial (x) position of the cursor at the time of the window
                event.
            side (:obj:`int`):
                Edge side.  Either -1 for left or 1 for right.
        """
        # Append to the book-keeping vectors
        self.ntrace += 1
        self.good = np.append(self.good, [True])
        self.side = np.append(self.side, [side])
        self.add = np.append(self.add, [True])
        self.remove = np.append(self.remove, [False])
        self.offset = np.append(self.offset, [0.])

        # Predict the spatial location of the trace
        if self.has_pca:
            new_trace = self.edges.predict_traces(spat, side=side)
        else:
            i = self.nearest_trace(spat)
            offset = spat - self.trace_cen[self.reference_row,i]
            new_trace = self.trace_cen[:,i] + offset
        # Add it
        self.trace_cen = np.hstack((self.trace_cen, new_trace[:,None]))
        # Plot it
        color = self._trace_color(side)
        self.trace_plot += [self.image_plot.axes.plot(new_trace, self.spec_pix,
                                                      color=color, lw=2)[0]]
        # Draw the plot
        pyplot.draw()
        # Report
        msgs.info(f'Added {"right" if side > 0 else "left"} trace passing through '
                  f'({spat:.1f}, {self.reference_row:.2f}).')

    def add_left(self, pos):
        """
        Add a new left trace nearest the cursor spatial (x) position.

        This is an event call-back function that must accept an array-like
        object giving the pointer coordinates.

        Args:
            pos (array-like):
                List with x and y position of the cursor at the time of the
                window event.
        """
        self.add_trace(pos[0], -1)

    def add_right(self, pos):
        """
        Add a new right trace nearest the cursor spatial (x) position.

        This is an event call-back function that must accept an array-like
        object giving the pointer coordinates.

        Args:
            pos (array-like):
                List with x and y position of the cursor at the time of the
                window event.
        """
        self.add_trace(pos[0], 1)



