"""
Implements a matplotlib GUI to inspect and interact with the slit edge tracing.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pathlib import Path

from IPython import embed

import numpy as np
from matplotlib import pyplot, rc, ticker, widgets
from matplotlib.backend_bases import MouseButton


#-----------------------------------------------------------------------
# Pointer class
class Pointer(widgets.AxesWidget):
    """
    A pointer widget

    Args:
        ax (matplotlib.image.AxesImage):
            Returned object after running ``matplotlib.pyplot.imshow`` plotting
            an image to point within.
        kwargs (dict):
            Passed directly to ``matplotlib.widgets.Cursor``.
    """
    def __init__(self, ax, **kwargs):
        super().__init__(ax)

        if 'name' in kwargs:
            self.name = kwargs['name']
            kwargs.pop('name')
        else:
            self.name = None

        self.cursor = widgets.Cursor(ax, useblit=True, **kwargs)

        self.connect_event('button_press_event', self._button_update)
        self.connect_event('button_release_event', self._button_update)
        self.connect_event('key_press_event', self._key_update)
        self.connect_event('key_release_event', self._key_update)
        self.observers = {}
        self.observers_descr = {}
        self.drag_active = False
        self.pos = None
        self.action = None

    def _event_update(self, event, event_type):
        """update the pointer position"""
        if self.ignore(event):
            return

        if event_type not in ['button', 'key']:
            raise ValueError(f'Event type must be button or key, not {event_type}.')

        if event.name == f'{event_type}_press_event' and event.inaxes is self.ax:
            self.drag_active = True
            event.canvas.grab_mouse(self.ax)

        if not self.drag_active:
            return

        if event.name == f'{event_type}_release_event' \
                or (event.name == f'{event_type}_press_event' and event.inaxes is not self.ax):
            self.drag_active = False
            event.canvas.release_mouse(self.ax)
            return

        self._set_event(getattr(event, event_type), event.xdata, event.ydata)

    def _button_update(self, event):
        self._event_update(event, 'button')

    def _key_update(self, event):
        self._event_update(event, 'key')

    def _set_event(self, action, x, y):
        self.action = action
        self.pos = (x, y)
        if not self.eventson:
            return
        if self.action not in self.observers:
            print(f'No action: {action} ({self.name}, {x}, {y})')
            return
        if self.action in self.observers:
            self.observers[self.action](self.pos)

    def register(self, action, func, descr=None):
        """
        Register a function to associate with a specific button or key press.
        """
        self.observers[action] = func
        self.observers_descr[action] = descr

    def disconnect(self, action):
        """remove the observer with connection id *cid*"""
        try:
            del self.observers[action]
            del self.observers_descr[action]
        except KeyError:
            pass

    def build_help(self):
        if '?' in self.observers:
            msgs.warn('Help key (?) is already registered and will be overwritten')
            self.disconnect('?')
        self.register('?', self.print_help, descr='Print this list of key bindings')

    def print_help(self):
        print('-'*50)
        print('Key bindings')
        for action in self.observers.keys():
            descr = 'No description given' if self.observers_descr[action] is None \
                        else self.observers_descr[action]
            print(f'    {action}: {descr}')
        print('-'*50)


class UpdateableRangeSlider(widgets.RangeSlider):
    def __init__(self, ax, label, valmin, valmax, valinit=None, valfmt=None, closedmin=True,
                 closedmax=True, dragging=True, valstep=None, orientation="horizontal",
                 track_color='lightgrey', handle_style=None, **kwargs):

        super().__init__(ax, label, valmin, valmax, valinit=valinit, valfmt=valfmt,
                         closedmin=closedmin, closedmax=closedmax, dragging=dragging,
                         valstep=valstep, orientation=orientation, track_color=track_color,
                         handle_style=handle_style, **kwargs)
        self.label.set_position((0.5, 0.8))
        self.label.set_verticalalignment('bottom')
        self.label.set_horizontalalignment('center')
        self.label.set_weight('bold')
        # "Removes" the labels showing the current range
        self.valtext.set_visible(False)

    def update_range(self, rng, label=None):
        self._active_handle = None
        xy = self.poly.get_xy()
        xy[:,0] = np.roll(np.repeat(rng, 3)[:-1],-1)
        self.poly.set_xy(xy)
        self.valmin, self.valmax = rng
        self.valinit = np.array(rng)
        self._handles[0].set_xdata(np.array([rng[0]]))
        self._handles[1].set_xdata(np.array([rng[1]]))
        self.ax.set_xlim(rng)
        if label is not None:
            self.label.set_text(label)
        self.set_val(rng)


class UpdateableImage:
    def __init__(self, images, image_plot, slider):
        self.showing = 0
        self.images = images if isinstance(images, list) else [images]
        self.nimages = len(self.images)
        self.image_plot = image_plot
        self.slider = slider
        self.slider.on_changed(self.change_range)

    def change_range(self, val):
        self.image_plot.set_clim(*val)
        pyplot.draw()

    def next_image(self, *args, **kwargs):
        self.showing += 1
        if self.showing >= self.nimages:
            self.showing = 0
        self.image_plot.set_data(self.images[self.showing])
        rng = [np.amin(self.images[self.showing]), np.amax(self.images[self.showing])]
        self.slider.update_range(rng)
        self.change_range(rng)


def clean_pyplot_keymap():
    """
    Remove all default key bindings for matplotlib plot window.
    """
    for key in pyplot.rcParams.keys():
        if 'keymap' in key:
            pyplot.rcParams[key] = []

class EdgeInspectorGUI:
    """
    Pointer used when analyzing a trace/sobel image.

    Args:
        edges (:class:`~pypeit.edgetrace.EdgeTraceSet`):
            Edge tracing to inspect/edit.
    """
    def __init__(self, edges):

        # Remove all matplotlib plotting window key bindings to make sure they
        # don't overlap with any defined by the function.
        clean_pyplot_keymap()

        # Layout 
        sx = 0.07
        sy = 0.15
        dx = 0.78
        dy = 0.78

        # Figure
        w,h = pyplot.figaspect(1)
        fig = pyplot.figure(figsize=(1.5*w,1.5*h))

        # Axes
        image_ax = fig.add_axes([sx, sy, dx, dy])
        image_ax.tick_params(which='both', top=True, right=True)
        image_ax.minorticks_on()
        image_cax = fig.add_axes([sx + dx + 0.02, sy, 0.01, dy])
        image_slider_ax = fig.add_axes([sx, sy-0.08, dx/1.2, 0.02])
        image_change_button_ax = fig.add_axes([sx + dx - 0.1, sy - 0.1, 0.1, 0.05])

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
        # Add the colorbar
        self.image_cb = fig.colorbar(self.image_plot, cax=image_cax)
        # And the colorbar slider
        self.image_slider = UpdateableRangeSlider(image_slider_ax, 'Data Range', lim[0], lim[1],
                                                  valinit=lim)
        # And the ability to update the image
        self.image_updater = UpdateableImage(self.images, self.image_plot, self.image_slider)
        # And the button that will toggle the images
        self.image_change_button = widgets.Button(image_change_button_ax, 'Toggle', color='0.7',
                                                  hovercolor='0.9')
        self.image_change_button.on_clicked(self.image_updater.next_image)

        # Set the pointer
        self.img_pointer = Pointer(self.image_plot.axes, name='image', color='C1', lw=0.5)
        self.img_pointer.register('m', self.move)
        self.img_pointer.build_help()

    def close(self):
        pyplot.rcdefaults()
        pyplot.close()

    def move(self, pos):
        print('move', pos)

"""
        # Build the image pointer
        self.img_pointer = Pointer(self.img_plt.axes, name='image', color='C1', lw=0.5)
        self.img_pointer.register(MouseButton.LEFT, self.set_level)
        self.img_pointer.register('b', self.set_bkg_lower_image)
        self.img_pointer.register('B', self.set_bkg_upper_image)
        self.img_pointer.register('D', self.rm_bkg)
        self.img_pointer.register('C', self.print_contour)

        # NOTE: Have to use ee_plt here because it is the "top" axis, overlaying
        # the flux axis.  This is okay because we only care about the radius
        # position.
#        self.flx_pointer = Pointer(self.flux_sc.axes, name='scatter')
        self.flx_pointer = Pointer(self.ee_plt.axes, name='scatter', color='C1', lw=0.5)
        self.flx_pointer.register('b', self.set_bkg_lower_scatter)
        self.flx_pointer.register('B', self.set_bkg_upper_scatter)
        self.flx_pointer.register('D', self.rm_bkg)

    def set_level(self, pos):
        self.level = self.ee.img[int(pos[1]), int(pos[0])]
        self.update()

    def get_image_radius(self, pos):
        return numpy.sqrt((self.ee.circ_p[0] - pos[0])**2 + (self.ee.circ_p[1]-pos[1])**2) \
                / self.ee.circ_p[2]

    def set_bkg_lower_image(self, pos):
        self.set_bkg_lower(self.get_image_radius(pos))

    def set_bkg_lower_scatter(self, pos):
        r = contour.convert_radius(pos[0], pixelsize=self.pixelsize, distance=self.distance,
                                   inverse=True)[1] / self.ee.circ_p[2]
        self.set_bkg_lower(r)

    def set_bkg_lower(self, r):
        if self.bkg_lim is None:
            self.bkg_lim = [r, None]
        else:
            self.bkg_lim[0] = r
        self.update()

    def set_bkg_upper_image(self, pos):
        self.set_bkg_upper(self.get_image_radius(pos))

    def set_bkg_upper_scatter(self, pos):
        r = contour.convert_radius(pos[0], pixelsize=self.pixelsize, distance=self.distance,
                                   inverse=True)[1] / self.ee.circ_p[2]
        self.set_bkg_upper(r)

    def set_bkg_upper(self, r):
        if self.bkg_lim is None:
            warnings.warn('Ignoring input.  Set the lower limit first!')
            return
        self.bkg_lim[1] = r
        self.update()

    def rm_bkg(self, pos):
        self.bkg_lim = None
        if self.bkg_lo_line is not None:
            self.bkg_lo_line.remove()
            self.bkg_lo_line = None
        if self.bkg_hi_line is not None:
            self.bkg_hi_line.remove()
            self.bkg_hi_line = None
        self.update()

    def print_contour(self, pos):
        filename = input('File for contour data: ')
        if len(filename) == 0:
            for t in self.ee.trace:
                print(f'{t[0]:7.2f} {t[1]:7.2f}')
            return
        numpy.savetxt(filename, self.ee.trace, fmt=' %7.2f %7.2f')

    def update(self):
        if self.bkg_lim is not None and self.bkg_lim[1] is not None \
                and self.bkg_lim[0] > self.bkg_lim[1]:
            warnings.warn('Lower background boundary at larger radius than upper boundary.  '
                          'Swapping order.')
            self.bkg_lim = [self.bkg_lim[1], self.bkg_lim[0]]
        _img = numpy.ma.MaskedArray(self.ee.img, mask=numpy.logical_not(self.ee.inp_gpm))
        ee, xc, yc, rc, radius, flux, model_radius, model_flux, normalized_ee, ee90, rms \
                = reset_ee(_img, self.level, self.bkg_lim, self.pixelsize, self.distance)
        
        self.ee = ee
        
        self.img_plt.set_data(self.ee.img)
        self.img_cen.set_offsets((xc,yc))
        self.img_contour.set_data((self.ee.trace[:,0], self.ee.trace[:,1]))

        self.mod_plt.set_data(self.ee.model_img)
        self.mod_cen.set_offsets((xc,yc))
        self.mod_contour.set_data((self.ee.trace[:,0], self.ee.trace[:,1]))

        self.res_plt.set_data(self.ee.img - self.ee.model_img)

        self.flux_sc.set_offsets(numpy.column_stack((radius, flux)))
        self.mod_flux_plt.set_data((model_radius, model_flux))

        self.circ_line.set_data(([rc, rc], [0,1]))
        if self.bkg_lim is not None:
            if self.bkg_lim[1] is None:
                lo = self.bkg_lim[0]*rc
                hi = None
            else:
                lo = self.bkg_lim[0]*rc
                hi = self.bkg_lim[1]*rc

            if self.bkg_lo_line is None:
                self.bkg_lo_line \
                        = self.flux_sc.axes.axvline(lo, color='0.5', lw=1, ls='--', zorder=6) 
            else:
                self.bkg_lo_line.set_data(([lo,lo],[0,1]))

            if hi is None and self.bkg_hi_line is not None:
                self.bkg_hi_line.remove()
                self.bkg_hi_line = None
            elif hi is not None:
                if self.bkg_hi_line is None:
                    self.bkg_hi_line \
                            = self.flux_sc.axes.axvline(hi, color='0.5', lw=1, ls='--', zorder=6)
                else:
                    self.bkg_hi_line.set_data(([hi,hi],[0,1]))
        self.ee_plt.set_data((radius, normalized_ee))

        self.text['bkg'].set_text(f'{self.ee.bkg:.4e}')
        self.text['ee90'].set_text(f'{ee90:.4e}')
        self.text['flux'].set_text(f'{self.ee.ee_norm:.4e}')
        self.text['rms_full'].set_text(f'{rms[0]:.4e}')
        self.text['rms_90'].set_text(f'{rms[1]:.4e}')

        pyplot.draw()
"""


