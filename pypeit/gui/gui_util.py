"""
GUI utilities

.. include:: ../include/links.rst
"""

import warnings

import numpy as np
from matplotlib import pyplot, widgets


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
        """
        Update the pointer position

        Args:
            event (`matplotlib.backend_bases.Event`_):
                Key (or mouse) or button event instance.
            event_type (:obj:`str`):
                Type of event.  Must be 'key' or 'button'.
        """
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
        """Execute an event created by a button press."""
        self._event_update(event, 'button')

    def _key_update(self, event):
        """Execute an event created by a key or mouse press."""
        self._event_update(event, 'key')

    def _set_event(self, action, x, y):
        """
        Execute an event.

        Args:
            action (:obj:`str`):
                The keyword for the event action to perform
            x (:obj:`float`):
                X position in the pyplot window where the event took place
            y (:obj:`float`):
                Y position in the pyplot window where the event took place
        """
        self.action = action
        self.pos = (x, y)
        if not self.eventson:
            return
        if self.action not in self.observers:
            print(f'No action: {action} ({self.name}, {x}, {y})')
            return
        self.observers[self.action](self.pos)

    def register(self, action, func, descr=None):
        """
        Register a function to associate with a specific button or key press.

        *All* functions must have the same calling sequence (see
        :func:`_set_event`), which is that they only accept a tuple with the
        coordinates of the cursor when the event occurred.

        Args:
            action (:obj:`str`):
                The keyword for the event action to perform
            func (callable):
                The function to call when the corresponding event is triggered
            descr (:obj:`str`, optional):
                The description of the action being taken.  Used to build a help
                dialog.
        """
        self.observers[action] = func
        self.observers_descr[action] = descr

    def disconnect(self, action):
        """
        Remove the action from the register.

        Args:
            action (:obj:`str`):
                The keyword for the event action to remove
        """
        try:
            del self.observers[action]
        except KeyError:
            pass
        try:
            del self.observers_descr[action]
        except KeyError:
            pass

    def build_help(self):
        """
        Register the help dialog.  Any action already assigned to the '?' key
        will be removed!
        """
        if '?' in self.observers:
            warnings.warn('Help key (?) is already registered and will be overwritten')
            self.disconnect('?')
        self.register('?', self.print_help, descr='Print this list of key bindings')

    def print_help(self, pos):
        """
        Print the help dialog.

        This is an event call-back function that must accept an array-like
        object giving the pointer coordinates.

        Args:
            pos (array-like):
                List with x and y position of the cursor at the time of the
                window event.
        """
        print('-'*50)
        print('Key bindings')
        for action in self.observers.keys():
            descr = 'No description given' if self.observers_descr[action] is None \
                        else self.observers_descr[action]
            print(f'    {action}: {descr}')
        print('-'*50)


class UpdateableRangeSlider(widgets.RangeSlider):
    """
    A range slider.

    This is virtually identical to the base class, but with a few customizations
    regarding where labels are placed (or removed).

    See `matplotlib.widgets.RangeSlider`_ for the argument descriptions.
    """
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
        """
        Update the slider to cover a different range.

        Args:
            rng (array-like):
                Two-element array like object setting the new limits for the
                slider.
            label (:obj:`str`, optional):
                New label for the updated slider.  If None, the label remains
                unchanged.
        """
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
    """
    Provides an interface to change the Z range of an image and toggle between a
    set of images.

    Args:
        images (`numpy.ndarray`_, :obj:`list`):
            One or more images to plot.  The images should have the same shape.
        image_plot (`matplotlib.image.AxesImage`_):
            Object returned by `matplotlib.pyplot.imshow`_.
        slider (:class:`UpdateableRangeSlider`):
            Slider used to adjust the image plot limits.
    """
    def __init__(self, images, image_plot, slider):
        self.showing = 0
        self.images = images if isinstance(images, list) else [images]
        self.nimages = len(self.images)
        self.image_plot = image_plot
        self.slider = slider
        self.slider.on_changed(self.change_range)

    def change_range(self, val):
        """
        Adjust the plotted range.

        Args:
            val (:obj:`list`):
                New range to plot
        """
        self.image_plot.set_clim(*val)
        pyplot.draw()

    def next_image(self, *args):
        """
        Got to the next image in the list.

        All arguments to this function are accepted but ignored.  The reason is
        because this is the function passed to
        `matplotlib.widgets.Button.on_clicked`_, but this function does not need
        any of the input from the ``Button`` event.
        """
        self.showing += 1
        if self.showing >= self.nimages:
            self.showing = 0
        self.image_plot.set_data(self.images[self.showing])
        rng = [np.amin(self.images[self.showing]), np.amax(self.images[self.showing])]
        self.slider.update_range(rng)
        self.change_range(rng)


def clean_pyplot_keymap():
    """
    Remove all default key bindings for a matplotlib plot window.
    """
    for key in pyplot.rcParams.keys():
        if 'keymap' in key:
            pyplot.rcParams[key] = []


