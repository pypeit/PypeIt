"""
View plugins that control the displays made by the PypeIt ginga scripts.

Plugin Type: Local
==================

``Spec2dView``, ``EdgesView`` and ``FlatsView`` are local plugins.  The
``PypeItRemote`` plugin opens one on the first channel of a display made by
``pypeit_show_2dspec``, ``pypeit_chk_edges`` or ``pypeit_chk_flats``,
respectively.  They are not meant to be opened by hand.

Usage
-----

The plugin panel has four parts:

- **File**: the file being displayed.  **Reload** reads it again, e.g. after
  re-running a reduction step.

- **Overlays**: check boxes that show or hide the slit edges, slit IDs,
  object traces and object names.  Each overlay is shared by all channels of
  the display, so a check box affects every channel at once.

- **Options**: options of the script that made the display, e.g. the
  detector shown by ``pypeit_show_2dspec`` or the flats shown by
  ``pypeit_chk_flats``.  **Show** re-displays the file with the new options.

- **Slits**: the slits in the display.  Selecting a slit pans all channels to
  it and zooms so that it fills about a third of the window width.

Errors from **Reload** and **Show** are reported in the ginga Errors tab.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import argparse
import json
from pathlib import Path

from ginga import GingaPlugin
from ginga.gw import Widgets

from pypeit import __version__
from pypeit.display import ginga_draw

__all__ = ['Spec2dView', 'EdgesView', 'FlatsView']


class PypeItView(GingaPlugin.LocalPlugin):
    """
    Base class for the plugins that control a PypeIt display.

    The plugin receives a :class:`~pypeit.display.remote_handlers.ViewState`
    from the command that made the display, through :meth:`set_state`, before
    it is started.

    Parameters
    ----------
    fv : ``ginga.rv.Control.GingaShell``
        The ginga reference viewer shell.
    fitsimage : ``ginga.ImageView.ImageViewBase``
        The viewer of the channel the plugin is opened on.
    """
    script = None
    """Name of the script that makes the displays controlled by the plugin."""

    def __init__(self, fv, fitsimage):
        super().__init__(fv, fitsimage)
        self.state = None
        self.w = None
        self.gui_up = False

    def set_state(self, state):
        """
        Set what the display shows; called before the plugin is started.

        Parameters
        ----------
        state : :class:`~pypeit.display.remote_handlers.ViewState`
            What was displayed.
        """
        self.state = state

    def build_options(self, vbox):
        """
        Add the widgets of the script options to the Options frame.

        Subclasses with options override this, and :meth:`options_args`.

        Parameters
        ----------
        vbox : ``ginga.gw.Widgets.VBox``
            Container for the widgets.

        Returns
        -------
        :obj:`bool`
            True if any widgets were added.
        """
        return False

    def options_args(self, args):
        """
        Return the command arguments updated with the options set in the panel.

        Parameters
        ----------
        args : `argparse.Namespace`_
            A copy of the current arguments, which may be changed.

        Returns
        -------
        `argparse.Namespace`_
            The arguments to re-run the command with.
        """
        return args

    def build_gui(self, container):
        """
        Construct the plugin panel.

        Parameters
        ----------
        container : ``ginga.gw.Widgets.VBox``
            Container for the panel.
        """
        top = Widgets.VBox()
        top.set_border_width(4)
        vbox = Widgets.VBox()
        vbox.set_spacing(4)
        self.w = {}

        if self.state is None:
            vbox.add_widget(
                Widgets.Label(f'Nothing to control: use {self.script} to display data.'),
                stretch=0
            )
        else:
            vbox.add_widget(self._build_file_frame(), stretch=0)
            if len(self.state.overlays) > 0:
                vbox.add_widget(self._build_overlay_frame(), stretch=0)
            fr = Widgets.Frame('Options')
            opts = Widgets.VBox()
            opts.set_spacing(2)
            if self.build_options(opts):
                btn = Widgets.Button('Show')
                btn.set_tooltip('Display the file again with these options')
                btn.add_callback('activated', lambda w: self.show_cb())
                opts.add_widget(btn, stretch=0)
                fr.set_widget(opts)
                vbox.add_widget(fr, stretch=0)
            if len(self.state.slits) > 0:
                vbox.add_widget(self._build_slit_frame(), stretch=1)
            self.w['status'] = Widgets.Label('')
            vbox.add_widget(self.w['status'], stretch=0)

        top.add_widget(vbox, stretch=1)
        btns = Widgets.HBox()
        btns.set_spacing(3)
        btn = Widgets.Button('Close')
        btn.add_callback('activated', lambda w: self.close())
        btns.add_widget(btn, stretch=0)
        btn = Widgets.Button('Help')
        btn.add_callback('activated', lambda w: self.help())
        btns.add_widget(btn, stretch=0)
        btns.add_widget(Widgets.Label(''), stretch=1)
        top.add_widget(btns, stretch=0)
        container.add_widget(top, stretch=1)
        self.gui_up = True

    def _build_file_frame(self):
        """
        Build the File frame.

        Returns
        -------
        ``ginga.gw.Widgets.Frame``
            The frame.
        """
        fr = Widgets.Frame('File')
        hbox = Widgets.HBox()
        hbox.set_spacing(4)
        name = Widgets.Label(Path(self.state.path).name)
        name.set_tooltip(self.state.path)
        hbox.add_widget(name, stretch=1)
        btn = Widgets.Button('Reload')
        btn.set_tooltip('Read the file again and redisplay it')
        btn.add_callback('activated', lambda w: self.reload_cb())
        hbox.add_widget(btn, stretch=0)
        fr.set_widget(hbox)
        return fr

    def _build_overlay_frame(self):
        """
        Build the Overlays frame.

        Returns
        -------
        ``ginga.gw.Widgets.Frame``
            The frame.
        """
        fr = Widgets.Frame('Overlays')
        vbox = Widgets.VBox()
        for overlay in self.state.overlays:
            cb = Widgets.CheckBox(overlay.name)
            cb.set_state(overlay.visible)
            cb.add_callback('activated', self.toggle_overlay_cb, overlay)
            vbox.add_widget(cb, stretch=0)
        fr.set_widget(vbox)
        return fr

    def _build_slit_frame(self):
        """
        Build the Slits frame with the slit navigator.

        Returns
        -------
        ``ginga.gw.Widgets.Frame``
            The frame.
        """
        fr = Widgets.Frame('Slits')
        table = Widgets.TreeView(selection='single', sortable=True)
        columns = [('Slit', 'Slit'), ('Mask ID', 'Mask ID'), ('Center', 'Center'),
                   ('Width', 'Width')]
        table.setup_table(columns, 1, 'Slit')
        tree = {}
        for row in self.state.slits:
            tree[str(row['id'])] = {
                'Slit': row['id'], 'Mask ID': '' if row['maskdef_id'] is None
                                              else row['maskdef_id'],
                'Center': f"{row['x']:.1f}", 'Width': f"{row['width']:.1f}"
            }
        table.set_tree(tree)
        table.add_callback('selected', self.select_slit_cb)
        self.w['slits'] = table
        fr.set_widget(table)
        return fr

    def toggle_overlay_cb(self, w, state, overlay):
        """
        Callback of an overlay check box; see :meth:`set_overlay`.

        Parameters
        ----------
        w : ``ginga.gw.Widgets.CheckBox``
            The check box.
        state : :obj:`bool`
            True to show the overlay.
        overlay : :class:`~pypeit.display.remote_handlers.Overlay`
            The overlay.
        """
        self.set_overlay(overlay.name, state)

    def set_overlay(self, name, visible):
        """
        Show or hide an overlay on all of its channels.

        Parameters
        ----------
        name : :obj:`str`
            Name of the overlay, e.g. ``'Slit edges'``.
        visible : :obj:`bool`
            True to show the overlay.
        """
        overlay = [o for o in self.state.overlays if o.name == name][0]
        overlay.visible = visible
        for chname in overlay.chnames:
            if not self.fv.has_channel(chname):
                continue
            channel = self.fv.get_channel(chname)
            if visible:
                ginga_draw.add_overlay(channel, overlay.canvas, overlay.tag)
            else:
                ginga_draw.remove_overlay(channel, overlay.tag)

    def select_slit_cb(self, w, selected):
        """
        Callback of the slit table; see :meth:`goto_slit`.

        Parameters
        ----------
        w : ``ginga.gw.Widgets.TreeView``
            The slit table.
        selected : :obj:`dict`
            The selected rows, keyed by slit ID.
        """
        if len(selected) > 0:
            self.goto_slit(int(list(selected.keys())[0]))

    def goto_slit(self, slit_id):
        """
        Pan all channels to a slit and zoom so that it fills about a third of
        the window width.

        Parameters
        ----------
        slit_id : :obj:`int`
            ID of the slit.
        """
        row = [r for r in self.state.slits if r['id'] == slit_id][0]
        # The channels are linked, so moving the reference moves them all
        viewer = self.fv.get_channel(self.state.chnames[0]).fitsimage
        width, _ = viewer.get_window_size()
        scale = min(width / (3*max(row['width'], 1.)), 8.)
        viewer.scale_to(scale, scale)
        viewer.set_pan(row['x'], row['y'])

    def rerun(self, args):
        """
        Display the file again, by sending the command to ``PypeItRemote``.

        Parameters
        ----------
        args : `argparse.Namespace`_
            Arguments for the command.
        """
        status = self.fv.call_global_plugin_method(
            'PypeItRemote', 'submit', [self.state.command, json.dumps(vars(args)), __version__],
            {'notify': True}
        )
        if status['state'] == 'error':
            self.fv.show_error(f'{self.state.command} failed: {status["message"]}')
        elif self.gui_up:
            self.w['status'].set_text('Loading...')

    def reload_cb(self):
        """
        Read the file again and redisplay it with the same options.
        """
        self.rerun(argparse.Namespace(**vars(self.state.args)))

    def show_cb(self):
        """
        Redisplay the file with the options set in the panel.
        """
        self.rerun(self.options_args(argparse.Namespace(**vars(self.state.args))))

    def close(self):
        """
        Close the plugin.
        """
        self.fv.stop_local_plugin(self.chname, str(self))

    def stop(self):
        """
        Called when the plugin is closed.
        """
        self.gui_up = False

    def __str__(self):
        return type(self).__name__.lower()


class Spec2dView(PypeItView):
    """
    Plugin that controls a ``pypeit_show_2dspec`` display.

    Its options select the detector, the channels to show, and the mask,
    wavelength-RMS and extraction-mask options of ``pypeit_show_2dspec``.
    """
    script = 'pypeit_show_2dspec'
    channel_names = ('sciimg', 'skysub', 'sky_resid', 'resid')
    """Names of the channels selected by the ``--channels`` option, in order."""

    def build_options(self, vbox):
        """
        Add the detector, channel, and image options.

        Parameters
        ----------
        vbox : ``ginga.gw.Widgets.VBox``
            Container for the widgets.

        Returns
        -------
        :obj:`bool`
            True, because options were added.
        """
        args = self.state.args
        detectors = self.state.options['detectors']
        hbox = Widgets.HBox()
        hbox.add_widget(Widgets.Label('Detector:'), stretch=0)
        combo = Widgets.ComboBox()
        for det in detectors:
            combo.append_text(det)
        combo.set_index(detectors.index(self.state.options['detname']))
        hbox.add_widget(combo, stretch=1)
        vbox.add_widget(hbox, stretch=0)
        self.w['det'] = combo

        shown = [0, 1, 2, 3] if args.channels is None \
                    else [int(c) for c in args.channels.split(',')]
        hbox = Widgets.HBox()
        self.w['channels'] = []
        for i, name in enumerate(self.channel_names):
            cb = Widgets.CheckBox(name)
            cb.set_state(i in shown)
            hbox.add_widget(cb, stretch=0)
            self.w['channels'].append(cb)
        vbox.add_widget(hbox, stretch=0)

        for key, label, value, tip in (
                ('showmask', 'Show mask', args.showmask is not None,
                 'Add a channel with the mask bits'),
                ('show_wave_rms', 'Wavelength RMS', args.show_wave_rms,
                 'Show the RMS of the wavelength solution in the cursor readout'),
                ('ignore_extract_mask', 'Ignore extraction mask', args.ignore_extract_mask,
                 'Ignore the extraction mask in the sky-subtracted images')):
            cb = Widgets.CheckBox(label)
            cb.set_state(value)
            cb.set_tooltip(tip)
            vbox.add_widget(cb, stretch=0)
            self.w[key] = cb
        return True

    def options_args(self, args):
        """
        Return the arguments with the detector, channels and image options set.

        Parameters
        ----------
        args : `argparse.Namespace`_
            A copy of the current arguments.

        Returns
        -------
        `argparse.Namespace`_
            The arguments to re-run ``show_2dspec`` with.
        """
        args.det = self.state.options['detectors'][self.w['det'].get_index()]
        shown = [str(i) for i, cb in enumerate(self.w['channels']) if cb.get_state()]
        args.channels = None if len(shown) == len(self.channel_names) else ','.join(shown)
        if not self.w['showmask'].get_state():
            args.showmask = None
        elif args.showmask is None:
            args.showmask = []
        args.show_wave_rms = self.w['show_wave_rms'].get_state()
        args.ignore_extract_mask = self.w['ignore_extract_mask'].get_state()
        return args


class FlatsView(PypeItView):
    """
    Plugin that controls a ``pypeit_chk_flats`` display.

    Its option selects which flats are shown.
    """
    script = 'pypeit_chk_flats'
    flat_types = ('all', 'pixel', 'illum')
    """Values of the ``--type`` option."""

    def build_options(self, vbox):
        """
        Add the flat-type option.

        Parameters
        ----------
        vbox : ``ginga.gw.Widgets.VBox``
            Container for the widgets.

        Returns
        -------
        :obj:`bool`
            True, because options were added.
        """
        hbox = Widgets.HBox()
        hbox.add_widget(Widgets.Label('Flats:'), stretch=0)
        combo = Widgets.ComboBox()
        for name in self.flat_types:
            combo.append_text(name)
        combo.set_index(self.flat_types.index(self.state.args.type))
        hbox.add_widget(combo, stretch=1)
        vbox.add_widget(hbox, stretch=0)
        self.w['type'] = combo
        return True

    def options_args(self, args):
        """
        Return the arguments with the flat type set.

        Parameters
        ----------
        args : `argparse.Namespace`_
            A copy of the current arguments.

        Returns
        -------
        `argparse.Namespace`_
            The arguments to re-run ``chk_flats`` with.
        """
        args.type = self.flat_types[self.w['type'].get_index()]
        return args


class EdgesView(PypeItView):
    """
    Plugin that controls a ``pypeit_chk_edges`` display.
    """
    script = 'pypeit_chk_edges'
