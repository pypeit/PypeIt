"""
Spec2dView is a plugin for the Ginga viewer that displays PypeIt spec2d files.

Plugin Type: Global
===================

Spec2dView is a global plugin.  It creates one channel for each image shown
(science image, sky-subtracted image, sky residuals, model residuals, and,
optionally, the mask), overplots the slit edges and the extracted object
traces (from the spec1d file, if it exists), and matches the pan and zoom of
all channels.

Usage
-----

The plugin is used by ``pypeit_show_2dspec``.  It can also be started from
the "PypeIt" category of the ginga "Plugins" menu, after which a spec2d file
can be entered in the "File" box.

UI
--

- File: the spec2d file to show.

- Detector: the detector or mosaic to show.  Changing the detector replaces
  the channels shown for the previous detector.

- Channels: toggle the images (channels) shown.  Only the channel that is
  toggled is added or removed; the others are not changed.

- Overlays: toggle the slit edges and object traces in all channels.  For
  files with many traces (e.g., fiber-fed IFUs), the traces are drawn only
  when the view is zoomed in enough to distinguish them; hover over a trace
  to see its name.

Buttons
-------

- Reload: read the file again (e.g., after it has been re-reduced) and show
  it with the current settings.

"""
from pathlib import Path

from ginga import GingaPlugin
from ginga.gw import Widgets

from pypeit.display import display
from pypeit.display import spec2d_data

#: Tag prefixes of the overlays, used to remove them from the channels
OVERLAY_TAGS = dict(slits='pypeit_slits', traces='pypeit_traces')


class Spec2dView(GingaPlugin.GlobalPlugin):
    """
    ginga global plugin that displays PypeIt spec2d files.

    See the documentation for `ginga GlobalPlugin`_.

    Parameters
    ----------
    fv : ``ginga.rv.Control.GingaShell``
        The ginga reference viewer shell.
    """
    def __init__(self, fv):
        super().__init__(fv)
        self.w = None
        self.gui_up = False
        self.data = None
        # Channels shown, keyed by their index in spec2d_data.CHANNEL_NAMES
        # (or 'mask'), in the order they were created
        self.shown = {}
        self.params = dict(file=None, det=None, channels=None, prefix='', showmask=None,
                           removetrace=False, ignore_extract_mask=False, spat_id=None,
                           maskID=None, try_old=False, show_slits=True, show_traces=True)
        fv.add_callback('delete-channel', self.delete_channel_cb)

    def build_gui(self, container):
        """
        Construct the plugin GUI.

        Parameters
        ----------
        container : ``ginga.gw.Widgets.VBox``
            Container provided by ginga for the plugin widgets.
        """
        top = Widgets.VBox()
        top.set_border_width(4)
        top.set_spacing(4)

        fr = Widgets.Frame('File')
        captions = (('File:', 'label', 'file', 'entryset'),
                    ('Detector:', 'label', 'det', 'combobox'))
        w, b = Widgets.build_info(captions, orientation='vertical')
        self.w = b
        b.file.set_tooltip('spec2d file to show')
        b.file.add_callback('activated', self.set_file_cb)
        b.det.set_tooltip('Detector or mosaic to show')
        b.det.add_callback('activated', self.set_det_cb)
        fr.set_widget(w)
        top.add_widget(fr, stretch=0)

        fr = Widgets.Frame('Channels')
        vbox = Widgets.VBox()
        self.w.channels = {}
        for key, name in list(enumerate(spec2d_data.CHANNEL_NAMES)) + [('mask', 'mask')]:
            cb = Widgets.CheckBox(name)
            cb.set_state(key != 'mask')
            cb.set_tooltip('Show the mask bit values' if key == 'mask'
                           else f'Show the {name} image')
            cb.add_callback('activated', self.toggle_channel_cb, key)
            vbox.add_widget(cb, stretch=0)
            self.w.channels[key] = cb
        fr.set_widget(vbox)
        top.add_widget(fr, stretch=0)

        fr = Widgets.Frame('Overlays')
        vbox = Widgets.VBox()
        self.w.overlays = {}
        for which, label in [('slits', 'Slit edges'), ('traces', 'Object traces')]:
            cb = Widgets.CheckBox(label)
            cb.set_state(True)
            cb.add_callback('activated', self.toggle_overlay_cb, which)
            vbox.add_widget(cb, stretch=0)
            self.w.overlays[which] = cb
        fr.set_widget(vbox)
        top.add_widget(fr, stretch=0)

        btns = Widgets.HBox()
        btns.set_spacing(3)
        btn = Widgets.Button('Reload')
        btn.set_tooltip('Read the file again and show it with the current settings')
        btn.add_callback('activated', lambda w: self.reload_cb())
        btns.add_widget(btn, stretch=0)
        btns.add_widget(Widgets.Label(''), stretch=1)
        btn = Widgets.Button('Close')
        btn.add_callback('activated', lambda w: self.close())
        btns.add_widget(btn, stretch=0)
        btn = Widgets.Button('Help')
        btn.add_callback('activated', lambda w: self.help())
        btns.add_widget(btn, stretch=0)
        top.add_widget(btns, stretch=0)
        top.add_widget(Widgets.Label(''), stretch=1)

        container.add_widget(top, stretch=1)
        self.gui_up = True
        self._params_to_gui()

    # ------------------------------------------------------------------
    # API

    def load_spec2d(self, file, det=None, channels=None, prefix='', showmask=None,
                    removetrace=False, ignore_extract_mask=False, spat_id=None, maskID=None,
                    try_old=False, show_slits=True, show_traces=True, clear=True, data=None):
        """
        Show a spec2d file.

        See :func:`pypeit.display.spec2d_data.load_spec2d_display` for the
        description of most of the parameters.

        Parameters
        ----------
        file : :obj:`str`
            spec2d file.
        show_slits : :obj:`bool`, optional
            Overplot the slit edges.
        show_traces : :obj:`bool`, optional
            Overplot the object traces.
        clear : :obj:`bool`, optional
            Delete all existing channels before showing the file.  If False,
            the channels are added to (or replace channels with the same
            name in) the existing ones.
        data : :class:`~pypeit.display.spec2d_data.Spec2DDisplay`, optional
            The data to show, already read from ``file`` with the other
            parameters (e.g., by ``pypeit_show_2dspec``, before opening the
            viewer).  If None, the file is read.
        """
        self.params = dict(file=str(Path(file).absolute()), det=det, channels=channels,
                           prefix=prefix, showmask=showmask, removetrace=removetrace,
                           ignore_extract_mask=ignore_extract_mask, spat_id=spat_id,
                           maskID=maskID, try_old=try_old, show_slits=show_slits,
                           show_traces=show_traces)
        self.shown = {}
        self.reload_file(clear=clear, data=data)

    # NOTE: Do *not* name this method "redo".  ginga calls the redo(channel,
    # image) method of active global plugins every time an image is shown in
    # any channel.
    def reload_file(self, clear=False, data=None):
        """
        Read the file and show it with the current parameters.

        Parameters
        ----------
        clear : :obj:`bool`, optional
            Delete all existing channels first.  Otherwise, only the channels
            shown the last time the file was shown are deleted.  In either
            case, existing channels with the same names as the channels to
            show are reused.
        data : :class:`~pypeit.display.spec2d_data.Spec2DDisplay`, optional
            The data already read from the file.  If None, the file is read.
        """
        viewer = self._viewer()
        if self.params['file'] is None:
            return
        if data is None:
            pars = {k: v for k, v in self.params.items()
                        if k not in ['show_slits', 'show_traces']}
            data = spec2d_data.load_spec2d_display(**pars)
        self.data = data
        self.params['det'] = self.data.detname
        self._params_to_gui()

        # Delete the channels that are not reused.  Creating and deleting
        # channels is slow, so channels with the same name are reused; their
        # image and overlays are replaced.
        keep = [chan.name for chan in self.data.channels]
        # NOTE: Reset the shown channels before deleting any, so that
        # delete_channel_cb does not change the GUI.
        previous = viewer.channel_names() if clear else list(self.shown.values())
        self.shown = {}
        for chname in previous:
            if chname not in keep:
                viewer.delete_channel(chname)

        for chan in self.data.channels:
            self._show_channel(viewer, chan)

        if len(self.shown) > 0:
            chnames = list(self.shown.values())
            # Sync the pan and zoom of all channels
            viewer.set_wcs_match(chnames[-1])
            self.fv.change_channel(chnames[0])
            # Fit the images of the new channels to the window
            viewer.fit_channels(chnames)

    # ------------------------------------------------------------------
    # Internals

    def _viewer(self):
        """Return the viewer wrapper used to display the data."""
        from pypeit.display import viewer as pypeit_viewer
        if pypeit_viewer.viewer_exists() and pypeit_viewer.get_viewer().shell is self.fv:
            return pypeit_viewer.get_viewer()
        return pypeit_viewer.GingaShellViewer(self.fv)

    def _show_channel(self, viewer, chan):
        """Show one channel image and its overlays."""
        ch = viewer.channel(chan.name)
        # Only redraw the channel once everything has been added
        with ch.fitsimage.suppress_redraw:
            display.show_image(chan.image, chname=chan.name, waveimg=self.data.waveimg,
                               cuts=chan.cuts, viewer=viewer)
            for which in OVERLAY_TAGS.keys():
                if self.params[f'show_{which}']:
                    self._add_overlay(viewer, ch, chan, which)
        self.shown[chan.key] = chan.name

    def _add_overlay(self, viewer, ch, chan, which):
        """Add the slit edges or the object traces to a channel."""
        data = self.data
        if which == 'traces' and chan.show_traces and data.traces is not None:
            display.show_trace(viewer, ch, data.traces, data.trace_names,
                               maskdef_extr=data.maskdef_extract,
                               manual_extr=data.manual_extract, tagpfx=OVERLAY_TAGS[which])
        elif which == 'slits' and chan.show_slits and data.left is not None:
            display.show_slits(viewer, ch, data.left, data.right, slit_ids=data.slit_ids,
                               maskdef_ids=data.maskdef_ids, tagpfx=OVERLAY_TAGS[which])

    def _params_to_gui(self):
        """Set the GUI controls from the current parameters."""
        if not self.gui_up:
            return
        file = self.params['file']
        self.w.file.set_text('' if file is None else file)
        self.w.det.clear()
        if self.data is not None:
            for name in self.data.detnames:
                self.w.det.append_text(name)
            if self.data.detname in self.data.detnames:
                self.w.det.set_index(self.data.detnames.index(self.data.detname))
        channels = self.params['channels']
        for key, cb in self.w.channels.items():
            if key == 'mask':
                cb.set_state(self.params['showmask'] is not None)
            else:
                cb.set_state(channels is None or key in [int(c) for c in channels])
        for which, cb in self.w.overlays.items():
            cb.set_state(self.params[f'show_{which}'])

    def _selected_channels(self):
        """Return the indices of the selected (non-mask) channels."""
        return [key for key, cb in self.w.channels.items() if key != 'mask' and cb.get_state()]

    def _report(self, msg):
        """Show an error message."""
        self.fv.show_error(msg)

    # ------------------------------------------------------------------
    # Callbacks

    def set_file_cb(self, w):
        """Callback for the file entry: show a new file."""
        self.params['file'] = self.w.file.get_text().strip() or None
        self.params['det'] = None
        self.data = None
        try:
            self.reload_file(clear=False)
        except Exception as e:
            self._report(f'Could not show {self.params["file"]}: {e}')

    def set_det_cb(self, w, idx):
        """Callback for the detector selection: show a different detector."""
        if self.data is None or idx < 0 or idx >= len(self.data.detnames):
            return
        self.params['det'] = self.data.detnames[idx]
        try:
            self.reload_file(clear=False)
        except Exception as e:
            self._report(f'Could not show {self.params["det"]}: {e}')

    def reload_cb(self):
        """Callback for the Reload button."""
        try:
            self.reload_file(clear=False)
        except Exception as e:
            self._report(f'Could not show {self.params["file"]}: {e}')

    def toggle_channel_cb(self, w, state, key):
        """
        Callback for the channel check boxes: add or remove one channel.

        The other channels are not changed.
        """
        if key == 'mask':
            if state and self.params['showmask'] is None:
                self.params['showmask'] = []
            elif not state:
                self.params['showmask'] = None
        else:
            self.params['channels'] = self._selected_channels()
        if self.data is None:
            return

        viewer = self._viewer()
        if not state:
            if key in self.shown:
                viewer.delete_channel(self.shown.pop(key))
            return

        if key == 'mask':
            self.data.showmask = self.params['showmask']
        chan = self.data.make_channel(key)
        if chan is None:
            w.set_state(False)
            self._report(f'The data needed to show this channel are not available.')
            return
        linked = list(self.shown.values())
        self._show_channel(viewer, chan)
        if len(linked) > 0:
            # Match the pan and zoom of the channels already shown
            viewer.set_wcs_match(linked[0])
            viewer.fit_channels([chan.name], fit=False)
        else:
            viewer.fit_channels([chan.name])
        self.fv.change_channel(chan.name)

    def toggle_overlay_cb(self, w, state, which):
        """
        Callback for the overlay check boxes: add or remove the slit edges
        or the object traces in all channels.

        The images are not changed.
        """
        self.params[f'show_{which}'] = state
        if self.data is None:
            return
        viewer = self._viewer()
        for key, chname in self.shown.items():
            ch = viewer.channel(chname)
            with ch.fitsimage.suppress_redraw:
                display.clear_overlays(ch, OVERLAY_TAGS[which])
                if state:
                    self._add_overlay(viewer, ch, self._channel_info(key), which)

    def _channel_info(self, key):
        """Return the overlay settings of a shown channel without its image."""
        return spec2d_data.ChannelImage(name=self.shown[key], image=None, key=key,
                                        show_slits=key != 'mask',
                                        show_traces=key not in ['mask']
                                            and (key == 0 or not self.data.removetrace))

    def delete_channel_cb(self, fv, channel):
        """Callback when a channel is deleted (e.g., by closing its tab)."""
        for key, chname in list(self.shown.items()):
            if chname == channel.name:
                del self.shown[key]
                if self.gui_up:
                    self.w.channels[key].set_state(False)

    def start(self):
        """Plugin start (nothing to do)."""
        pass

    def stop(self):
        """Plugin stop."""
        self.gui_up = False

    def close(self):
        """Close the plugin."""
        self.fv.stop_global_plugin(str(self))
        return True

    def __str__(self):
        return 'spec2dview'
