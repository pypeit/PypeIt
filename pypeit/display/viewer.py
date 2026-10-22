"""
In-process ginga viewer used by all PypeIt displays.

PypeIt builds a single ginga reference viewer per python process the first
time something is displayed and keeps it for the lifetime of the process.
Images and overlays are handed to ginga as python objects; nothing is
serialized or sent over a socket.

Displays are *blocking*: after the images and overlays are added, calling
:func:`show` runs the ginga event loop until the user closes the window (which
only hides it).  The window title says what closing it does.  Control then
returns to the caller, and the next display reuses the same window.

A viewer can optionally listen on a local socket (see :meth:`PypeItViewer.listen`)
so that a second PypeIt process can ask it to load a file (e.g.,
``pypeit_show_2dspec --no_clear``).  Only a small JSON command with the file
path and options is sent; the receiving viewer reads the file itself.

.. include:: ../include/links.rst
"""
import argparse
import atexit
import contextlib
import json
import os
import select
import socket
import sys
import threading
import time
from pathlib import Path

from pypeit import log
from pypeit import PypeItError

#: Default local port used by a viewer to accept commands from other PypeIt
#: processes.  Can be overridden by the ``PYPEIT_GINGA_PORT`` environment
#: variable.
DEFAULT_PORT = 9919

_viewer = None


def get_port():
    """
    Return the local port used to send commands to an open PypeIt viewer.

    Returns
    -------
    :obj:`int`
        Port number; ``PYPEIT_GINGA_PORT`` if defined, otherwise
        :data:`DEFAULT_PORT`.
    """
    return int(os.environ.get('PYPEIT_GINGA_PORT', DEFAULT_PORT))


def get_ginga_home():
    """
    Return the directory with the ginga preferences used by PypeIt's viewer.

    PypeIt's viewer keeps its own ginga preferences (layout, plugin
    settings, etc.), separate from those of a standalone ginga (in
    ``~/.ginga``).  The directory is ``ginga`` in PypeIt's astropy
    configuration directory (see ``astropy.config.paths.get_config_dir``;
    e.g., ``~/.pypeit/config/ginga`` or ``~/.config/pypeit/ginga``) unless the
    ``PYPEIT_GINGA_HOME`` environment variable is defined.

    Returns
    -------
    `Path`_
        Directory with the ginga preferences.
    """
    if 'PYPEIT_GINGA_HOME' in os.environ:
        return Path(os.environ['PYPEIT_GINGA_HOME']).expanduser()
    from astropy.config.paths import get_config_dir
    return Path(get_config_dir('pypeit')) / 'ginga'


def _patch_ginga_hidpi_plots():
    """
    Work around a ginga bug that draws plots at half size on HiDPI displays.

    On HiDPI (e.g., Retina) displays, matplotlib multiplies the figure dpi by
    the device pixel ratio of the Qt canvas, but
    ``ginga.plot.PlotView.PlotViewBase.set_window_size`` (ginga 7.5.0) sets
    the figure size by dividing the window size in *logical* pixels by that
    dpi.  The figure is then half the size of the window, so plots (e.g., the
    ``Spec1dView`` spectrum) fill only a quarter of their panel.  This
    replaces the method with one that includes the device pixel ratio.

    .. todo::
        Remove this workaround once the bug is fixed in ginga.

    The patch is only applied if ginga's implementation has the expected
    internals; otherwise nothing is done.
    """
    try:
        from ginga.plot.PlotView import PlotViewBase
    except ImportError:
        return
    if getattr(PlotViewBase, '_pypeit_hidpi_patch', False):
        return
    if not all(hasattr(PlotViewBase, attr)
               for attr in ['set_window_size', 'configure', 'configure_window', 'get_figure']):
        return

    def set_window_size(self, wd_px, ht_px):
        with self._resize_lock:
            self.time_last_resize = time.time()
            fig = self.get_figure()
            dpr = getattr(fig.canvas, 'device_pixel_ratio', 1.) or 1.
            fig.set_size_inches(float(wd_px) * dpr / fig.dpi, float(ht_px) * dpr / fig.dpi)
            self.logger.debug("figure resized to %dx%d" % (wd_px, ht_px))
            self.make_callback('configure', wd_px, ht_px)

    PlotViewBase.set_window_size = set_window_size
    PlotViewBase.configure = set_window_size
    PlotViewBase.configure_window = set_window_size
    PlotViewBase._pypeit_hidpi_patch = True


def get_viewer():
    """
    Return the process-wide PypeIt viewer, creating it if necessary.

    Returns
    -------
    :class:`PypeItViewer`
        The viewer.
    """
    global _viewer
    if _viewer is None or _viewer.closed:
        _viewer = PypeItViewer()
    return _viewer


def viewer_exists():
    """
    Check if a PypeIt viewer has already been created in this process.

    Returns
    -------
    :obj:`bool`
        True if the viewer exists and has not been shut down.
    """
    return _viewer is not None and not _viewer.closed


def show(block=True, message=None):
    """
    Show the PypeIt viewer and, by default, wait for the user.

    This is a convenience wrapper for :meth:`PypeItViewer.show` that does
    nothing if nothing has been displayed in this process.

    Parameters
    ----------
    block : :obj:`bool`, optional
        If True, run the ginga event loop until the user closes the window.
        If False, only process pending GUI events and return.
    message : :obj:`str`, optional
        Message shown in the window title while waiting.
    """
    if not viewer_exists():
        return
    _viewer.show(block=block, message=message)


def embed(header='', local_ns=None):
    """
    Start an embedded IPython shell that keeps the PypeIt viewer responsive.

    While the IPython prompt waits for input, the ginga event loop keeps
    running, so the images can be inspected (panned, zoomed, etc.) at the same
    time as the variables of the calling code.

    Parameters
    ----------
    header : :obj:`str`, optional
        Header printed when the shell starts.
    local_ns : :obj:`dict`, optional
        Namespace for the shell.  If None, use the local namespace of the
        calling function.
    """
    import inspect
    from IPython.terminal import pt_inputhooks
    from IPython.terminal.embed import InteractiveShellEmbed

    if local_ns is None:
        local_ns = inspect.currentframe().f_back.f_locals
    if viewer_exists():
        def inputhook(context):
            while not context.input_is_ready() and not _viewer.closed:
                _viewer._update(0.01)
        pt_inputhooks.register('pypeit_ginga', inputhook)
    shell = InteractiveShellEmbed.instance()
    if viewer_exists():
        shell.enable_gui('pypeit_ginga')
    shell(header=header, local_ns=local_ns)


def _builtin_command(name):
    """
    Return the function that executes one of the commands that all PypeIt
    viewers accept from other processes.

    Parameters
    ----------
    name : :obj:`str`
        Command name: ``'spec2d'`` (show a spec2d file; see
        :meth:`~pypeit.display.spec2dview.Spec2dView.load_spec2d`) or
        ``'view_fits'`` (show an image; see
        :meth:`~pypeit.scripts.view_fits.ViewFits.show`).

    Returns
    -------
    callable
        The function, or None if the command is not known.
    """
    if name == 'spec2d':
        return lambda **kwargs: get_viewer().get_global_plugin('Spec2dView').load_spec2d(**kwargs)
    if name == 'view_fits':
        from pypeit.scripts.view_fits import ViewFits
        return lambda **kwargs: ViewFits.show(argparse.Namespace(**kwargs))
    return None


def is_listening(port=None):
    """
    Check if a PypeIt viewer is accepting commands on a local port.

    Parameters
    ----------
    port : :obj:`int`, optional
        Port to check.  If None, use :func:`get_port`.

    Returns
    -------
    :obj:`bool`
        True if a connection to the port could be made.
    """
    try:
        with socket.create_connection(('127.0.0.1', get_port() if port is None else port),
                                      timeout=1.):
            return True
    except OSError:
        return False


def send_command(cmd, port=None, timeout=10., **kwargs):
    """
    Send a command to a PypeIt viewer running in another process.

    Parameters
    ----------
    cmd : :obj:`str`
        Name of the command; must have been registered by the receiving
        viewer using :meth:`PypeItViewer.register_command`.
    port : :obj:`int`, optional
        Port to connect to.  If None, use :func:`get_port`.
    timeout : :obj:`float`, optional
        Time in seconds to wait for the receiving viewer to execute the
        command.
    **kwargs
        Keyword arguments passed to the command.  Must be JSON serializable.

    Returns
    -------
    :obj:`bool`
        True if a viewer accepted the command, False if no viewer is listening
        on the port.  The command is executed by the receiving viewer after it
        is accepted; any errors are reported by that viewer.

    Raises
    ------
    PypeItError
        Raised if the receiving viewer rejects the command.
    """
    _port = get_port() if port is None else port
    try:
        sock = socket.create_connection(('127.0.0.1', _port), timeout=1.)
    except OSError:
        return False
    with sock:
        sock.settimeout(timeout)
        try:
            sock.sendall((json.dumps(dict(cmd=cmd, kwargs=kwargs)) + '\n').encode())
            reply = b''
            while not reply.endswith(b'\n'):
                chunk = sock.recv(4096)
                if not chunk:
                    break
                reply += chunk
        except OSError:
            # Includes timeouts
            log.warning(f'The PypeIt viewer on port {_port} did not respond.')
            return False
    try:
        reply = json.loads(reply.decode())
    except json.JSONDecodeError:
        raise PypeItError(f'Unexpected reply from the PypeIt viewer on port {_port}.')
    if reply.get('status') != 'ok':
        raise PypeItError(f'PypeIt viewer rejected {cmd}: {reply.get("error")}')
    return True


class GingaShellViewer:
    """
    Thin wrapper around a ginga reference-viewer shell used by the PypeIt
    display functions (see :mod:`pypeit.display.display`).

    The wrapper can be used with any ginga shell, including that of a
    standalone ginga that has loaded the PypeIt plugins.  PypeIt's own
    in-process viewer is :class:`PypeItViewer`.

    Parameters
    ----------
    shell : ``ginga.rv.Control.GingaShell``
        The ginga reference-viewer shell.
    """
    def __init__(self, shell):
        self.shell = shell
        #: Set when something has been displayed but not yet shown; see
        #: :meth:`PypeItViewer.show`.
        self.pending = False
        self.dc = None
        # Channels already fit to the window; see fit_channels
        self._fit_channels = set()

    def _update(self, duration=0.):
        """Process pending GUI events; nothing to do in the base class."""
        pass

    def fit_channels(self, chnames=None, fit=True):
        """
        Fit the images of channels to the window.

        Images loaded before the window has its final size are not fit to the
        window, so this should be done once the window is shown.  Each
        channel is fit only once, so that the pan and zoom chosen by the user
        are kept.  Channels with linked pan and zoom (see
        :meth:`set_wcs_match`) are fit together.

        Parameters
        ----------
        chnames : :obj:`list`, optional
            Names of the channels to fit.  If None, fit all channels that have
            not been fit before.
        fit : :obj:`bool`, optional
            If False, only record that the channels do not need to be fit
            (e.g., because they were linked to channels already shown).
        """
        names = [name for name in (self.channel_names() if chnames is None else chnames)
                    if name not in self._fit_channels]
        if len(names) == 0:
            return
        if not fit:
            self._fit_channels.update(names)
            return
        self._update(0.05)
        fitted = []
        for name in names:
            fi = self.channel(name).fitsimage
            self._fit_channels.add(name)
            # Fitting one of a set of linked channels fits all of them
            scale = fi.get_settings().group['scale']
            if fi.get_image() is None or any(scale is s for s in fitted):
                continue
            fi.zoom_fit()
            fitted.append(scale)

    # ------------------------------------------------------------------
    # Channels and plugins

    def canvas_types(self):
        """
        Return the ginga canvas types (see ``ginga.canvas.CanvasObject``).

        Returns
        -------
        ``ginga.misc.Bunch.Bunch``
            The canvas types; e.g., ``dc.Path``, ``dc.Text``.
        """
        if self.dc is None:
            from ginga.canvas.CanvasObject import get_canvas_types
            self.dc = get_canvas_types()
        return self.dc

    def channel(self, chname):
        """
        Return a ginga channel, creating it if it does not exist.

        Parameters
        ----------
        chname : :obj:`str`
            Channel name.

        Returns
        -------
        ``ginga.rv.Channel.Channel``
            The channel.
        """
        return self.shell.get_channel_on_demand(chname)

    def channel_names(self):
        """
        Return the names of the existing channels.

        Returns
        -------
        :obj:`list`
            Channel names.
        """
        return list(self.shell.get_channel_names())

    def delete_channel(self, chname):
        """
        Delete a channel, if it exists.

        Parameters
        ----------
        chname : :obj:`str`
            Channel name.
        """
        if self.shell.has_channel(chname):
            self.shell.delete_channel(chname)
        self._fit_channels.discard(chname)

    def clear_all(self):
        """Delete all channels."""
        for chname in self.channel_names():
            self.shell.delete_channel(chname)
        self._fit_channels.clear()
        self._update()

    def get_global_plugin(self, name):
        """
        Start (if necessary) and return a global ginga plugin.

        Parameters
        ----------
        name : :obj:`str`
            Plugin name.

        Returns
        -------
        ``ginga.GingaPlugin.GlobalPlugin``
            The plugin object, with its GUI built.
        """
        gpmon = self.shell.gpmon
        if not gpmon.is_active(name):
            self.shell.start_global_plugin(name)
            self._wait_for_plugin(name)
        return gpmon.get_plugin(name)

    def _wait_for_plugin(self, name):
        """Wait for a global plugin to start; nothing to do in the base class."""
        pass

    #: Viewer settings shared by channels linked with :meth:`set_wcs_match`
    linked_settings = ('scale', 'pan', 'flip_x', 'flip_y', 'swap_xy', 'rot_deg')

    def set_wcs_match(self, chname):
        """
        Match the pan, zoom, and orientation of all channels to a reference.

        PypeIt images do not have a sky WCS, so the channels are matched in
        pixel coordinates: the channel viewers share the same pan, scale,
        flip, swap, and rotation settings.  Panning or zooming in any one of
        the channels changes all of them.

        Parameters
        ----------
        chname : :obj:`str`
            Name of the reference channel.
        """
        ref = self.channel(chname).fitsimage.get_settings()
        # Hold the redraws of all channels until all settings are shared, so
        # that each channel is only redrawn once.
        with contextlib.ExitStack() as stack:
            for name in self.channel_names():
                stack.enter_context(self.channel(name).fitsimage.suppress_redraw)
            for name in self.channel_names():
                if name == chname:
                    continue
                other = self.channel(name).fitsimage.get_settings()
                keys = [key for key in self.linked_settings
                            if other.group.get(key) is not ref.group[key]]
                if len(keys) > 0:
                    ref.share_settings(other, keylist=keys)
        self._update()


class PypeItViewer(GingaShellViewer):
    """
    Process-wide, in-process ginga reference viewer.

    Do not instantiate this class directly; use :func:`get_viewer`.

    The ginga event loop only runs while :meth:`show` waits for the user (and
    briefly while displays are being constructed).
    """

    def __init__(self):
        # NOTE: ginga must choose its widget toolkit before any of its widget
        # or canvas-type modules are imported, so all ginga imports that touch
        # the GUI happen here, not at the top of the module.
        from ginga.rv.main import ReferenceViewer, default_layout, plugins

        from pypeit.display import pypeit_plugin_specs

        class _ReferenceViewer(ReferenceViewer):
            def setup(self):
                # The in-process viewer does not need ginga's HTTP server
                self.settings.set(http_server=False)
                return super().setup()

            def run(self):
                # Do not enter ginga's own (blocking) event loop or call
                # sys.exit when it ends; PypeIt drives the loop in show().
                pass

        self.closed = False
        self._continue = False
        self._server = None
        self._commands = {}

        t = time.perf_counter()
        basedir = get_ginga_home()
        basedir.mkdir(parents=True, exist_ok=True)
        self.rv = _ReferenceViewer(layout=default_layout, plugins=plugins, appname='ginga',
                                   basedir=str(basedir), channels=[])
        # The remote-control and SAMP plugins are not needed
        self.rv.add_default_plugins(except_global=['RC', 'SAMP'])
        for spec in pypeit_plugin_specs():
            self.rv.add_plugin_spec(spec)
        parser = argparse.ArgumentParser()
        self.rv.add_default_options(parser)
        self.rv.main(parser.parse_args(['--nosplash', f'--basedir={basedir}']), [])

        super().__init__(self.rv.ginga_shell)
        # TODO: Remove when ginga's HiDPI plot bug is fixed
        _patch_ginga_hidpi_plots()
        self._set_title()
        # PypeIt drives the event loop from this thread
        self.shell.gui_thread_id = threading.get_ident()

        # Closing the window only hides it and returns control to the caller
        for win in self.shell.ds.toplevels:
            win.clear_callback('close')
            win.add_callback('close', self._window_close_cb)

        self._update(0.2)
        log.info(f'Started the ginga viewer in {time.perf_counter() - t:.1f} s')
        # NOTE: The ginga threads are not daemon threads, and python waits for
        # all such threads to finish *before* calling the functions registered
        # with atexit.  Use the (private) threading hook, which is called
        # before waiting for the threads, so that the viewer is shut down and
        # the program can exit.
        if hasattr(threading, '_register_atexit'):
            threading._register_atexit(self._atexit)
        else:
            atexit.register(self._atexit)

    # ------------------------------------------------------------------
    # Event loop

    def _update(self, duration=0.):
        """Process pending GUI events for (at least) the given duration."""
        tstart = time.perf_counter()
        while True:
            self.shell.update_pending(timeout=0.001)
            if time.perf_counter() - tstart >= duration:
                break

    def _wait_for_plugin(self, name):
        """Process GUI events until a global plugin has started."""
        tstart = time.perf_counter()
        while not self.shell.gpmon.is_active(name):
            self.shell.update_pending(timeout=0.001)
            if time.perf_counter() - tstart > 10.:
                raise PypeItError(f'Timeout starting ginga plugin {name}.')

    def _window_close_cb(self, w):
        """Window-close callback: hide the window and stop waiting."""
        w.hide()
        self._continue = True

    def continue_(self):
        """Stop waiting in :meth:`show`, as if the window had been closed."""
        self._continue = True

    #: Title of the viewer window
    title = 'PypeIt'

    def _set_title(self, message=None):
        """Set the window title, including a message if provided."""
        title = self.title if message is None else f'{self.title}: {message}'
        for win in self.shell.ds.toplevels:
            win.set_title(title)

    def show(self, block=True, message=None):
        """
        Raise the viewer window and, by default, wait for the user to close it.

        Closing the window only hides it; it is shown again, with the same
        channels, by the next call.

        Parameters
        ----------
        block : :obj:`bool`, optional
            If True, run the ginga event loop until the user closes the
            window.  If False, only process pending GUI events and return.
        message : :obj:`str`, optional
            Message shown in the window title (and logged) while waiting.  If
            None, the message is "close the window to continue".
        """
        if self.closed:
            raise PypeItError('The PypeIt viewer has been shut down.')
        for win in self.shell.ds.toplevels:
            win.show()
            win.raise_()
        self.pending = False
        # Fit the images of new channels to the window
        self.fit_channels()
        if not block:
            self._update()
            return

        self._continue = False
        if message is None:
            message = 'close the window to continue'
        self._set_title(message)
        log.info(f'PypeIt viewer: {message}.')
        try:
            next_poll = 0.
            while not self._continue and not self.shell.ev_quit.is_set():
                self.shell.update_pending(timeout=0.001)
                if self._server is not None and time.perf_counter() > next_poll:
                    self._poll_server()
                    next_poll = time.perf_counter() + 0.05
        finally:
            # Anything displayed while waiting (e.g., files sent by other
            # processes) has been seen by the user.
            self.pending = False
            if not self.shell.ev_quit.is_set():
                self._set_title()
                self._update()

    def _atexit(self):
        """
        At exit, show anything displayed but never shown, then shut down.

        This makes sure that a display is not silently lost if the calling
        code does not call :meth:`show`.
        """
        # NOTE: sys.last_value is set if the program is ending because of an
        # uncaught exception; never wait for the user in that case.
        if not self.closed and self.pending and not hasattr(sys, 'last_value'):
            try:
                self.show(message='close the window to exit')
            except Exception:
                pass
        self.shutdown()

    def shutdown(self):
        """Stop the viewer and release its resources."""
        if self.closed:
            return
        self.closed = True
        if self._server is not None:
            self._server.close()
            self._server = None
        try:
            self.shell.stop()
        except Exception:
            pass
        listener = getattr(self.rv, '_gui_log_listener', None)
        if listener is not None:
            listener.stop()
            self.rv._gui_log_listener = None

    # ------------------------------------------------------------------
    # Commands from other processes

    def register_command(self, name, func):
        """
        Register a function that other processes can call via
        :func:`send_command`.

        All viewers accept the commands listed in :func:`_builtin_command`;
        use this to add other commands.

        Parameters
        ----------
        name : :obj:`str`
            Name of the command.
        func : callable
            Function called with the keyword arguments sent by the client.
        """
        self._commands[name] = func

    def listen(self, port=None):
        """
        Accept commands from other PypeIt processes on a local port.

        Only connections from the local host are accepted, and only commands
        registered with :meth:`register_command` are executed.  Commands are
        processed while the viewer is waiting in :meth:`show`.

        Parameters
        ----------
        port : :obj:`int`, optional
            Port to listen on.  If None, use :func:`get_port`.

        Returns
        -------
        :obj:`bool`
            True if the viewer is listening; False if the port is already in
            use (e.g., by another viewer).
        """
        if self._server is not None:
            return True
        _port = get_port() if port is None else port
        server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        try:
            server.bind(('127.0.0.1', _port))
        except OSError:
            server.close()
            log.warning(f'Port {_port} is in use (perhaps by another PypeIt viewer); this viewer '
                        'will not accept files sent with --no_clear.')
            return False
        server.listen(4)
        server.setblocking(False)
        self._server = server
        return True

    def _poll_server(self):
        """Execute any commands sent by other processes."""
        while True:
            try:
                ready, _, _ = select.select([self._server], [], [], 0)
            except (OSError, ValueError):
                return
            if len(ready) == 0:
                return
            try:
                conn, _ = self._server.accept()
            except BlockingIOError:
                return
            with conn:
                conn.settimeout(2.)
                try:
                    self._handle_connection(conn)
                except Exception as e:
                    log.warning(f'Error handling command sent to the PypeIt viewer: {e}')

    def _handle_connection(self, conn):
        """
        Read and execute one command.

        The reply is sent as soon as the command is accepted, so that the
        sending process does not wait for the command to be executed.  Errors
        executing the command are reported by this viewer.
        """
        buf = b''
        while not buf.endswith(b'\n'):
            chunk = conn.recv(4096)
            if not chunk:
                break
            buf += chunk
        if len(buf) == 0:
            # Connection only used to check if the viewer is listening
            return
        try:
            request = json.loads(buf.decode())
            func = self._commands.get(request['cmd'], _builtin_command(request['cmd']))
            kwargs = request.get('kwargs', {})
            if func is None or not isinstance(kwargs, dict):
                raise TypeError
        except (json.JSONDecodeError, KeyError, TypeError):
            conn.sendall((json.dumps(dict(status='error', error='Unknown or malformed command.'))
                          + '\n').encode())
            return
        conn.sendall((json.dumps(dict(status='ok')) + '\n').encode())
        conn.close()

        try:
            func(**kwargs)
        except Exception as e:
            msg = f'Error executing {request["cmd"]} sent by another process: {e}'
            log.error(msg)
            self.shell.show_error(msg)
        for win in self.shell.ds.toplevels:
            win.show()
            win.raise_()
