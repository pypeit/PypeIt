"""
Client side of the :class:`~pypeit.display.pypeit_remote.PypeItRemote` ginga
plugin.

The PypeIt display scripts use :func:`send_command` to run their work inside a
long-lived ginga viewer instead of in the script process.  The script parses
its command line, sends the parsed arguments to ``PypeItRemote`` over ginga's
remote control (RC) interface, and waits for the command to finish.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import json
from pathlib import Path
import time
import xmlrpc.client

from ginga.util import grc

from pypeit import __version__
from pypeit import log
from pypeit import PypeItError
from pypeit.display import display


def add_timeout_argument(parser):
    """
    Add the ``--timeout`` option shared by the scripts that use
    :func:`send_command`.

    Parameters
    ----------
    parser : `argparse.ArgumentParser`_
        Parser to add the option to.
    """
    parser.add_argument(
        '--timeout', type=float, default=120.,
        help='Number of seconds to wait for ginga to display the data.  If the display is '
             'not done by then, the script stops waiting, but ginga finishes the display.'
    )


def args_to_json(args, path_args=None):
    """
    Encode parsed command-line arguments for :func:`send_command`.

    Parameters
    ----------
    args : `argparse.Namespace`_
        Parsed command-line arguments.  All values must be JSON-serializable
        after the path arguments are converted to strings.
    path_args : :obj:`list`, optional
        Names of the arguments that are file or directory paths.  Each is
        made absolute relative to the current working directory, because the
        ginga process may have a different working directory.  Arguments that
        are None are left as None.

    Returns
    -------
    :obj:`str`
        JSON-encoded dictionary of the arguments.
    """
    _args = dict(vars(args))
    if path_args is not None:
        for key in path_args:
            if _args[key] is not None:
                _args[key] = str(Path(_args[key]).absolute())
    return json.dumps(_args)


class _TimeoutTransport(xmlrpc.client.Transport):
    """
    XML-RPC transport with a socket timeout.

    Parameters
    ----------
    timeout : :obj:`float`
        Socket timeout in seconds.
    """
    def __init__(self, timeout):
        super().__init__()
        self._timeout = timeout

    def make_connection(self, host):
        """
        Return the HTTP connection to the server, with the timeout set.

        Parameters
        ----------
        host : :obj:`str`
            Server host.

        Returns
        -------
        `http.client.HTTPConnection`_
            The connection.
        """
        connection = super().make_connection(host)
        connection.timeout = self._timeout
        return connection


class RemoteViewer:
    """
    Minimal RC client for the ``PypeItRemote`` plugin.

    Unlike `ginga.util.grc.RemoteClient`, each call has a socket timeout, so a
    viewer that has stopped responding cannot hang the calling script.

    Parameters
    ----------
    host : :obj:`str`, optional
        Host running ginga.
    port : :obj:`int`, optional
        Port of ginga's RC server.
    rpc_timeout : :obj:`float`, optional
        Socket timeout in seconds for each call.
    """
    def __init__(self, host='localhost', port=grc.default_rc_port, rpc_timeout=10.):
        self._proxy = xmlrpc.client.ServerProxy(
            f'http://{host}:{port}', allow_none=True, transport=_TimeoutTransport(rpc_timeout)
        )

    def call(self, method_name, *args):
        """
        Call a method of the ``PypeItRemote`` plugin.

        Parameters
        ----------
        method_name : :obj:`str`
            Name of the method.
        *args :
            Positional arguments for the method.  Must be strings, numbers,
            booleans, lists or dictionaries of these (not None).

        Returns
        -------
        object
            The value returned by the method.
        """
        p_args = grc.marshall(
            ('call_global_plugin_method', 'PypeItRemote', method_name, list(args), {})
        )
        return grc.unmarshall(self._proxy.dispatch_call('ginga', p_args, grc.marshall({})))


def connect(viewer=None):
    """
    Connect to a ginga viewer running the ``PypeItRemote`` plugin.

    A new ginga viewer is started if none is running.

    Parameters
    ----------
    viewer : :class:`RemoteViewer`, optional
        Viewer to check.  If None, a new connection to the default host and
        port is made.

    Returns
    -------
    :class:`RemoteViewer`
        The connection.

    Raises
    ------
    :class:`~pypeit.PypeItError`
        Raised if ginga cannot be reached or does not have the
        ``PypeItRemote`` plugin.
    """
    if viewer is None:
        display.connect_to_ginga(raise_err=True, allow_new=True)
        viewer = RemoteViewer()
    try:
        viewer.call('ping')
    except xmlrpc.client.Fault as e:
        raise PypeItError(
            'The running ginga viewer does not have the PypeItRemote plugin.  Close ginga and '
            'restart it with "ginga --modules=RC" from the environment where PypeIt is '
            'installed.'
        ) from e
    except OSError as e:
        raise PypeItError(f'Could not communicate with the ginga viewer: {e}') from e
    return viewer


def send_command(command, args_json, timeout=120., poll_interval=0.25, viewer=None):
    """
    Run a PypeIt display command in ginga and wait for it to finish.

    Warnings raised while the command runs are logged here, and an error in
    the command is raised here, so that both appear in the user's terminal.

    Parameters
    ----------
    command : :obj:`str`
        Name of the command; see
        :func:`pypeit.display.remote_handlers.get_handlers`.
    args_json : :obj:`str`
        JSON-encoded command arguments; see :func:`args_to_json`.
    timeout : :obj:`float`, optional
        Number of seconds to wait for the command to finish.  If it is still
        running after this time, a warning is logged and the function returns;
        the command keeps running in ginga.
    poll_interval : :obj:`float`, optional
        Number of seconds between status checks.
    viewer : :class:`RemoteViewer`, optional
        Connection to use.  If None, :func:`connect` is called.

    Returns
    -------
    :obj:`dict`
        Last status of the job; see
        :meth:`pypeit.display.remote_jobs.RemoteJob.status`.

    Raises
    ------
    :class:`~pypeit.PypeItError`
        Raised if the command fails, or if ginga cannot be reached.
    """
    viewer = connect(viewer=viewer)
    status = None
    try:
        status = viewer.call('submit', command, args_json, __version__)
        deadline = time.monotonic() + timeout
        while status['state'] in ('queued', 'running'):
            if time.monotonic() > deadline:
                break
            time.sleep(poll_interval)
            status = viewer.call('status', status['id'])
    except KeyboardInterrupt:
        if status is not None:
            viewer.call('cancel', status['id'])
        raise
    except OSError as e:
        raise PypeItError(f'Lost communication with the ginga viewer: {e}') from e

    for warning in status['warnings']:
        log.warning(warning)
    if status['state'] == 'error':
        raise PypeItError(status['message'])
    if status['state'] == 'cancelled':
        log.warning(f'{command} was cancelled in ginga: {status["message"]}')
    elif status['state'] != 'done':
        log.warning(
            f'{command} is still running in ginga after {timeout} s; it will finish there.  '
            'Use --timeout to wait longer.'
        )
    return status
