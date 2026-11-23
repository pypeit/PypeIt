"""
Command handlers run by the :class:`~pypeit.display.pypeit_remote.PypeItRemote`
ginga plugin.

Each handler runs one PypeIt display command inside the ginga process.  It is
called on a ginga worker thread as ``handler(fv, args, job)``, where ``fv`` is
the ginga reference viewer shell, ``args`` is the `argparse.Namespace`_ parsed
by the script that sent the command, and ``job`` is the
:class:`~pypeit.display.remote_jobs.RemoteJob` running it.

A handler must:

    - do file I/O and computation directly, on the worker thread;
    - make every change to the viewer through ``fv.gui_call``;
    - call ``job.check_cancel()`` between steps so that a newer command can
      supersede it;
    - raise :class:`~pypeit.PypeItError` for errors the user should see in the
      terminal.

Warnings logged with the PypeIt logger by the worker thread are reported to
the script that sent the command.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""


def get_handlers():
    """
    Return the command handlers known to the ``PypeItRemote`` plugin.

    The command names match the names used by the scripts that send them
    (e.g., ``'show_2dspec'`` for ``pypeit_show_2dspec``).

    Returns
    -------
    :obj:`dict`
        Handler functions keyed by command name.
    """
    return {}
