"""
PypeItRemote is a hidden ginga global plugin that runs PypeIt display commands
inside the ginga reference viewer.

Plugin Type: Global
===================

PypeItRemote has no user interface.  It is loaded whenever ginga starts in an
environment where PypeIt is installed, and is driven over ginga's remote
control (RC) interface by the PypeIt display scripts (e.g.,
``pypeit_show_2dspec``).

Usage
-----

A script sends a command name and its parsed command-line arguments with
:meth:`PypeItRemote.submit`, which returns at once with a job ID.  The command
then runs on a ginga worker thread, so the viewer stays responsive while files
are read.  The script polls :meth:`PypeItRemote.status` until the job is done
and reports any warnings or errors in the terminal.

Only one command runs at a time.  Submitting a new command cancels any
command that is still queued or running ("last command wins").

The commands that can be run are defined in
:func:`pypeit.display.remote_handlers.get_handlers`.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
from ginga import GingaPlugin

from pypeit import __version__
from pypeit.display import remote_handlers
from pypeit.display import remote_jobs

__all__ = ['PypeItRemote']


class PypeItRemote(GingaPlugin.GlobalPlugin):
    """
    Hidden ginga global plugin that runs PypeIt display commands.

    The public methods are called over ginga's RC interface with
    ``call_global_plugin_method``; see :mod:`pypeit.display.remote` for the
    client side.

    Parameters
    ----------
    fv : ``ginga.rv.Control.GingaShell``
        The ginga reference viewer shell.
    """
    def __init__(self, fv):
        super().__init__(fv)
        self.jobs = remote_jobs.JobManager(fv, remote_handlers.get_handlers(), logger=self.logger)

    def ping(self):
        """
        Return the PypeIt version loaded in ginga.

        Returns
        -------
        :obj:`str`
            PypeIt version.
        """
        return __version__

    def submit(self, command, args_json, version):
        """
        Queue a command; see :meth:`~pypeit.display.remote_jobs.JobManager.submit`.
        """
        return self.jobs.submit(command, args_json, version)

    def status(self, job_id):
        """
        Return the status of a job; see :meth:`~pypeit.display.remote_jobs.JobManager.status`.
        """
        return self.jobs.status(job_id)

    def cancel(self, job_id):
        """
        Ask a job to stop; see :meth:`~pypeit.display.remote_jobs.JobManager.cancel`.
        """
        return self.jobs.cancel(job_id)

    def __str__(self):
        return 'pypeitremote'
