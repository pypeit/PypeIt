"""
Job handling for the :class:`~pypeit.display.pypeit_remote.PypeItRemote` ginga
plugin.

This module holds the logic that queues and runs PypeIt display commands on
ginga worker threads.  It is kept separate from the plugin class so that it
does not import ginga's widget toolkit and can be tested without a GUI.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import argparse
import json
import logging
import threading

from pypeit import __version__
from pypeit import log
from pypeit import PypeItError


class JobCancelled(Exception):
    """
    Raised inside a command handler when its job has been cancelled.
    """
    pass


class RemoteJob:
    """
    State of a single command submitted to :class:`PypeItRemote`.

    A job moves through the states ``'queued'``, ``'running'`` and then one of
    ``'done'``, ``'error'`` or ``'cancelled'``.

    Parameters
    ----------
    job_id : :obj:`int`
        Unique ID of the job.
    command : :obj:`str`
        Name of the command to run.
    args : `argparse.Namespace`_
        Parsed command-line arguments passed to the command handler.
    """
    finished_states = ('done', 'error', 'cancelled')
    """
    States of a job that is no longer running.
    """

    def __init__(self, job_id, command, args):
        self.id = job_id
        self.command = command
        self.args = args
        self.state = 'queued'
        self.message = ''
        self.warnings = []
        self._cancel_event = threading.Event()
        self._lock = threading.Lock()

    @property
    def cancelled(self):
        """
        True if the job has been asked to stop.
        """
        return self._cancel_event.is_set()

    @property
    def finished(self):
        """
        True if the job is no longer queued or running.
        """
        return self.state in self.finished_states

    def cancel(self):
        """
        Ask the job to stop.

        Cancellation is cooperative: a running handler stops the next time it
        calls :meth:`check_cancel`.
        """
        self._cancel_event.set()

    def check_cancel(self):
        """
        Raise :class:`JobCancelled` if the job has been asked to stop.

        Command handlers call this between steps, e.g., before displaying each
        channel.
        """
        if self.cancelled:
            raise JobCancelled(f'Command {self.command} (job {self.id}) was cancelled.')

    def warn(self, message):
        """
        Record a warning to report to the script that submitted the job.

        Parameters
        ----------
        message : :obj:`str`
            Warning message.
        """
        with self._lock:
            self.warnings.append(str(message))

    def set_state(self, state, message=''):
        """
        Set the state of the job.

        Parameters
        ----------
        state : :obj:`str`
            New state.
        message : :obj:`str`, optional
            Message describing the state, e.g. an error message.
        """
        with self._lock:
            self.state = state
            self.message = str(message)

    def status(self):
        """
        Return the job status in a form that can be sent over RC.

        Returns
        -------
        :obj:`dict`
            Dictionary with the job ``id``, ``command``, ``state``,
            ``message`` and ``warnings``.
        """
        with self._lock:
            return dict(
                id=self.id, command=self.command, state=self.state, message=self.message,
                warnings=list(self.warnings)
            )


class _JobLogHandler(logging.Handler):
    """
    Collect the warnings logged by the thread running a job.

    Parameters
    ----------
    job : :class:`RemoteJob`
        Job that receives the warnings.
    thread_id : :obj:`int`
        Identifier of the thread running the job; records logged by other
        threads are ignored.
    """
    def __init__(self, job, thread_id):
        super().__init__(level=logging.WARNING)
        self.job = job
        self.thread_id = thread_id

    def emit(self, record):
        """
        Add a log record to the job warnings if it came from the job thread.

        Parameters
        ----------
        record : `logging.LogRecord`_
            The record to handle.
        """
        if record.thread == self.thread_id:
            self.job.warn(record.getMessage())


class JobManager:
    """
    Queue and run PypeIt display commands on ginga worker threads.

    This holds all of the logic of :class:`PypeItRemote`, separated from the
    plugin so that it can be tested without a ginga GUI.

    Parameters
    ----------
    fv : ``ginga.rv.Control.GingaShell``
        The ginga reference viewer shell.  Only ``nongui_do`` is used
        directly; it is also passed to the command handlers.
    handlers : :obj:`dict`
        Command handlers, keyed by command name.  Each handler is called as
        ``handler(fv, args, job)``; see
        :func:`pypeit.display.remote_handlers.get_handlers`.
    logger : `logging.Logger`_, optional
        Logger for errors.  If None, the PypeIt logger is used.
    max_jobs : :obj:`int`, optional
        Maximum number of jobs to keep.  The oldest finished jobs are
        forgotten beyond this number.
    """
    def __init__(self, fv, handlers, logger=None, max_jobs=20):
        self.fv = fv
        self.handlers = handlers
        self.logger = log if logger is None else logger
        self.max_jobs = max_jobs
        self.jobs = {}
        self._next_id = 1
        self._lock = threading.Lock()
        self._run_lock = threading.Lock()

    def submit(self, command, args_json, version):
        """
        Queue a command and return at once.

        Any job that is still queued or running is cancelled.

        Parameters
        ----------
        command : :obj:`str`
            Name of the command to run.
        args_json : :obj:`str`
            JSON-encoded dictionary of the parsed command-line arguments.
        version : :obj:`str`
            PypeIt version of the script submitting the command.

        Returns
        -------
        :obj:`dict`
            Status of the new job; see :meth:`RemoteJob.status`.  If the
            command is unknown or the arguments cannot be decoded, the state
            is ``'error'``.
        """
        with self._lock:
            job_id = self._next_id
            self._next_id += 1
            for job in self.jobs.values():
                if not job.finished:
                    job.cancel()

        try:
            args = argparse.Namespace(**json.loads(args_json))
        except (TypeError, ValueError) as e:
            job = RemoteJob(job_id, command, None)
            job.set_state('error', message=f'Could not decode the arguments for {command}: {e}')
            return self._add_job(job).status()

        job = self._add_job(RemoteJob(job_id, command, args))
        if version != __version__:
            job.warn(
                f'Ginga is running PypeIt {__version__}, but the script is from PypeIt '
                f'{version}.  Restart ginga to use the same version.'
            )
        if command not in self.handlers:
            job.set_state(
                'error',
                message=f'Unknown command: {command}.  Available commands are: '
                        f'{", ".join(sorted(self.handlers.keys()))}.'
            )
            return job.status()

        self.fv.nongui_do(self.run_job, job)
        return job.status()

    def _add_job(self, job):
        """
        Add a job to the table, forgetting the oldest finished jobs if needed.

        Parameters
        ----------
        job : :class:`RemoteJob`
            Job to add.

        Returns
        -------
        :class:`RemoteJob`
            The job that was added.
        """
        with self._lock:
            self.jobs[job.id] = job
            finished = [j.id for j in self.jobs.values() if j.finished and j.id != job.id]
            nremove = len(self.jobs) - self.max_jobs
            for job_id in finished[:max(nremove, 0)]:
                del self.jobs[job_id]
        return job

    def run_job(self, job):
        """
        Run a job's command handler.

        This is called on a ginga worker thread.  Jobs run one at a time; a
        job that is cancelled while it waits for the previous one is not run.

        Parameters
        ----------
        job : :class:`RemoteJob`
            Job to run.
        """
        with self._run_lock:
            if job.cancelled:
                job.set_state('cancelled', message='Superseded by a newer command.')
                return
            job.set_state('running')
            log_handler = _JobLogHandler(job, threading.get_ident())
            log.addHandler(log_handler)
            try:
                self.handlers[job.command](self.fv, job.args, job)
            except JobCancelled:
                job.set_state('cancelled', message='Superseded by a newer command.')
            except PypeItError as e:
                job.set_state('error', message=str(e))
            except Exception as e:
                self.logger.error(f'Error running {job.command}: {e}', exc_info=True)
                job.set_state('error', message=f'{type(e).__name__}: {e}')
            else:
                job.set_state('done')
            finally:
                log.removeHandler(log_handler)

    def status(self, job_id):
        """
        Return the status of a job.

        Parameters
        ----------
        job_id : :obj:`int`
            ID of the job.

        Returns
        -------
        :obj:`dict`
            Status of the job; see :meth:`RemoteJob.status`.  If the job is
            not known, the state is ``'error'``.
        """
        with self._lock:
            job = self.jobs.get(job_id)
        if job is None:
            return dict(
                id=job_id, command='', state='error', message=f'Unknown job ID: {job_id}.',
                warnings=[]
            )
        return job.status()

    def cancel(self, job_id):
        """
        Ask a job to stop.

        Parameters
        ----------
        job_id : :obj:`int`
            ID of the job.

        Returns
        -------
        :obj:`dict`
            Status of the job; see :meth:`status`.
        """
        with self._lock:
            job = self.jobs.get(job_id)
        if job is not None:
            job.cancel()
        return self.status(job_id)
