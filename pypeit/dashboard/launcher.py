"""
Launch PypeIt's external tools as subprocesses.

Uses Qt's :obj:`QProcess` so launches are **non-blocking** and integrate
with the event loop: stdout/stderr/return code are captured and the
:class:`~pypeit.dashboard.view.activity.ActivityBar` is updated on start and
finish, so every launch is observable in the GUI.
"""

from pathlib import Path

from qtpy.QtCore import QProcess

from pypeit import log


class Launcher:
    """
    Launches external tools as non-blocking subprocesses and reports
    progress/outcome to the activity bar.

    Two kinds of subprocess are distinguished:

    * :meth:`launch` — a fire-and-forget **inspection viewer** (e.g.
      ``ginga``, ``pypeit_chk_*``).  It does not touch the run lock, and
      its progress is reported on the activity bar's Inspection channel.
    * :meth:`run` — a **PypeIt reduction run** (e.g. ``run_pypeit``,
      ``pypeit_run_to_calibstep``).  It holds the single-run lock while
      active (so a second run cannot be launched), reports on the Build
      channel, and invokes a completion callback so the views can refresh
      the reduction state.

    Parameters
    ----------
    activity : :class:`~pypeit.dashboard.view.activity.ActivityBar`, optional
        The status/activity bar to report to (may be ``None``).
    run_lock : :class:`~pypeit.dashboard.runlock.RunLock`, optional
        The single-run lock; a run launched via :meth:`run` holds it while
        active.  May be ``None``.
    """

    def __init__(self, activity=None, run_lock=None):
        self.activity = activity
        # The single-run lock; set so a launched run holds the lock while
        # it is active.  May be None (e.g. headless tests).
        self.run_lock = run_lock
        # Keep references to running processes so they are not garbage
        # collected while still running.
        self._procs = []

    def _report(self, message, busy=False, build=True):
        """
        Report a message to the activity bar (and the log), on the **Build**
        channel (reduction runs) or the **Inspection** channel (viewer
        launches) per ``build``.

        Parameters
        ----------
        message : :obj:`str`
            The message.
        busy : :obj:`bool`, optional
            Show the busy indicator.
        build : :obj:`bool`, optional
            Route to the Build channel (True, default) or the Inspection
            channel (False).
        """
        log.info(message)
        if self.activity is not None:
            if build:
                self.activity.set_build(message, busy=busy)
            else:
                self.activity.set_inspection(message, busy=busy)

    def launch(self, argv, description=None, hint='viewer window'):
        """
        Launch an **inspection viewer** as a fire-and-forget subprocess.

        Unlike :meth:`run`, this never engages the single-run lock: any
        number of viewers may be open at once, and their outcome is
        reported on the activity bar's Inspection channel.

        Parameters
        ----------
        argv : :obj:`list`
            The command and its arguments.
        description : :obj:`str`, optional
            Human description for the activity bar (defaults to the
            program name).
        hint : :obj:`str`, optional
            Where the result appears (e.g. ``'Ginga window'``); shown in
            the finished message so the user knows where to look.

        Returns
        -------
        :obj:`bool`
            True if the process was started, False if it could not be
            launched (e.g. the target file is missing).
        """
        if not argv:
            return False
        program, args = argv[0], [str(a) for a in argv[1:]]
        # The exact command line, shown in quotes so the user can copy and
        # run it themselves.
        cmd = '"' + ' '.join([program] + args) + '"'

        # If the command references a file that does not exist, report and
        # bail rather than spawning a doomed process.
        for a in args:
            if a.startswith('-'):
                continue
            if ('/' in a or a.endswith('.fits') or a.endswith('.gz')) \
                    and not Path(a).exists():
                self._report(f'Cannot launch {cmd}: file not found ({a})',
                             busy=False, build=False)
                return False

        proc = QProcess()
        proc.setProcessChannelMode(QProcess.MergedChannels)
        self._procs.append(proc)

        # Show the running command (with the busy indicator) once the
        # process has actually started.
        proc.started.connect(
            lambda: self._report(f'Running {cmd}', busy=True, build=False))
        # Report a process that failed to start (e.g. command not found).
        proc.errorOccurred.connect(
            lambda _err, p=proc, c=cmd: self._on_error(p, c))
        # Report the outcome (exit code + captured output) when it ends.
        proc.finished.connect(
            lambda code, _status, p=proc, c=cmd, h=hint:
            self._on_finished(p, c, code, h))

        proc.start(program, args)
        return True

    def run(self, argv, description='Rebuilding', on_finished=None):
        """
        Launch a PypeIt **reduction run** (e.g. ``pypeit_run_to_calibstep``).

        Unlike :meth:`launch` (inspection viewers), a run engages the
        :attr:`run_lock` for its duration — so the launch controls disable
        while it executes and a second run is refused — and calls
        ``on_finished`` when it ends so the views can re-read the reduction
        state.  Progress is reported on the activity bar's Build channel.

        Parameters
        ----------
        argv : :obj:`list`
            The run command and its arguments.
        description : :obj:`str`, optional
            Human label for the activity bar (e.g. ``'(Re)building
            wv_calib'``).
        on_finished : callable, optional
            Called as ``on_finished(code)`` when the run ends (success or
            failure), for the refresh-on-completion.

        Returns
        -------
        :obj:`bool`
            True if the process was started, False otherwise.
        """
        if not argv:
            return False
        program, args = argv[0], [str(a) for a in argv[1:]]
        cmd = '"' + ' '.join([program] + args) + '"'

        proc = QProcess()
        proc.setProcessChannelMode(QProcess.MergedChannels)
        self._procs.append(proc)

        # Engage the lock up front so a second launch is refused immediately.
        if self.run_lock is not None:
            self.run_lock.set_dashboard_running(True)

        # Show the running command (with the busy indicator) once started.
        proc.started.connect(
            lambda d=description, c=cmd:
            self._report(f'{d}: running {c}', busy=True))
        # A run that failed to start: release the lock and report.
        proc.errorOccurred.connect(
            lambda _err, p=proc, c=cmd, cb=on_finished:
            self._on_run_error(p, c, cb))
        # A finished run: release the lock, report, and trigger the
        # refresh-on-completion callback.
        proc.finished.connect(
            lambda code, _status, p=proc, c=cmd, cb=on_finished:
            self._on_run_finished(p, c, code, cb))

        proc.start(program, args)
        return True

    def _on_run_finished(self, proc, cmd, code, on_finished):
        """
        Handle a finished run: release the lock, report, refresh.

        Parameters
        ----------
        proc : :obj:`QProcess`
            The finished process.
        cmd : :obj:`str`
            The quoted command line.
        code : :obj:`int`
            The exit code.
        on_finished : callable
            Called as ``on_finished(code)`` (or ``None``).
        """
        output = bytes(proc.readAll()).decode('utf-8', errors='replace')
        if code == 0:
            self._report(f'Ran {cmd} — refreshing state.', busy=False)
        else:
            self._report(f'Ran {cmd} — exited with code {code}.', busy=False)
            log.warning(f'{cmd} output:\n{output}')
        if self.run_lock is not None:
            self.run_lock.set_dashboard_running(False)
        if on_finished is not None:
            on_finished(code)
        self._drop(proc)

    def _on_run_error(self, proc, cmd, on_finished):
        """
        Handle a run that failed to start: release the lock and report.

        Parameters
        ----------
        proc : :obj:`QProcess`
            The failed process.
        cmd : :obj:`str`
            The quoted command line.
        on_finished : callable
            Called as ``on_finished(-1)`` (or ``None``).
        """
        self._report(f'Failed to run {cmd}: {proc.errorString()}',
                     busy=False)
        if self.run_lock is not None:
            self.run_lock.set_dashboard_running(False)
        if on_finished is not None:
            on_finished(-1)
        self._drop(proc)

    def _on_finished(self, proc, cmd, code, hint='viewer window'):
        """
        Handle inspection-process completion: report the outcome and
        capture output.  The exact command stays shown after finishing so
        the user can copy and rerun it.

        Parameters
        ----------
        proc : :obj:`QProcess`
            The finished process.
        cmd : :obj:`str`
            The quoted command line.
        code : :obj:`int`
            The exit code.
        hint : :obj:`str`, optional
            Where the result appears.
        """
        output = bytes(proc.readAll()).decode('utf-8', errors='replace')
        if code == 0:
            self._report(f'Ran {cmd} — see the {hint}.', busy=False,
                         build=False)
        else:
            self._report(f'Ran {cmd} — exited with code {code}.',
                         busy=False, build=False)
            log.warning(f'{cmd} output:\n{output}')
        self._drop(proc)

    def _on_error(self, proc, cmd):
        """
        Handle an inspection process that failed to start (e.g. command not
        found).

        Parameters
        ----------
        proc : :obj:`QProcess`
            The failed process.
        cmd : :obj:`str`
            The quoted command line.
        """
        self._report(f'Failed to run {cmd}: {proc.errorString()}',
                     busy=False, build=False)
        self._drop(proc)

    def _drop(self, proc):
        """
        Release a finished process reference.

        Parameters
        ----------
        proc : :obj:`QProcess`
            The process to release.
        """
        if proc in self._procs:
            self._procs.remove(proc)
