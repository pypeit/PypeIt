"""
Tests for the PypeItRemote ginga plugin job logic and its client.
"""
import argparse
import json
import logging
from pathlib import Path
import threading
import xmlrpc.client

import numpy as np
import pytest

from pypeit import __version__
from pypeit import log
from pypeit import PypeItError
from pypeit import slittrace
from pypeit.display import remote_handlers
from pypeit.display import remote_jobs
from pypeit.display import remote


class FakeShell:
    """
    Stand-in for the ginga shell: GUI calls run inline, and worker-thread
    calls run inline or on real threads.
    """
    def __init__(self, threaded=False):
        self.threaded = threaded
        self.threads = []
        self.logger = logging.getLogger('test_ginga_remote')

    def gui_call(self, method, *args, **kwargs):
        return method(*args, **kwargs)

    def nongui_do(self, method, *args, **kwargs):
        if not self.threaded:
            return method(*args, **kwargs)
        thread = threading.Thread(target=method, args=args, kwargs=kwargs)
        self.threads.append(thread)
        thread.start()

    def join(self):
        for thread in self.threads:
            thread.join(timeout=10)


def _args_json(**kwargs):
    return json.dumps(kwargs)


def test_submit_runs_handler():
    received = {}

    def handler(fv, args, job):
        received['args'] = args
        received['fv'] = fv

    fv = FakeShell()
    manager = remote_jobs.JobManager(fv, {'show': handler})
    status = manager.submit('show', _args_json(file='/a/b.fits', det=None), __version__)
    assert status['state'] == 'done', 'A handler that returns normally should finish the job'
    assert status['warnings'] == [], 'Matching versions should not produce a warning'
    assert isinstance(received['args'], argparse.Namespace), \
        'The handler should receive an argparse.Namespace'
    assert received['args'].file == '/a/b.fits', 'String arguments should be passed through'
    assert received['args'].det is None, 'None arguments should survive the JSON round trip'
    assert received['fv'] is fv, 'The handler should receive the ginga shell'
    assert manager.status(status['id'])['state'] == 'done', \
        'status() should report the same state as submit()'


def test_submit_unknown_command():
    manager = remote_jobs.JobManager(FakeShell(), {'show': lambda fv, args, job: None})
    status = manager.submit('nope', _args_json(), __version__)
    assert status['state'] == 'error', 'An unknown command should be an error'
    assert 'nope' in status['message'] and 'show' in status['message'], \
        'The error should name the unknown command and the available commands'


def test_submit_bad_json():
    manager = remote_jobs.JobManager(FakeShell(), {'show': lambda fv, args, job: None})
    status = manager.submit('show', 'not json', __version__)
    assert status['state'] == 'error', 'Arguments that are not valid JSON should be an error'


def test_handler_errors():
    def pypeit_error(fv, args, job):
        raise PypeItError('bad file')

    def other_error(fv, args, job):
        raise ValueError('bad value')

    manager = remote_jobs.JobManager(
        FakeShell(), {'a': pypeit_error, 'b': other_error}, logger=logging.getLogger('quiet')
    )
    status = manager.submit('a', _args_json(), __version__)
    assert status['state'] == 'error', 'A PypeItError should end the job with an error'
    assert status['message'] == 'bad file', 'A PypeItError message should be passed unchanged'
    status = manager.submit('b', _args_json(), __version__)
    assert status['state'] == 'error', 'Any other exception should end the job with an error'
    assert status['message'] == 'ValueError: bad value', \
        'Other exceptions should be reported with their type'


def test_version_mismatch_warning():
    manager = remote_jobs.JobManager(FakeShell(), {'show': lambda fv, args, job: None})
    status = manager.submit('show', _args_json(), '0.0.0')
    assert status['state'] == 'done', 'A version mismatch should not stop the command'
    assert len(status['warnings']) == 1 and '0.0.0' in status['warnings'][0], \
        'A version mismatch should produce one warning naming the script version'


def test_log_warnings_are_collected():
    def handler(fv, args, job):
        log.warning('slit file missing')
        log.info('just information')

    manager = remote_jobs.JobManager(FakeShell(), {'show': handler})
    status = manager.submit('show', _args_json(), __version__)
    assert status['warnings'] == ['slit file missing'], \
        'Warnings logged by the job thread, and only warnings, should be reported'
    assert all(not isinstance(h, remote_jobs._JobLogHandler) for h in log.handlers), \
        'The job log handler should be removed when the job finishes'


def test_last_command_wins():
    started = threading.Event()
    release = threading.Event()

    def slow(fv, args, job):
        started.set()
        release.wait(timeout=10)
        job.check_cancel()

    fv = FakeShell(threaded=True)
    manager = remote_jobs.JobManager(fv, {'slow': slow, 'fast': lambda fv, args, job: None})
    first = manager.submit('slow', _args_json(), __version__)
    assert started.wait(timeout=10), 'The first job should start'
    second = manager.submit('fast', _args_json(), __version__)
    release.set()
    fv.join()
    assert manager.status(first['id'])['state'] == 'cancelled', \
        'A running job should be cancelled when a new command is submitted'
    assert manager.status(second['id'])['state'] == 'done', \
        'The newest command should run after the cancelled one stops'


def test_cancel():
    def handler(fv, args, job):
        job.cancel()
        job.check_cancel()

    manager = remote_jobs.JobManager(FakeShell(), {'show': handler})
    status = manager.submit('show', _args_json(), __version__)
    assert status['state'] == 'cancelled', 'check_cancel() should stop a cancelled job'
    assert manager.cancel(999)['state'] == 'error', 'Cancelling an unknown job should be an error'


def test_old_jobs_are_forgotten():
    manager = remote_jobs.JobManager(
        FakeShell(), {'show': lambda fv, args, job: None}, max_jobs=5
    )
    ids = [manager.submit('show', _args_json(), __version__)['id'] for i in range(12)]
    assert len(manager.jobs) == 5, 'Only max_jobs jobs should be kept'
    assert ids[-1] in manager.jobs, 'The newest job should be kept'
    assert manager.status(ids[0])['state'] == 'error', 'A forgotten job should be unknown'


def test_args_to_json():
    args = argparse.Namespace(file='spec2d.fits', slits=None, det='DET01', flag=True)
    decoded = json.loads(remote.args_to_json(args, path_args=['file', 'slits']))
    assert decoded['file'] == str(Path('spec2d.fits').absolute()), \
        'Path arguments should be made absolute'
    assert decoded['slits'] is None, 'Path arguments that are None should stay None'
    assert decoded['det'] == 'DET01' and decoded['flag'] is True, \
        'Other arguments should be unchanged'


class FakeViewer:
    """
    Stand-in for RemoteViewer that talks to a JobManager directly.
    """
    def __init__(self, manager, ping_error=None):
        self.manager = manager
        self.ping_error = ping_error
        self.calls = []

    def call(self, method_name, *args):
        self.calls.append(method_name)
        if method_name == 'ping':
            if self.ping_error is not None:
                raise self.ping_error
            return __version__
        return getattr(self.manager, method_name)(*args)


def test_send_command_done(caplog):
    def handler(fv, args, job):
        log.warning('watch out')

    viewer = FakeViewer(remote_jobs.JobManager(FakeShell(), {'show': handler}))
    with caplog.at_level(logging.WARNING, logger='pypeit'):
        status = remote.send_command('show', _args_json(), viewer=viewer)
    assert status['state'] == 'done', 'send_command should return the finished status'
    assert 'watch out' in caplog.text, 'Job warnings should be logged by the client'


def test_send_command_error():
    def handler(fv, args, job):
        raise PypeItError('no such detector')

    viewer = FakeViewer(remote_jobs.JobManager(FakeShell(), {'show': handler}))
    with pytest.raises(PypeItError, match='no such detector'):
        remote.send_command('show', _args_json(), viewer=viewer)


def test_send_command_timeout():
    release = threading.Event()

    def handler(fv, args, job):
        release.wait(timeout=10)

    fv = FakeShell(threaded=True)
    viewer = FakeViewer(remote_jobs.JobManager(fv, {'show': handler}))
    status = remote.send_command(
        'show', _args_json(), timeout=0.2, poll_interval=0.05, viewer=viewer
    )
    release.set()
    fv.join()
    assert status['state'] in ('queued', 'running'), \
        'On timeout, send_command should return without waiting for the job'
    assert 'cancel' not in viewer.calls, 'A timeout should not cancel the job in ginga'


def test_send_command_without_plugin():
    viewer = FakeViewer(None, ping_error=xmlrpc.client.Fault(1, 'KeyError'))
    with pytest.raises(PypeItError, match='PypeItRemote'):
        remote.send_command('show', _args_json(), viewer=viewer)


class FakeSpecObjs:
    """
    Stand-in for SpecObjs with just what spec1d_extension uses.
    """
    def __init__(self, names, opt):
        self.NAME = np.array(names)
        self._opt = opt

    def __len__(self):
        return len(self.NAME)

    def __getitem__(self, i):
        return {'OPT_WAVE': np.arange(3) if self._opt[i] else None}


def test_spec1d_extension():
    sobjs = FakeSpecObjs(['SPAT0100-SLIT0001-DET01', 'SPAT0200-SLIT0002-DET01'], [True, False])
    assert remote_handlers.spec1d_extension(sobjs) == 0, \
        'The default extension 1 should select the first spectrum'
    assert remote_handlers.spec1d_extension(sobjs, obj='SPAT0100-SLIT0001-DET01') == 0, \
        'A spectrum should be selectable by name'
    assert remote_handlers.spec1d_extension(sobjs, exten=2, extract='BOX') == 1, \
        'A spectrum should be selectable by 1-based extension number'
    with pytest.raises(PypeItError, match='2 spectra'):
        remote_handlers.spec1d_extension(sobjs, exten=3)
    with pytest.raises(PypeItError, match='2 spectra'):
        remote_handlers.spec1d_extension(sobjs, exten=0)
    with pytest.raises(PypeItError, match='--list'):
        remote_handlers.spec1d_extension(sobjs, obj='nope')
    with pytest.raises(PypeItError, match='BOX'):
        remote_handlers.spec1d_extension(sobjs, exten=2)


def test_handlers_are_registered():
    handlers = remote_handlers.get_handlers()
    assert 'show_1dspec' in handlers, 'pypeit_show_1dspec should have a handler'
    assert all(callable(h) for h in handlers.values()), 'Every handler should be callable'


def test_find_slits_file(tmp_path):
    name = slittrace.SlitTraceSet.construct_file_name('A_0_DET01')
    recorded = tmp_path / 'recorded'
    moved = tmp_path / 'moved'
    recorded.mkdir()
    moved.mkdir()
    assert remote_handlers.find_slits_file('A_0_DET01', str(recorded), moved) is None, \
        'No slits file should be found when none exists'
    (moved / name).touch()
    assert remote_handlers.find_slits_file('A_0_DET01', str(recorded), moved) == moved / name, \
        'The slits file next to the calibration file should be found if the recorded one is missing'
    (recorded / name).touch()
    assert remote_handlers.find_slits_file('A_0_DET01', str(recorded), moved) \
                == recorded / name, 'The slits file in the recorded directory should come first'
    assert remote_handlers.find_slits_file(None, str(recorded), moved) is None, \
        'Without a calibration key, no slits file can be found'


def test_read_slits_unreadable(tmp_path):
    bad = tmp_path / 'Slits_A_0_DET01.fits.gz'
    bad.write_text('not a fits file')
    assert remote_handlers.read_slits(None) is None, 'No file should give no slits'
    assert remote_handlers.read_slits(bad) is None, \
        'An unreadable slits file should give no slits rather than an error'
