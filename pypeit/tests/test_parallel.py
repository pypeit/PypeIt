"""
Tests of the per-detector parallel mapping helper.
"""
import io
import logging
import multiprocessing as mp
import os
import sys
import warnings

import numpy as np
import pytest

from pypeit import log
from pypeit import PypeItError
from pypeit import parallel


# Module-level, so that they can be pickled by reference to a worker.
def _work(det, scale=1.0):
    """
    Deterministic float work: a per-detector seed, some linear algebra, a log
    record, and a scalar return.

    NOTE: this particular workload was checked and found *not* to be sensitive
    to the BLAS thread count at this problem size, i.e., it gives the same bits
    under ``threadpool_limits(1)`` and ``threadpool_limits(24)``.  The
    serial/parallel comparison below therefore verifies the plumbing -- result
    assembly order, payload routing and marshalling -- and not the thread
    pinning, which ``test_worker_threads_are_pinned`` covers directly.
    """
    rng = np.random.default_rng(1000 + (det if isinstance(det, int) else sum(det)))
    a = rng.normal(size=(200, 200))
    u, s, vt = np.linalg.svd(a @ a.T)
    log.info(f'worked on {det}')
    return float(scale) * s.sum()


def _fail(det):
    if det == 2:
        raise RuntimeError('boom')
    return det


def _env(det):
    from threadpoolctl import threadpool_info
    # Touch the BLAS so that a runtime that has not been loaded yet is loaded
    # and therefore reported by threadpool_info; otherwise the assertions in
    # test_worker_threads_are_pinned could pass vacuously.
    np.ones((8, 8)) @ np.ones((8, 8))
    return ({var: os.environ.get(var) for var in parallel.THREAD_ENV_VARS},
            [i['num_threads'] for i in threadpool_info()])


def _config(det):
    import matplotlib
    from pypeit.core import fitting
    return fitting.POLYFIT2D_CHUNK, matplotlib.get_backend()


def _shout(det):
    log.info(f'unique-token-{det}')
    return det


def _warn(det):
    warnings.warn(f'worker-warning-{det}', UserWarning)
    return det


DETS = [1, 2, 3, 4]

# The start methods the default can choose (see parallel.start_method).  Every
# pool test runs under each of them that this platform offers.
METHODS = [m for m in ('forkserver', 'spawn') if m in mp.get_all_start_methods()]


@pytest.fixture(params=METHODS)
def start_method(request, monkeypatch):
    """Force the pool start method for the duration of a test."""
    monkeypatch.setenv(parallel.START_METHOD_ENV_VAR, request.param)
    return request.param


DETS = [1, 2, 3, 4]


def test_nworkers():
    # NOTE: the expectations are derived from os.cpu_count() rather than
    # hard-coded.  A CI runner can have as few as 3 cores, so e.g.
    # nworkers(8, 3) is 3 on a workstation but 2 there, and asserting the
    # literal 3 fails on the small machine for the wrong reason.
    ncores = max(1, (os.cpu_count() or 1) - 1)
    assert parallel.nworkers(1, 4) == 1, 'ncpu=1 must collapse to a single worker'
    assert parallel.nworkers(None, 4) == 1, 'ncpu=None must collapse to a single worker'
    assert parallel.nworkers(8, 3) == min(3, ncores), \
            'number of workers must be capped by the item count'
    assert parallel.nworkers(1000, 1) == 1, 'a single item needs a single worker'
    assert parallel.nworkers(1000, 100) == ncores, \
            'number of workers must be capped by the available CPUs'
    assert parallel.nworkers(4, 0) == 1, 'the worker count must never drop below 1'


def test_start_method(monkeypatch):
    monkeypatch.delenv(parallel.START_METHOD_ENV_VAR, raising=False)
    expected = 'forkserver' if sys.platform.startswith('linux') else 'spawn'
    assert parallel.start_method() == expected, \
            'default must be forkserver on Linux and spawn elsewhere'
    monkeypatch.setenv(parallel.START_METHOD_ENV_VAR, 'spawn')
    assert parallel.start_method() == 'spawn', f'{parallel.START_METHOD_ENV_VAR} must override'
    monkeypatch.setenv(parallel.START_METHOD_ENV_VAR, 'not-a-method')
    with pytest.raises(PypeItError, match='not-a-method'):
        parallel.start_method()


def test_serial_and_parallel_are_bit_identical(start_method):
    """
    The serial and parallel paths must return the same bits, in the same order.

    See the note in :func:`_work`: this pins down the mapping machinery.  The
    end-to-end ``--ncpu`` regression in ``test_runpypeit`` is what exercises the
    requirement on real reduction products.
    """
    serial = parallel.map_over_detectors(_work, DETS, ncpu=1, kwargs=dict(scale=2.0))
    par = parallel.map_over_detectors(_work, DETS, ncpu=3, kwargs=dict(scale=2.0))
    assert list(par.keys()) == DETS, 'results must be assembled in detector order'
    for det in DETS:
        assert serial[det] == par[det], f'det {det} differs bit-for-bit'


def test_parallel_is_deterministic(start_method):
    a = parallel.map_over_detectors(_work, DETS, ncpu=3)
    b = parallel.map_over_detectors(_work, DETS, ncpu=3)
    assert all(a[d] == b[d] for d in DETS), 'repeated parallel runs must agree bit-for-bit'


def test_payload_is_per_detector(start_method):
    payload = {d: dict(scale=float(d)) for d in DETS}
    out = parallel.map_over_detectors(_work, DETS, ncpu=2, payload=payload)
    ref = {d: _work(d, scale=float(d)) for d in DETS}
    assert all(out[d] == ref[d] for d in DETS), \
            'each worker must receive only its own payload entry'


def test_worker_threads_are_pinned(start_method):
    """
    The workers must start with the thread variables set, and with every BLAS
    runtime threadpoolctl can see limited to one thread.

    The environment check holds even where threadpoolctl sees nothing (e.g.
    Accelerate on macOS), and it is what pins such a runtime.
    """
    from threadpoolctl import threadpool_info
    before = {var: os.environ.get(var) for var in parallel.THREAD_ENV_VARS}
    out = parallel.map_over_detectors(_env, [1, 2], ncpu=2)
    for det, (env, nthreads) in out.items():
        assert all(v == '1' for v in env.values()), \
                f'thread variables not set in the worker for detector {det}: {env}'
        assert all(n == 1 for n in nthreads), \
                f'BLAS/OpenMP threads not pinned for detector {det}: {nthreads}'
    assert {var: os.environ.get(var) for var in parallel.THREAD_ENV_VARS} == before, \
            "the parent's environment must be restored after the pool exits"
    if len(threadpool_info()) == 0:
        pytest.skip('no runtime visible to threadpoolctl; only the environment was checked')


def test_parent_config_reaches_worker(start_method):
    """
    Process-global configuration set in the parent after import must be
    re-applied in the worker, which under spawn/forkserver imports PypeIt afresh.
    """
    import matplotlib
    from pypeit.core import fitting
    chunk = fitting.POLYFIT2D_CHUNK
    try:
        fitting.set_polyfit2d_chunk(1234)
        out = parallel.map_over_detectors(_config, [1, 2], ncpu=2)
    finally:
        fitting.set_polyfit2d_chunk(chunk)
    for det, (wchunk, backend) in out.items():
        assert wchunk == 1234, f'fit2d_chunk did not reach the worker for detector {det}'
        assert backend.lower() == matplotlib.get_backend().lower(), \
                f'matplotlib backend differs in the worker for detector {det}'


def test_worker_logs_appear_once(start_method, monkeypatch, capfd, caplog):
    """
    A worker's log lines must reach the output only through the replay in the
    parent.  A freshly imported pypeit has its own stderr handler, which the
    worker must remove, or every line would also be printed directly.

    The worker's stderr is the inherited file descriptor, which capfd captures.
    The parent's own stream handlers are redirected for the test, so that
    anything capfd sees came from a worker.
    """
    for logger in (log, logging.getLogger('py.warnings')):
        for handler in logger.handlers:
            if isinstance(handler, logging.StreamHandler):
                monkeypatch.setattr(handler, 'stream', io.StringIO())
    with caplog.at_level('INFO', logger='pypeit'):
        parallel.map_over_detectors(_shout, DETS, ncpu=2)
    err = capfd.readouterr().err
    for det in DETS:
        token = f'unique-token-{det}'
        assert token not in err, f'a worker wrote the log line for detector {det} directly'
        assert sum(token in m for m in caplog.messages) == 1, \
                f'log line for detector {det} must be replayed exactly once'


def test_worker_warnings_replayed(start_method, caplog):
    """Warnings raised in a worker are buffered and replayed in the parent's log."""
    with caplog.at_level('WARNING', logger='pypeit'):
        parallel.map_over_detectors(_warn, [1, 2], ncpu=2)
    for det in [1, 2]:
        assert any(f'worker-warning-{det}' in m for m in caplog.messages), \
                f'warning from detector {det} was not replayed'


def test_worker_logs_replayed_in_detector_order(start_method, caplog):
    with caplog.at_level('INFO', logger='pypeit'):
        parallel.map_over_detectors(_work, DETS, ncpu=3)
    order = [int(m.split()[-1]) for m in caplog.messages if m.startswith('worked on ')]
    assert order == DETS, 'worker log records must be replayed grouped in detector order'


def test_worker_exception_names_the_detector(start_method):
    with pytest.raises(PypeItError, match='detector 2'):
        parallel.map_over_detectors(_fail, DETS, ncpu=3)


def test_serial_exception_propagates():
    """Without a worker to attribute it to, the original exception is raised."""
    with pytest.raises(RuntimeError, match='boom'):
        parallel.map_over_detectors(_fail, DETS, ncpu=1)


def test_surviving_detectors_rebuild():
    """Regression for the detectors.remove()-during-iteration bug."""
    status = {1: True, 2: False, 3: True, 4: True}
    assert parallel.surviving_detectors([1, 2, 3, 4], status) == [1, 3, 4], \
            'the surviving list must keep every successful detector, in order'
    # The old in-loop remove() would have skipped detector 3 entirely.
