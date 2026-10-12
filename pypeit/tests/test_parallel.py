"""
Tests of the per-detector parallel mapping helper.
"""
import os

import numpy as np
import pytest

from pypeit import log
from pypeit import PypeItError
from pypeit import parallel


# Module-level so it is importable in a forked worker.
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
    return os.environ.get('OMP_NUM_THREADS'), [i['num_threads'] for i in threadpool_info()]


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


def test_serial_and_parallel_are_bit_identical():
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


def test_parallel_is_deterministic():
    a = parallel.map_over_detectors(_work, DETS, ncpu=3)
    b = parallel.map_over_detectors(_work, DETS, ncpu=3)
    assert all(a[d] == b[d] for d in DETS), 'repeated parallel runs must agree bit-for-bit'


def test_payload_is_per_detector():
    payload = {d: dict(scale=float(d)) for d in DETS}
    out = parallel.map_over_detectors(_work, DETS, ncpu=2, payload=payload)
    ref = {d: _work(d, scale=float(d)) for d in DETS}
    assert all(out[d] == ref[d] for d in DETS), \
            'each worker must receive only its own payload entry'


@pytest.mark.skipif(not parallel.can_fork(),
                    reason='no fork: map_over_detectors runs serially, with no worker to pin')
def test_worker_threads_are_pinned():
    from threadpoolctl import threadpool_info
    if len(threadpool_info()) == 0:
        pytest.skip('no BLAS/OpenMP runtime to pin in this environment')
    out = parallel.map_over_detectors(_env, [1, 2], ncpu=2)
    for det, (omp, nthreads) in out.items():
        assert omp == '1', f'OMP_NUM_THREADS not pinned in the worker for detector {det}'
        assert len(nthreads) > 0, f'no runtime reported by the worker for detector {det}'
        assert all(n == 1 for n in nthreads), \
                f'BLAS/OpenMP threads not pinned for detector {det}: {nthreads}'


def test_worker_logs_replayed_in_detector_order(caplog):
    with caplog.at_level('INFO', logger='pypeit'):
        parallel.map_over_detectors(_work, DETS, ncpu=3)
    order = [int(m.split()[-1]) for m in caplog.messages if m.startswith('worked on ')]
    assert order == DETS, 'worker log records must be replayed grouped in detector order'


@pytest.mark.skipif(not parallel.can_fork(),
                    reason='no fork: the serial path propagates the original exception, '
                           'rather than wrapping it per detector')
def test_worker_exception_names_the_detector():
    with pytest.raises(PypeItError, match='detector 2'):
        parallel.map_over_detectors(_fail, DETS, ncpu=3)


def test_serial_exception_propagates():
    """Without a worker to attribute it to, the original exception is raised."""
    with pytest.raises(RuntimeError, match='boom'):
        parallel.map_over_detectors(_fail, DETS, ncpu=1)


def test_fork_unavailable_falls_back_to_serial(monkeypatch, caplog):
    """
    On a platform without usable ``fork``, ``ncpu>1`` must warn and run serially.

    On Linux the guard is otherwise dead code, so it is forced here.
    """
    monkeypatch.setattr(parallel, 'can_fork', lambda: False)
    with caplog.at_level('WARNING', logger='pypeit'):
        out = parallel.map_over_detectors(_work, DETS, ncpu=3, kwargs=dict(scale=2.0))
    assert any('fork' in m for m in caplog.messages), 'the fallback must warn about "fork"'
    ref = {d: _work(d, scale=2.0) for d in DETS}
    assert all(out[d] == ref[d] for d in DETS), 'the fallback must still return every result'


def test_surviving_detectors_rebuild():
    """Regression for the detectors.remove()-during-iteration bug."""
    status = {1: True, 2: False, 3: True, 4: True}
    assert parallel.surviving_detectors([1, 2, 3, 4], status) == [1, 3, 4], \
            'the surviving list must keep every successful detector, in order'
    # The old in-loop remove() would have skipped detector 3 entirely.
