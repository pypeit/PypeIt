"""
Helpers for distributing per-detector work over multiple processes.

The reduction is already factored into module-level functions that operate on a
single detector/mosaic and return a per-detector result
(:func:`~pypeit.pypeit_steps.calib_one`,
:func:`~pypeit.pypeit_steps.process_one_det`,
:func:`~pypeit.pypeit_steps.extract_det`).  :func:`map_over_detectors` maps such
a function over a list of detectors, either serially (``ncpu<=1``; the literal
pre-existing code path) or over a ``fork``-based process pool (``ncpu>1``).

The two paths are built to produce *bit-identical* products.  Two things make
that true.  First, results are always assembled in the order of the input
detector list, never in completion order.  Second, **both** paths pin the
BLAS/OpenMP thread count to one: the number of threads changes the order in
which a floating-point reduction is accumulated, and the sigma-clipped iterative
fits downstream amplify that last-ulp difference into a visibly different sky
model, so pinning only the workers would make the products depend on ``ncpu``.

.. note::

    Pinning the serial path means a default ``ncpu=1`` reduction no longer uses
    multi-threaded BLAS.  That is a deliberate trade: it costs a few percent of
    wall-clock at ``ncpu=1``, and in exchange the products stop depending on
    ``ncpu`` -- and, as a side effect, on the number of cores the machine has.

.. note::

    This is version 1 of the detector parallelism and it requires the ``fork``
    start method, i.e., Linux.  On platforms where ``fork`` is unavailable
    (Windows) or unsafe (macOS), :func:`map_over_detectors` warns and falls back
    to the serial loop.  Supporting ``spawn`` is a future goal, which is why
    ``func`` must be a module-level function: under ``fork`` a closure happens to
    work, but it must not be relied on.

.. include:: ../include/links.rst
"""
import logging
import multiprocessing as mp
import os
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed

from pypeit import log
from pypeit import PypeItError

__all__ = ['nworkers', 'map_over_detectors', 'surviving_detectors']


THREAD_ENV_VARS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                   'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS')
"""
Environment variables honoured by the BLAS/OpenMP runtimes numpy and scipy may
be linked against.  Pinned to 1 in each worker so that (a) N workers do not each
spawn N threads, and (b) the floating-point reduction order inside a worker
matches the serial run bit-for-bit.
"""

_FORK_PAYLOAD = None
"""
Module-global set by :func:`map_over_detectors` immediately *before* the pool is
created.  Because the workers are forked, they inherit this copy-on-write and
only the (tiny) detector key has to be pickled on submit.  Never read outside a
worker.
"""


def nworkers(ncpu, nitem):
    """
    Number of worker processes to use.

    Parameters
    ----------
    ncpu : :obj:`int`
        Requested number of CPUs (``par['rdx']['ncpu']``).
    nitem : :obj:`int`
        Number of independent work items (detectors).

    Returns
    -------
    :obj:`int`
        ``min(ncpu, nitem, os.cpu_count()-1)``, at least 1.
    """
    if ncpu is None or ncpu <= 1:
        return 1
    return max(1, min(int(ncpu), int(nitem), max(1, (os.cpu_count() or 1) - 1)))


class _RecordBuffer(logging.Handler):
    """
    Logging handler that accumulates picklable copies of the emitted records.

    The record message is rendered eagerly and ``args``/``exc_info`` are dropped
    (the approach used by `logging.handlers.QueueHandler`_), so the buffer can be
    shipped back to the parent process.  ``record.created`` is preserved, so the
    replayed lines carry the worker's real timestamps, as do ``filename``,
    ``funcName`` and ``lineno``, which the stream and file formatters include.

    Attributes
    ----------
    records : :obj:`list`
        The accumulated :class:`logging.LogRecord` objects, in emission order.
    """
    def __init__(self):
        super().__init__(level=logging.DEBUG)
        self.records = []

    def emit(self, record):
        """
        Append a picklable copy of ``record`` to :attr:`records`.

        Parameters
        ----------
        record : :class:`logging.LogRecord`
            The record to buffer.
        """
        rec = logging.makeLogRecord(record.__dict__)
        rec.msg = record.getMessage()
        rec.args = None
        rec.exc_info = None
        rec.exc_text = None
        self.records.append(rec)


def _worker_init():
    """
    Run once in each forked worker process.

    Pins the BLAS/OpenMP thread count, detaches the log handlers inherited from
    the parent (they point at the parent's log file and stderr), and disables the
    inherited QA thread pool (its threads do not survive the fork; see
    :func:`~pypeit.pkg.qawriter.QAWriter.init`).
    """
    from pypeit import qaWriter
    for var in THREAD_ENV_VARS:
        os.environ[var] = '1'
    for handler in log.handlers[:]:
        log.removeHandler(handler)
    warn_log = logging.getLogger('py.warnings')
    for handler in warn_log.handlers[:]:
        warn_log.removeHandler(handler)
    # QA figures are written in-line inside a worker: the parent's encoding
    # threads do not exist in the forked child, and a pool per worker would
    # oversubscribe the machine.  QA threading remains active only in the
    # single-process serial stages of the reduction.
    qaWriter.init(ncpu=1)


def _run_one(key):
    """
    Execute one work item in a worker process.

    The function to call and its arguments are read from :attr:`_FORK_PAYLOAD`,
    which the worker inherited from the parent when it was forked, so only
    ``key`` crosses the process boundary on submit.

    Parameters
    ----------
    key : :obj:`int`, :obj:`tuple`
        The detector or mosaic to work on.

    Returns
    -------
    :obj:`tuple`
        ``(key, result, traceback_or_None, log_records)``, where ``result`` is
        None if the call raised, ``traceback_or_None`` is the formatted
        traceback (a :obj:`str`) if it did and None otherwise, and
        ``log_records`` is the :obj:`list` of :class:`logging.LogRecord` objects
        the call emitted.
    """
    from threadpoolctl import threadpool_limits
    from pypeit import qaWriter

    func, args, kwargs, payload = _FORK_PAYLOAD

    buf = _RecordBuffer()
    log.addHandler(buf)
    try:
        # threadpool_limits handles runtimes that were already initialised in
        # the parent, which the environment variables set by _worker_init
        # cannot.
        with threadpool_limits(limits=1):
            result = func(key, *args, **kwargs, **payload.get(key, {}))
            # Make sure this worker's QA figures are on disk before the process
            # is handed the next work item.
            qaWriter.flush()
        tb = None
    except BaseException as e:          # reported to the parent, not swallowed
        result = None
        tb = ''.join(traceback.format_exception(type(e), e, e.__traceback__))
    finally:
        log.removeHandler(buf)
    return key, result, tb, buf.records


def map_over_detectors(func, detectors, ncpu=1, args=(), kwargs=None, payload=None,
                       label='detector'):
    """
    Map a per-detector function over a list of detectors.

    ``func`` is called as ``func(det, *args, **kwargs, **payload[det])`` and must
    be a module-level function (a requirement for the future ``spawn`` support;
    with ``fork`` a closure happens to work, but do not rely on it).

    When ``ncpu <= 1`` this is a plain ``for`` loop executed in the calling
    process: no pool, no pickling, no log capture.  That is the default.  The
    loop does pin the BLAS/OpenMP thread count to one, exactly as the workers
    do, so that the products do not depend on ``ncpu``; see :func:`_run_one`.

    When ``ncpu > 1`` the pool is entered *even for a single detector*, so that
    the parallel code path is exercised by the single-detector CI reduction.  If
    the platform cannot ``fork``, a warning is issued and the serial loop is used
    instead.

    Parameters
    ----------
    func : callable
        Per-detector function.  First argument is the detector.
    detectors : :obj:`list`
        Detectors/mosaics.  Ints for single detectors, tuples for mosaics.
        Defines the deterministic result-assembly order.
    ncpu : :obj:`int`, optional
        ``par['rdx']['ncpu']``.
    args : :obj:`tuple`, optional
        Positional arguments common to every detector.
    kwargs : :obj:`dict`, optional
        Keyword arguments common to every detector.
    payload : :obj:`dict`, optional
        ``{det: {kwarg: value}}`` of *per-detector* keyword arguments.  Keeping
        the large per-detector objects here (instead of in ``args``) means each
        worker only ever touches its own slice.
    label : :obj:`str`, optional
        Noun used in the progress log lines.

    Returns
    -------
    :obj:`dict`
        ``{det: result}``, iterated in the order of ``detectors``.

    Raises
    ------
    PypeItError
        Raised, in detector order, if any worker raised.
    """
    global _FORK_PAYLOAD
    _kwargs = {} if kwargs is None else kwargs
    _payload = {} if payload is None else payload

    def _serial():
        """Run every detector in the calling process, with threads pinned."""
        # NOTE: threadpool_limits is applied here as well as in the workers
        # (see _run_one).  It has to be: pinning only the workers would make the
        # result depend on ncpu, because BLAS computes a different -- equally
        # valid -- answer at 1 thread than at N, and the sigma-clipped iterative
        # fits downstream amplify that last-ulp difference into a visibly
        # different sky model.  Pinning both paths is what makes the ncpu=1 and
        # ncpu>1 products bit-identical.
        from threadpoolctl import threadpool_limits
        with threadpool_limits(limits=1):
            return {det: func(det, *args, **_kwargs, **_payload.get(det, {}))
                    for det in detectors}

    # --------------------------- serial loop -------------------------------
    if ncpu is None or ncpu <= 1:
        return _serial()

    if 'fork' not in mp.get_all_start_methods():
        # Version 1 of the detector parallelism is fork-only; see the module
        # docstring.
        log.warning(f'The "fork" start method is not available on this platform, so ncpu={ncpu} '
                    f'cannot be used to process {label}s concurrently.  Falling back to '
                    f'processing the {label}s serially.')
        return _serial()

    # ---------------------------- parallel ---------------------------------
    nw = nworkers(ncpu, len(detectors))
    log.info(f'Processing {len(detectors)} {label}(s) using {nw} worker process(es).')
    _FORK_PAYLOAD = (func, tuple(args), _kwargs, _payload)
    results, errors, records = {}, {}, {}
    try:
        # NOTE: a fresh executor per call; the workers are forked here and
        # inherit _FORK_PAYLOAD copy-on-write.
        with ProcessPoolExecutor(max_workers=nw, mp_context=mp.get_context('fork'),
                                 initializer=_worker_init) as ex:
            futures = {}
            for det in detectors:
                futures[ex.submit(_run_one, det)] = det
                log.info(f'Started {label} {det}')
            for fut in as_completed(futures):
                det, result, tb, recs = fut.result()
                results[det], errors[det], records[det] = result, tb, recs
                log.info(f'Finished {label} {det}')
    finally:
        _FORK_PAYLOAD = None

    # Replay the worker logs grouped in detector order, so each detector's block
    # stays contiguous in the run log.  Formatter timestamps come from
    # record.created, i.e. the worker's real clock.
    for det in detectors:
        recs = records.get(det, [])
        if len(recs) == 0:
            continue
        log.info(f'----- begin log for {label} {det} -----')
        for rec in recs:
            log.handle(rec)
        log.info(f'----- end log for {label} {det} -----')

    # Re-raise the first failure, in detector order (deterministic).
    for det in detectors:
        if errors.get(det) is not None:
            raise PypeItError(f'{label} {det} failed in a worker process:\n{errors[det]}')

    return {det: results[det] for det in detectors}


def surviving_detectors(detectors, success):
    """
    Rebuild the list of detectors that should continue to be reduced.

    Parameters
    ----------
    detectors : :obj:`list`
        The detectors that were attempted, in order.
    success : :obj:`dict`
        ``{det: bool}``.

    Returns
    -------
    :obj:`list`
        The subset of ``detectors``, in the original order, for which
        ``success[det]`` is True.
    """
    return [det for det in detectors if success.get(det, False)]
