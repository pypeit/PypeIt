"""
Helpers for distributing per-detector work over multiple processes.

The reduction is already factored into module-level functions that operate on a
single detector/mosaic and return a per-detector result
(:func:`~pypeit.pypeit_steps.calib_one`,
:func:`~pypeit.pypeit_steps.process_one_det`,
:func:`~pypeit.pypeit_steps.extract_det`).  :func:`map_over_detectors` maps such
a function over a list of detectors, either serially (``ncpu<=1``; the literal
pre-existing code path) or over a process pool (``ncpu>1``).

The pool does not rely on ``fork``.  The arguments common to every detector are
pickled once per worker and handed to :func:`_worker_init`; each detector's own
slice of the ``payload`` is pickled with its submission.  Any start method
therefore works.  By default PypeIt uses ``forkserver`` on Linux and ``spawn``
elsewhere; see :func:`start_method`.

The serial and parallel paths are built to give the same products.  Two things
make that possible.  First, results are always assembled in the order of the
input detector list, never in completion order.  Second, **both** paths pin the
BLAS/OpenMP thread count to one: the number of threads changes the order in
which a floating-point reduction is accumulated, and the sigma-clipped iterative
fits downstream amplify that last-ulp difference into a visibly different sky
model, so pinning only the workers would make the products depend on ``ncpu``.
The products are bit-identical only where the pinning actually takes effect
(:mod:`threadpoolctl` cannot see every BLAS, e.g. Apple's Accelerate, which is
then pinned in the workers only, through :data:`THREAD_ENV_VARS`) and where the
reduction is reproducible run-to-run in the first place (see `Issue 2208
<https://github.com/pypeit/PypeIt/issues/2208>`__).

.. note::

    Pinning the serial path means a default ``ncpu=1`` reduction no longer uses
    multi-threaded BLAS.  That is a deliberate trade: it costs a few percent of
    wall-clock at ``ncpu=1``, and in exchange the products stop depending on
    ``ncpu`` -- and, as a side effect, on the number of cores the machine has.

.. note::

    A worker does not see module-level state that was configured in the parent
    after import (under ``spawn`` and ``forkserver`` it imports PypeIt afresh).
    Such state must be captured by :func:`_parent_config` and re-applied by
    :func:`_worker_init`.  For the same reason ``func`` must be a module-level
    function: it is pickled by reference.

.. include:: ../include/links.rst
"""
import contextlib
import logging
import multiprocessing as mp
import os
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed

from pypeit import log
from pypeit import PypeItError

__all__ = ['nworkers', 'start_method', 'map_over_detectors', 'surviving_detectors']


THREAD_ENV_VARS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                   'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS')
"""
Environment variables honoured by the BLAS/OpenMP runtimes numpy and scipy may
be linked against.  Set to 1 in the parent while the pool exists, so that each
worker starts with them set *before* numpy loads its BLAS.  That (a) stops N
workers each spawning N threads, (b) makes the floating-point reduction order
inside a worker match the serial run, and (c) is the only way to pin a runtime
that :mod:`threadpoolctl` cannot control, such as Accelerate.
"""

START_METHOD_ENV_VAR = 'PYPEIT_MP_START'
"""
Environment variable that, if set, overrides the default start method of the
worker pool; see :func:`start_method`.  Intended for debugging.
"""

_WORKER_STATE = None
"""
``(func, args, kwargs)`` common to every work item, set in each worker by
:func:`_worker_init`.  Never set in the parent.
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


def start_method():
    """
    The :mod:`multiprocessing` start method used for the worker pool.

    The default is ``forkserver`` on Linux and ``spawn`` elsewhere.  Plain
    ``fork`` is avoided: the parent of a reduction already has threads running
    (the QA encoder pool and the BLAS runtime's threads), and forking a
    multi-threaded process can deadlock the child on a lock one of those threads
    held.  ``forkserver`` sidesteps that by forking from a single-threaded
    server process, which has the pipeline preloaded (so each worker does not
    pay to import it); ``spawn`` starts a fresh interpreter.

    The default can be overridden by setting :data:`START_METHOD_ENV_VAR`
    (``PYPEIT_MP_START``) to any method in
    :func:`multiprocessing.get_all_start_methods`.

    Returns
    -------
    :obj:`str`
        The start method.

    Raises
    ------
    PypeItError
        Raised if :data:`START_METHOD_ENV_VAR` names a method that is not
        available on this platform.
    """
    available = mp.get_all_start_methods()
    method = os.environ.get(START_METHOD_ENV_VAR)
    if method:
        if method not in available:
            raise PypeItError(f'{START_METHOD_ENV_VAR}={method} is not a start method available '
                              f'on this platform; options are: {", ".join(available)}.')
        return method
    return 'forkserver' if sys.platform.startswith('linux') and 'forkserver' in available \
                else 'spawn'


def _get_context(method):
    """
    Return the :mod:`multiprocessing` context for ``method``.

    For ``forkserver`` the pipeline is preloaded into the server.  The server is
    started once per process, on first use, so the preload and the environment
    it inherits are fixed by the first pool to use it.
    """
    ctx = mp.get_context(method)
    if method == 'forkserver':
        ctx.set_forkserver_preload(['pypeit.pypeit_steps'])
    return ctx


@contextlib.contextmanager
def _worker_environment(cfg):
    """
    Temporarily set the environment that new worker processes must inherit.

    Sets every :data:`THREAD_ENV_VARS` to ``'1'`` and, if a matplotlib backend
    has been selected in the parent, ``MPLBACKEND``, so that both take effect in
    a worker before it imports numpy or matplotlib.  The parent's environment is
    restored on exit.  Changing these in the parent is harmless: its BLAS and
    matplotlib are already loaded and do not re-read them.
    """
    env = {var: '1' for var in THREAD_ENV_VARS}
    if cfg.get('mpl_backend') is not None:
        env['MPLBACKEND'] = cfg['mpl_backend']
    saved = {var: os.environ.get(var) for var in env}
    os.environ.update(env)
    try:
        yield
    finally:
        for var, value in saved.items():
            if value is None:
                os.environ.pop(var, None)
            else:
                os.environ[var] = value


def _parent_config():
    """
    Capture the process-global configuration a worker must reproduce.

    Returns
    -------
    :obj:`dict`
        ``polyfit2d_chunk``: :data:`pypeit.core.fitting.POLYFIT2D_CHUNK`, set by
        :class:`~pypeit.pypeit.PypeIt` from ``[rdx] fit2d_chunk``.
        ``mpl_backend``: the matplotlib backend, which the scripts force to
        ``Agg`` unless figures are shown, or None if matplotlib has not been
        imported.
    """
    from pypeit.core import fitting
    cfg = dict(polyfit2d_chunk=fitting.POLYFIT2D_CHUNK, mpl_backend=None)
    if 'matplotlib' in sys.modules:
        cfg['mpl_backend'] = sys.modules['matplotlib'].get_backend()
    return cfg


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


def _worker_init(func, args, kwargs, cfg):
    """
    Run once in each worker process.

    Stores the arguments common to every work item, re-applies the parent's
    process-global configuration (see :func:`_parent_config`), detaches the log
    handlers the worker starts with, and makes QA figures be written in line.

    The handlers have to go whatever the start method: a forked worker inherits
    the parent's (pointing at the parent's log file and stderr), and a freshly
    imported ``pypeit`` attaches its own stderr handler to both ``pypeit`` and
    ``py.warnings``.  Left in place, every line would reach the terminal
    directly, unordered, as well as through the replay in the parent.

    Parameters
    ----------
    func : callable
        Per-detector function.
    args : :obj:`tuple`
        Positional arguments common to every detector.
    kwargs : :obj:`dict`
        Keyword arguments common to every detector.
    cfg : :obj:`dict`
        Configuration returned by :func:`_parent_config`.
    """
    global _WORKER_STATE
    from pypeit import qaWriter
    from pypeit.core import fitting
    _WORKER_STATE = (func, args, kwargs)

    fitting.set_polyfit2d_chunk(cfg['polyfit2d_chunk'])
    if cfg.get('mpl_backend') is not None:
        import matplotlib
        matplotlib.use(cfg['mpl_backend'], force=True)

    for logger in (log, logging.getLogger('py.warnings')):
        for handler in logger.handlers[:]:
            logger.removeHandler(handler)
    # QA figures are written in-line inside a worker: a pool per worker would
    # oversubscribe the machine, and under fork the parent's encoding threads do
    # not exist in the child.  QA threading remains active only in the
    # single-process serial stages of the reduction.
    qaWriter.init(ncpu=1)


def _run_one(key, det_kwargs):
    """
    Execute one work item in a worker process.

    The function to call and its common arguments are read from
    :data:`_WORKER_STATE`, set by :func:`_worker_init`.

    Parameters
    ----------
    key : :obj:`int`, :obj:`tuple`
        The detector or mosaic to work on.
    det_kwargs : :obj:`dict`
        The keyword arguments specific to this detector (its slice of the
        ``payload`` passed to :func:`map_over_detectors`).

    Returns
    -------
    :obj:`tuple`
        ``(key, result, traceback_or_None, log_records)``, where ``result`` is
        None if the call raised, ``traceback_or_None`` is the formatted
        traceback (a :obj:`str`) if it did and None otherwise, and
        ``log_records`` is the :obj:`list` of :class:`logging.LogRecord` objects
        the call emitted, including those of captured Python warnings.
    """
    from threadpoolctl import threadpool_limits
    from pypeit import qaWriter

    func, args, kwargs = _WORKER_STATE

    buf = _RecordBuffer()
    warn_log = logging.getLogger('py.warnings')
    log.addHandler(buf)
    warn_log.addHandler(buf)
    try:
        # The environment variables set by the parent (see
        # _worker_environment) are the primary pinning; threadpool_limits
        # additionally covers any runtime that was initialised before they
        # could take effect, e.g. under fork.
        with threadpool_limits(limits=1):
            result = func(key, *args, **kwargs, **det_kwargs)
            # Make sure this worker's QA figures are on disk before the process
            # is handed the next work item.
            qaWriter.flush()
        tb = None
    except BaseException as e:          # reported to the parent, not swallowed
        result = None
        tb = ''.join(traceback.format_exception(type(e), e, e.__traceback__))
    finally:
        log.removeHandler(buf)
        warn_log.removeHandler(buf)
    return key, result, tb, buf.records


def map_over_detectors(func, detectors, ncpu=1, args=(), kwargs=None, payload=None,
                       label='detector'):
    """
    Map a per-detector function over a list of detectors.

    ``func`` is called as ``func(det, *args, **kwargs, **payload[det])`` and must
    be a module-level function, because it is pickled by reference.  ``args``,
    ``kwargs`` and every ``payload`` entry must be picklable.

    When ``ncpu <= 1`` this is a plain ``for`` loop executed in the calling
    process: no pool, no pickling, no log capture.  That is the default.  The
    loop does pin the BLAS/OpenMP thread count to one, exactly as the workers
    do, so that the products do not depend on ``ncpu``; see :func:`_run_one`.

    When ``ncpu > 1`` the pool is entered *even for a single detector*, so that
    the parallel code path is exercised by the single-detector CI reduction.
    The workers are started with :func:`start_method`.  ``args`` and ``kwargs``
    are pickled once per worker; each ``payload`` entry is pickled only for the
    worker that processes that detector.

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

    # ---------------------------- parallel ---------------------------------
    nw = nworkers(ncpu, len(detectors))
    method = start_method()
    log.info(f'Processing {len(detectors)} {label}(s) using {nw} worker process(es) '
             f'(start method: {method}).')
    cfg = _parent_config()
    results, errors, records = {}, {}, {}
    # NOTE: a fresh executor per call.  The environment has to stay set for the
    # life of the executor, not just its construction, because the workers are
    # started on demand as work is submitted.
    with _worker_environment(cfg), \
            ProcessPoolExecutor(max_workers=nw, mp_context=_get_context(method),
                                initializer=_worker_init,
                                initargs=(func, tuple(args), _kwargs, cfg)) as ex:
        futures = {}
        for det in detectors:
            futures[ex.submit(_run_one, det, _payload.get(det, {}))] = det
            log.info(f'Started {label} {det}')
        for fut in as_completed(futures):
            det, result, tb, recs = fut.result()
            results[det], errors[det], records[det] = result, tb, recs
            log.info(f'Finished {label} {det}')

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
