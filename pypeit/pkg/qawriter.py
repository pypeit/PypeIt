"""
Deferred writing of the QA figures produced by a reduction.

.. include:: ../include/links.rst
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from PIL import Image


class QAWriter:
    """
    Write QA figures to disk, optionally encoding them in background threads.

    A single instance is created when the package is imported
    (:attr:`pypeit.qaWriter`) and configured once the reduction parameters are
    known, in the same way as :attr:`pypeit.log`::

        from pypeit import qaWriter
        qaWriter.init(ncpu=par['rdx']['ncpu'])
        ...
        qaWriter.save_figure(fig, outfile, dpi=400)
        ...
        qaWriter.flush()

    Rendering and encoding a QA PNG (matplotlib text metrics and draw, then the
    PNG encode) is a significant fraction of the run time of reductions that
    produce many hundreds of QA files.  When ``ncpu>1``, :func:`save_figure`
    rasterizes the figure on the calling thread and hands only the PNG encoding
    to a small thread pool.

    .. warning::

        The rendering *must* stay on the calling thread.  Matplotlib's mathtext
        parser and text-metric machinery are process-global and not
        thread-safe: rendering in a worker races the figure layout being done
        on the main thread and corrupts the parser (e.g., "ParseFatalException:
        Unknown symbol: \\mathdefault" raised for a log-scaled axis).  Encoding
        the rasterized buffer is thread-safe and releases the GIL, so it is the
        part that is deferred.

    The default-constructed instance is fully serial, meaning the QA output of
    a default reduction is written exactly as it was before this class existed.

    **When to call** :func:`flush`: once :func:`save_figure` returns, a figure
    may not be on disk yet.  Call :func:`flush` wherever subsequent code
    expects the files to exist, or wherever the process (or a worker process)
    may end.  In the reduction, that is at the end of each unit of work: the
    end of :func:`~pypeit.pypeit_steps.calib_one` (one detector's
    calibrations), :func:`~pypeit.exposure.reduce_exposure` (one exposure),
    :meth:`~pypeit.pypeit.PypeIt.calib_all` and
    :meth:`~pypeit.pypeit.PypeIt.reduce_all`, and immediately before the QA
    HTML pages that link to the PNGs are built.  Independently of these,
    :func:`save_figure` flushes automatically once :attr:`max_pending` encodes
    are queued, which bounds the memory held by the queue.

    Parameters
    ----------
    ncpu : :obj:`int`, optional
        Number of figure-encoding threads, capped at :attr:`max_threads`.
        Values less than 2 write the figures in-line, with no thread pool.

    Attributes
    ----------
    pool : :class:`concurrent.futures.ThreadPoolExecutor`
        Pool used to encode the figures, or None when writing serially.
    pending : :obj:`list`
        The :class:`concurrent.futures.Future` objects for the figures handed
        to :attr:`pool`, one per figure.  A future is added by
        :func:`save_figure`, and the list is emptied by :func:`flush`, which
        waits for each future to finish.  Always empty when writing serially.
    max_pending : :obj:`int`
        Number of queued encodes that triggers an automatic :func:`flush`.
        This bounds the memory held by the queued rasters; a large QA figure
        rasterizes to more than 100 MB.
    """

    #: Maximum number of encoder threads, regardless of ``ncpu``
    max_threads = 8

    def __init__(self, ncpu:int=1):
        self.pool = None
        self.pending = []
        self.max_pending = 4
        self.init(ncpu=ncpu)

    def init(self, ncpu:int=1):
        """
        (Re)initialize the encoding thread pool.

        Call this once the reduction parameters are known.  Calling it with
        ``ncpu<2`` restores fully serial, in-line figure writing.

        This is safe to call more than once.  In particular, a process forked
        from one that was using a pool inherits an instance whose threads do
        not exist in the child, so the child must call this to build its own
        pool (or to go serial).  Any queued encodes are dropped rather than
        waited on, since they belong to the threads of the calling process;
        call :func:`flush` first if they still need to be written.

        Parameters
        ----------
        ncpu : :obj:`int`, optional
            Number of figure-encoding threads, capped at :attr:`max_threads`.
            Values less than 2 disable the pool.
        """
        if self.pool is not None:
            self.pool.shutdown(wait=False)
        self.pool = None
        self.pending = []
        if ncpu is not None and ncpu > 1:
            self.pool = ThreadPoolExecutor(max_workers=min(int(ncpu), self.max_threads),
                                           thread_name_prefix='pypeit-qa')

    @property
    def parallel(self):
        """Flag that figures are being encoded by a thread pool."""
        return self.pool is not None

    def save_figure(self, fig, outfile=None, dpi=None, show:bool=False):
        """
        Write a matplotlib figure to disk and close it.

        When the pool is active (see :func:`init`), the figure is rasterized
        here and only its PNG encoding is deferred to a background thread; the
        file written is pixel-identical to the one written by
        :meth:`matplotlib.figure.Figure.savefig`.  Use :func:`flush` to ensure
        the deferred writes have completed.

        Parameters
        ----------
        fig : :class:`matplotlib.figure.Figure`
            Figure to write.  It is closed by this call and must not be used
            afterwards.
        outfile : :obj:`str`, `Path`_, optional
            Output file.  If None, the figure is not written.
        dpi : :obj:`float`, optional
            Resolution of the written figure.  If None, use the matplotlib
            default (``rcParams['savefig.dpi']``).
        show : :obj:`bool`, optional
            Show the figure interactively.  This forces the synchronous path.
        """
        if outfile is None or not self._can_defer(outfile, show):
            # Write (and/or show) the figure synchronously
            if outfile is not None:
                fig.savefig(outfile, **({} if dpi is None else {'dpi': dpi}))
            if show:
                plt.show()
            plt.close(fig)
            return

        # Rasterize here (see the class warning), then queue only the encode.
        if dpi is None:
            dpi = matplotlib.rcParams['savefig.dpi']
        if dpi != 'figure':
            fig.set_dpi(dpi)
        fig.canvas.draw()
        rgba = np.asarray(fig.canvas.buffer_rgba()).copy()
        # Deregister the figure from pyplot now, so that pyplot-state plotting
        # done by a subsequent QA function cannot draw into a figure that is
        # still queued.
        plt.close(fig)
        self.pending.append(self.pool.submit(self._encode_png, rgba, outfile, fig.dpi))
        if len(self.pending) >= self.max_pending:
            self.flush()

    def _can_defer(self, outfile, show):
        """
        Check if a figure can be written by the deferred path.

        The deferred path writes the raster that ``fig.canvas.draw()`` produces
        directly to a PNG.  That is the same image
        :meth:`~matplotlib.figure.Figure.savefig` writes under matplotlib's
        default ``savefig`` settings, but not if the output format is not PNG,
        or if the ``savefig.*`` rcParams request a change to the image (e.g., a
        ``'tight'`` bounding box, a transparent background, or a face color
        that differs from the figure's).  Those cases, and interactive display,
        are written synchronously instead.  This only affects how quickly the
        file is written, not its contents, so it is not relevant to users.

        Parameters
        ----------
        outfile : :obj:`str`, `Path`_
            Output file.
        show : :obj:`bool`
            Whether the figure is also to be shown interactively.

        Returns
        -------
        :obj:`bool`
            True if the figure can be written by the deferred path.
        """
        rc = matplotlib.rcParams
        return self.pool is not None and not show \
                and Path(outfile).suffix.lower() == '.png' \
                and rc['savefig.bbox'] is None and not rc['savefig.transparent'] \
                and rc['savefig.facecolor'] == 'auto' and rc['savefig.edgecolor'] == 'auto'

    def flush(self):
        """
        Block until every deferred figure has been written.

        Exceptions raised while encoding are re-raised here, on the calling
        thread.  This is a no-op when writing serially.
        """
        try:
            for future in self.pending:
                future.result()
        finally:
            # Empty the queue even if an encode failed, so that a subsequent
            # flush does not wait on (and re-raise from) the same futures.
            self.pending = []

    @staticmethod
    def _encode_png(rgba, outfile, dpi):
        """
        Write a rasterized figure to disk as a PNG.

        This is the unit of work handed to the encoding threads.

        Parameters
        ----------
        rgba : `numpy.ndarray`_
            The rasterized figure; shape is ``(nrows, ncols, 4)``, type is
            ``uint8``.
        outfile : :obj:`str`, `Path`_
            Output PNG file.
        dpi : :obj:`float`
            Resolution recorded in the PNG metadata.
        """
        Image.fromarray(rgba).save(outfile, dpi=(dpi, dpi))

    def __repr__(self):
        if self.pool is None:
            return f'<{self.__class__.__name__}: serial>'
        return f'<{self.__class__.__name__}: {self.pool._max_workers} encoding threads, ' \
               f'{len(self.pending)} pending>'
