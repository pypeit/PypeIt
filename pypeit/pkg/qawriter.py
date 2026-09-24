"""
Deferred writing of the QA figures produced by a reduction.

.. include:: ../include/links.rst
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

import numpy as np


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

    Parameters
    ----------
    ncpu : :obj:`int`, optional
        Number of figure-encoding threads.  Values less than 2 write the
        figures in-line, with no thread pool.

    Attributes
    ----------
    pool : `concurrent.futures.ThreadPoolExecutor`_
        Pool used to encode the figures, or None when writing serially.
    pending : :obj:`list`
        Futures for the encodes that have not yet been reaped.
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
            Number of figure-encoding threads.  Values less than 2 disable the
            pool.
        """
        old = self.pool
        self.pool = None
        self.pending = []
        if old is not None:
            old.shutdown(wait=False)
        if ncpu is not None and ncpu > 1:
            self.pool = ThreadPoolExecutor(max_workers=min(int(ncpu), self.max_threads),
                                           thread_name_prefix='pypeit-qa')

    @property
    def parallel(self):
        """Flag that figures are being encoded by a thread pool."""
        return self.pool is not None

    def save_figure(self, fig, outfile, show:bool=False, close:bool=True, **kwargs):
        """
        Write a matplotlib figure to disk.

        When the pool is active (see :func:`init`), the figure is rasterized
        here and only its PNG encoding is deferred to a background thread; the
        file written is identical to the one written by
        `matplotlib.figure.Figure.savefig`_.  Use :func:`flush` to ensure the
        deferred writes have completed.

        Parameters
        ----------
        fig : `matplotlib.figure.Figure`_
            Figure to write.  It must not be modified after this call.
        outfile : :obj:`str`, `Path`_, optional
            Output file.  If None, the figure is not written.
        show : :obj:`bool`, optional
            Show the figure interactively.  This forces the synchronous path.
        close : :obj:`bool`, optional
            Close the figure once it has been written.
        **kwargs
            Passed to `matplotlib.figure.Figure.savefig`_ (e.g., ``dpi``).
            Only a plain ``dpi`` is compatible with the deferred path; any
            other keyword, or an ``outfile`` that is not a PNG, is written
            synchronously.
        """
        from matplotlib import pyplot as plt

        if show or self.pool is None:
            if outfile is not None:
                fig.savefig(outfile, **kwargs)
            if show:
                plt.show()
            if close:
                plt.close(fig)
            return
        if outfile is None:
            if close:
                plt.close(fig)
            return
        if Path(outfile).suffix.lower() != '.png' or not close \
                or any(k != 'dpi' for k in kwargs):
            # Only the plain PNG + close case is handled by the deferred path
            fig.savefig(outfile, **kwargs)
            if close:
                plt.close(fig)
            return

        # Rasterize here (see the class warning), then queue only the encode.
        fig.set_dpi(kwargs.get('dpi', fig.dpi))
        fig.canvas.draw()
        rgba = np.asarray(fig.canvas.buffer_rgba()).copy()
        # Deregister the figure from pyplot now, so that pyplot-state plotting
        # done by a subsequent QA function cannot draw into a figure that is
        # still queued.
        plt.close(fig)
        self.pending.append(self.pool.submit(self._encode_png, rgba, outfile, fig.dpi))
        if len(self.pending) >= self.max_pending:
            self.flush()

    def flush(self):
        """
        Block until every deferred figure has been written.

        Exceptions raised while encoding are re-raised here, on the calling
        thread.  This is a no-op when writing serially.
        """
        pending, self.pending = self.pending, []
        for future in pending:
            future.result()

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
        from PIL import Image
        Image.fromarray(rgba).save(outfile, dpi=(dpi, dpi))

    def __repr__(self):
        if self.pool is None:
            return f'<{self.__class__.__name__}: serial>'
        return f'<{self.__class__.__name__}: {self.pool._max_workers} encoding threads, ' \
               f'{len(self.pending)} pending>'
