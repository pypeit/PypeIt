"""
Instrument-independent up-the-ramp fitting.

Detectors read out non-destructively "up the ramp" record a series of reads
of the accumulating charge rather than a single integration.  Fitting the
per-pixel count rate to those reads (with cosmic-ray/jump detection) recovers
the flux while rejecting jumps and beating down read noise.  This module holds
the parts of that machinery that do not depend on any particular instrument:
the count-rate fit, the single-read-noise calibration, the effective-read-noise
formula, the worker-thread dispatch, and the on-disk sidecar convention for
caching a fitted 2D count-rate image next to a reduction.

The numerical fitting is delegated to the ``fitramp`` algorithm of Brandt
(2024, https://arxiv.org/abs/2404.01326; reference implementation:
https://github.com/t-brandt/fitramp), vendored in :mod:`pypeit.ext.fitramp`.

Instruments plug in through the ramp hooks on
:class:`~pypeit.spectrographs.spectrograph.Spectrograph` (``_load_ramp`` and
``_count_reads``), which read an instrument's raw cube into the reads,
covariance, and header this module operates on; MMT/MMIRS
(:mod:`pypeit.spectrographs.mmt_mmirs`) is the first such instrument.

.. include:: ../include/links.rst
"""
import os
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from astropy.io import fits

from pypeit.ext.fitramp import fitramp


def ramp_diffs(reads, covar):
    """
    Compute scaled resultant differences from ramp reads.

    Parameters
    ----------
    reads : `numpy.ndarray`_
        Reads in time order, shape ``(ngroups, ny, nx)``, in electrons.
    covar : :class:`~pypeit.ext.fitramp.fitramp.Covar`
        Covariance object providing the time intervals ``delta_t``.

    Returns
    -------
    diffs : `numpy.ndarray`_
        ``(reads[i+1] - reads[i]) / covar.delta_t[i]``, shape
        ``(ngroups-1, ny, nx)``, in e-/s.
    """
    return np.diff(reads, axis=0) / np.asarray(covar.delta_t)[:, None, None]


def effective_ronoise(sig, ngroups):
    """
    Effective read noise of an up-the-ramp-fitted image.

    For ``N`` uniformly spaced reads with single-read noise ``sig``, the
    read-noise contribution to the total-count uncertainty of the fitted
    slope is ``sig * sqrt(12 (N-1) / (N (N+1)))`` (Brandt 2024a).

    Parameters
    ----------
    sig : :obj:`float`
        Single-read noise in electrons.
    ngroups : :obj:`int`
        Number of reads in the ramp.

    Returns
    -------
    :obj:`float`
        Effective read noise in electrons.
    """
    return sig * np.sqrt(12. * (ngroups - 1) / (ngroups * (ngroups + 1)))


def _resolve_workers(workers):
    """
    Resolve the requested number of ramp-fit worker threads.

    ``None`` selects ``min(6, os.cpu_count())``; any explicit value is passed
    through (clamped to at least 1).  The per-pixel fit is memory-bandwidth
    bound, so more than ~6 threads does not help and can hurt.
    """
    if workers is None:
        return max(1, min(6, os.cpu_count() or 1))
    return max(1, int(workers))


def _fit_rows(diffs, nb, workers, worker):
    """
    Dispatch ``worker((row_start, row_stop))`` over blocks of detector rows.

    The per-pixel ramp fit is independent, so the ``ny`` rows of ``diffs`` are
    split into contiguous blocks of ``nb`` rows and each block is handed to
    ``worker``, which is expected to write its results into a preallocated
    output array.  With ``workers > 1`` the blocks are fit concurrently in a
    thread pool; NumPy releases the GIL during the element-wise arithmetic that
    dominates :func:`~pypeit.ext.fitramp.fitramp.fit_ramps`, so threads scale
    until the memory bus saturates (empirically ~3x at 6 threads).

    Parameters
    ----------
    diffs : `numpy.ndarray`_
        Scaled resultant differences, shape ``(ndiffs, ny, nx)``.
    nb : :obj:`int`
        Number of rows per block.
    workers : :obj:`int`
        Number of worker threads (already resolved; ``1`` runs serially).
    worker : callable
        Called with a ``(row_start, row_stop)`` tuple for each block.
    """
    ny = diffs.shape[1]
    ranges = [(r0, min(r0 + nb, ny)) for r0 in range(0, ny, nb)]
    if workers <= 1 or len(ranges) == 1:
        for rng in ranges:
            worker(rng)
        return
    with ThreadPoolExecutor(max_workers=workers) as ex:
        # Consume the iterator so exceptions in workers propagate.
        list(ex.map(worker, ranges))


def calibrate_sigma(diffs, covar, sig_guess=9.0, nrows=200, workers=None,
                    nb=16, return_err=False, n_boot=200, seed=1234):
    """
    Calibrate the single-read noise from ramp differences.

    Fits a subsample of rows without jump detection (two-pass, with the
    count-rate guess clamped to non-negative values to debias the second
    pass) and rescales ``sig_guess`` so that the median chi-squared matches
    the expected degrees of freedom (``ngroups - 2``).

    Parameters
    ----------
    diffs : `numpy.ndarray`_
        Scaled resultant differences, shape ``(ndiffs, ny, nx)``, electrons.
    covar : :class:`~pypeit.ext.fitramp.fitramp.Covar`
        Covariance object matching ``diffs``.
    sig_guess : :obj:`float`, optional
        Initial guess for the single-read noise in electrons.
    nrows : :obj:`int`, optional
        Number of evenly spaced rows (from the central 80% of the detector)
        to include in the calibration.
    workers : :obj:`int`, optional
        Number of worker threads; ``None`` selects ``min(6, os.cpu_count())``
        and ``1`` disables threading.
    nb : :obj:`int`, optional
        Number of subsampled rows fit per :func:`~pypeit.ext.fitramp.fitramp.fit_ramps`
        call.
    return_err : :obj:`bool`, optional
        If True, also return a bootstrap estimate of the uncertainty on the
        calibrated noise (see ``n_boot``), for inverse-variance weighting when
        combining multiple darks.
    n_boot : :obj:`int`, optional
        Number of bootstrap resamples of the per-pixel chi-squared ensemble
        used to estimate the uncertainty when ``return_err`` is True.
    seed : :obj:`int`, optional
        Seed for the bootstrap resampling, so the uncertainty is deterministic.

    Returns
    -------
    :obj:`float` or :obj:`tuple`
        Calibrated single-read noise in electrons (unclamped).  If
        ``return_err`` is True, a ``(sigma, sigma_err)`` tuple is returned
        instead, with ``sigma_err`` the bootstrap standard deviation.
    """
    ndiffs, ny, nx = diffs.shape
    margin = int(ny * 0.10)
    row_candidates = np.arange(margin, ny - margin)
    nrows = min(nrows, len(row_candidates))
    indices = np.linspace(0, len(row_candidates) - 1, nrows, dtype=int)
    # Gather the subsampled rows; the per-pixel fit is independent, so they can
    # be fit in blocks with rows folded into the pixel axis.
    sub = np.ascontiguousarray(diffs[:, row_candidates[indices], :])
    chisq = np.empty((nrows, nx), dtype=np.float64)

    def worker(rng):
        r0, r1 = rng
        flat = sub[:, r0:r1, :].reshape(ndiffs, (r1 - r0) * nx)
        sig_row = np.full(flat.shape[1], sig_guess, dtype=np.float64)
        with np.errstate(divide='ignore', invalid='ignore'):
            result = fitramp.fit_ramps(flat, covar, sig_row)
            guess = result.countrate * (result.countrate > 0)
            result = fitramp.fit_ramps(flat, covar, sig_row,
                                       countrateguess=guess)
        chisq[r0:r1] = result.chisq.reshape(r1 - r0, nx)

    _fit_rows(sub, nb, _resolve_workers(workers), worker)
    expected_chisq = float(ndiffs - 1)
    flat_chisq = chisq.ravel()
    sigma = sig_guess * np.sqrt(float(np.median(flat_chisq)) / expected_chisq)
    if not return_err:
        return sigma
    # Bootstrap the per-pixel chi-squared ensemble to estimate the uncertainty
    # on the calibrated noise.  Looped (rather than a (n_boot, Npix) array) to
    # avoid a large allocation for full-detector calibrations.
    npix = flat_chisq.size
    rng = np.random.default_rng(seed)
    boot = np.empty(n_boot, dtype=np.float64)
    for b in range(n_boot):
        m_b = np.median(flat_chisq[rng.integers(0, npix, npix)])
        boot[b] = sig_guess * np.sqrt(m_b / expected_chisq)
    return sigma, float(np.std(boot))


def fit_ramp(diffs, covar, sig, workers=None, nb=16):
    """
    Fit all pixels of a ramp, in blocks of rows, with jump detection.

    The per-pixel fit is independent, so the detector is fit in blocks of
    ``nb`` rows (rows folded into the pixel axis of
    :func:`~pypeit.ext.fitramp.fitramp.fit_ramps`) and, for ``workers > 1``, the
    blocks are fit concurrently.  Results are numerically identical to a
    row-by-row fit.

    Parameters
    ----------
    diffs : `numpy.ndarray`_
        Scaled resultant differences, shape ``(ndiffs, ny, nx)``, electrons.
    covar : :class:`~pypeit.ext.fitramp.fitramp.Covar`
        Covariance object matching ``diffs``.
    sig : :obj:`float`
        Single-read noise in electrons.
    workers : :obj:`int`, optional
        Number of worker threads; ``None`` selects ``min(6, os.cpu_count())``
        and ``1`` disables threading.
    nb : :obj:`int`, optional
        Number of rows fit per :func:`~pypeit.ext.fitramp.fitramp.fit_ramps` call.

    Returns
    -------
    countrate : `numpy.ndarray`_
        Fitted count rates in e-/s, shape ``(ny, nx)``.
    """
    ndiffs, ny, nx = diffs.shape
    countrate = np.empty((ny, nx), dtype=np.float64)

    def worker(rng):
        r0, r1 = rng
        flat = np.ascontiguousarray(diffs[:, r0:r1, :]).reshape(ndiffs,
                                                                (r1 - r0) * nx)
        sig_row = np.full(flat.shape[1], sig, dtype=np.float64)
        with np.errstate(divide='ignore', invalid='ignore'):
            diffs2use, guess = fitramp.mask_jumps(flat, covar, sig_row)
            result = fitramp.fit_ramps(flat, covar, sig_row,
                                       diffs2use=diffs2use,
                                       countrateguess=guess * (guess > 0))
        countrate[r0:r1] = result.countrate.reshape(r1 - r0, nx)

    _fit_rows(diffs, nb, _resolve_workers(workers), worker)
    return countrate


def rampfit_path(raw_file, redux_path, rampfit_dir='RampFit'):
    """
    Return the preprocessed-image path for a raw up-the-ramp cube.

    Preprocessed 2D count-rate images live in the ramp-fit directory
    inside the reduction directory (alongside ``Calibrations``,
    ``Science``, etc.), named after the raw cube with a ``_rampfit`` marker
    inserted before the extension so the processed product is distinct from
    the raw file.  The directory name is set by the ``[rdx] rampfit_dir``
    parameter.

    Parameters
    ----------
    raw_file : :obj:`str`, `Path`_
        Path to the raw cube.
    redux_path : :obj:`str`, `Path`_
        Path to the reduction directory.
    rampfit_dir : :obj:`str`, optional
        Name of the ramp-fit subdirectory, relative to ``redux_path``.

    Returns
    -------
    `Path`_
        ``<redux_path>/<rampfit_dir>/<raw stem>_rampfit<ext>``
    """
    raw = Path(raw_file)
    return Path(redux_path) / rampfit_dir / f'{raw.stem}_rampfit{raw.suffix}'


def rampfit_fresh(rampfit_file, raw_file):
    """
    Check whether a preprocessed image exists and is up to date.

    A preprocessed image is fresh when its ``RAWMTIME`` header card matches
    the raw cube's current modification time to within 1 second.  Missing
    or unreadable files (or files without the card) are not fresh.

    Parameters
    ----------
    rampfit_file : :obj:`str`, `Path`_
        Path to the candidate preprocessed image.
    raw_file : :obj:`str`, `Path`_
        Path to the source raw cube.

    Returns
    -------
    :obj:`bool`
        True if the preprocessed image can be used in place of the cube.
    """
    rampfit_file = Path(rampfit_file)
    if not rampfit_file.exists():
        return False
    # Reading the sidecar header is the only I/O here: a corrupt/unreadable
    # file is not fresh.
    try:
        header = fits.getheader(rampfit_file)
    except OSError:
        return False
    # A file without the freshness card is not a usable preprocessed image.
    if 'RAWMTIME' not in header:
        return False
    # If the sidecar records which raw cube it came from, require it to match:
    # a same-named raw file from a different directory (raw cube names are only
    # unique within a program for some instruments) must never reuse this
    # image.  Sidecars written by older versions carry no RAWPATH and fall back
    # to the mtime check.
    raw_path = header.get('RAWPATH')
    if raw_path is not None \
            and str(raw_path) != str(Path(raw_file).resolve()):
        return False
    # A sidecar that outlives its raw source cannot be checked, so it is not
    # fresh.
    try:
        raw_mtime = Path(raw_file).stat().st_mtime
    except OSError:
        return False
    return abs(float(header['RAWMTIME']) - raw_mtime) < 1.


def write_rampfit(rampfit_file, rate, hdu, sig, eff_ronoise, ngroups, raw_mtime,
                  raw_file=None):
    """
    Write a preprocessed 2D count-rate image.

    The output carries a copy of the raw primary header (``hdu[0]``) plus the
    cards ``RAMPFIT`` (marker), ``RAMPSIG``, ``RAMPRON``, ``NGROUPS``,
    ``RAWMTIME``, and ``RAWPATH`` (the resolved path of the source raw cube),
    and a single image extension holding the fitted count rate in e-/s
    (float32) under a copy of the raw first-extension header (``hdu[1]``), so
    all metadata used by ``pypeit_setup`` is preserved.  ``RAWPATH`` lets the
    freshness check reject a sidecar that a same-named raw cube from a
    *different* directory would otherwise map onto (some instruments' raw file
    names are only unique within a program, not globally).

    The file is written atomically: the FITS data are first written to a
    temporary file in the same directory, which is then renamed onto the
    final path.  This ensures that a crash or full disk mid-write can never
    leave a truncated sidecar whose header (and hence its ``RAWMTIME``
    freshness check) is already flushed, which would otherwise be treated
    as fresh forever while being unreadable.

    Parameters
    ----------
    rampfit_file : :obj:`str`, `Path`_
        Output path; its parent directory is created if needed.
    rate : `numpy.ndarray`_
        Fitted count rate in e-/s, trimmed to the data section.
    hdu : `astropy.io.fits.HDUList`_
        Opened source raw cube (headers are copied from ``hdu[0]`` and
        ``hdu[1]``).
    sig : :obj:`float`
        Single-read noise used in the fit (electrons).
    eff_ronoise : :obj:`float`
        Effective read noise of the fitted image (electrons).
    ngroups : :obj:`int`
        Number of reads in the source ramp (recorded in the ``NGROUPS`` card).
    raw_mtime : :obj:`float`
        Modification time of the source raw cube.
    raw_file : :obj:`str`, `Path`_, optional
        Source raw cube; when given, its resolved path is recorded in the
        ``RAWPATH`` card for the freshness check.

    Raises
    ------
    OSError
        If the output directory cannot be created or the file cannot be
        written.
    """
    rampfit_file = Path(rampfit_file)
    prihead = hdu[0].header.copy()
    prihead['RAMPFIT'] = (True, 'PypeIt up-the-ramp preprocessed image')
    prihead['RAMPSIG'] = (float(sig), 'Single-read noise used in the fit (e-)')
    prihead['RAMPRON'] = (float(eff_ronoise), 'Effective read noise (e-)')
    prihead['NGROUPS'] = (int(ngroups), 'Number of reads in the source ramp')
    prihead['RAWMTIME'] = (float(raw_mtime),
                           'Modification time of the source raw cube')
    if raw_file is not None:
        prihead['RAWPATH'] = (str(Path(raw_file).resolve()),
                              'Resolved path of the source raw cube')
    head1 = hdu[1].header.copy()
    head1['DATASEC'] = f'[1:{rate.shape[0]},1:{rate.shape[1]}]'
    head1['BUNIT'] = 'e-/s'
    out = fits.HDUList([fits.PrimaryHDU(header=prihead),
                        fits.ImageHDU(data=rate.astype(np.float32),
                                      header=head1)])
    rampfit_file.parent.mkdir(parents=True, exist_ok=True)
    tmp_file = rampfit_file.with_name(rampfit_file.name + f'.tmp{os.getpid()}')
    try:
        out.writeto(tmp_file, overwrite=True)
        tmp_file.replace(rampfit_file)
    finally:
        tmp_file.unlink(missing_ok=True)
