"""
Mixin adding up-the-ramp fitting to a spectrograph.

Detectors read out non-destructively "up the ramp" record a series of reads of
the accumulating charge rather than a single integration.  A spectrograph read
out this way mixes in :class:`RampSpectrograph` alongside
:class:`~pypeit.spectrographs.spectrograph.Spectrograph`::

    class MyRampSpectrograph(RampSpectrograph, Spectrograph):
        ...

and implements the two instrument hooks :func:`RampSpectrograph._load_ramp` and
:func:`RampSpectrograph._count_reads`; the mixin then provides the generic
orchestration (noise calibration, count-rate fit, sidecar caching), delegating
the instrument-independent numerics to :mod:`pypeit.core.ramp`.  MMT/MMIRS
(:mod:`pypeit.spectrographs.mmt_mmirs`) is the first such instrument.

The base :class:`~pypeit.spectrographs.spectrograph.Spectrograph` carries
``is_up_the_ramp = False`` and a :func:`preprocess_ramp_file` stub that raises;
mixing in this class flips ``is_up_the_ramp`` to True and supplies the working
implementation, so ``pypeit_fit_ramp`` and the reduction detect ramp support by
the ``is_up_the_ramp`` flag alone -- there is no registry to keep in sync.

.. include:: ../include/links.rst
"""
from pathlib import Path

import numpy as np

from pypeit import log
from pypeit import PypeItError
from pypeit import io
from pypeit.core import ramp


class RampSpectrograph:
    """
    Mixin providing up-the-ramp fitting to a
    :class:`~pypeit.spectrographs.spectrograph.Spectrograph`.

    Not a spectrograph itself (it is not a ``Spectrograph`` subclass, so it is
    never discovered as one); it must be combined with ``Spectrograph`` and
    relies on the host for ``name`` and ``get_detector_par``.  Subclasses must
    implement :func:`_load_ramp` and :func:`_count_reads`.
    """

    #: Flags a spectrograph read out up-the-ramp; see
    #: :func:`~pypeit.spectrographs.spectrograph.Spectrograph.is_up_the_ramp`.
    is_up_the_ramp = True

    # Up-the-ramp fitting configuration.  Instruments override the calibrated
    # values (e.g. ramp_sig_guess/ramp_sig_range) as needed.  See
    # :mod:`pypeit.core.ramp`.
    ramp_sig_guess = 9.0
    """Initial guess for the *effective* single-read noise in electrons (the
    instantaneous read noise plus accumulated dark-current/flux shot noise),
    used to seed the chi-square rescaling and as the fallback when a frame has
    too few reads to calibrate.  Refined from darks or self-calibrated per
    frame; the header ``RDNOISE`` is used as a physical floor (see
    :func:`get_ramp_sigma`)."""
    ramp_sig_range = (3.0, 50.0)
    """Absolute sanity range for the calibrated single-read noise (electrons).
    The lower bound is raised to the header ``RDNOISE`` when available, since a
    derived noise below the instantaneous read noise is unphysical."""
    ramp_min_reads = 5
    """Minimum number of reads for up-the-ramp fitting.  Frames with fewer
    reads do not sample the ramp well enough to fit reliably and fall back to
    correlated double sampling."""
    ramp_min_cal_groups = 10
    """Minimum number of reads for a frame (dark or science) to calibrate the
    single-read noise.  A dark with fewer reads is not used; a science frame
    with fewer reads falls back to the published guess
    (:attr:`ramp_sig_guess`) instead of self-calibrating."""
    ramp_fit_workers = None
    """Number of worker threads for up-the-ramp fitting.  ``None`` selects
    ``min(6, os.cpu_count())``; set to ``1`` to disable threading.  The
    per-pixel fit is memory-bandwidth bound, so throughput plateaus at roughly
    6 threads (measured ~3x over the serial fit on a 10-core machine)."""
    ramp_fit_chunk_rows = 16
    """Number of detector rows fit per
    :func:`~pypeit.ext.fitramp.fitramp.fit_ramps` call.  Small chunks keep each
    thread's working set cache-resident; 16 was the empirical sweet spot."""
    _ramp_fitstbl = None
    _ramp_sigma = None
    _ramp_sigma_cache = None
    _ramp_output_dir = None
    _rampfit_dir = 'RampFit'
    """
    str: Name of the subdirectory (relative to the reduction directory) where
    preprocessed up-the-ramp count-rate images are written and reused.  The
    default is overridden from the ``[rdx] rampfit_dir`` parameter in
    :func:`cache_metadata`.
    """

    def cache_metadata(self, fitstbl):
        """
        Record the reduction directory and dark frames for later ramp fitting.

        Overrides the base
        :func:`~pypeit.spectrographs.spectrograph.Spectrograph.cache_metadata`
        no-op.  The reduction directory determines where preprocessed ramp
        images are written (its ``[rdx] rampfit_dir`` subdirectory).  The
        metadata table is kept so the dark frames used to calibrate the
        single-read noise can be looked up lazily by :func:`_ramp_dark_sigmas`
        (via :func:`~pypeit.metadata.PypeItMetaData.find_frame_files`) once
        frame types are assigned -- they are not yet set when this hook runs
        during ``PypeItMetaData`` construction.  Cheap and idempotent, as
        required of the hook.

        Args:
            fitstbl (:class:`~pypeit.metadata.PypeItMetaData`):
                The class holding the metadata for all the frames.
        """
        self._ramp_output_dir = Path(fitstbl.par['rdx']['redux_path'])
        self._rampfit_dir = fitstbl.par['rdx']['rampfit_dir']
        # Let the user override the ramp-fit threading/chunking from the
        # [rdx] block of the pypeit file; unset (None) keeps the class default.
        if fitstbl.par['rdx']['ramp_fit_cores'] is not None:
            self.ramp_fit_workers = fitstbl.par['rdx']['ramp_fit_cores']
        if fitstbl.par['rdx']['ramp_fit_chunk_rows'] is not None:
            self.ramp_fit_chunk_rows = fitstbl.par['rdx']['ramp_fit_chunk_rows']
        self._ramp_fitstbl = fitstbl

    def rampfit_path(self, raw_file):
        """
        Path to the cached ramp-fit image for a raw cube during a reduction.

        Resolves the ``[rdx] rampfit_dir`` subdirectory of the reduction
        directory recorded by :func:`cache_metadata`.  When that hook never
        fired (e.g. direct API use of :func:`get_rawimage` outside a
        reduction), the current working directory is used instead.  This is
        the reduction-time counterpart of :func:`preprocess_ramp_file`, which
        takes the reduction directory explicitly from the ``pypeit_fit_ramp``
        script.

        Args:
            raw_file (:obj:`str`, `Path`_):
                Path to the raw up-the-ramp cube.

        Returns:
            `Path`_: Path to the (possibly not-yet-written) preprocessed image.
        """
        redux_path = self._ramp_output_dir if self._ramp_output_dir is not None \
            else Path.cwd()
        return ramp.rampfit_path(raw_file, redux_path, self._rampfit_dir)

    def _load_ramp(self, hdu, detector_par):
        """
        Load the non-destructive reads of a raw up-the-ramp cube.

        Instrument-specific hook backing the ramp orchestration
        (:func:`_ramp_fit_image`, :func:`_ramp_dark_sigmas`); subclasses must
        implement it.  Implementations read the raw cube into the reads (in
        electrons, so any gain and reference-pixel correction is applied here),
        build the covariance from the read times, and return the header
        carrying the frame metadata (``EXPTIME``, ``RDNOISE``, ...).

        Args:
            hdu (`astropy.io.fits.HDUList`_):
                Opened raw up-the-ramp cube.
            detector_par (:class:`~pypeit.images.detector_container.DetectorContainer`):
                Detector parameters (e.g. the gain applied to convert the
                reads to electrons).

        Returns:
            :obj:`tuple`: The reads in electrons (`numpy.ndarray`_, shape
            ``(ngroups, ny, nx)``), the matching
            :class:`~pypeit.ext.fitramp.fitramp.Covar`, and the metadata
            `astropy.io.fits.Header`_.
        """
        raise NotImplementedError(
            f'{self.name} must implement _load_ramp to support up-the-ramp '
            'fitting.')

    def _count_reads(self, hdu):
        """
        Count the non-destructive reads in a raw up-the-ramp cube.

        Instrument-specific hook; subclasses must implement it.

        Args:
            hdu (`astropy.io.fits.HDUList`_):
                Opened raw up-the-ramp cube.

        Returns:
            :obj:`int`: Number of reads in the ramp.
        """
        raise NotImplementedError(
            f'{self.name} must implement _count_reads to support up-the-ramp '
            'fitting.')

    def _ramp_dark_sigmas(self, exptime=None):
        """
        Calibrate the single-read noise and its uncertainty from every
        recorded dark frame that matches the science ramp depth.

        The darks are the ``dark`` frames of the cached metadata table (via
        :func:`~pypeit.metadata.PypeItMetaData.find_frame_files`); each with at
        least :attr:`ramp_min_cal_groups` reads is fit independently.  When
        ``exptime`` is provided, only darks whose own ``EXPTIME`` matches it
        (within ``rtol=1e-3``) are used, because the calibrated noise is the
        effective per-read noise -- instantaneous read noise plus accumulated
        dark-current/flux shot noise -- which grows with ramp length, so it
        must be measured at the science exposure time.  Darks that cannot be
        opened, or that yield a non-finite noise or a non-positive uncertainty,
        are skipped.

        Args:
            exptime (:obj:`float`, optional):
                Science-frame exposure time in seconds.  If given, darks with a
                different ``EXPTIME`` are excluded from the calibration.

        Returns:
            :obj:`list`: List of ``(name, sigma, sigma_err)`` tuples (in
            electrons) for the qualifying darks; empty if none qualify.
        """
        # find_frame_files needs the frame types, which are assigned after the
        # cache_metadata hook runs; before then (e.g. during setup) there is
        # nothing to calibrate against.
        fitstbl = self._ramp_fitstbl
        if fitstbl is None or 'framebit' not in fitstbl.keys():
            return []
        results = []
        for f in fitstbl.find_frame_files('dark'):
            try:
                with io.fits_open(f) as dhdu:
                    if self._count_reads(dhdu) < self.ramp_min_cal_groups:
                        continue
                    detpar = self.get_detector_par(1, hdu=dhdu)
                    dreads, dcovar, dhead = self._load_ramp(dhdu, detpar)
            except (OSError, PypeItError):
                log.warning(f'Could not open recorded dark frame {f}; skipping '
                            'it for read-noise calibration.')
                continue
            if exptime is not None and not np.isclose(dhead.get('EXPTIME', np.nan),
                                                      exptime, rtol=1e-3):
                continue
            ddiffs = ramp.ramp_diffs(dreads, dcovar)
            sig, err = ramp.calibrate_sigma(ddiffs, dcovar,
                                            sig_guess=self.ramp_sig_guess,
                                            workers=self.ramp_fit_workers,
                                            nb=self.ramp_fit_chunk_rows,
                                            return_err=True)
            if np.isfinite(sig) and np.isfinite(err) and err > 0:
                results.append((Path(f).name, float(sig), float(err)))
            else:
                log.warning(f'Dark {Path(f).name} gave an unusable read-noise '
                            f'calibration (sigma={sig}, err={err}); skipping it.')
        return results

    def get_ramp_sigma(self, diffs, covar, exptime=None, ron_floor=None):
        """
        Determine the single-read noise for up-the-ramp fitting.

        Preferentially calibrates the noise from the dark frames recorded by
        :func:`cache_metadata` that match the science ramp depth: every dark
        with the same ``EXPTIME`` as the science frame (see ``exptime``) and at
        least :attr:`ramp_min_cal_groups` reads is calibrated independently and
        the results are combined as an inverse-variance weighted mean, weighting
        each dark by the (bootstrap) uncertainty on its own calibrated noise
        (the result is cached for the rest of the run).  Matching the exposure
        time matters because the calibrated value is the effective per-read
        noise, which includes accumulated dark-current/flux shot noise and so
        grows with ramp length.  If no matching dark is available and the frame
        has at least :attr:`ramp_min_cal_groups` reads, self-calibrates from the
        provided ramp differences (not cached); with fewer reads the fit is too
        poorly constrained, so the published guess (:attr:`ramp_sig_guess`) is
        used.

        A calibrated result is clamped to :attr:`ramp_sig_range`, with the lower
        bound raised to ``ron_floor`` (the header ``RDNOISE``, the instantaneous
        read noise) when provided: a derived noise below the instantaneous read
        noise is unphysical, so such values are floored.

        Args:
            diffs (`numpy.ndarray`_):
                Scaled resultant differences of the frame being processed,
                shape ``(ndiffs, ny, nx)``, in electrons.
            covar (:class:`~pypeit.ext.fitramp.fitramp.Covar`):
                Covariance object matching ``diffs``.
            exptime (:obj:`float`, optional):
                Science-frame exposure time in seconds, used to select darks of
                the same ramp depth.  If ``None``, darks are not filtered by
                exposure time.
            ron_floor (:obj:`float`, optional):
                Instantaneous read noise in electrons (header ``RDNOISE``) used
                as the physical lower bound on the derived noise.  If ``None``,
                the lower bound of :attr:`ramp_sig_range` is used.

        Returns:
            :obj:`float`: Single-read noise in electrons.
        """
        # An explicitly forced value (``_ramp_sigma`` set on the instance) is
        # global and always wins.
        if self._ramp_sigma is not None:
            return self._ramp_sigma
        # Automatic dark calibration is cached per (exptime, ron_floor): the
        # effective per-read noise grows with ramp length, so science and
        # standard frames of different EXPTIME (or with different RDNOISE
        # floors) must not share a single cached value.
        if self._ramp_sigma_cache is None:
            self._ramp_sigma_cache = {}
        cache_key = (exptime, ron_floor)
        if cache_key in self._ramp_sigma_cache:
            return self._ramp_sigma_cache[cache_key]
        lo = self.ramp_sig_range[0] if ron_floor is None \
            else max(self.ramp_sig_range[0], float(ron_floor))
        hi = self.ramp_sig_range[1]
        darks = self._ramp_dark_sigmas(exptime=exptime)
        if len(darks) > 0:
            names = [d[0] for d in darks]
            sigs = np.array([d[1] for d in darks])
            errs = np.array([d[2] for d in darks])
            log.info(f'Calibrating single-read noise from {len(darks)} '
                     f'dark(s): {", ".join(names)}')
            weights = 1.0 / errs ** 2
            sig = float(np.sum(weights * sigs) / np.sum(weights))
            # Report the larger of the inverse-variance (within-dark) error and
            # the between-dark standard error of the mean, so real dark-to-dark
            # scatter is not hidden by tiny per-dark bootstrap uncertainties.
            ivar_err = float(np.sqrt(1.0 / np.sum(weights)))
            sem = float(np.std(sigs, ddof=1) / np.sqrt(len(sigs))) \
                if len(sigs) > 1 else 0.0
            comb_err = max(ivar_err, sem)
            result = float(np.clip(sig, lo, hi))
            self._ramp_sigma_cache[cache_key] = result
            log.info(f'Calibrated single-read noise: {result:.2f} '
                     f'+/- {comb_err:.2f} e-')
            return result
        ngroups = diffs.shape[0] + 1
        if ngroups < self.ramp_min_cal_groups:
            guess = float(max(lo, self.ramp_sig_guess))
            log.info(f'No suitable dark listed and only {ngroups} reads '
                     f'(< {self.ramp_min_cal_groups}); using the guess '
                     f'single-read noise of {guess:.2f} e-')
            return guess
        log.info('No suitable dark listed; self-calibrating single-read '
                 'noise from the frame itself')
        sig = ramp.calibrate_sigma(diffs, covar, sig_guess=self.ramp_sig_guess,
                                   workers=self.ramp_fit_workers,
                                   nb=self.ramp_fit_chunk_rows)
        sig = float(np.clip(sig, lo, hi))
        log.info(f'Self-calibrated single-read noise: {sig:.2f} e-')
        return sig

    def _ramp_fit_image(self, hdu, detector_par):
        """
        Perform up-the-ramp fitting of a multi-read frame.

        Instrument-independent orchestration: the instrument hook
        :func:`_load_ramp` provides the reads (in electrons), the covariance,
        and the metadata header; the count rate, single-read noise, and
        effective read noise are then computed with :mod:`pypeit.core.ramp`.

        Args:
            hdu (`astropy.io.fits.HDUList`_):
                Opened raw cube with at least :attr:`ramp_min_reads`
                non-destructive reads.
            detector_par (:class:`~pypeit.images.detector_container.DetectorContainer`):
                Detector parameters, passed through to :func:`_load_ramp`.

        Returns:
            :obj:`tuple`: The fitted count-rate image in e-/s
            (`numpy.ndarray`_, shape of the trimmed data section), the
            single-read noise in electrons (:obj:`float`), and the
            effective read noise in electrons (:obj:`float`).
        """
        reads, covar, header = self._load_ramp(hdu, detector_par)
        ngroups = reads.shape[0]
        diffs = ramp.ramp_diffs(reads, covar)
        del reads
        sig = self.get_ramp_sigma(diffs, covar, exptime=header.get('EXPTIME'),
                                  ron_floor=header.get('RDNOISE'))
        log.info(f'Up-the-ramp fitting {ngroups} reads '
                 f'(single-read noise {sig:.2f} e-)')
        countrate = ramp.fit_ramp(diffs, covar, sig,
                                  workers=self.ramp_fit_workers,
                                  nb=self.ramp_fit_chunk_rows)
        eff_ronoise = ramp.effective_ronoise(sig, ngroups)
        log.info(f'Effective read noise: {eff_ronoise:.2f} e-')
        return countrate, sig, eff_ronoise

    def preprocess_ramp_file(self, raw_file, redux_path, rampfit_dir='RampFit',
                             force=False):
        """
        Fit one raw up-the-ramp cube and cache its 2D count-rate image.

        Overrides the base
        :func:`~pypeit.spectrographs.spectrograph.Spectrograph.preprocess_ramp_file`
        stub.  Backs the ``pypeit_fit_ramp`` script and mirrors the fit that
        ``get_rawimage`` performs during a reduction, writing the preprocessed
        image where the reduction will find and reuse it.  Frames that are
        already up to date, already preprocessed, or too short to fit are
        skipped.  The orchestration is instrument-independent; the raw cube is
        read through the :func:`_load_ramp` / :func:`_count_reads` hooks.

        Args:
            raw_file (:obj:`str`, `Path`_):
                Path to the raw up-the-ramp cube to fit.
            redux_path (:obj:`str`, `Path`_):
                The reduction directory holding the ramp-fit subdirectory.
            rampfit_dir (:obj:`str`, optional):
                Name of the ramp-fit subdirectory, relative to ``redux_path``
                (the ``[rdx] rampfit_dir`` parameter).
            force (:obj:`bool`, optional):
                Re-fit and overwrite an existing, up-to-date preprocessed
                image instead of skipping it.

        Returns:
            `Path`_: The path to the preprocessed image, or None if the frame
            was skipped (already up to date, already preprocessed, or too few
            reads).
        """
        raw = Path(raw_file)
        rampfit_file = ramp.rampfit_path(raw, redux_path, rampfit_dir)
        if not force and ramp.rampfit_fresh(rampfit_file, raw):
            log.info(f'{raw.name}: up-to-date preprocessed image exists; '
                     'skipping (use force=True to re-fit)')
            return None
        with io.fits_open(raw) as hdu:
            if hdu[0].header.get('RAMPFIT') is not None:
                log.warning(f'{raw.name} is already a preprocessed image; '
                            'skipping')
                return None
            n_reads = self._count_reads(hdu)
            if n_reads < self.ramp_min_reads:
                log.info(f'{raw.name}: only {n_reads} read(s); up-the-ramp '
                         f'fitting requires at least {self.ramp_min_reads} '
                         '(the reduction uses correlated double sampling). '
                         'Skipping.')
                return None
            log.info(f'{raw.name}: fitting {n_reads} reads')
            detector_par = self.get_detector_par(1, hdu=hdu)
            rate, sig, eff_ronoise = self._ramp_fit_image(hdu, detector_par)
            ramp.write_rampfit(rampfit_file, rate, hdu, sig, eff_ronoise,
                               n_reads, raw.stat().st_mtime, raw_file=raw)
        log.info(f'{raw.name}: single-read noise {sig:.2f} e-, effective '
                 f'read noise {eff_ronoise:.2f} e- -> {rampfit_file}')
        return rampfit_file
