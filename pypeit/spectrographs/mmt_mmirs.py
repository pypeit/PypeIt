"""
Module for MMT MMIRS

.. include:: ../include/links.rst
"""
import os
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from scipy.signal import savgol_filter

from astropy.table import Table
from astropy.time import Time
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.coordinates import Angle
from astropy import units

from pypeit import log
from pypeit import PypeItError
from pypeit import telescopes
from pypeit import utils
from pypeit import io
from pypeit.core import parse
from pypeit.core import framematch
from pypeit.ext.fitramp import fitramp
from pypeit.images import detector_container
from pypeit.spectrographs import spectrograph
from pypeit.spectrographs.slitmask import SlitMask
from pypeit.spectrographs import mmirs_maskfile
from pypeit.par import parset


class MMTMMIRSSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle MMT/MMIRS specific code
    """
    ndet = 1
    name = 'mmt_mmirs'
    telescope = telescopes.MMTTelescopePar()
    camera = 'MMIRS'
    url = 'https://lweb.cfa.harvard.edu/mmti/mmirs.html'
    header_name = 'mmirs'
    supported = True

    # Up-the-ramp fitting configuration
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
    ``min(6, os.cpu_count())``; set to ``1`` to disable threading.  The per-pixel
    fit is memory-bandwidth bound, so throughput plateaus at roughly 6 threads
    (measured ~3x over the serial fit on a 10-core machine)."""
    ramp_fit_chunk_rows = 16
    """Number of detector rows fit per :func:`~pypeit.ext.fitramp.fitramp.fit_ramps`
    call.  Small chunks keep each thread's working set cache-resident; 16 was
    the empirical sweet spot."""
    _ramp_dark_files = None
    _ramp_sigma = None
    _ramp_output_dir = None
    _ramp_match_dark_exptime = True
    """
    bool: Restrict read-noise-calibration darks to those matching the science
    exposure time.  True for the reduction (darks are auto-discovered and may
    span several exposure times); set False when a dark is supplied explicitly
    (e.g. ``pypeit_fit_ramp --dark``), where the user's choice should be used
    as given.
    """

    nod_min_offset = 1.0
    """
    float: Minimum peak-to-peak along-slit dither offset (arcsec) for a
    sequence to be treated as nodded.  Below this, frames are assumed to be a
    stare and are not paired for A-B subtraction.  Comfortably above typical
    pointing jitter (< 0.5") and below any real MMIRS nod throw (several ").
    """

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=1, card='RA')
        self.meta['dec'] = dict(ext=1, card='DEC')
        self.meta['target'] = dict(ext=1, card='OBJECT')
        self.meta['decker'] = dict(ext=1, card='APERTURE')
        self.meta['dichroic'] = dict(ext=1, card='FILTER')
        self.meta['binning'] = dict(ext=1, card=None, default='1,1')

        self.meta['mjd'] = dict(ext=0, card=None, compound=True)
        self.meta['exptime'] = dict(ext=1, card='EXPTIME')
        self.meta['airmass'] = dict(ext=1, card='AIRMASS')
        # Extras for config and frametyping
        self.meta['dispname'] = dict(ext=1, card='DISPERSE')
        self.meta['idname'] = dict(ext=1, card='IMAGETYP')
        self.meta['instrument'] = dict(ext=1, card='INSTRUME')

        # Dither metadata for automatic A-B nod pairing.  MMIRS has no dither
        # header card, so the along-slit offset is derived (see compound_meta).
        self.meta['dithoff'] = dict(ext=1, card=None, compound=True)
        self.meta['frameno'] = dict(ext=1, card=None, compound=True)
        self.meta['posang'] = dict(ext=1, card='POSANGLE')
        # Labels filled in by get_comb_group; default keeps the columns present
        # even when the user has pre-set comb_id (get_comb_group is skipped).
        self.meta['dithpat'] = dict(ext=1, card=None, default='None')
        self.meta['dithpos'] = dict(ext=1, card=None, default='None')

    def compound_meta(self, headarr, meta_key):
        """
        Methods to generate metadata requiring interpretation of the header
        data, instead of simply reading the value of a header card.

        Args:
            headarr (:obj:`list`):
                List of `astropy.io.fits.Header`_ objects.
            meta_key (:obj:`str`):
                Metadata keyword to construct.

        Returns:
            object: Metadata value read from the header(s).
        """
        # TODO: This should be how we always deal with timeunit = 'isot'. Are
        # we doing that for all the relevant spectrographs?
        if meta_key == 'mjd':
            time = headarr[1]['DATE-OBS']
            ttime = Time(time, format='isot')
            return ttime.mjd
        if meta_key == 'dithoff':
            # Along-slit dither offset in arcsec: projection of the telescope
            # pointing (RA in hours, DEC in deg) minus the catalog target
            # (CAT-RA/CAT-DEC, sexagesimal deg) onto the slit PA (POSANGLE).
            hdr = headarr[1]
            try:
                ra = float(hdr['RA']) * 15.0
                dec = float(hdr['DEC'])
                catra = Angle(hdr['CAT-RA'], unit=units.deg).deg
                catdec = Angle(hdr['CAT-DEC'], unit=units.deg).deg
                pa = float(hdr['POSANGLE'])
            except (KeyError, TypeError, ValueError):
                return 0.0
            dra = (ra - catra) * np.cos(np.radians(dec)) * 3600.0
            ddec = (dec - catdec) * 3600.0
            return dra * np.sin(np.radians(pa)) + ddec * np.cos(np.radians(pa))
        if meta_key == 'frameno':
            # Frame number is the trailing token of the ext-1 FILENAME card
            # (e.g. 'MMIRS/2019.0913/nep.as1_mos.1822' -> 1822).
            fname = headarr[1].get('FILENAME')
            if fname is None:
                return -1
            try:
                return int(str(fname).split('.')[-1])
            except (ValueError, IndexError):
                return -1
        raise PypeItError("Not ready for this compound meta")

    def cache_metadata(self, fitstbl):
        """
        Record the reduction directory and the dark frames in the metadata
        table for later use when up-the-ramp fitting raw frames.

        The reduction directory determines where preprocessed ramp images
        are written (its ``RampFit`` subdirectory); the darks are used to
        calibrate the single-read noise.  Cheap and idempotent; the darks
        are only opened (lazily) by :func:`get_ramp_sigma`.

        Args:
            fitstbl (:class:`~pypeit.metadata.PypeItMetaData`):
                The class holding the metadata for all the frames.
        """
        self._ramp_output_dir = Path(fitstbl.par['rdx']['redux_path'])
        # Let the user override the ramp-fit threading/chunking from the
        # [rdx] block of the pypeit file; unset (None) keeps the class default.
        if fitstbl.par['rdx']['ramp_fit_cores'] is not None:
            self.ramp_fit_workers = fitstbl.par['rdx']['ramp_fit_cores']
        if fitstbl.par['rdx']['ramp_fit_chunk_rows'] is not None:
            self.ramp_fit_chunk_rows = fitstbl.par['rdx']['ramp_fit_chunk_rows']
        tbl = fitstbl.table
        if 'directory' not in tbl.colnames or 'filename' not in tbl.colnames:
            return
        if 'frametype' in tbl.colnames:
            indx = np.array([ft is not None and 'dark' in str(ft)
                             for ft in tbl['frametype']])
        elif 'idname' in tbl.colnames:
            indx = np.array([str(idn).strip().lower() == 'dark'
                             for idn in tbl['idname']])
        else:
            return
        self._ramp_dark_files = [Path(str(d)) / str(f)
                                 for d, f in zip(tbl['directory'][indx],
                                                 tbl['filename'][indx])
                                 if not str(f).startswith('#')]

    def _ramp_dark_sigmas(self, exptime=None):
        """
        Calibrate the single-read noise and its uncertainty from every
        recorded dark frame that matches the science ramp depth.

        Each dark with at least :attr:`ramp_min_cal_groups` reads is fit
        independently.  When ``exptime`` is provided, only darks whose own
        ``EXPTIME`` matches it (within ``rtol=1e-3``) are used, because the
        calibrated noise is the effective per-read noise -- instantaneous read
        noise plus accumulated dark-current/flux shot noise -- which grows with
        ramp length, so it must be measured at the science exposure time.
        Darks that cannot be opened, or that yield a non-finite noise or a
        non-positive uncertainty, are skipped.

        Args:
            exptime (:obj:`float`, optional):
                Science-frame exposure time in seconds.  If given, darks with a
                different ``EXPTIME`` are excluded from the calibration.

        Returns:
            :obj:`list`: List of ``(name, sigma, sigma_err)`` tuples (in
            electrons) for the qualifying darks; empty if none qualify.
        """
        results = []
        for f in self._ramp_dark_files or []:
            try:
                with io.fits_open(f) as dhdu:
                    if mmirs_count_reads(dhdu) < self.ramp_min_cal_groups:
                        continue
                    dreads, dhead = mmirs_load_ramp(dhdu)
            except (OSError, PypeItError):
                log.warning(f'Could not open recorded dark frame {f}; skipping it '
                            'for read-noise calibration.')
                continue
            if exptime is not None and not np.isclose(dhead.get('EXPTIME', np.nan),
                                                      exptime, rtol=1e-3):
                continue
            ngroups = dreads.shape[0]
            dcovar = fitramp.Covar([dhead['GRPTIME'] * (i + 1)
                                    for i in range(ngroups)])
            ddiffs = mmirs_ramp_diffs(dreads * dhead['GAIN'], dcovar)
            sig, err = mmirs_calibrate_sigma(ddiffs, dcovar,
                                             sig_guess=self.ramp_sig_guess,
                                             workers=self.ramp_fit_workers,
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
        if self._ramp_sigma is not None:
            return self._ramp_sigma
        lo = self.ramp_sig_range[0] if ron_floor is None \
            else max(self.ramp_sig_range[0], float(ron_floor))
        hi = self.ramp_sig_range[1]
        match_exptime = exptime if self._ramp_match_dark_exptime else None
        darks = self._ramp_dark_sigmas(exptime=match_exptime)
        if darks:
            names = [d[0] for d in darks]
            sigs = np.array([d[1] for d in darks])
            errs = np.array([d[2] for d in darks])
            log.info(f'Calibrating MMIRS single-read noise from {len(darks)} '
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
            self._ramp_sigma = float(np.clip(sig, lo, hi))
            log.info(f'Calibrated single-read noise: {self._ramp_sigma:.2f} '
                     f'+/- {comb_err:.2f} e-')
            return self._ramp_sigma
        ngroups = diffs.shape[0] + 1
        if ngroups < self.ramp_min_cal_groups:
            guess = float(max(lo, self.ramp_sig_guess))
            log.info(f'No suitable dark listed and only {ngroups} reads '
                     f'(< {self.ramp_min_cal_groups}); using the guess '
                     f'single-read noise of {guess:.2f} e-')
            return guess
        log.info('No suitable dark listed; self-calibrating MMIRS single-read '
                 'noise from the frame itself')
        sig = mmirs_calibrate_sigma(diffs, covar, sig_guess=self.ramp_sig_guess,
                                    workers=self.ramp_fit_workers)
        sig = float(np.clip(sig, lo, hi))
        log.info(f'Self-calibrated single-read noise: {sig:.2f} e-')
        return sig

    def raw_header_cards(self):
        """
        Return additional raw header cards to be propagated in
        downstream output files for configuration identification.

        The list of raw data FITS keywords should be those used to populate
        the :meth:`~pypeit.spectrographs.spectrograph.Spectrograph.configuration_keys`
        or are used in :meth:`~pypeit.spectrographs.spectrograph.Spectrograph.config_specific_par`
        for a particular spectrograph, if different from the name of the
        PypeIt metadata keyword.

        This list is used by :meth:`~pypeit.spectrographs.spectrograph.Spectrograph.subheader_for_spec`
        to include additional FITS keywords in downstream output files.

        Returns:
            :obj:`list`: List of keywords from the raw data files that should
            be propagated in output files.
        """
        return ['DISPERSE']

    def get_detector_par(self, det, hdu=None):
        """
        Return metadata for the selected detector.

        Args:
            det (:obj:`int`):
                1-indexed detector number.
            hdu (`astropy.io.fits.HDUList`_, optional):
                The open fits file with the raw image of interest.  If not
                provided, frame-dependent parameters are set to a default.

        Returns:
            :class:`~pypeit.images.detector_container.DetectorContainer`:
            Object with the detector metadata.
        """
        # Read the instantaneous read noise from the header (RDNOISE, ext 1);
        # it comes from the instrument and is the ground truth.  Fall back to
        # the long-standing value only when no header is available.
        ronoise = 3.14
        if hdu is not None and len(hdu) > 1:
            ronoise = float(hdu[1].header.get('RDNOISE', ronoise))
        # Detector 1
        detector_dict = dict(
            binning='1,1',
            det             = 1,
            dataext         = 1,
            specaxis        = 0,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.2012,
            darkcurr        = 36.0,  # e-/pixel/hour  (=0.01 e-/pixel/s)
            saturation      = 700000., #155400.,
            nonlinear       = 1.0,
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d(0.95),
            ronoise         = np.atleast_1d(ronoise),
            datasec         = np.atleast_1d('[:,:]'),
            oscansec        = None, #np.atleast_1d('[:,:]')
            )
        return detector_container.DetectorContainer(**detector_dict)

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.
        
        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        # Image processing steps
        turn_off = dict(use_illumflat=False, use_biasimage=False, use_overscan=False,
                        use_darkimage=False)
        par.reset_all_processimages_par(**turn_off)
        #par['calibrations']['traceframe']['process']['use_darkimage'] = True
        #par['calibrations']['pixelflatframe']['process']['use_darkimage'] = True
        #par['calibrations']['illumflatframe']['process']['use_darkimage'] = True
        #par['scienceframe']['process']['use_darkimage'] = True
        par['scienceframe']['process']['use_illumflat'] = True

        # Wavelengths
        # 1D wavelength solution with arc lines
        par['calibrations']['wavelengths']['rms_thresh_frac_fwhm'] = 0.125
        par['calibrations']['wavelengths']['sigdetect']=5
        par['calibrations']['wavelengths']['fwhm'] = 4.
        par['calibrations']['wavelengths']['n_first']=2
        par['calibrations']['wavelengths']['n_final']=4
        par['calibrations']['wavelengths']['lamps'] = ['OH_NIRES']
        par['calibrations']['wavelengths']['match_toler']=5.0

        # Set slits and tilts parameters
        par['calibrations']['tilts']['tracethresh'] = 5
        par['calibrations']['tilts']['spat_order'] = 7
        par['calibrations']['tilts']['spec_order'] = 5
        par['calibrations']['slitedges']['trace_thresh'] = 10.
        par['calibrations']['slitedges']['edge_thresh'] = 100.
        par['calibrations']['slitedges']['fit_min_spec_length'] = 0.4
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'
        par['calibrations']['slitedges']['bound_detector'] = True

        # Set the default exposure time ranges for the frame typing
        par['calibrations']['standardframe']['exprng'] = [None, 60]
        par['calibrations']['tiltframe']['exprng'] = [60, None]
        par['calibrations']['arcframe']['exprng'] = [60, None]
        par['calibrations']['darkframe']['exprng'] = [30, None]
        par['scienceframe']['exprng'] = [30, None]

        # dark
        # TODO: This is now the default.
        par['calibrations']['darkframe']['process']['apply_gain'] = True

        # cosmic ray rejection
        par['scienceframe']['process']['sigclip'] = 5.0
        par['scienceframe']['process']['objlim'] = 2.0
        par['scienceframe']['process']['grow'] = 0.5

        # Science reduction
        par['reduce']['findobj']['snr_thresh'] = 5.0
        par['reduce']['skysub']['sky_sigrej'] = 5.0
        par['reduce']['findobj']['find_trim_edge'] = [5,5]
        # Object tracing: bound_detector=True gives synthetic straight slit
        # edges, so the initial object trace is a straight line while the real
        # trace is inclined; a linear fit with loose rejection lets the
        # centroids reach and follow the real trace.
        par['reduce']['findobj']['trace_npoly'] = 1
        par['reduce']['findobj']['trace_maxdev'] = 50.
        par['reduce']['findobj']['trace_maxshift'] = 20.
        # Do not correct for flexure
        par['flexure']['spec_method'] = 'skip'

        # Sensitivity function parameters
        par['sensfunc']['algorithm'] = 'IR'
        par['sensfunc']['polyorder'] = 8
        # ToDo: replace the telluric grid file for MMT site.
        par['sensfunc']['IR']['telgridfile'] = 'TellPCA_3000_26000_R10000.fits'

        return par

    def config_specific_par(
            self,
            inp:str|list|Path|fits.Header|Table,
            inp_par:parset.ParSet|None=None
        ) -> parset.ParSet:
        """
        Modify the PypeIt parameters to hard-wired values used for
        specific instrument configurations.

        Args:
            inp (:obj:`str`, :obj:`list`, `Path`_, `astropy.io.fits.Header`_, `astropy.table.Table`_):
                Input filename, an `astropy.io.fits.Header`_ object, or a list
                of `astropy.io.fits.Header`_ objects.  Or a row from the
                metadata table.
            inp_par (:class:`~pypeit.par.parset.ParSet`, optional):
                Parameter set used for the full run of PypeIt.  If None,
                use :func:`default_pypeit_par`.

        Returns:
            :class:`~pypeit.par.parset.ParSet`: The PypeIt parameter set
            adjusted for configuration specific parameter values.
        """
        # Start with instrument-wide parameters
        par = super().config_specific_par(inp, inp_par=inp_par)

       # Adjust parameters based on grating & dichroic used
        grating = self.get_meta_value(inp, 'dispname')
        dichroic = self.get_meta_value(inp, 'dichroic')

        if (grating=='HK') and (dichroic=='zJ'):
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = 'mmt_mmirs_HK_zJ.fits'
        elif (grating=='K3000') and (dichroic=='Kspec'):
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = 'mmt_mmirs_K3000_Kspec.fits'
        elif (grating=='J') and (dichroic=='zJ'):
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = 'mmt_mmirs_J_zJ.fits'

        # Auto-enable slitmask design if a <decker>.msk file sits next to the
        # data.  Never fatal: any problem falls back to generic tracing and the
        # SLIT-id keyed coadd workflow.
        maskfile = self._find_maskfile(inp)
        if maskfile is not None:
            par['calibrations']['slitedges']['use_maskdesign'] = True
            par['calibrations']['slitedges']['maskdesign_filename'] = str(maskfile)
            par['reduce']['slitmask']['assign_obj'] = True
            par['reduce']['slitmask']['extract_missing_objs'] = True
            par['reduce']['slitmask']['use_alignbox'] = True

        return par

    def _find_maskfile(self, inp):
        """
        Locate a ``<decker>.msk`` mask-design file next to the input frame.

        The mask label is the ``decker`` metadata value (header ``APERTURE`` ==
        ``MOSID``); the discovery file is ``<decker>.msk`` in the same directory
        as the frame.  Any failure to resolve the label or directory returns
        ``None`` so the caller falls back to generic slit tracing.

        Parameters
        ----------
        inp : :obj:`str`, `Path`_, `astropy.io.fits.Header`_, or `astropy.table.Table`_
            The input passed to :func:`config_specific_par`.

        Returns
        -------
        `Path`_ or :obj:`None`
            Path to the mask file if found, else ``None``.
        """
        try:
            label = self.get_meta_value(inp, 'decker')
        except Exception:
            return None
        if label is None:
            return None
        # Resolve the directory of the frame.  The full run_pypeit path sends a
        # single-row metadata Table whose 'filename' column holds the full path
        # (there is no 'directory' column at this stage), while scripts may send
        # a filename string/Path directly.
        directory = None
        if isinstance(inp, (str, Path)):
            directory = Path(inp).parent
        elif isinstance(inp, Table):
            if 'directory' in inp.colnames:
                directory = Path(str(inp['directory'][0]))
            elif 'filename' in inp.colnames:
                directory = Path(str(inp['filename'][0])).parent
        if directory is None:
            return None
        maskfile = directory / f'{label}.msk'
        return maskfile if maskfile.exists() else None

    def get_slitmask(self, filename, det=1):
        """
        Parse an MMIRS ``.msk`` mask-design file into :attr:`slitmask`.

        The mask ``y`` coordinate (mm) maps to the detector spatial direction
        and ``x`` to the spectral direction; slit corners are built as
        rectangles in on-sky arcseconds with the spatial axis anti-aligned to
        mask ``y`` (see :func:`get_maskdef_slitedges`).  ``BOX`` slits are
        flagged as alignment slits and ``TARGET`` slits as science slits.

        Parameters
        ----------
        filename : :obj:`str` or `Path`_
            Path to the ``.msk`` file.
        det : :obj:`int`, optional
            1-indexed detector number.  Ignored (MMIRS is single-detector),
            retained for API compatibility.

        Returns
        -------
        :class:`~pypeit.spectrographs.slitmask.SlitMask`
            The slitmask, also stored in :attr:`slitmask`.
        """
        header, slits = mmirs_maskfile.read_mmirs_maskfile(filename)

        arcsec_per_mm = 1.0 / header['arc2mm']
        n = len(slits)

        # SlitMask corners: shape (N, 4, 2), (x=spatial, y=spectral).
        # Spatial axis is -y_mm, spectral axis is x_mm (see spec geometry).
        # Build a rectangle centered on each slit from height (spatial) and
        # width (spectral), in arcsec.
        half_len = 0.5 * np.asarray(slits['height_mm'], dtype=float) * arcsec_per_mm
        half_wid = 0.5 * np.asarray(slits['width_mm'], dtype=float) * arcsec_per_mm
        cx = -np.asarray(slits['y_mm'], dtype=float) * arcsec_per_mm   # spatial center
        cy = np.asarray(slits['x_mm'], dtype=float) * arcsec_per_mm    # spectral center
        corners = np.zeros((n, 4, 2), dtype=float)
        # top-right, bottom-right, bottom-left, top-left (clockwise): the input
        # order SlitMask expects (high x/low y, low x/low y, low x/high y,
        # high x/high y).
        corners[:, 0, :] = np.column_stack([cx + half_len, cy - half_wid])
        corners[:, 1, :] = np.column_stack([cx - half_len, cy - half_wid])
        corners[:, 2, :] = np.column_stack([cx - half_len, cy + half_wid])
        corners[:, 3, :] = np.column_stack([cx + half_len, cy + half_wid])

        is_target = np.array([t == 'TARGET' for t in slits['type']])
        # top/bottom object distances from the slit edges (arcsec); the target
        # sits at the slit center, so both are the half-length.
        top = half_len.copy()
        bot = half_len.copy()
        objects = np.array([
            np.array(slits['slit'], dtype=int),
            np.array(slits['target'], dtype=int),
            np.array(slits['ra_deg'], dtype=float),
            np.array(slits['dec_deg'], dtype=float),
            np.array(slits['object'], dtype=object),
            np.zeros(n, dtype=float),               # no magnitude
            np.array(['None'] * n, dtype=object),   # no band
            top, bot], dtype=object).T

        onsky = np.column_stack([
            np.array(slits['ra_deg'], dtype=float),
            np.array(slits['dec_deg'], dtype=float),
            np.asarray(slits['height_mm'], dtype=float) * arcsec_per_mm,  # length arcsec
            np.asarray(slits['width_mm'], dtype=float) * arcsec_per_mm,   # width arcsec
            np.full(n, header['pa'], dtype=float)])

        self.slitmask = SlitMask(corners,
                                 slitid=np.array(slits['slit'], dtype=int),
                                 align=~is_target,
                                 science=is_target,
                                 onsky=onsky,
                                 objects=objects,
                                 object_names=np.array(slits['object'], dtype=object),
                                 mask_radec=(header['ra_deg'], header['dec_deg']),
                                 posx_pa=header['pa'])
        return self.slitmask

    def get_maskdef_slitedges(self, filename=None, det=1, debug=None,
                              binning=None, trc_path=None):
        """
        Predict slit-edge spatial-pixel positions from the ``.msk`` design.

        The mask ``y`` coordinate (mm) maps linearly and anti-aligned to the
        detector spatial pixel; the global offset/registration is fit
        downstream by
        :func:`~pypeit.edgetrace.EdgeTraceSet.maskdesign_matching`, so only the
        relative scale must be correct here.  The scale is
        ``(1/arc2mm)/platescale/bin_spat`` pixels per mm.  The mask ``x``
        coordinate is the dispersion axis and does not enter the spatial edge
        prediction.

        Parameters
        ----------
        filename : :obj:`str` or :obj:`list`
            Path to the ``.msk`` file (a list uses the first element).
        det : :obj:`int`, optional
            1-indexed detector number.
        debug : :obj:`bool`, optional
            Unused; retained for API compatibility.
        binning : :obj:`str`, optional
            ``'spec,spat'`` binning of the trace image.
        trc_path : :obj:`str`, optional
            Directory of the trace image, used to resolve a relative
            ``filename``.

        Returns
        -------
        left_edges : `numpy.ndarray`_
            Predicted left slit edges in spatial pixels, ordered to match
            ``slitmask.slitid``.
        right_edges : `numpy.ndarray`_
            Predicted right slit edges in spatial pixels, ordered to match
            ``slitmask.slitid``.
        sortindx : `numpy.ndarray`_
            Indices ordering the slits left to right.
        slitmask : :class:`~pypeit.spectrographs.slitmask.SlitMask`
            The mask design, also stored in :attr:`slitmask`.
        """
        _fname = filename[0] if isinstance(filename, (list, tuple)) else filename
        _fname = str(_fname)
        if trc_path is not None and not Path(_fname).exists():
            _fname = str(Path(trc_path) / Path(_fname).name)
        if not Path(_fname).exists():
            raise PypeItError(f'The mask design file {_fname} does not exist.')

        bin_spat = 1
        if binning is not None:
            _, bin_spat = parse.parse_binning(binning)
        platescale = self.get_detector_par(det=det)['platescale']

        header, slits = mmirs_maskfile.read_mmirs_maskfile(_fname)
        self.get_slitmask(_fname)

        arcsec_per_mm = 1.0 / header['arc2mm']
        scale = arcsec_per_mm / platescale / bin_spat        # px/mm
        y = np.asarray(slits['y_mm'], dtype=float)
        # arbitrary constant; absolute registration is fit by the matcher
        const = 1024.0 - scale * np.median(y)
        centers = const - scale * y
        half = 0.5 * np.asarray(slits['height_mm'], dtype=float) * arcsec_per_mm \
            / platescale / bin_spat
        left_edges = centers - half
        right_edges = centers + half

        # order edges to match slitmask.slitid, as GMOS does
        idx = utils.index_of_x_eq_y(np.asarray(slits['slit'], dtype=int),
                                    self.slitmask.slitid, strict=True)
        left_edges = left_edges[idx]
        right_edges = right_edges[idx]
        sortindx = np.argsort(left_edges)
        return left_edges, right_edges, sortindx, self.slitmask

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.

        Args:
            ftype (:obj:`str`):
                Type of frame to check. Must be a valid frame type; see
                frame-type :ref:`frame_type_defs`.
            fitstbl (`astropy.table.Table`_):
                The table with the metadata for one or more frames to check.
            exprng (:obj:`list`, optional):
                Range in the allowed exposure time for a frame of type
                ``ftype``. See
                :func:`pypeit.core.framematch.check_frame_exptime`.

        Returns:
            `numpy.ndarray`_: Boolean array with the flags selecting the
            exposures in ``fitstbl`` that are ``ftype`` type frames.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype in ['pinhole', 'bias']:
            # No pinhole or bias frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            return good_exp & (fitstbl['idname'] == 'flat')
        if ftype == 'standard':
            return good_exp & (fitstbl['idname'] == 'object')
        if ftype == 'science':
            return good_exp & (fitstbl['idname'] == 'object')
        if ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['idname'] == 'object')
        if ftype == 'dark':
            return good_exp & (fitstbl['idname'] == 'dark')
        log.debug('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

    def get_comb_group(self, fitstbl):
        """
        Automatically assign A-B nod combination/background groups from the
        derived along-slit dither offsets.

        Called by
        :func:`~pypeit.metadata.PypeItMetaData.set_combination_groups` after a
        unique ``comb_id`` has been assigned to every science/standard frame.
        For each instrument configuration, frames are split into two nods by the
        midpoint of the observed ``dithoff`` range, and each frame is paired --
        by greedy walk in time order -- with the temporally-adjacent still-
        unpaired frame on the opposite nod side.  The pair's ``bkg_id`` values
        are cross-linked so PypeIt subtracts one from the other.  A sequence
        whose peak-to-peak ``dithoff`` range is below :attr:`nod_min_offset` is
        treated as a stare and left unpaired.  Works identically for long-slit
        and MOS (it does not use any mask/decker metadata).

        Args:
            fitstbl (`astropy.table.Table`_):
                Metadata table for all frames.  Modified in place.

        Returns:
            `astropy.table.Table`_: The modified table.
        """
        sci_std = np.array(['science' in ft or 'standard' in ft
                            for ft in fitstbl['frametype']])
        if not np.any(sci_std):
            return fitstbl
        # Make sure the label columns exist and can hold arbitrary-length
        # strings.  init_meta seeds them from the 4-char default 'None', giving
        # a width-4 unicode column that would silently truncate longer labels
        # (e.g. "ABA'B'" -> "ABA'"), so coerce any fixed-width string column to
        # object dtype before writing.
        for col in ['dithpat', 'dithpos']:
            if col not in fitstbl.colnames:
                fitstbl[col] = np.full(len(fitstbl), 'None', dtype=object)
            elif fitstbl[col].dtype.kind in ('U', 'S'):
                fitstbl[col] = np.asarray(fitstbl[col], dtype=object)

        setups = np.unique(fitstbl['setup'][sci_std]) \
            if 'setup' in fitstbl.colnames else np.array(['A'])
        for setup in setups:
            in_cfg = sci_std & (fitstbl['setup'] == setup) \
                if 'setup' in fitstbl.colnames else sci_std
            idx = np.where(in_cfg)[0]
            if idx.size < 2:
                continue
            dithoff = np.asarray(fitstbl['dithoff'][idx], dtype=float)
            if dithoff.max() - dithoff.min() < self.nod_min_offset:
                # Stare / not nodded: leave bkg_id = -1.
                continue

            midpoint = 0.5 * (dithoff.max() + dithoff.min())
            side = dithoff > midpoint                      # True = "A" side
            order = np.argsort(np.asarray(fitstbl['mjd'][idx], dtype=float))
            combid = np.asarray(fitstbl['comb_id'][idx])
            bkgid = np.asarray(fitstbl['bkg_id'][idx])
            paired = np.zeros(idx.size, dtype=bool)

            # Greedy sequential pairing across the nod split, in time order.
            for a in range(idx.size):
                i = order[a]
                if paired[i]:
                    continue
                for b in range(a + 1, idx.size):
                    j = order[b]
                    if paired[j] or side[j] == side[i]:
                        continue
                    bkgid[i] = combid[j]
                    bkgid[j] = combid[i]
                    paired[i] = paired[j] = True
                    break
            fitstbl['bkg_id'][idx] = bkgid

            # Informational A/B (+prime) labels: within each nod side, distinct
            # dithoff values (rounded to 0.1") get a prime suffix by first
            # appearance in time.
            dithpos = np.array(['None'] * idx.size, dtype=object)
            for is_A, base in [(True, 'A'), (False, 'B')]:
                grp = np.where(side == is_A)[0]
                grp = grp[np.argsort(order.argsort()[grp])]   # time order
                seen = []
                for g in grp:
                    val = round(float(dithoff[g]), 1)
                    if val not in seen:
                        seen.append(val)
                    dithpos[g] = base + "'" * seen.index(val)
            fitstbl['dithpos'][idx] = np.asarray(dithpos, dtype=str)
            # Pattern string = unique labels in time order, e.g. "ABA'B'".
            seq = list(np.asarray(dithpos)[order])
            fitstbl['dithpat'][idx] = ''.join(dict.fromkeys(seq))
        return fitstbl

    def pypeit_file_keys(self):
        """
        Define the list of columns written to the pypeit file, adding the
        derived dither columns so the A-B nod grouping is visible/editable.

        Returns:
            :obj:`list`: Column keywords for the pypeit file.
        """
        return super().pypeit_file_keys() + ['dithpat', 'dithpos', 'dithoff',
                                             'frameno']

    def bpm(self, filename, det, shape=None, msbias=None):
        """
        Generate a default bad-pixel mask.

        Even though they are both optional, either the precise shape for
        the image (``shape``) or an example file that can be read to get
        the shape (``filename`` using :func:`get_image_shape`) *must* be
        provided.

        Args:
            filename (:obj:`str` or None):
                An example file to use to get the image shape.
            det (:obj:`int`):
                1-indexed detector number to use when getting the image
                shape from the example file.
            shape (tuple, optional):
                Processed image shape
                Required if filename is None
                Ignored if filename is not None
            msbias (`numpy.ndarray`_, optional):
                Processed bias frame used to identify bad pixels

        Returns:
            `numpy.ndarray`_: An integer array with a masked value set
            to 1 and an unmasked value set to 0.  All values are set to
            0.
        """
        # Call the base-class method to generate the empty bpm
        bpm_img = super().bpm(filename, det, shape=shape, msbias=msbias)

        log.info("Using hard-coded BPM for det=1 on MMIRS")

        # Get the binning
        hdu = io.fits_open(filename)
        binning = hdu[1].header['CCDSUM']
        hdu.close()

        # Apply the mask
        xbin, ybin = int(binning.split(' ')[0]), int(binning.split(' ')[1])
        bpm_img[:, 187 // ybin] = 1

        return bpm_img

    def get_rawimage(self, raw_file, det):
        """
        Read raw images and generate a few other bits and pieces
        that are key for image processing.

        Parameters
        ----------
        raw_file : :obj:`str`
            File to read
        det : :obj:`int`
            1-indexed detector to read

        Returns
        -------
        detector_par : :class:`pypeit.images.detector_container.DetectorContainer`
            Detector metadata parameters.
        raw_img : `numpy.ndarray`_
            Raw image for this detector.
        hdu : `astropy.io.fits.HDUList`_
            Opened fits file
        exptime : :obj:`float`
            Exposure time read from the file header
        rawdatasec_img : `numpy.ndarray`_
            Data (Science) section of the detector as provided by setting the
            (1-indexed) number of the amplifier used to read each detector
            pixel. Pixels unassociated with any amplifier are set to 0.
        oscansec_img : `numpy.ndarray`_
            Overscan section of the detector as provided by setting the
            (1-indexed) number of the amplifier used to read each detector
            pixel. Pixels unassociated with any amplifier are set to 0.

        Frames with at least :attr:`ramp_min_reads` non-destructive reads are
        combined using up-the-ramp fitting (see :func:`_ramp_fit_image`);
        frames with fewer reads use correlated double sampling, as before.

        Fitted images are persisted as 2D count-rate files in the ``RampFit``
        directory inside the reduction directory (written on first load,
        reused on subsequent loads while the raw file is unchanged).  The
        reduction directory is recorded by :func:`cache_metadata`; when that
        hook never fired (e.g. direct API use), the current working
        directory is used instead.  Preprocessed files — created here or by
        ``pypeit_fit_ramp`` — are identified by the ``RAMPFIT`` header
        card and loaded directly.
        """
        fil = utils.find_single_file(f'{raw_file}*', required=True)

        # Read
        log.info(f'Reading MMIRS file: {fil}')
        hdu = io.fits_open(fil)

        redux_path = self._ramp_output_dir if self._ramp_output_dir is not None \
                else Path.cwd()

        if hdu[0].header.get('RAMPFIT') is None \
                and mmirs_count_reads(hdu) >= self.ramp_min_reads:
            # Multi-read cube: swap in a fresh preprocessed 2D image if one
            # exists in the reduction directory
            rampfit_file = mmirs_rampfit_path(fil, redux_path)
            if mmirs_rampfit_fresh(rampfit_file, fil):
                log.info(f'Loading preprocessed ramp image: {rampfit_file}')
                hdu.close()
                hdu = io.fits_open(rampfit_file)

        head1 = hdu[1].header

        detector_par = self.get_detector_par(det if det is not None else 1, hdu=hdu)

        # get the x and y binning factors...
        binning = head1['CCDSUM']
        xbin, ybin = [int(ibin) for ibin in binning.split(' ')]

        # Need the exposure time
        exptime = hdu[self.meta['exptime']['ext']].header[self.meta['exptime']['card']]
        gain = detector_par['gain'][0]

        if hdu[0].header.get('RAMPFIT') is not None:
            # Preprocessed 2D count-rate image (e-/s): convert to ADU
            array = hdu[1].data.astype(np.float64) * exptime / gain
            detector_par['ronoise'] = np.atleast_1d(hdu[0].header['RAMPRON'])
        elif mmirs_count_reads(hdu) >= self.ramp_min_reads:
            # Up-the-ramp fitting with jump detection
            rate, sig, eff_ronoise = self._ramp_fit_image(hdu, detector_par)
            detector_par['ronoise'] = np.atleast_1d(eff_ronoise)
            array = rate * exptime / gain
            # Persist the fit so later loads (and other scripts) reuse it
            rampfit_file = mmirs_rampfit_path(fil, redux_path)
            try:
                mmirs_write_rampfit(rampfit_file, rate, hdu, sig, eff_ronoise,
                                    Path(fil).stat().st_mtime)
                log.info(f'Wrote preprocessed ramp image: {rampfit_file}')
            except OSError as e:
                log.warning(f'Could not write preprocessed ramp image '
                            f'{rampfit_file} ({e}); continuing without it.')
        else:
            # Correlated double sampling (first minus last read)
            datasec = head1['DATASEC']
            x1, x2, y1, y2 = np.array(parse.load_sections(datasec,
                                                          fmt_iraf=False)).flatten()
            if len(hdu) > 2:
                data = mmirs_read_amp(hdu[1].data.astype('float64')) \
                        - mmirs_read_amp(hdu[2].data.astype('float64'))
            else:
                data = mmirs_read_amp(hdu[1].data.astype('float64'))
            array = data[x1-1:x2, y1-1:y2]

        ## ToDo: This is a hack. Need to solve this issue. I cut at 998 due to the HK zero order contaminating
        ## the blue part of the zJ+HK spectrum. For other setup, you do not need to cut the detector.
        if (head1['FILTER']=='zJ') and (head1['DISPERSE']=='HK'):
            array = array[:int(998/ybin),:]
        rawdatasec_img = np.ones_like(array,dtype='int')
        # NOTE: If there is no overscan, must be set to 0s
        oscansec_img = np.zeros_like(array,dtype='int')

        # Return, transposing array back to orient the overscan properly
        return detector_par, np.flipud(array), hdu, exptime, np.flipud(rawdatasec_img),\
               np.flipud(np.flipud(oscansec_img))

    def _ramp_fit_image(self, hdu, detector_par):
        """
        Perform up-the-ramp fitting of a multi-read MMIRS frame.

        This uses the ``fitramp`` algorithm of Brandt (2024,
        https://arxiv.org/abs/2404.01326) and was inspired by the prototype
        at https://github.com/zhechenghu/mmt-mmirs-up-the-ramp-pypeit.

        Args:
            hdu (`astropy.io.fits.HDUList`_):
                Opened raw file with at least :attr:`ramp_min_reads`
                non-destructive reads.
            detector_par (:class:`~pypeit.images.detector_container.DetectorContainer`):
                Detector parameters; provides the gain.

        Returns:
            :obj:`tuple`: The fitted count-rate image in e-/s
            (`numpy.ndarray`_, shape of the trimmed data section), the
            single-read noise in electrons (:obj:`float`), and the
            effective read noise in electrons (:obj:`float`).
        """
        reads, head1 = mmirs_load_ramp(hdu)
        ngroups = reads.shape[0]
        gain = detector_par['gain'][0]
        reads *= gain      # ADU -> electrons
        covar = fitramp.Covar([head1['GRPTIME'] * (i + 1)
                               for i in range(ngroups)])
        diffs = mmirs_ramp_diffs(reads, covar)
        del reads
        sig = self.get_ramp_sigma(diffs, covar, exptime=head1.get('EXPTIME'),
                                  ron_floor=head1.get('RDNOISE'))
        log.info(f'Up-the-ramp fitting {ngroups} reads '
                 f'(single-read noise {sig:.2f} e-)')
        countrate = mmirs_fit_ramp(diffs, covar, sig,
                                   workers=self.ramp_fit_workers)
        eff_ronoise = mmirs_effective_ronoise(sig, ngroups)
        log.info(f'Effective read noise: {eff_ronoise:.2f} e-')
        return countrate, sig, eff_ronoise

def mmirs_read_amp(img, namps=32):
    """
    MMIRS has 32 reading out channels. Need to deal with this issue a little
    bit. The pypeit overscan subtraction is not used; reference-pixel correction
    is done here instead.

    Imported from MMIRS IDL pipeline refpix.pro
    """

    # number of channels for reading out
    if namps is None:
        namps = 32

    data_shape = np.shape(img)
    ampsize = int(data_shape[0] / namps)

    refpix1 = np.array([1, 2, 3])
    refpix2 = np.arange(4) + data_shape[0] - 4
    refpix_all = np.hstack([[0, 1, 2, 3], np.arange(4) + data_shape[0] - 4])
    refvec = np.sum(img[:, refpix_all], axis=1) / np.size(refpix_all)
    svec = savgol_filter(refvec, 11, polyorder=5)

    refvec_2d = np.reshape(np.repeat(svec, data_shape[0], axis=0), data_shape)
    img_out = img - refvec_2d

    for amp in range(namps):
        img_out_ref = img_out[np.hstack([refpix1, refpix2]), :]
        ref1, _, _ = sigma_clipped_stats(
            img_out_ref[:, amp * ampsize + 2 * np.arange(int(ampsize / 2))], sigma=3
        )
        ref2, _, _ = sigma_clipped_stats(
            img_out_ref[:, amp * ampsize + 2 * np.arange(int(ampsize / 2)) + 1], sigma=3
        )
        ref12 = (ref1 + ref2) / 2.
        img_out[:, amp * ampsize:(amp + 1) * ampsize] -= ref12

    return img_out


def mmirs_load_ramp(hdu, namps=32):
    """
    Load the non-destructive reads of an MMIRS ramp in time order.

    Image extensions are sorted by ``EXTVER`` (in raw MMIRS files, ext 1 is
    the *final* read and holds the complete metadata).  Each read is
    reference-pixel corrected with :func:`mmirs_read_amp` on the full frame
    and then trimmed to ``DATASEC``.

    Parameters
    ----------
    hdu : `astropy.io.fits.HDUList`_
        Opened raw MMIRS file.
    namps : :obj:`int`, optional
        Number of readout amplifiers passed to :func:`mmirs_read_amp`.

    Returns
    -------
    reads : `numpy.ndarray`_
        Float64 array with shape ``(ngroups, ny, nx)`` holding the
        reference-pixel-corrected reads in ADU, trimmed to ``DATASEC``.
    head1 : `astropy.io.fits.Header`_
        Header of extension 1 (the metadata-complete final read).
    """
    head1 = hdu[1].header
    img_hdus = sorted([h for h in hdu if h.header.get('NAXIS') == 2
                       and h.header.get('NAXIS1', 0) > 0],
                      key=lambda h: h.header['EXTVER'])
    x1, x2, y1, y2 = np.array(parse.load_sections(head1['DATASEC'],
                                                  fmt_iraf=False)).flatten()
    ngroups = len(img_hdus)
    reads = np.empty((ngroups, x2 - x1 + 1, y2 - y1 + 1), dtype=np.float64)
    for i, h in enumerate(img_hdus):
        frame = mmirs_read_amp(h.data.astype(np.float64), namps=namps)
        reads[i] = frame[x1-1:x2, y1-1:y2]
    return reads, head1


def mmirs_count_reads(hdu):
    """
    Count the non-destructive reads (non-empty 2D image extensions) in an
    MMIRS file.

    Parameters
    ----------
    hdu : `astropy.io.fits.HDUList`_
        Opened MMIRS file.

    Returns
    -------
    :obj:`int`
        Number of non-empty 2D image extensions.
    """
    return sum(1 for h in hdu if h.header.get('NAXIS') == 2
               and h.header.get('NAXIS1', 0) > 0)


def mmirs_rampfit_path(raw_file, redux_path):
    """
    Return the preprocessed-image path for a raw MMIRS cube.

    Preprocessed 2D count-rate images live in a ``RampFit`` directory
    inside the reduction directory (alongside ``Calibrations``,
    ``Science``, etc.), with the same file name as the raw cube.

    Parameters
    ----------
    raw_file : :obj:`str`, `Path`_
        Path to the raw MMIRS cube.
    redux_path : :obj:`str`, `Path`_
        Path to the reduction directory.

    Returns
    -------
    `Path`_
        ``<redux_path>/RampFit/<raw filename>``
    """
    return Path(redux_path) / 'RampFit' / Path(raw_file).name


def mmirs_rampfit_fresh(rampfit_file, raw_file):
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
    try:
        mtime = float(fits.getval(rampfit_file, 'RAWMTIME'))
        raw_mtime = Path(raw_file).stat().st_mtime
    except (KeyError, OSError):
        return False
    return abs(mtime - raw_mtime) < 1.


def mmirs_write_rampfit(rampfit_file, rate, hdu, sig, eff_ronoise, raw_mtime):
    """
    Write a preprocessed MMIRS 2D count-rate image.

    The output carries a copy of the raw primary header plus the cards
    ``RAMPFIT`` (marker), ``RAMPSIG``, ``RAMPRON``, ``NGROUPS`` and
    ``RAWMTIME``, and a single image extension holding the fitted count
    rate in e-/s (float32) under a copy of the raw final-read header, so
    all metadata used by ``pypeit_setup`` is preserved.

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
        Opened source raw cube (headers are copied from it).
    sig : :obj:`float`
        Single-read noise used in the fit (electrons).
    eff_ronoise : :obj:`float`
        Effective read noise of the fitted image (electrons).
    raw_mtime : :obj:`float`
        Modification time of the source raw cube.

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
    prihead['NGROUPS'] = (mmirs_count_reads(hdu),
                          'Number of reads in the source ramp')
    prihead['RAWMTIME'] = (float(raw_mtime),
                           'Modification time of the source raw cube')
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


def mmirs_ramp_diffs(reads, covar):
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
    diffs = np.diff(reads, axis=0)
    diffs /= np.asarray(covar.delta_t)[:, None, None]
    return diffs


def _resolve_ramp_workers(workers):
    """
    Resolve the requested number of ramp-fit worker threads.

    ``None`` selects ``min(6, os.cpu_count())``; any explicit value is passed
    through (clamped to at least 1).  The per-pixel fit is memory-bandwidth
    bound, so more than ~6 threads does not help and can hurt.
    """
    if workers is None:
        return max(1, min(6, os.cpu_count() or 1))
    return max(1, int(workers))


def _fit_ramp_rows(diffs, nb, workers, worker):
    """
    Dispatch ``worker((row_start, row_stop))`` over blocks of detector rows.

    The per-pixel ramp fit is independent, so the ``ny`` rows of ``diffs`` are
    split into contiguous blocks of ``nb`` rows and each block is handed to
    ``worker``, which is expected to write its results into a preallocated
    output array.  With ``workers > 1`` the blocks are fit concurrently in a
    thread pool; NumPy releases the GIL during the element-wise arithmetic that
    dominates :func:`~pypeit.ext.fitramp.fitramp.fit_ramps`, so threads scale until the
    memory bus saturates (empirically ~3x at 6 threads).

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


def mmirs_calibrate_sigma(diffs, covar, sig_guess=9.0, nrows=200, workers=None,
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

    _fit_ramp_rows(sub, nb, _resolve_ramp_workers(workers), worker)
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


def mmirs_fit_ramp(diffs, covar, sig, workers=None, nb=16):
    """
    Fit all pixels of a ramp, in blocks of rows, with jump detection.

    The per-pixel fit is independent, so the detector is fit in blocks of
    ``nb`` rows (rows folded into the pixel axis of
    :func:`~pypeit.ext.fitramp.fitramp.fit_ramps`) and, for ``workers > 1``, the blocks
    are fit concurrently.  Results are numerically identical to a row-by-row fit.

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

    _fit_ramp_rows(diffs, nb, _resolve_ramp_workers(workers), worker)
    return countrate


def mmirs_effective_ronoise(sig, ngroups):
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
