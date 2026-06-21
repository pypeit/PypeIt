"""
Module for P200/NGPS specific methods.

NGPS raw FITS files store each channel ('U', 'G', 'R', 'I') in its own
image extension labelled by ``EXTNAME``. The order of the image
extensions is *not* stable across files, and frames where only some
channels were read out (e.g. U+G dome flats) drop the unused
extensions entirely. Per-channel subclasses below resolve their image
extension by scanning ``EXTNAME``, never by hardcoded HDU index.

Per-channel CCD readout orientation is expressed through PypeIt's stock
``specaxis``/``specflip``/``spatflip`` detector parameters, so the raw
image is reoriented by
:meth:`~pypeit.spectrographs.spectrograph.Spectrograph.orient_image`
exactly as for every other spectrograph (and the
:ref:`instr_par-detectors` table reports the true orientation):

* R, I: dispersion along NAXIS1 -> ``specaxis=1``.
* G:    dispersion along NAXIS1 but flipped left-right relative to R/I
        -> ``specaxis=1, specflip=True``.
* U:    dispersion along NAXIS2 (read out rotated 90 deg) and flipped
        -> ``specaxis=0, specflip=True``.

``DATASEC``/``OSCANSEC`` are read from each frame's own image-extension
header in the native readout frame and reoriented alongside the image.
The U channel ships malformed ``DATASEC`` cards, so it falls back to the
full native frame.

For operational use at the telescope (whole-night ingestion, per-channel
setup/reduction orchestration, and real-time quicklook products), see
the separately maintained wrapper around these classes:
https://github.com/cfremling/ngps_pipeline

.. include:: ../include/links.rst
"""
from pathlib import Path

import numpy as np

from astropy.io import fits
from astropy.time import Time

from pypeit import log
from pypeit import PypeItError
from pypeit import io
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.spectrographs import spectrograph
from pypeit.core import parse
from pypeit.images import detector_container


class P200NGPSSpectrograph(spectrograph.Spectrograph):
    """
    Common base for P200/NGPS per-channel subclasses.

    Each per-channel subclass sets ``_extname`` to the EXTNAME of its
    image extension in the raw FITS file ('U', 'G', 'R', or 'I') and,
    where its CCD readout orientation differs from R/I, sets
    ``_specaxis``/``_specflip`` so that PypeIt's stock
    :meth:`~pypeit.spectrographs.spectrograph.Spectrograph.orient_image`
    reorients the raw frame the same way for every channel (dispersion
    along the spectral axis, slices in a uniform spatial order).
    """
    ndet = 1
    telescope = telescopes.P200TelescopePar()

    # -- Per-channel overrides -----------------------------------------
    #: EXTNAME of this channel in the raw FITS HDU list.
    _extname = None
    #: Native CCD readout orientation, applied by PypeIt's stock
    #: :meth:`~pypeit.spectrographs.spectrograph.Spectrograph.orient_image`.
    #: R/I read out dispersion-along-NAXIS1 (``specaxis=1``); G is the
    #: same but left-right flipped (``specflip=True``); U reads out
    #: dispersion-along-NAXIS2 (``specaxis=0``) and flipped.  Subclasses
    #: override these as needed.
    _specaxis = 1
    _specflip = False
    _spatflip = False
    #: Whether to use the bias overscan region in raw images.
    _use_overscan = True
    #: Hardcoded fallback ``DATASEC`` (FITS-style ``[NAXIS1, NAXIS2]``)
    #: used only when the per-frame extension header has no usable
    #: ``DATASEC`` card (e.g. the U channel ships malformed
    #: ``[-1:1097,...]``).  When ``None``, the per-frame raw header
    #: ``DATASEC`` card is used as-is (native readout frame).
    _datasec_fallback = None
    #: Hardcoded fallback ``OSCANSEC`` (bias) for the same case.
    _oscansec_fallback = None
    #: Whether to use the master bias to build the bad-pixel mask.
    _bpm_usebias = False
    #: Per-channel arc-line wavelength reidentification template
    #: (``reid_arxiv``).  Lives in
    #: ``pypeit/data/arc_lines/reid_arxiv/``.
    _wvarxiv = None
    #: Default wavelength-calibration lamps and raw IMGTYPE values used
    #: for arc/tilt frame typing.  All four NGPS channels use ThAr;
    #: FeAr was investigated and rejected because pypeit's FeAr line
    #: list is too sparse on g/r/i (only 7-31 cataloged lines vs
    #: 2000-7000 ThAr lines), and pypeit's reidentify second-pass
    #: rejection makes the sparse-list fits unreliable.
    _wave_lamps = ['ThAr']
    _arc_idnames = ('THAR',)
    #: PypeIt wavelength-solution method.  'full_template' uses the
    #: single-slice central-slit ThAr archive and a global
    #: cross-correlation shift to anchor the per-slit polynomial fit.
    _wavecal_method = 'full_template'
    #: Minimum unmasked spectral length of a slit edge as a fraction of
    #: the detector size.  PypeIt's default 0.6 is too aggressive for
    #: the U channel where dome flat illumination only spans part of
    #: the spectral range.
    _fit_min_spec_length = 0.6
    #: Per-channel overrides for ``par['calibrations']['flatfield']``.
    #: Each subclass can tune b-spline knot spacing / median-filter
    #: width independently, since each NGPS channel has different
    #: optics, slit-edge sharpness and dome-flat illumination
    #: characteristics.  Empty by default; subclasses fill in.
    _flatfield_overrides = {}
    #: Per-channel overrides for ``par['reduce']['skysub']``.
    #:
    #:   - ``bspline_spacing = 0.5`` (vs PypeIt default 0.6): tighter
    #:     knot spacing along the spectral axis.
    #:   - ``no_local_sky = True``: skip the local-window sky re-fit
    #:     so the spec2d's sky-subtracted image is continuous.
    _skysub_overrides = {
        'bspline_spacing': 0.5,
        'no_local_sky': True,
    }
    #: Per-channel overrides for ``par['flexure']``; default
    #: ``spec_method='skip'``.
    _flexure_overrides = {}
    #: Per-channel overrides for ``par['sensfunc']``.  Top-level keys go
    #: in ``par['sensfunc']``; nested ``UVIS``/``IR`` dicts go in
    #: ``par['sensfunc']['UVIS']``/``['IR']``.
    #:
    #: Why ``resolution = 4000`` is constant rather than slit-aware:
    _sensfunc_overrides = dict(
        extr='OPT',
        UVIS=dict(polycorrect=False, resolution=4000),
    )

    #: PypeIt's ``full_template`` drops detected arc peaks whose
    #: nearest catalog line is more than ``match_toler`` *pixels* away.
    #: 0.3 binned pixels at NGPS' binspec=3 corresponds to a
    #: ~0.5-0.6 A window across u/g/r/i.
    _match_toler_pixels = 0.3
    #: Per-channel detector gain (electrons / ADU); set by each subclass.
    _gain = 1.0
    #: Per-channel readnoise (electrons).
    _ronoise = 8.5
    #: Per-channel arc-line detection SNR threshold for the wavecal.
    #: PypeIt's default is 5.0. We use 5.0 to keep only clean,
    #: well-isolated peaks -- the previous 1.0 was letting noise
    #: spikes / blended low-SNR features into the reidentify pool
    #: which then mis-matched against the catalog and biased the
    #: polynomial.
    _sigdetect = 5.0

    #: Detector saturation level (ADU).  Set conservatively to 40k
    #: for NGPS -- the CCDs nominally hold ~57 kADU but flats start
    #: to deviate from linearity well before that, and dome lamps
    #: occasionally drive bright slit regions past the conservative
    #: threshold.  Pypeit masks pixels above this value as saturated;
    #: ``nonlinear`` then sets the soft warning threshold as a fraction
    #: of saturation.
    _saturation_adu = 40000.0

    # -- Helpers -------------------------------------------------------
    def _extname_index(self, hdulist_or_headarr):
        """
        Return the HDU/header index whose ``EXTNAME`` equals
        ``self._extname``, or ``None`` if the channel is not present in
        this file (e.g. a U+G-only calibration drops the R and I
        extensions).
        """
        if hdulist_or_headarr is None:
            return None
        for i, item in enumerate(hdulist_or_headarr):
            hdr = item.header if hasattr(item, 'header') else item
            if hdr is None:
                continue
            extname = (hdr.get('EXTNAME') or '').strip()
            if extname == self._extname:
                return i
        return None

    def _datasec_string(self, ext_hdr):
        """Resolve the native-frame ``DATASEC`` for this extension.

        Default (R/I/G): use the per-frame raw ``DATASEC`` card (this
        picks up the binning-correct section the instrument writes out);
        fall back to :attr:`_datasec_fallback` if the raw section is
        missing or malformed.  ``orient_image`` reorients the resulting
        pixel map alongside the image, so the section is given in the
        native readout frame.  Subclasses (e.g. U, where every raw
        ``DATASEC`` is malformed) can override this method.
        """
        raw = ext_hdr.get('DATASEC')
        if self._parse_section(raw) is not None:
            return raw
        return self._datasec_fallback

    def _oscansec_string(self, ext_hdr):
        """Same as :meth:`_datasec_string` but for OSCANSEC/BIASSEC."""
        raw = ext_hdr.get('BIASSEC')
        if self._parse_section(raw) is not None:
            return raw
        return self._oscansec_fallback

    @staticmethod
    def _parse_section(sec_str):
        """Parse a FITS section string into (x_lo, x_hi, y_lo, y_hi),
        1-indexed inclusive.  Returns ``None`` if the section is
        malformed (e.g. the U channel's ``[-1:1097,...]``)."""
        if not sec_str or not isinstance(sec_str, str):
            return None
        try:
            x, y = sec_str.strip('[]').split(',')
            x_lo, x_hi = (int(v) for v in x.split(':'))
            y_lo, y_hi = (int(v) for v in y.split(':'))
            if x_lo < 1 or x_hi < 1 or y_lo < 1 or y_hi < 1:
                return None
            return x_lo, x_hi, y_lo, y_hi
        except (ValueError, AttributeError):
            return None

    # -- Standard Spectrograph hooks -----------------------------------
    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.
        """
        self.meta = {}
        # Required (core).  ra/dec are compound so that frames with a
        # missing TELRA/TELDEC card (e.g. truncated 2880-byte raw files
        # written by the instrument host on aborted exposures) yield
        # ``np.nan`` instead of ``None``, which would otherwise corrupt
        # the dtype of the metadata table column and break downstream
        # ``np.isnan`` checks in :meth:`vet_assigned_ftypes`.
        self.meta['ra'] = dict(card=None, compound=True,
                               required_ftypes=['science', 'standard'])
        self.meta['dec'] = dict(card=None, compound=True,
                                required_ftypes=['science', 'standard'])
        self.meta['target'] = dict(ext=0, card='NAME', compound=True,
                                   required_ftypes=['science', 'standard'])

        self.meta['dispname'] = dict(card=None, compound=True, default='VPH')
        self.meta['decker'] = dict(ext=0, card='SLITW', rtol=1e-2)
        self.meta['binning'] = dict(card=None, compound=True)

        self.meta['mjd'] = dict(ext=0, card='MJD')
        # Exposure-time units are fixed: ``SHUTTIME`` (the measured
        # shutter-open time) is always seconds, ``EXPTIME`` (the
        # requested time) is always milliseconds.  Prefer SHUTTIME and
        # fall back to EXPTIME/1000 when it is absent.  See
        # :meth:`compound_meta`.
        self.meta['exptime'] = dict(card=None, compound=True)
        self.meta['airmass'] = dict(ext=0, card='AIRMASS',
                                    required_ftypes=['science', 'standard'])

        self.meta['dichroic'] = dict(card=None, compound=True)
        self.meta['dispangle'] = dict(card=None, rtol=1e-2, compound=True)
        self.meta['slitwid'] = dict(ext=0, card='SLITW', rtol=1e-2)
        self.meta['idname'] = dict(ext=0, card='IMGTYPE')
        self.meta['instrument'] = dict(ext=0, card='INSTRUME')

        self.meta['lampstat01'] = dict(ext=0, card='LAMPBLUC')
        self.meta['lampstat02'] = dict(ext=0, card='LAMPFEAR')
        self.meta['lampstat03'] = dict(ext=0, card='LAMPREDC')
        self.meta['lampstat04'] = dict(ext=0, card='LAMPTHAR')

    def configuration_keys(self):
        """
        Return the metadata keys that define a unique instrument
        configuration.

        This list is used by :class:`~pypeit.metadata.PypeItMetaData` to
        identify the unique configurations among the list of frames read
        for a given reduction.

        Returns:
            :obj:`list`: List of keywords of data pulled from file headers
            and used to construct the :class:`~pypeit.metadata.PypeItMetaData`
            object.
        """
        return ['binning']

    def raw_header_cards(self):
        """
        Return additional raw header cards to be propagated in
        downstream output files for configuration identification.

        This list is used by
        :meth:`~pypeit.spectrographs.spectrograph.Spectrograph.subheader_for_spec`
        to include additional FITS keywords in downstream output files.

        Returns:
            :obj:`list`: List of keywords from the raw data files that
            should be propagated in output files.
        """
        return ['GRATING', 'ANGLE', 'APERTURE']

    def pypeit_file_keys(self):
        """
        Define the list of keys to be output into a standard PypeIt file.

        Returns:
            :obj:`list`: The list of keywords in the relevant
            :class:`~pypeit.metadata.PypeItMetaData` instance to print
            to the :ref:`pypeit_file`.
        """
        return super().pypeit_file_keys()

    _extname_cache: dict = {}

    @classmethod
    def _file_extnames(cls, path) -> frozenset[str]:
        """Cached EXTNAME set of the image HDUs in this raw FITS."""
        key = str(path)
        cached = cls._extname_cache.get(key)
        if cached is None:
            try:
                with fits.open(path, memmap=True) as h:
                    cached = frozenset(
                        str(hdu.header.get('EXTNAME', '')).strip().upper()
                        for hdu in h[1:])
            except Exception:
                cached = frozenset()
            cls._extname_cache[key] = cached
        return cached

    # True on u and g subclasses: their arcs and dome flats must come
    # from the long ``getcalib_ug`` exposures (image extensions exactly
    # {G, U}, with R and I detectors off) -- the short 4-channel cal
    # frames have ~zero flux at u/g wavelengths.
    _require_ug_only_cals = False

    @classmethod
    def find_raw_files(cls, root, extension=None):
        """Same as base class, but only returns FITS files that carry
        this channel's image extension (``_extname``).  An NGPS raw
        directory contains a mix of long ``getcalib_ug`` exposures
        (only G+U image HDUs) and short 4-channel exposures
        (G+I+R+U), and only one subset is meaningful for any given
        per-channel spectrograph.  Filtering here means pypeit never
        attempts to read this channel's BINSPEC/BINSPAT from a file
        where that extension doesn't exist."""
        files = super().find_raw_files(root, extension=extension)
        if not cls._extname:
            return files
        return [p for p in files if cls._extname.upper() in cls._file_extnames(p)]

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Frame-typing rule, shared by all four channels.  IMGTYPE values
        in the raw header drive the baseline classification: ``SCI`` is
        science / standard, ``BIAS`` is bias, ``DOMEFLAT`` and ``CONT``
        are flats, and the channel-specific arc IMGTYPE values define
        the arc/tilt frames.  Discrimination between science and
        standard is left to PypeIt's downstream coordinate-matching
        against the archived flux-standard registry.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)

        # Type by IMGTYPE first.
        if ftype in ['science', 'standard']:
            mask = good_exp & (fitstbl['idname'] == 'SCI')
        elif ftype == 'bias':
            mask = good_exp & (fitstbl['idname'] == 'BIAS')
        elif ftype in ['pinhole', 'dark']:
            return np.zeros(len(fitstbl), dtype=bool)
        elif ftype in ['arc', 'tilt']:
            mask = good_exp & np.isin(fitstbl['idname'], self._arc_idnames)
        elif ftype in ['pixelflat', 'trace', 'illumflat']:
            mask = ((good_exp & (fitstbl['idname'] == 'DOMEFLAT'))
                    | (good_exp & (fitstbl['idname'] == 'CONT')))
        else:
            log.debug('Cannot determine if frames are of type {0}.'.format(ftype))
            return np.zeros(len(fitstbl), dtype=bool)

        # On u/g, restrict arc/flat frames to the long ``getcalib_ug``
        # exposures (image extensions exactly {G, U}); the short
        # 4-channel cal frames have ~zero flux at u/g wavelengths.
        if (self._require_ug_only_cals
                and ftype in ('arc', 'tilt', 'pixelflat',
                              'trace', 'illumflat')):
            for i in range(len(fitstbl)):
                if not mask[i]:
                    continue
                fpath = Path(fitstbl['directory'][i]) / fitstbl['filename'][i]
                if self._file_extnames(fpath) != frozenset({'G', 'U'}):
                    mask[i] = False
        return mask

    def compound_meta(self, headarr, meta_key):
        """
        Compound metadata.  Shared across all four channels; binning is
        read from the channel's own image extension (resolved by
        EXTNAME), so frames where the channel is absent return ``None``
        and are filtered out of the configuration.
        """
        retval = super().compound_meta(headarr, meta_key)
        if retval is not None:
            return retval

        ph = headarr[0] if headarr else None

        if meta_key == 'ra':
            v = ph.get('TELRA') if ph is not None else None
            if v is None or not str(v).strip():
                return np.nan
            return str(v).strip()

        if meta_key == 'dec':
            v = ph.get('TELDEC') if ph is not None else None
            if v is None or not str(v).strip():
                return np.nan
            return str(v).strip()

        if meta_key == 'mjd':
            for k in ('UTSHUT', 'EXPSTART', 'DATE-OBS'):
                v = ph.get(k) if ph is not None else None
                if v:
                    try:
                        return Time(v).mjd
                    except Exception:
                        continue
            return None

        if meta_key == 'exptime':
            # Fixed units (no value-based guessing): SHUTTIME is in
            # seconds, EXPTIME is in milliseconds.  Prefer SHUTTIME;
            # fall back to EXPTIME/1000 only when SHUTTIME is absent.
            if ph is not None:
                shut = ph.get('SHUTTIME')
                if shut is not None:
                    try:
                        return float(shut)
                    except (TypeError, ValueError):
                        pass
                ex = ph.get('EXPTIME')
                if ex is not None:
                    try:
                        return float(ex) / 1000.0
                    except (TypeError, ValueError):
                        pass
            return None

        if meta_key == 'dispangle':
            return 0

        if meta_key == 'binning':
            # Read BINSPEC/BINSPAT from this channel's image extension.
            # ``find_raw_files`` already filtered out any file that
            # doesn't carry this channel's image extension, so we know
            # ``_extname_index`` will return a valid index here.
            idx = self._extname_index(headarr)
            return parse.binning2string(int(headarr[idx]['BINSPEC']),
                                        int(headarr[idx]['BINSPAT']))

        if meta_key == 'target':
            for k in ('TARGET', 'NAME'):
                v = ph.get(k) if ph is not None else None
                if v is not None and str(v).strip():
                    return str(v).strip()
            return ph.get('IMGTYPE') if ph is not None else None

        if meta_key == 'dichroic':
            return None

        if meta_key == 'dispname':
            return 'VPH'

        raise PypeItError(f"Not ready for this compound meta: {meta_key}")

    # -- Detector + raw-image read ------------------------------------
    def get_detector_par(self, det, hdu=None):
        """
        Build the detector container.

        When ``hdu`` is supplied, ``dataext`` is set to the index of
        the per-channel extension in this file (looked up by
        ``EXTNAME``), and ``DATASEC``/``OSCANSEC`` are read from the
        per-frame raw extension header (so they match the actual
        on-chip binning for this exposure).  The sections are given in
        the native readout frame and reoriented alongside the image by
        :meth:`~pypeit.spectrographs.spectrograph.Spectrograph.orient_image`.
        If the raw section is malformed (e.g. U ships ``[-1:1097,...]``),
        :attr:`_datasec_fallback` / :attr:`_oscansec_fallback` are used
        instead.
        """
        if hdu is None:
            binning = '1,1'
            dataext = 1
            datasec = None
            oscansec = None
        else:
            idx = self._extname_index(hdu)
            if idx is None:
                msg = (f"Channel {self._extname!r} not present in raw file."
                       f"  Available EXTNAMEs: "
                       f"{[hdu[i].header.get('EXTNAME','') for i in range(1, len(hdu))]}.")
                raise PypeItError(msg)
            binning = self.get_meta_value(self.get_headarr(hdu), 'binning')
            dataext = idx
            ext_hdr = hdu[idx].header

            datasec_str = self._datasec_string(ext_hdr)
            datasec = (np.atleast_1d(parse.flip_fits_slice(datasec_str))
                       if datasec_str else None)
            oscansec_str = self._oscansec_string(ext_hdr)
            oscansec = (np.atleast_1d(parse.flip_fits_slice(oscansec_str))
                        if oscansec_str else None)

        det_dict = dict(
            binning=binning,
            det=1,
            dataext=dataext,
            specaxis=self._specaxis,
            specflip=self._specflip,
            spatflip=self._spatflip,
            # arcsec / *unbinned* pixel (= 50" slice / 134 binned-px,
            # binspat=2 → 0.373"/binned-px = 0.187"/unbinned-px).
            platescale=0.187,
            darkcurr=0.0,
            saturation=self._saturation_adu,
            nonlinear=40./45.,
            mincounts=-1e10,
            numamplifiers=1,
            gain=np.atleast_1d(self._gain),
            ronoise=np.atleast_1d(self._ronoise),
            datasec=datasec,
            oscansec=oscansec,
        )
        return detector_container.DetectorContainer(**det_dict)

    def get_rawimage(self, raw_file, det):
        """
        Read a raw NGPS image, find the channel by ``EXTNAME``, and
        return the native-frame image plus the DATASEC/OSCANSEC pixel
        maps.  The geometric reorientation into the PypeIt frame is left
        to the stock
        :meth:`~pypeit.spectrographs.spectrograph.Spectrograph.orient_image`
        via the channel's ``specaxis``/``specflip``/``spatflip``.

        This is a custom drop-in replacement for the base class
        :meth:`pypeit.spectrographs.spectrograph.Spectrograph.get_rawimage`
        method, needed only because the channel is located by
        ``EXTNAME`` (not a fixed HDU index).  The ``DATASEC``/``OSCANSEC``
        strings carried by the
        :class:`~pypeit.images.detector_container.DetectorContainer`
        already include the on-chip binning, so we set the per-amplifier
        section binning to ``None`` when converting the strings to
        slices (mirroring ``sec_includes_binning=True`` in the base
        flow).
        """
        self._check_extensions(raw_file)
        hdu = io.fits_open(raw_file, ignore_missing_end=True,
                           output_verify='ignore', ignore_blank=True)

        idx = self._extname_index(hdu)
        if idx is None:
            extnames = [hdu[i].header.get('EXTNAME', '')
                        for i in range(1, len(hdu))]
            raise PypeItError(
                f"Channel {self._extname!r} not present in {raw_file}; "
                f"available EXTNAMEs: {extnames}"
            )

        # Read the native-frame image (orient_image handles geometry).
        raw_img = hdu[idx].data.astype(float)
        if raw_img.ndim != 2:
            raw_img = np.squeeze(raw_img)
            if raw_img.ndim != 2:
                raise PypeItError(
                    f"Raw image extension {idx} of {raw_file} is not 2D."
                )

        # Detector container (with native-frame DATASEC/OSCANSEC)
        detector = self.get_detector_par(det=1, hdu=hdu)

        # Header metadata
        headarr = self.get_headarr(hdu)
        exptime = self.get_meta_value(headarr, 'exptime')

        # Build datasec / oscansec pixel maps
        rawdatasec_img = np.zeros(raw_img.shape, dtype=int)
        oscansec_img = np.zeros(raw_img.shape, dtype=int)
        for sec_name, target in (('datasec', rawdatasec_img),
                                 ('oscansec', oscansec_img)):
            sections = detector[sec_name]
            if sections is None:
                continue
            for j, secstr in enumerate(np.atleast_1d(sections)):
                try:
                    sl = parse.sec2slice(secstr, one_indexed=True,
                                         include_end=True, require_dim=2,
                                         binning=None)
                    target[sl] = j + 1
                except Exception as e:
                    log.warn(f"Failed to parse {sec_name}={secstr!r} for "
                             f"{self.name}: {e}")

        return detector, raw_img, hdu, exptime, rawdatasec_img, oscansec_img

    # -- Default parameters -------------------------------------------
    @classmethod
    def default_pypeit_par(cls):
        """
        Default reduction parameters.  Per-channel overrides
        (``_use_overscan``, ``_bpm_usebias``, ``_wvarxiv``) are pulled
        from class attributes so subclasses only need to set those.
        """
        par = super().default_pypeit_par()

        par['calibrations']['slitedges']['sync_predict'] = 'nearest'
        # Lower than the upstream 50 because the NGPS image-slicer
        # produces relatively low-contrast edges between slices and the
        # inter-slice gaps; 20 (PypeIt's default) finds all 6 edges
        # (left+right of each of the 3 image-slicer slices) reliably.
        par['calibrations']['slitedges']['edge_thresh'] = 20.
        # Median-filter the trace image 3x before the Sobel gradient to
        # wash out single-pixel hot/cold artefacts that occasionally
        # chop a slice into sub-fragments failing the slit-length gate.
        par['calibrations']['slitedges']['filt_iter'] = 3
        # Slicer slices are ~65" long; 30" filters spurious sub-slit
        # detections from inter-slice gaps without dropping real ones.
        par['calibrations']['slitedges']['minimum_slit_length'] = 30
        par['calibrations']['slitedges']['min_edge_side_sep'] = 1.0
        par['calibrations']['slitedges']['fit_min_spec_length'] = cls._fit_min_spec_length
        # NGPS slits are physical image-slicer mirror edges -- straight to
        # sub-pixel.  PypeIt's default fit_order=5 polynomial absorbs ~5-12
        # px of edge-detection noise across spec as fake curvature, then
        # spat_coo = (spat - left)/slitwidth propagates that into the
        # spatial bspline evaluation, painting wavy bands in the illumflat
        # where there should be straight vertical stripes.  Linear (order 1)
        # is enough to capture any small alignment tilt of the slits relative
        # to the detector axes without picking up noise.
        par['calibrations']['slitedges']['fit_order'] = 3
        # Mask low-SNR spec rows from the edge polynomial fit.  Without
        # this, PypeIt's per-row Sobel centroid measurement at G's blue
        # end (where the dome lamp is essentially dead) returns
        # noise-driven centroids that are 5-10 px from the real edge,
        # then the polynomial fits dutifully chase those bad points and
        # produce wavy slit edges that bleed into the illumflat.  With
        # ``trace_thresh`` set, spec rows where the rectified + median-
        # filtered Sobel signal falls below this threshold are dropped
        # from the edge fit entirely.
        par['calibrations']['slitedges']['trace_thresh'] = 10.0
        # exclude_regions intentionally not set: empirically the actual
        # bright bands on the dome flat extend a few pixels past the
        # slice_data.m bounds, so the natural exclude-the-gaps strategy
        # ends up cutting off the slice-edge gradient transitions
        # themselves.  PypeIt's edge_thresh=20 default reliably finds
        # the 3 image-slicer slice edges on the dome flat for binspec=3
        # data without needing exclude_regions.
        # ``fwhm_fromlines=False`` should make wave-cal use the
        # constant ``fwhm`` parameter instead of fitting per-slit FWHM
        # from arc lines.  Upstream currently calls ``map_fwhm``
        # unconditionally, so this switch alone is not sufficient if
        # any slit has no detectable lines; we rely on
        # ``exclude_regions`` above to keep the slits clean.
        par['calibrations']['wavelengths']['fwhm_fromlines'] = False

        par['scienceframe']['process']['combine'] = 'median'
        par['calibrations']['standardframe']['process']['combine'] = 'median'

        # Apply per-channel overscan policy to every frame-type processing
        # block.  U has no usable overscan region (raw header DATASEC is
        # malformed and the strip is only ~2 binned pixels wide), so
        # ``_use_overscan=False`` must take effect for biasframe /
        # arcframe / pixelflatframe processing too, not just science.
        # Spatial flexure correction is enabled for every science /
        # standard frame: PypeIt cross-correlates the dome-flat slit
        # pattern with the slit pattern in each science exposure to
        # measure the spatial shift in pixels and propagates it
        # through extraction.  This is essential at airmass >1.2 where
        # atmospheric refraction shifts the slit relative to the
        # dome-flat reference position.
        for ftype in ('scienceframe', 'biasframe', 'arcframe',
                      'pixelflatframe', 'tiltframe', 'illumflatframe',
                      'traceframe', 'darkframe'):
            try:
                par['calibrations'][ftype]['process']['use_overscan'] = cls._use_overscan
            except KeyError:
                pass
            try:
                par[ftype]['process']['use_overscan'] = cls._use_overscan
            except KeyError:
                pass
        par['calibrations']['standardframe']['process']['use_overscan'] = cls._use_overscan
        par['scienceframe']['process']['spat_flexure_method'] = 'detector'
        par['calibrations']['standardframe']['process']['spat_flexure_method'] = 'detector'

        # Disable LACosmic for standard-star frames only.  Empirically
        # the bright sharp trace of a flux standard (Feige34 at
        # ~30k ADU peak, well below the 57 kADU saturation threshold)
        # gets false-flagged as cosmic-ray hits across ~100 px in the
        # central column on R / I, which makes the displayed spec2d
        # look like the trace centre is "saturated" and feeds bad
        # rows into the optimal-extraction profile fit.  Science
        # frames are typically much fainter and need real CR rejection,
        # so this override is scoped to ``standardframe`` only.  Std
        # exposures are typically short and almost never accumulate
        # true CRs in numbers that matter.
        for key in ('use_lacosmic', 'mask_cr'):
            try:
                par['calibrations']['standardframe']['process'][key] = False
            except KeyError:
                pass

        # Object finding: lower S/N threshold so faint U/G sources are
        # detected, and explicitly use the standard star's trace as
        # the crutch (this is PypeIt's default but we set it
        # explicitly).  ``use_std_trace=True`` means: when an object
        # is detected on a science slit, its spatial trace shape is
        # taken from the standard star's trace, with only the
        # zero-point offset re-fit on the science S/N.  This makes
        # faint-source extraction robust even when the science S/N is
        # too low to fit the trace shape independently.
        par['reduce']['findobj']['use_std_trace'] = True
        par['reduce']['findobj']['snr_thresh'] = 3.0
        par['reduce']['findobj']['maxnumber_sci'] = 1
        par['reduce']['findobj']['maxnumber_std'] = 1

        # ±5 binned-px BOX aperture = 1.87" at NGPS platescale.
        par['reduce']['extraction']['boxcar_radius'] = 1.87
        par['scienceframe']['process']['sigclip'] = 4.0
        par['scienceframe']['process']['objlim'] = 5.0

        par['calibrations']['bpm_usebias'] = cls._bpm_usebias

        par['calibrations']['pixelflatframe']['process']['combine'] = 'median'

        # Push the slicer / slit-edge pattern fully into the
        # illumination model so ``pixelflat_norm`` becomes pure PRNU.
        # Defaults median-filter the slit illumination function across
        # 5 spatial pixels and use 50-pixel b-spline knot spacing in
        # spec, both of which smooth out the high-frequency
        # slicer/slit-edge wiggles and leak them into
        # ``pixelflat_norm`` as visible stripes.  Each channel needs
        # its own tuning -- the per-channel ``_flatfield_overrides``
        # dict provides that.  These values are reasonable starting
        # points; subclasses can override.
        par['calibrations']['flatfield']['spat_samp'] = 1.0
        par['calibrations']['flatfield']['spec_samp_coarse'] = 5.0
        par['calibrations']['flatfield']['spec_samp_fine'] = 1.0
        par['calibrations']['flatfield']['slit_illum_pad'] = 1
        par['calibrations']['flatfield']['slit_illum_smooth_npix'] = 1
        par['calibrations']['flatfield']['illum_iter'] = 3
        par['calibrations']['flatfield']['illum_rej'] = 4.0
        # slit_illum_finecorr OFF across the board (2026-05-10):
        # the bspline-with-sigma-rejection inside the fine-correction
        # 2D Legendre fit is not seeded, so each run lands on a
        # slightly different solution. The drift propagates into
        # PIXELFLAT_NORM (~10 % pixel-level swing between fresh runs
        # observed on 20260501 U) and then into the local sky model,
        # producing the "unpredictable gradients + sharp artefacts"
        # reported across all 4 channels. PypeIt's coarse illumination
        # correction still runs and is deterministic.
        par['calibrations']['flatfield']['slit_illum_finecorr'] = False
        par['calibrations']['flatfield']['tweak_slits'] = True
        # Per-channel overrides win.
        for k, v in (cls._flatfield_overrides or {}).items():
            par['calibrations']['flatfield'][k] = v

        # Per-channel skysub tuning.  Empty by default; subclasses opt
        # in via ``_skysub_overrides`` if a channel benefits from a
        # specific tweak.
        for k, v in (cls._skysub_overrides or {}).items():
            par['reduce']['skysub'][k] = v

        # Per-channel spectral-flexure tuning.  Pypeit defaults
        # ``spec_method='skip'`` (no correction).  Channels with
        # reliable sky-line cross-correlation against a sky template
        # opt in via ``_flexure_overrides`` (typically R and I; U
        # has no sky lines, G has only one strong line).
        for k, v in (cls._flexure_overrides or {}).items():
            par['flexure'][k] = v

        # Per-channel sensfunc tuning.  Top-level keys (``polyorder``,
        # ``samp_fact``, etc.) go in ``par['sensfunc']``; ``UVIS`` and
        # ``IR`` keys go in the matching nested subsection.  Default
        # empty.  pypeit_sensfunc reads the spectrograph's
        # ``config_specific_par`` when invoked without ``--sens_file``,
        # so these defaults flow through automatically.
        for k, v in (cls._sensfunc_overrides or {}).items():
            if isinstance(v, dict) and k in ('UVIS', 'IR'):
                for sk, sv in v.items():
                    par['sensfunc'][k][sk] = sv
            else:
                par['sensfunc'][k] = v

        par['calibrations']['wavelengths']['lamps'] = list(cls._wave_lamps)
        # Force pypeit to use every arc line it can detect, including
        # the weak ones near the spectral edges of each channel.
        # PypeIt's default sigdetect=5 misses edge lines; the
        # polynomial fit then extrapolates wildly past the last
        # matched line.  Setting sigdetect=1 keeps every peak with
        # SNR > 1, including the few faint Th lines at the chip
        # edges, which anchors the polynomial fit out to the data
        # range.
        par['calibrations']['wavelengths']['sigdetect'] = cls._sigdetect
        # Single quartic polynomial fit (no order escalation); 5th-order
        # risks over-fitting slits with sparse anchor lines.
        par['calibrations']['wavelengths']['n_first'] = 4
        par['calibrations']['wavelengths']['n_final'] = 4
        par['calibrations']['wavelengths']['method'] = cls._wavecal_method
        if cls._wvarxiv is not None:
            par['calibrations']['wavelengths']['reid_arxiv'] = cls._wvarxiv
        # Single-template archive => one match per line is sufficient.
        par['calibrations']['wavelengths']['nreid_min'] = 1
        # Loosen the cross-correlation thresholds vs. the PypeIt
        # default 0.8 because (a) the U channel arc has fewer strong
        # lines and lower CC peak by construction, and (b) the
        # NGPS_DRP-derived templates have ~5% flux-amplitude
        # differences vs the on-sky arcs (different lamp current,
        # slightly different optics).
        par['calibrations']['wavelengths']['cc_thresh'] = 0.6
        par['calibrations']['wavelengths']['cc_local_thresh'] = 0.6
        par['calibrations']['wavelengths']['rms_thresh_frac_fwhm'] = 1.0
        # ``match_toler`` in binned pixels; see _match_toler_pixels docs.
        par['calibrations']['wavelengths']['match_toler'] = cls._match_toler_pixels

        par['sensfunc']['algorithm'] = 'UVIS'

        # No exposure-time bounds on any frame type: IMGTYPE +
        # (for cals) the {G,I,R,U} vs {G,U} extension-set check are
        # the physically meaningful discriminators.  Exposure time
        # legitimately varies with binning, slit width, lamp current,
        # and target brightness; any hardcoded numeric range will
        # bite at some configuration sooner or later.

        return par


# ---------------------------------------------------------------------
# Per-channel subclasses
# ---------------------------------------------------------------------


class P200NGPSSpectrograph_r(P200NGPSSpectrograph):
    """P200/NGPS r-channel."""
    name = 'p200_ngps_r'
    camera = 'NGPS_r'
    header_name = 'NGPS_r'
    supported = True
    comment = 'r-Channel'

    _extname = 'R'
    # Disable overscan correction: NGPS overscan strips are too narrow
    # for the default ``savgol(5,65)``; bias master handles it cleanly.
    _use_overscan = False
    _bpm_usebias = True
    _wvarxiv = 'wvarxiv_p200_ngps_r_thar_central.fits'
    _gain = 0.9   # e-/ADU
    _flatfield_overrides = dict(
        # Dome-lamp signal usable between ~5850 and ~7800 A on r;
        # NORM=1 outside.
        pixelflat_min_wave=5850.0,
        pixelflat_max_wave=7800.0,
    )


class P200NGPSSpectrograph_i(P200NGPSSpectrograph):
    """P200/NGPS i-channel."""
    name = 'p200_ngps_i'
    camera = 'NGPS_i'
    header_name = 'NGPS_i'
    supported = True
    comment = 'i-Channel'

    _extname = 'I'
    _use_overscan = False
    _bpm_usebias = False
    _wvarxiv = 'wvarxiv_p200_ngps_i_thar_central.fits'
    _gain = 0.86  # e-/ADU
    _flatfield_overrides = dict(
        # Dome-lamp signal usable between ~7800 and ~10100 A on i;
        # NORM=1 outside.
        pixelflat_min_wave=7800.0,
        pixelflat_max_wave=10100.0,
    )


class P200NGPSSpectrograph_g(P200NGPSSpectrograph):
    """P200/NGPS g-channel.  Read out dispersion-along-NAXIS1 like R/I
    but left-right flipped, i.e. ``specaxis=1, specflip=True``.
    """
    name = 'p200_ngps_g'
    camera = 'NGPS_g'
    header_name = 'NGPS_g'
    supported = True
    comment = 'g-Channel'

    _extname = 'G'
    _specflip = True
    _require_ug_only_cals = True
    # Disable overscan correction for G.  PypeIt's default ``savgol(5,65)``
    # introduces a smooth row-dependent offset (~-540 ADU at the blue end,
    # decreasing to ~-22 ADU at the red end) -- visible as the discontinuity
    # at row ~150 of the rawflat where the dome lamp emerges from the noise
    # floor, and breaking PypeIt's spline-based flat fit.  With overscan
    # disabled, the full ~1009 ADU bias level is removed by the bias master
    # alone; manual tests confirm flat - bias gives ~0 ADU off-slit cleanly.
    _use_overscan = False
    _bpm_usebias = False
    _wvarxiv = 'wvarxiv_p200_ngps_g_thar_central.fits'
    _gain = 2.88  # e-/ADU, from NGPS_DRP/slice_data.m
    _flatfield_overrides = dict(
        # Finecorr disabled across the board (2026-05-10) -- the
        # bspline-with-sigma-rejection is non-deterministic and the
        # ~10 % swing in PIXELFLAT_NORM was driving the sky-sub
        # gradients we were chasing.
        slit_illum_finecorr=False,
        spec_samp_coarse=1.0,
    )


class P200NGPSSpectrograph_u(P200NGPSSpectrograph):
    """P200/NGPS u-channel.  Read out with the dispersion axis along
    NAXIS2 (rotated 90 deg relative to the other channels) and flipped,
    i.e. ``specaxis=0, specflip=True``.  The raw FITS header ships
    malformed ``DATASEC``/``BIASSEC`` cards (e.g. ``[-1:1097,1:4302]``),
    so we use the full native frame as the data section and leave
    ``_use_overscan=False`` since the bias strip is only 2 pixels wide.
    """
    name = 'p200_ngps_u'
    camera = 'NGPS_u'
    header_name = 'NGPS_u'
    supported = True
    comment = 'u-Channel (read out rotated 90 deg)'

    _extname = 'U'
    _specaxis = 0
    _specflip = True
    _require_ug_only_cals = True
    _use_overscan = False
    # U dome-flat illumination drops off toward the blue end, so slit
    # edge traces don't span 60% of the spectral length.
    _fit_min_spec_length = 0.2
    # Use the full native image as the data section; bias is unused
    # (master bias handles it).
    _datasec_fallback = None
    _oscansec_fallback = None
    _bpm_usebias = False
    _wvarxiv = 'wvarxiv_p200_ngps_u_thar_central.fits'
    _wave_lamps = ['ThAr']
    _arc_idnames = ('THAR',)
    _gain = 0.755  # e-/ADU, from NGPS_DRP/slice_data.m
    _flatfield_overrides = dict(
        pixelflat_min_wave=3400.0,
        spec_samp_coarse=2.0,
        slit_illum_finecorr=False,
    )
    def _datasec_string(self, ext_hdr):
        """U's raw DATASEC is malformed at every binning (negative low
        bounds), so we ignore it and use the entire native image
        (``[1:NAXIS1,1:NAXIS2]``); ``orient_image`` reorients it."""
        n1 = ext_hdr.get('NAXIS1')
        n2 = ext_hdr.get('NAXIS2')
        if not n1 or not n2:
            return self._datasec_fallback
        return f'[1:{n1},1:{n2}]'

    def _oscansec_string(self, ext_hdr):
        """No usable overscan strip for U (the bias region is only
        2 binned pixels wide and the raw card is malformed).  We rely
        on the master bias frame instead."""
        return self._oscansec_fallback




