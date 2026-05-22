"""
Module for P200/NGPS specific methods.

NGPS raw FITS files store each channel ('U', 'G', 'R', 'I') in its own
image extension labelled by ``EXTNAME``. The order of the image
extensions is *not* stable across files, and frames where only some
channels were read out (e.g. U+G dome flats) drop the unused
extensions entirely. Per-channel subclasses below resolve their image
extension by scanning ``EXTNAME``, never by hardcoded HDU index.

Per-channel orientation conventions (same as the NGPS MATLAB DRP's
``slice_data.m``):

* R, I: native orientation (dispersion along NAXIS1).
* G:    flipped left-right (``np.fliplr``) so that wavelengths increase with
        column number, matching R/I.
* U:    rotated 90 degrees clockwise (``np.rot90(k=-1)``) so that
        dispersion runs along NAXIS1, matching R/I.

The transforms are applied inside :meth:`get_rawimage` so that the rest
of the PypeIt pipeline only ever sees data in the canonical
"dispersion-along-axis-1, spatial-along-axis-0" orientation, with
hardcoded ``DATASEC``/``OSCANSEC`` regions describing the post-transform
image.

.. include:: ../include/links.rst
"""
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
    image extension in the raw FITS file ('U', 'G', 'R', or 'I') and
    optionally overrides ``_canonicalize`` to apply a per-channel
    geometric transform (no-op for R/I, ``fliplr`` for G, ``rot90(-1)``
    for U).

    All subclasses present their data to PypeIt with dispersion along
    array axis 1 (i.e. ``specaxis=1``) regardless of native CCD
    orientation, so that downstream reduction (slit edge tracing,
    wavelength calibration, sky subtraction) sees a uniform layout.
    """
    ndet = 1
    telescope = telescopes.P200TelescopePar()

    # -- Per-channel overrides -----------------------------------------
    #: EXTNAME of this channel in the raw FITS HDU list.
    _extname = None
    #: Whether to use the bias overscan region in raw images.
    _use_overscan = True
    #: Hardcoded fallback ``DATASEC`` (FITS-style ``[NAXIS1, NAXIS2]``)
    #: used only when the per-frame extension header has no usable
    #: ``DATASEC`` card (e.g. the U channel ships malformed
    #: ``[-1:1097,...]``).  When ``None``, the per-frame raw header
    #: values are used and adjusted for the canonical-orientation
    #: transform applied in :meth:`_canonicalize`.
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
    #: single-slice central-slit archive (built from a known-good
    #: NGPS_Pypeit WaveCalib) and a global cross-correlation shift to
    #: anchor the per-slit polynomial fit.  Holds the central slice
    #: (which drives the merged 1D output via ngps_pipeline.finalize)
    #: to sub-0.3 A RMS / <0.6 A peak vs the source wave solution on
    #: 20260501 -- no framework patches required.
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
    #: Per-channel overrides for ``par['reduce']['skysub']``.  Each NGPS
    #: channel has different sky-line density (i is bright with OH/H2O,
    #: r is moderate, g and u have nearly no sky lines), so the optimal
    #: bspline knot spacing and rejection threshold differ per channel.
    #: Determined empirically via diagnostics/sweep_extraction.py with
    #: SNR + 2D residual-quality metrics on UT 20260501 / iter.  Empty
    #: by default; subclasses fill in.
    #:
    #:   - ``bspline_spacing = 0.5`` (vs PypeIt default 0.6): tighter
    #:     knot spacing along the spectral axis.  Validated by the
    #:     i and r 18-/9-variant sweeps in
    #:     ``experiments/aajunnc_{i,r}_iter/clean_sweep_table.md``.
    #:   - ``no_local_sky = True``: skip the local-window sky re-fit
    #:     so the spec2d's sky-subtracted image is continuous (no
    #:     discontinuities at the local-window edges).
    _skysub_overrides = {
        'bspline_spacing': 0.5,
        'no_local_sky': True,
    }
    #: Per-channel overrides for ``par['flexure']``.  Empty by default
    #: (= ``spec_method='skip'``).  Channels with a rich enough sky-line
    #: forest to make a paranal-template cross-correlation reliable
    #: opt in to ``spec_method='boxcar'`` with the air-wavelength
    #: template.  u and g leave the default skip because u has no
    #: significant atmospheric emission below 4400 A and g has only
    #: 5577 OI as a strong line, so cross-correlation is noise-driven
    #: or single-line-locked respectively.
    _flexure_overrides = {}
    #: Per-channel overrides for ``par['sensfunc']``.  Top-level keys go
    #: in ``par['sensfunc']``; nested ``UVIS`` / ``IR`` dicts go in
    #: ``par['sensfunc']['UVIS']`` / ``['IR']``.
    #:
    #: ``polycorrect = False`` and ``resolution = 4000`` were chosen by
    #: the 20260501 i-band sensfunc parameter sweep
    #: (``diagnostics/sweep_sensfunc_i.md``); they drop the
    #: telluric-corrected Fritz-match RMS from 17.3% to 11.5% on i, with
    #: no regression on g/r/u.
    #:
    #: Why ``resolution = 4000`` is constant rather than slit-aware:
    #: the bspline knot spacing in
    #: :func:`pypeit.core.flux_calib.standard_zeropoint` is
    #: ``bkspace ~ (lambda / R) * nresln``; this is a knot-density
    #: target, NOT the physical spectrograph R despite the name.  Lower
    #: R -> sparser knots -> bspline averages through telluric structure
    #: in the std/CALSPEC ratio and biases the flux-cal in those
    #: windows.  Empirically the i-band Fritz residual is monotonically
    #: better as we go from R=3000 (PypeIt default) to 4000 then
    #: saturates by 5000-6000.  All four NGPS channels' actual measured
    #: R is in [3900, 4300] (project update slide 9, 26 Feb 2026), so
    #: 4000 is a sensible single all-channel anchor.
    _sensfunc_overrides = dict(
        extr='OPT',
        UVIS=dict(polycorrect=False, resolution=4000),
    )

    #: PypeIt's ``reidentify``/``full_template`` drops detected arc peaks
    #: whose nearest catalog line is more than ``match_toler`` *pixels*
    #: away.  0.3 binned pixels at NGPS' binspec=3 corresponds to a
    #: ~0.5-0.6 A window across u/g/r/i -- tight enough that the
    #: order-4 polynomial can't absorb mis-identifications, loose
    #: enough that real edge lines stay anchored.  Selected by the
    #: 60-point sweep (see ``n_first/n_final`` comment in
    #: ``default_pypeit_par``).
    _match_toler_pixels = 0.3
    #: Per-channel detector gain (electrons / ADU).  PypeIt multiplies
    #: this into the science ADU image at spec2d build time so the
    #: sensfunc and downstream products are in electrons.  Values come
    #: from the matlab DRP's ``slice_data.m``.  Each subclass sets this.
    _gain = 1.0
    #: Per-channel readnoise (electrons).  Currently the same upstream
    #: default for all four; revisit when proper per-channel
    #: characterizations are available.
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

    def _canonicalize(self, raw):
        """
        Apply the per-channel geometric transform that puts the raw
        image into the canonical "dispersion along axis 1" orientation.
        Default: identity (R, I).  G overrides with ``np.fliplr``; U
        overrides with ``np.rot90(k=-1)``.
        """
        return raw

    def _canonicalize_section(self, sec_str, ext_hdr):
        """
        Adjust a FITS section string ``[NAXIS1_lo:NAXIS1_hi, NAXIS2_lo:NAXIS2_hi]``
        from the *raw* extension header into the canonical post-transform
        coordinate system (matches what :meth:`_canonicalize` produces).

        Default: identity (R, I).  G and U override.
        """
        return sec_str

    def _datasec_string(self, ext_hdr):
        """Resolve the canonical-orientation DATASEC for this extension.

        Default (R/I/G): try the per-frame raw DATASEC card via
        ``_canonicalize_section`` (this picks up the binning-correct
        section the instrument writes out); fall back to
        ``_datasec_fallback`` if the raw section is missing or
        malformed.  Subclasses (e.g. U, where every raw DATASEC is
        malformed) can override this method to bypass the parse check.
        """
        raw = ext_hdr.get('DATASEC')
        if self._parse_section(raw) is not None:
            return self._canonicalize_section(raw, ext_hdr)
        return self._datasec_fallback

    def _oscansec_string(self, ext_hdr):
        """Same as :meth:`_datasec_string` but for OSCANSEC/BIASSEC."""
        raw = ext_hdr.get('BIASSEC')
        if self._parse_section(raw) is not None:
            return self._canonicalize_section(raw, ext_hdr)
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
        # ``SHUTTIME`` is the actual shutter-open time in seconds; the
        # primary-header ``EXPTIME`` is in *milliseconds* for some
        # frames (it's the requested integration time written by the
        # detector controller).  Use SHUTTIME when present, otherwise
        # fall back to EXPTIME with an automatic ms-to-s conversion
        # if the value looks like ms (>1e4).  Implemented in
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
        return ['binning']

    def raw_header_cards(self):
        return ['GRATING', 'ANGLE', 'APERTURE']

    def pypeit_file_keys(self):
        return super().pypeit_file_keys()

    _standard_targets = (
        'feige34', 'feige110',
        'bd+28', 'bd+284211', 'bd+28d4211',
        'bd+25', 'bd+33',
        'g191b2b', 'g191-b2b',
        'hz44',
    )

    def _standard_target_mask(self, fitstbl):
        """Return True for recognized spectrophotometric standard targets."""
        if 'target' not in fitstbl.keys():
            return np.zeros(len(fitstbl), dtype=bool)

        mask = np.zeros(len(fitstbl), dtype=bool)
        for i, target in enumerate(fitstbl['target']):
            if target is None:
                continue
            norm = str(target).strip().lower().replace(' ', '')
            if not norm or norm in ('none', '--'):
                continue
            mask[i] = any(std in norm for std in self._standard_targets)
        return mask

    # True on U and G so their arcs/flats are restricted to the long
    # getcalib_ug frames (image extensions == {G, U}).  R/I use the
    # short 4-channel arcs/flats, so this stays False.
    _require_ug_only_cals = False
    _ug_only_cache: dict = {}

    @classmethod
    def _is_ug_only_cal(cls, path) -> bool:
        """True iff this FITS's image extensions are exactly {G, U}."""
        cached = cls._ug_only_cache.get(str(path))
        if cached is not None:
            return cached
        from astropy.io import fits as _fits
        try:
            with _fits.open(path, memmap=True) as h:
                extnames = {str(hdu.header.get('EXTNAME', '')).strip().upper()
                            for hdu in h[1:]}
        except Exception:
            extnames = set()
        result = extnames == {'G', 'U'}
        cls._ug_only_cache[str(path)] = result
        return result

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Frame-typing rule, shared by all four channels.  IMGTYPE values
        in the raw header drive the baseline classification: ``SCI`` is
        science / standard, ``BIAS`` is bias, ``DOMEFLAT`` and ``CONT``
        are flats, and channel-specific arc IMGTYPE values define the
        arc/tilt frames.  Science and standard frames are separated by
        known standard-star target names; otherwise short science frames
        like the 60 s NGC 2440 validation exposure get mis-typed as
        standards.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)

        if ftype in ['science', 'standard']:
            sci = fitstbl['idname'] == 'SCI'
            std_target = self._standard_target_mask(fitstbl)
            return good_exp & sci & (std_target if ftype == 'standard' else ~std_target)
        if ftype == 'bias':
            return good_exp & (fitstbl['idname'] == 'BIAS')
        if ftype in ['pinhole', 'dark']:
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc', 'tilt']:
            mask = good_exp & np.isin(fitstbl['idname'], self._arc_idnames)
        elif ftype in ['pixelflat', 'trace', 'illumflat']:
            mask = ((good_exp & (fitstbl['idname'] == 'DOMEFLAT'))
                    | (good_exp & (fitstbl['idname'] == 'CONT')))
        else:
            log.debug('Cannot determine if frames are of type {0}.'.format(ftype))
            return np.zeros(len(fitstbl), dtype=bool)

        if self._require_ug_only_cals:
            import os as _os
            for i in range(len(fitstbl)):
                if not mask[i]:
                    continue
                fpath = _os.path.join(str(fitstbl['directory'][i]),
                                       str(fitstbl['filename'][i]))
                if not self._is_ug_only_cal(fpath):
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
            # Prefer SHUTTIME (seconds); fall back to EXPTIME with an
            # ms-to-s conversion if the value is large enough that it
            # can't be a plausible visible-band exposure in seconds
            # (ZTF/SN spectra are 60-3600 s; values >= 10000 are ms).
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
                        v = float(ex)
                        return v / 1000.0 if v >= 10000 else v
                    except (TypeError, ValueError):
                        pass
            return None

        if meta_key == 'dispangle':
            return 0

        if meta_key == 'binning':
            # Raw NGPS frames carry BINSPAT/BINSPEC on the per-channel
            # image extension (resolved by EXTNAME).  Pypeit-written
            # products (spec1d, sens) drop those extensions but write
            # the canonical ``BINNING`` string in the primary header.
            # Try the per-channel extension first (raw); if absent,
            # fall back to the primary BINNING card (pypeit products);
            # only return the '0,0' sentinel if both fail, signalling
            # a frame the configuration code should drop.
            idx = self._extname_index(headarr)
            if idx is not None:
                try:
                    binspat = int(headarr[idx]['BINSPAT'])
                    binspec = int(headarr[idx]['BINSPEC'])
                    return parse.binning2string(binspec, binspat)
                except (KeyError, IndexError, TypeError, ValueError):
                    pass
            if ph is not None:
                v = ph.get('BINNING')
                if v is not None and str(v).strip() not in ('', '0,0'):
                    return str(v).strip()
            return '0,0'

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
        on-chip binning for this exposure).  Per-channel subclasses
        can override :meth:`_canonicalize_section` to translate those
        raw-orientation sections into the canonical post-transform
        frame; if the raw section is malformed (e.g. U ships
        ``[-1:1097,...]``), :attr:`_datasec_fallback` /
        :attr:`_oscansec_fallback` are used instead.
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
            specaxis=1,
            specflip=False,
            spatflip=False,
            # ``platescale`` is in arcsec / *unbinned* pixel.  Each
            # NGPS image-slicer slice spans ~50" on sky and projects
            # onto ~134 binned-px (binspat=2) on the detector, giving
            # 50/134 = 0.373"/binned-px = 0.187"/unbinned-px.  PypeIt
            # uses this value to convert ``boxcar_radius`` (arcsec)
            # to pixels for BOX extraction; with the correct value,
            # the default ``boxcar_radius=1.5"`` lands at ±4 binned-px,
            # matching the matlab DRP's 8-binned-px box.  An older
            # incorrect 0.5"/unbinned was inherited from upstream and
            # made ``boxcar_radius`` 2.7x too small in pixels, so any
            # arcsec-based extraction tuning behaved unexpectedly.
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
        Read a raw NGPS image, find the channel by ``EXTNAME``, apply
        the per-channel canonical transform (no-op for R/I, ``fliplr``
        for G, ``rot90(-1)`` for U), and return the raw image plus the
        DATASEC/OSCANSEC pixel maps in the post-transform orientation.

        This is a custom drop-in replacement for the base class
        :meth:`pypeit.spectrographs.spectrograph.Spectrograph.get_rawimage`
        method.  The ``DATASEC``/``OSCANSEC`` strings carried by the
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

        # Read & canonicalize
        raw_img = hdu[idx].data.astype(float)
        if raw_img.ndim != 2:
            raw_img = np.squeeze(raw_img)
            if raw_img.ndim != 2:
                raise PypeItError(
                    f"Raw image extension {idx} of {raw_file} is not 2D."
                )
        raw_img = self._canonicalize(raw_img)

        # Detector container (with hardcoded post-transform DATASEC/OSCANSEC)
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
        # Install the per-slice rawflat exporter on first call
        # (deferred from module-load to avoid the
        # ``pypeit.flatfield`` circular import).
        _install_ngps_rawflat_exporter()

        par = super().default_pypeit_par()

        par['calibrations']['slitedges']['sync_predict'] = 'nearest'
        # Lower than the upstream 50 because the NGPS image-slicer
        # produces relatively low-contrast edges between slices and the
        # inter-slice gaps; 20 (PypeIt's default) finds all 6 edges
        # (left+right of each of the 3 image-slicer slices) reliably.
        par['calibrations']['slitedges']['edge_thresh'] = 20.
        # Median-filter the trace image 3x before the Sobel gradient is
        # computed.  Single-pixel hot columns and ~1-px cold/dim
        # features (seen empirically on the G detector at spat~92 and
        # spat~160 with 2x1 binning -- documented in
        # ``diagnostics/g_slit0_dropout.md``) were chopping slit 0
        # into <80-px sub-fragments that then failed the
        # ``minimum_slit_length`` gate, causing pypeit to drop the
        # leftmost slice entirely.  Smoothing 3 iterations washes
        # those out while shifting edge centroids on clean (2x3) data
        # by <0.15 binned px (verified on 20260501).  Used ONLY for
        # slit-edge detection in EdgeTraceSet.initial_trace; the
        # master Flat / wavecal / tilts / science extraction all see
        # the unsmoothed image.
        par['calibrations']['slitedges']['filt_iter'] = 3
        # Each of the 3 image-slicer slices is ~65 arcsec long along
        # the spatial direction (~130 binned pixels at platescale
        # 0.5 arcsec/binned-pix), so 30 arcsec is well under the real
        # value.  Empirically, PypeIt's gradient edge detector also
        # finds occasional 5-15-pixel-wide spurious slits sitting in
        # the inter-slice gaps; a 30-arcsec minimum filters those out
        # while always letting the 3 real slices through.
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
        par['scienceframe']['process']['spat_flexure_correct'] = True
        par['calibrations']['standardframe']['process']['spat_flexure_correct'] = True

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

        # Match the matlab DRP's standard-star BOX width: matlab uses
        # ±5 binned-px (= 11-px aperture, see
        # NGPS_DRP/extract_2dbox_trace_fast.m line 79 with
        # tracew=width_extract=5; called with ``(14, 5, ...)``).  At
        # platescale 0.187"/unbinned-px = 0.374"/binned-px, ±5 binned-px
        # corresponds to boxcar_radius = 5*0.374 = 1.87".  Pypeit only
        # supports a single global ``boxcar_radius``; this matches the
        # std side, while science extractions in matlab were chosen
        # case-by-case (1.5"/4 binned-px for the 20260501 reduction).
        # The width difference between std and sci in matlab introduces
        # a small (~5-10%) flux-scale bias that pypeit-with-1.87 cannot
        # exactly reproduce without per-frametype boxcar_radius support.
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

        # Extraction: pypeit defaults (use_user_fwhm=False, sn_gauss=4)
        # so the central-slit profile is measured from data and OPT
        # extraction uses the true wing shape.  The external
        # ``ngps_pipeline.finalize`` post-processor (shipped
        # alongside the NGPS_Pypeit ``reduce.py``, NOT inside the
        # pypeit/ tree) currently uses ONLY the central slit per
        # channel for the merged 1D, so side-slice profile-fit
        # quality does not affect the delivered spectrum.  Combining
        # all three slicer slices (matlab DRP behaviour) requires
        # custom code outside pypeit_coadd_1dspec.

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
        # n_first=4, n_final=4: a quartic polynomial fit with no
        # order escalation.  Selected by a 60-point (match_toler,
        # n_first, n_final, method) sweep on 20260501 cals (see
        # work/wavecal_sweep_results/summary_central_rms.csv) which
        # measured central-slice Δλ vs the NGPS_Pypeit-reference
        # WaveCalib: this combo gave the lowest u-band RMS
        # (0.07 A) and the lowest r-band RMS (0.17 A), with g/i
        # comparable to the other top combos (g=0.03, i=0.26 A).
        # Sticking to a single quartic order avoids 5th-order over-
        # fitting risk on slits with sparse anchor lines.
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
    # Disable overscan correction across all NGPS channels.  PypeIt's
    # default ``savgol(5,65)`` overscan introduces a row-dependent
    # offset given how narrow the NGPS overscan strips are; the
    # bias master alone (built without overscan) handles bias
    # subtraction cleanly.
    _use_overscan = False
    _bpm_usebias = True
    _wvarxiv = 'wvarxiv_p200_ngps_r_thar_central.fits'
    _gain = 0.9   # e-/ADU, from NGPS_DRP/slice_data.m
    _flatfield_overrides = dict(
        # Empirically: dome-lamp signal stays above ~30% of peak between
        # ~5850 and ~7800 A; outside that range the b-spline fit is
        # unreliable and produces visible structure in pixelflat_norm
        # at the spectral edges.  Force NORM=1 outside.
        pixelflat_min_wave=5850.0,
        pixelflat_max_wave=7800.0,
    )
    # R-channel sky-sub historical notes (kept for future tuning when
    # we revisit local vs no-local trade-offs; the NGPS-wide default
    # set on the base class is currently ``no_local_sky=True`` because
    # the GUI needs a discontinuity-free 2D for custom extractions):
    #
    #   bsp=0.5 + npoly=6 + no_local_sky=True was the validated combo.
    #   npoly=6 is required (npoly=auto=3 is too rigid for R's
    #   spat-direction sky illumination); the higher Legendre order
    #   lets the global b-spline track the smooth ~20-count spat
    #   variation visible in column-medians of off-object pixels.
    #
    #   no_local_sky=True: the local skysub re-fits in a window
    #   around the object using bspline, which fights the global model
    #   and reintroduces oscillations near the right slit edge
    #   (spat 290-310 on the 20260501 aajunnc_o1 commissioning frame).
    #
    #   bsp=0.5 is the sweet spot: at bsp<=0.4 with npoly=6 the
    #   bspline silently fails the 35-iter convergence budget and
    #   falls back to no_poly (grep "Maximum iterations reached in
    #   bspline_profile" in run_pypeit.log).  At bsp=0.5/n6 all 3
    #   slits converge cleanly.
    #
    #   Validated by the 9-variant r-channel sweep in
    #   ``experiments/aajunnc_r_iter/clean_sweep_table.md`` (best
    #   clean MAD=26.05 vs 27.57 at pypeit defaults).
    # R has many strong OI/OH sky lines (5577, 6300, 6363, 7340, 7794),
    # so spec-flexure cross-correlation against the air-wavelength
    # paranal template is reliable.  spec_maxshift=2 px guards against
    # any spuriously large measurement; if exceeded, the
    # excessive_shift='use_median' default falls back to the median of
    # other slits' shifts.
    _flexure_overrides = dict(
        spec_method='boxcar',
        spectrum='paranal_sky_air.fits',
        spec_maxshift=2,
    )


class P200NGPSSpectrograph_i(P200NGPSSpectrograph):
    """P200/NGPS i-channel."""
    name = 'p200_ngps_i'
    camera = 'NGPS_i'
    header_name = 'NGPS_i'
    supported = True
    comment = 'i-Channel'

    _extname = 'I'
    # Disable overscan correction (see G subclass for details).
    _use_overscan = False
    _bpm_usebias = False
    _wvarxiv = 'wvarxiv_p200_ngps_i_thar_central.fits'
    _gain = 0.86  # e-/ADU, from NGPS_DRP/slice_data.m
    _flatfield_overrides = dict(
        # I-band lamp signal is strong from ~7800 to ~10100 A; below
        # 7800 it's a steep climb in throughput and above 10100 the
        # lamp drops off rapidly.
        pixelflat_min_wave=7800.0,
        pixelflat_max_wave=10100.0,
    )
    # I-channel sky-sub historical notes:
    #   Same shape as R.  bsp=0.5 + npoly=6 minimises off-source 2D
    #   MAD around the dense 8200-10000 A OH forest; npoly=auto (=3)
    #   under-fits the spat-direction sky gradient and
    #   no_local_sky=True prevents the local re-fit from
    #   re-introducing OH-line residuals under the source trace.
    #   Validated by the 18-variant i-channel sweep in
    #   ``experiments/aajunnc_i_iter/clean_sweep_table.md`` (best
    #   clean MAD=25.47 at bsp=0.5/n6 vs ~27 at pypeit defaults;
    #   no_poly variants run at MAD~32, conclusively worse).
    # I has the densest sky-line forest (OH bands across 8200-10000 A
    # plus the H2O 9300 envelope).  Cross-correlation against the air
    # paranal template is reliable; on UT 20260501 the central-slit
    # trace contamination at sky lines drops 3.8sigma -> 1.1sigma when
    # this is enabled.
    _flexure_overrides = dict(
        spec_method='boxcar',
        spectrum='paranal_sky_air.fits',
        spec_maxshift=2,
    )


class P200NGPSSpectrograph_g(P200NGPSSpectrograph):
    """P200/NGPS g-channel.  Raw image is left-right flipped relative
    to R/I, so :meth:`_canonicalize` applies ``np.fliplr``.  After the
    flip, the bias overscan moves from the left edge to the right edge.
    """
    name = 'p200_ngps_g'
    camera = 'NGPS_g'
    header_name = 'NGPS_g'
    supported = True
    comment = 'g-Channel'

    _extname = 'G'
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
    # G-channel sky-sub historical notes:
    #   No per-channel deep sweep was run on G (the only test target
    #   so far has very few sky lines in the G window), but the
    #   bsp-only UGRI sweep in
    #   ``diagnostics/sweep_skysub_all_channels.md`` showed G prefers
    #   bsp<=0.6 (chi^2 around 4861/5577 climbs monotonically above
    #   0.6), and the slit geometry is identical to R/I.
    _require_ug_only_cals = True

    def _canonicalize(self, raw):
        return np.fliplr(raw)

    def _canonicalize_section(self, sec_str, ext_hdr):
        """G applies ``np.fliplr`` to the raw image, so the FITS
        ``DATASEC`` ``[x_lo:x_hi, y_lo:y_hi]`` (in raw NAXIS1 / NAXIS2
        coords) reflects to ``[(NAXIS1-x_hi+1):(NAXIS1-x_lo+1), y_lo:y_hi]``
        in the canonical post-flip frame."""
        parsed = self._parse_section(sec_str)
        if parsed is None:
            return sec_str
        x_lo, x_hi, y_lo, y_hi = parsed
        naxis1 = ext_hdr.get('NAXIS1') or (x_hi + 1)
        new_lo = naxis1 - x_hi + 1
        new_hi = naxis1 - x_lo + 1
        return f'[{new_lo}:{new_hi},{y_lo}:{y_hi}]'


class P200NGPSSpectrograph_u(P200NGPSSpectrograph):
    """P200/NGPS u-channel.  Raw image has the dispersion axis along
    NAXIS2 (rotated 90 deg relative to the other channels), so
    :meth:`_canonicalize` applies ``np.rot90(k=-1)``.  The raw FITS
    header ships malformed ``DATASEC``/``BIASSEC`` cards (e.g.
    ``[-1:1097,1:4302]``), so we hardcode the post-rotation regions
    and leave ``_use_overscan=False`` since the bias strip is only
    2 pixels wide.
    """
    name = 'p200_ngps_u'
    camera = 'NGPS_u'
    header_name = 'NGPS_u'
    supported = True
    comment = 'u-Channel (raw image rotated 90 deg)'

    _extname = 'U'
    _use_overscan = False
    # U dome-flat illumination drops off toward the blue end, so slit
    # edge traces don't span 60% of the spectral length.
    _fit_min_spec_length = 0.2
    # Post-rot90(-1) shape is (1099, 4302).  Use the full image as the
    # data section; bias is unused (master bias handles it).
    _datasec_fallback = '[1:4302,1:1099]'
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
    # U-channel sky-sub historical notes:
    #   No deep validation possible on U yet -- the test target is at
    #   the noise floor in the 3300-4200 A range, so the bsp-only
    #   UGRI sweep (``diagnostics/sweep_skysub_all_channels.md``)
    #   couldn't separate bsp values at any chi^2 level.
    _require_ug_only_cals = True

    def _canonicalize(self, raw):
        return np.rot90(raw, k=-1)

    def _datasec_string(self, ext_hdr):
        """U's raw DATASEC is malformed at every binning (negative low
        bounds), so we ignore it and return the entire post-rot90(-1)
        image based on NAXIS1 / NAXIS2 of the raw extension.  After
        ``np.rot90(k=-1)``, post-transform NAXIS1 = raw NAXIS2 (now
        spectral) and NAXIS2 = raw NAXIS1 (now spatial)."""
        n1 = ext_hdr.get('NAXIS1')
        n2 = ext_hdr.get('NAXIS2')
        if not n1 or not n2:
            return self._datasec_fallback
        return f'[1:{n2},1:{n1}]'

    def _oscansec_string(self, ext_hdr):
        """No usable overscan strip for U (the bias region is only
        2 binned pixels wide and the raw card is malformed).  We rely
        on the master bias frame instead."""
        return self._oscansec_fallback




# ---------------------------------------------------------------------
# NGPS per-slice rawflat exporter
# ---------------------------------------------------------------------
# Hooks ``FlatField.run`` to dump, BEFORE PypeIt does any flat fitting:
#   - the rawflat image
#   - the wavelength image
#   - the slit-id image (which slit each pixel belongs to)
#   - one extension per slit with the cropped per-slit data
# to a fits file in the Calibrations directory.  PypeIt's normal flat
# processing then runs unchanged.  The export gives the user a clean
# per-slice view of the raw flat data + slit shape, so they can
# design a method to fix the bad-blue rows externally.
#
# Activated lazily from ``P200NGPSSpectrograph.default_pypeit_par``.

#: Channels to export per-slice.  Restrict to the channels with
#: known bad-blue regions (G, U).  R, I are skipped.
_NGPS_EXPORT_CHANNELS = ('p200_ngps_g', 'p200_ngps_u')


def _ngps_export_rawflat_per_slice(self):
    """Save the rawflat / waveimg / slit info to a per-slice fits
    file.  Called from the patched ``FlatField.run`` BEFORE PypeIt's
    own flat fitting runs.
    """
    from pathlib import Path
    rawflat = self.rawflatimg.image if self.rawflatimg is not None else None
    waveimg = self.waveimg
    slits   = self.slits
    if rawflat is None or waveimg is None or slits is None:
        return
    if rawflat.shape != waveimg.shape:
        return

    # Output to the experiments/ directory at the NGPS_Pypeit project
    # root so the user can iterate on the patch design separately
    # from PypeIt's calibration outputs.
    qa_path = getattr(self, 'qa_path', None)
    if qa_path is None:
        return
    p = Path(qa_path).resolve()
    out_dir = None
    for parent in p.parents:
        cand = parent / 'experiments'
        if cand.is_dir():
            out_dir = cand
            break
    if out_dir is None:
        out_dir = Path(qa_path).parent.parent / 'experiments'
        out_dir.mkdir(parents=True, exist_ok=True)
    calib_key = getattr(self, 'calib_key', None) or 'export'
    spec_name_short = (getattr(self.spectrograph, 'name', 'p200_ngps')
                        .replace('p200_ngps_', ''))
    out = out_dir / f'rawflat_export_{spec_name_short}_{calib_key}.fits'

    spec_name = getattr(self.spectrograph, 'name', 'p200_ngps')
    nspec, nspat = rawflat.shape

    # Build the slit-id image: -1 off-slit, slit_spat_id on-slit.
    slitid_img = np.full(rawflat.shape, -1, dtype=np.int32)
    spat_ids = []
    for slit_idx in range(slits.nslits):
        if slits.mask[slit_idx] != 0:
            continue
        sid = int(slits.spat_id[slit_idx])
        m = slits.slit_img(slitidx=slit_idx, initial=False) == slits.spat_id[slit_idx]
        slitid_img[m] = sid
        spat_ids.append(sid)

    # Per-slit cropped HDUs
    hdul = fits.HDUList([fits.PrimaryHDU()])
    hdul[0].header['SPECTROG'] = spec_name
    hdul[0].header['CALIBKEY'] = calib_key
    hdul[0].header['NSPEC']    = nspec
    hdul[0].header['NSPAT']    = nspat
    hdul[0].header['NSLITS']   = len(spat_ids)
    hdul[0].header['COMMENT']  = ('NGPS per-slice rawflat export, dumped before '
                                  'PypeIt flat fitting.')

    hdul.append(fits.ImageHDU(rawflat.astype(np.float32),
                              name='RAWFLAT'))
    hdul.append(fits.ImageHDU(waveimg.astype(np.float32),
                              name='WAVEIMG'))
    hdul.append(fits.ImageHDU(slitid_img, name='SLITID_IMG'))
    hdul.append(fits.ImageHDU(np.asarray(spat_ids, dtype=np.int32),
                              name='SPAT_ID'))

    for sid in spat_ids:
        on_slit = slitid_img == sid
        if not on_slit.any():
            continue
        # Slit bbox: cols range, full row range (slit spans whole spec).
        cols = np.where(on_slit.any(axis=0))[0]
        col_lo, col_hi = int(cols.min()), int(cols.max())
        sub_raw  = rawflat[:, col_lo:col_hi+1].astype(np.float32)
        sub_wave = waveimg[:, col_lo:col_hi+1].astype(np.float32)
        sub_mask = on_slit[:, col_lo:col_hi+1].astype(np.uint8)
        # Replace off-slit pixels in the cropped raw with NaN so
        # downstream code can tell which pixels are inside this slit.
        sub_raw_masked = np.where(sub_mask, sub_raw, np.nan)
        sub_wave_masked = np.where(sub_mask, sub_wave, 0.0)

        hdr = fits.Header()
        hdr['SPAT_ID'] = sid
        hdr['COL_LO']  = col_lo
        hdr['COL_HI']  = col_hi
        hdr['NSPEC']   = nspec
        hdr['NSPAT_S'] = col_hi - col_lo + 1
        hdr['COMMENT'] = ('Per-slice rawflat: off-slit pixels are NaN '
                          'in DATA and 0 in WAVE.')
        hdul.append(fits.ImageHDU(sub_raw_masked,
                                  name=f'SLIT_{sid:04d}_DATA',
                                  header=hdr))
        hdul.append(fits.ImageHDU(sub_wave_masked,
                                  name=f'SLIT_{sid:04d}_WAVE',
                                  header=hdr))
        hdul.append(fits.ImageHDU(sub_mask,
                                  name=f'SLIT_{sid:04d}_MASK',
                                  header=hdr))

    hdul.writeto(out, overwrite=True)
    log.info(f'NGPS {spec_name}: per-slice rawflat exported to {out}')


def _ngps_clean_flat_model(self, spec_name):
    """Zero out PypeIt's PIXELFLAT_MODEL for (a) off-slit pixels and
    (b) on-slit pixels outside the channel's
    ``[pixelflat_min_wave, pixelflat_max_wave]`` band.

    Doesn't affect science correction (which uses NORM x runtime
    illum_flat, not MODEL directly).  Just makes
    ``flat_diagnostic_<date>.png`` show clean panels in the bad-wave
    regions instead of bspline-fit chaos.
    """
    bounds = _NGPS_ILLUMFLAT_WAVE_BOUNDS.get(spec_name)
    if bounds is None:
        return
    wave_lo, wave_hi = bounds
    waveimg = self.waveimg
    slits   = self.slits
    if waveimg is None or slits is None:
        return
    # Build on-slit mask (union of all slits).
    on_slit = np.zeros_like(self.flat_model, dtype=bool)
    for slit_idx in range(slits.nslits):
        if slits.mask[slit_idx] != 0:
            continue
        slitimg = slits.slit_img(slitidx=slit_idx, initial=False)
        on_slit |= (slitimg == slits.spat_id[slit_idx])
    # Off-slit -> 0
    self.flat_model[~on_slit] = 0.0
    # On-slit but outside wave bounds -> 0
    bad_blue = on_slit & (waveimg > 0) & (waveimg < wave_lo)
    self.flat_model[bad_blue] = 0.0
    if wave_hi is not None:
        bad_red = on_slit & (waveimg > 0) & (waveimg > wave_hi)
        self.flat_model[bad_red] = 0.0


#: Per-channel wavelength bounds outside which the illum_flat is
#: forced to 1.0 (no flat correction applied at science time).
#: Matches the existing ``pixelflat_min_wave`` / ``pixelflat_max_wave``
#: clip on PIXELFLAT_NORM.  G no longer needs this after disabling
#: overscan -- PypeIt's spline now fits the rawflat correctly at the
#: blue end.  U still benefits because its lamp truly dies below
#: ~3400 A so the spline has nothing real to fit there.
_NGPS_ILLUMFLAT_WAVE_BOUNDS = {
    'p200_ngps_u': (3400.0, None),
}


def _install_ngps_rawflat_exporter():
    """Monkey-patch ``FlatField.run`` to dump the rawflat + slit info
    per slice BEFORE the normal flat fitting runs (used for the
    user's MATLAB / external flat-design experiments).  PypeIt's
    normal flat fitting runs unchanged.  Idempotent.

    Also installs an override on ``FlatImages.fit2illumflat`` that
    returns 1.0 for pixels whose wavelength is outside the
    ``_NGPS_ILLUMFLAT_WAVE_BOUNDS`` band for the channel.  This
    matches PypeIt's existing ``pixelflat_min_wave`` clip on
    PIXELFLAT_NORM so the COMBINED ``total_flat`` is 1.0 in those
    regions (no flat correction applied), instead of being divided
    by an unreliable bspline-extrapolated illum_flat.
    """
    from pypeit.flatfield import FlatField, FlatImages
    if getattr(FlatField.run, '_ngps_exporter_installed', False):
        return
    _orig_run = FlatField.run

    def _patched_run(self, *args, **kwargs):
        spec_name = getattr(self.spectrograph, 'name', None)
        is_ngps = spec_name in _NGPS_EXPORT_CHANNELS
        if (is_ngps
                and self.rawflatimg is not None
                and self.slits is not None
                and self.waveimg is not None
                and not self.spat_illum_only):
            try:
                _ngps_export_rawflat_per_slice(self)
            except Exception as exc:
                log.warning(f'NGPS rawflat export raised '
                            f'{type(exc).__name__}: {exc}; continuing.')
        result = _orig_run(self, *args, **kwargs)
        # Post-fit clean-up of the saved Flat HDUs so the diagnostic
        # plots (panel 3 = MODEL/row_med, panel 4 = RAW/MODEL) don't
        # show bad-pixel chaos in the bad-wave regions.  Sets MODEL=0
        # off-slit and outside [pixelflat_min_wave, pixelflat_max_wave];
        # those pixels are then masked out by the diagnostic's
        # (model>0) and (|raw/model - 1| < 0.5) filters.  Science
        # correction is unaffected (it uses NORM x runtime_illumflat,
        # not MODEL directly).
        if (is_ngps and result is not None
                and self.flat_model is not None
                and self.slits is not None
                and self.waveimg is not None
                and not self.spat_illum_only):
            try:
                _ngps_clean_flat_model(self, spec_name)
            except Exception as exc:
                log.warning(f'NGPS flat-model cleanup raised '
                            f'{type(exc).__name__}: {exc}; continuing.')
            try:
                if getattr(result, 'pixelflat_model', None) is not None:
                    result.pixelflat_model = self.flat_model
            except Exception:
                pass
        return result

    _patched_run._ngps_exporter_installed = True
    FlatField.run = _patched_run

    # fit2illumflat override: return 1.0 outside the wavelength bounds
    # so the bad-blue / bad-red regions get no flat correction (rather
    # than being multiplied by PypeIt's wonky bspline-extrapolated
    # illumflat).
    if not getattr(FlatImages.fit2illumflat, '_ngps_waveclip_patched', False):
        _orig_fit2 = FlatImages.fit2illumflat
        def _patched_fit2(self, slits, *args, **kwargs):
            img = _orig_fit2(self, slits, *args, **kwargs)
            spec_name = getattr(self, 'PYP_SPEC', None)
            bounds = _NGPS_ILLUMFLAT_WAVE_BOUNDS.get(spec_name)
            if bounds is None:
                return img
            wave_lo, wave_hi = bounds
            waveimg = getattr(self, 'pixelflat_waveimg', None)
            if waveimg is None or waveimg.shape != img.shape:
                return img
            out = img.copy()
            bad = (waveimg > 0) & (waveimg < wave_lo)
            if wave_hi is not None:
                bad |= (waveimg > 0) & (waveimg > wave_hi)
            out[bad] = 1.0
            return out
        _patched_fit2._ngps_waveclip_patched = True
        FlatImages.fit2illumflat = _patched_fit2
