"""
Module for MMT MMIRS

.. include:: ../include/links.rst
"""
import re
from pathlib import Path

import numpy as np
from scipy.signal import savgol_filter

from astropy.table import Table
from astropy.time import Time
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.coordinates import Angle, SkyCoord
from astropy import units

from pypeit import log
from pypeit import PypeItError
from pypeit import telescopes
from pypeit import utils
from pypeit import io
from pypeit.core import parse
from pypeit.core import framematch
from pypeit.core import ramp
from pypeit.ext.fitramp import fitramp
from pypeit.images import detector_container
from pypeit.spectrographs import spectrograph
from pypeit.spectrographs.ramp_spectrograph import RampSpectrograph
from pypeit.spectrographs.slitmask import SlitMask
from pypeit.par import parset


class MMTMMIRSSpectrograph(RampSpectrograph, spectrograph.Spectrograph):
    """
    Child to handle MMT/MMIRS specific code.

    MMIRS is read out up-the-ramp, so it mixes in
    :class:`~pypeit.spectrographs.ramp_spectrograph.RampSpectrograph` (which
    provides the ramp-fitting machinery) and implements the two instrument
    hooks :func:`_load_ramp` and :func:`_count_reads`.
    """
    ndet = 1
    name = 'mmt_mmirs'
    telescope = telescopes.MMTTelescopePar()
    camera = 'MMIRS'
    url = 'https://lweb.cfa.harvard.edu/mmti/mmirs.html'
    header_name = 'mmirs'
    supported = True

    # Up-the-ramp fitting: MMIRS-calibrated single-read-noise values, overriding
    # the RampSpectrograph defaults.  The generic ramp machinery (orchestration,
    # dark handling, thread/chunk defaults, and the ``[rdx] rampfit_dir``
    # subdirectory) comes from the RampSpectrograph mixin; MMIRS implements only
    # the instrument hooks :func:`_load_ramp` and :func:`_count_reads`.  See
    # :mod:`pypeit.core.ramp` and
    # :class:`~pypeit.spectrographs.ramp_spectrograph.RampSpectrograph`.
    ramp_sig_guess = 9.0
    """Initial guess for the *effective* single-read noise in electrons (the
    instantaneous read noise plus accumulated dark-current/flux shot noise),
    used to seed the chi-square rescaling and as the fallback when a frame has
    too few reads to calibrate.  Refined from darks or self-calibrated per
    frame; the header ``RDNOISE`` is used as a physical floor (see
    :func:`~pypeit.spectrographs.ramp_spectrograph.RampSpectrograph.get_ramp_sigma`)."""
    ramp_sig_range = (3.0, 50.0)
    """Absolute sanity range for the calibrated single-read noise (electrons).
    The lower bound is raised to the header ``RDNOISE`` when available, since a
    derived noise below the instantaneous read noise is unphysical."""
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
            # pointing minus the catalog target (CAT-RA/CAT-DEC) onto the slit
            # PA (POSANGLE).
            hdr = headarr[1]
            # The RA card is always in hours (decimal, e.g. '13.70246500', in
            # 2017-2019 data; sexagesimal, e.g. '+02:59:16.62', in newer data)
            # and DEC always in degrees; Angle parses both forms.  A missing or
            # unparseable card means this is not an on-sky science frame, so
            # there is no nod offset.
            try:
                ra = Angle(str(hdr['RA']), unit=units.hourangle).deg
                dec = Angle(str(hdr['DEC']), unit=units.deg).deg
                catdec = Angle(str(hdr['CAT-DEC']), unit=units.deg).deg
                pa = float(hdr['POSANGLE'])
            except (KeyError, TypeError, ValueError):
                return 0.0
            # CAT-RA units are ambiguous: sexagesimal degrees ('+260:36:50.55',
            # old MOS mask tool), decimal hours ('13.70246500', old longslit),
            # or sexagesimal hours ('+02:59:16.80', newer data).  Parse it both
            # as degrees and as hours and keep whichever lands closest to the
            # actual pointing, so RA is never mis-scaled by 15x.
            cat_cands = []
            for unit in (units.deg, units.hourangle):
                try:
                    cat_cands.append(Angle(str(hdr['CAT-RA']), unit=unit).deg)
                except (KeyError, ValueError, TypeError):
                    pass
            if len(cat_cands) == 0:
                return 0.0
            catra = min(cat_cands,
                        key=lambda c: abs((c - ra + 180.0) % 360.0 - 180.0))
            # Off-sky calibrations (darks/flats) carry sentinel coordinates
            # (e.g. DEC = -100) that SkyCoord rejects as an invalid latitude;
            # caught here and returned as a zero offset -- they are never nodded.
            try:
                target = SkyCoord(catra * units.deg, catdec * units.deg)
                pointing = SkyCoord(ra * units.deg, dec * units.deg)
            except (ValueError, TypeError):
                return 0.0
            # Signed on-sky offset of the pointing from the catalog target;
            # spherical_offsets_to applies cos(dec) and wraps the RA difference,
            # so a sequence straddling 0h RA yields a small offset.
            dra, ddec = target.spherical_offsets_to(pointing)
            return dra.arcsec * np.sin(np.radians(pa)) \
                + ddec.arcsec * np.cos(np.radians(pa))
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

        # Name collated coadd outputs after the slitmask-design object name
        # (MASKDEF_OBJNAME) when a mask is used, so per-target coadds carry the
        # real catalog target names rather than sky coordinates.
        par['collate1d']['outfile_from'] = 'maskdef_objname'

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
        elif (grating=='H3000') and (dichroic=='H'):
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = 'mmt_mmirs_H3000_H.fits'
        elif (grating=='HK') and (dichroic=='HK3'):
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = 'mmt_mmirs_HK_HK3.fits'

        return par

    def config_specific_setup_lines(self, subtbl, paths):
        """
        Bake mask-design parameters into the PypeIt file when a ``<decker>.msk``
        sidecar sits next to the raw data.

        MMIRS stores its slitmask design in a separate xfitmask ``.msk`` file,
        unlike DEIMOS/MOSFIRE/Binospec which embed the design in the raw FITS
        frames.  The sidecar is named ``<decker>.msk`` (``decker`` is header
        ``APERTURE`` == ``MOSID``) and lives in the raw-data directory; because
        the MMIRS ``decker`` label carries no longslit/MOS naming convention,
        the *presence* of the file is also the signal that the observation is a
        MOS mask rather than a longslit.  When found, this enables mask-design
        slit tracing and object assignment by writing the parameters into the
        generated PypeIt file (:func:`config_specific_par` does not probe the
        filesystem).  If no ``.msk`` is found, no lines are added and the
        reduction falls back to generic slit tracing and the SLIT-id keyed
        coadd workflow.

        Parameters
        ----------
        subtbl : `astropy.table.Table`_
            The metadata rows for this setup; ``decker`` is constant within a
            setup (it is a configuration key).
        paths : :obj:`list`
            The unique raw-data directories for this setup.

        Returns
        -------
        :obj:`list`
            Configuration lines enabling mask design, or an empty list.
        """
        if 'decker' not in subtbl.colnames or len(subtbl) == 0:
            return []
        label = subtbl['decker'][0]
        if label is None:
            return []
        for directory in paths:
            maskfile = (Path(directory) / f'{label}.msk').absolute()
            if maskfile.exists():
                return [
                    '[calibrations]',
                    '  [[slitedges]]',
                    '    use_maskdesign = True',
                    f'    maskdesign_filename = {maskfile}',
                    '[reduce]',
                    '  [[slitmask]]',
                    '    assign_obj = True',
                    '    extract_missing_objs = True',
                    '    use_alignbox = True',
                ]
        return []

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
        header, slits = read_mmirs_maskfile(filename)
        return self._build_slitmask(header, slits)

    def _build_slitmask(self, header, slits):
        """
        Build :attr:`slitmask` from an already-parsed ``.msk`` design.

        Factored out of :func:`get_slitmask` so that
        :func:`get_maskdef_slitedges`, which also needs the raw ``header`` and
        ``slits``, can reuse a single parse of the file instead of reading it
        twice.

        Parameters
        ----------
        header : :obj:`dict`
            The mask-level header, as returned by :func:`read_mmirs_maskfile`.
        slits : `astropy.table.Table`_
            The per-slit design table, as returned by
            :func:`read_mmirs_maskfile`.

        Returns
        -------
        :class:`~pypeit.spectrographs.slitmask.SlitMask`
            The slitmask, also stored in :attr:`slitmask`.
        """
        arcsec_per_mm = 1.0 / header['arc2mm']
        n = len(slits)

        # SlitMask corners: shape (N, 4, 2), (x=spatial, y=spectral).
        # Spatial axis is -y_mm, spectral axis is x_mm (see spec geometry).
        # Build a rectangle centered on each slit from height (spatial) and
        # width (spectral), in arcsec, then rotate it by the design tilt theta.
        half_len = 0.5 * np.asarray(slits['height_mm'], dtype=float) * arcsec_per_mm
        half_wid = 0.5 * np.asarray(slits['width_mm'], dtype=float) * arcsec_per_mm
        cx = -np.asarray(slits['y_mm'], dtype=float) * arcsec_per_mm   # spatial center
        cy = np.asarray(slits['x_mm'], dtype=float) * arcsec_per_mm    # spectral center

        # Base (unrotated) corner offsets from the slit center, ordered
        # top-right, bottom-right, bottom-left, top-left (the input order
        # SlitMask expects): axis 0 is spatial (slit length, +/- half_len),
        # axis 1 is spectral (slit width, +/- half_wid).
        sign_len = np.array([+1., -1., -1., +1.])
        sign_wid = np.array([-1., -1., +1., +1.])
        doff_len = sign_len[None, :] * half_len[:, None]   # (n, 4) spatial
        doff_wid = sign_wid[None, :] * half_wid[:, None]   # (n, 4) spectral

        # Rotate each slit's rectangle by its design tilt theta (deg) about the
        # slit center, in the detector (spatial, spectral) plane.  theta = 0
        # (the usual MOS case) leaves the rectangle axis-aligned, reproducing
        # the previous behavior exactly.
        theta = np.radians(np.asarray(slits['theta_deg'], dtype=float))
        ct = np.cos(theta)[:, None]
        st = np.sin(theta)[:, None]
        rot_len = doff_len * ct - doff_wid * st
        rot_wid = doff_len * st + doff_wid * ct

        corners = np.zeros((n, 4, 2), dtype=float)
        corners[:, :, 0] = cx[:, None] + rot_len
        corners[:, :, 1] = cy[:, None] + rot_wid

        is_target = np.array([t == 'TARGET' for t in slits['type']])
        # top/bottom object distances from the slit edges (arcsec).  The design
        # `offset` (mm) displaces the target from the slit center along the slit
        # long axis; positive offset is toward +mask-y, which maps to smaller
        # spatial pixel (spatial = -y_mm), i.e. toward the left/top edge, so it
        # decreases OBJ_TOPDIST (distance from the left edge) and increases
        # OBJ_BOTDIST.  offset = 0 (the usual case) leaves the target centered
        # (top == bot == half_len).  (Sign convention follows the same
        # -y -> spatial mapping validated for the slit centers; it should be
        # re-checked against a real off-center mask when one is available.)
        off_as = np.asarray(slits['offset_mm'], dtype=float) * arcsec_per_mm
        top = np.clip(half_len - off_as, 0.0, 2.0 * half_len)
        bot = np.clip(half_len + off_as, 0.0, 2.0 * half_len)
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
                                 align=np.logical_not(is_target),
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

        header, slits = read_mmirs_maskfile(_fname)
        self._build_slitmask(header, slits)

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
        treated as a stare and left unpaired.  Works identically for longslit
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
        # init_meta seeds dithpat/dithpos (registered in the core meta data
        # model), so _build always creates them; they arrive as fixed-width
        # unicode columns (from their 'None' default) that would silently
        # truncate longer labels (e.g. "ABA'B'" -> "ABA'").  Widen to object
        # dtype before writing.
        for col in ['dithpat', 'dithpos']:
            fitstbl[col] = np.asarray(fitstbl[col], dtype=object)

        # Partition frames so nod pairing never crosses instrument
        # configurations (setup) *or* distinct targets.  A single longslit
        # setup can hold several science targets and standard stars that share
        # the same disperser/slit (hence the same setup); greedily pairing
        # across them could background-subtract an unrelated object.  MOS masks
        # are already separated by the decker (part of the setup); within a
        # setup we additionally split on the target name, which is constant
        # across a nod sequence.
        # set_combination_groups always runs after set_configurations, so the
        # 'setup' column must exist here; its absence means the metadata table
        # was built out of order.  ('target' is a plain header card, so treat a
        # missing target column as a single unnamed target.)
        if 'setup' not in fitstbl.colnames:
            raise PypeItError("'setup' column missing from the metadata table; "
                              'get_comb_group must run after set_configurations.')
        setup_col = np.asarray(fitstbl['setup'])
        target_col = np.asarray(fitstbl['target']) if 'target' in fitstbl.colnames \
            else np.full(len(fitstbl), '')
        part_key = np.array(['{0}\x00{1}'.format(s, t)
                             for s, t in zip(setup_col, target_col)])
        for gk in np.unique(part_key[sci_std]):
            idx = np.where(sci_std & (part_key == gk))[0]
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
            # appearance in time.  `side` (dithoff > midpoint) is the A/B split.
            dithpos = np.array(['None'] * idx.size, dtype=object)
            for base, on_side in [('A', side), ('B', np.logical_not(side))]:
                # local indices on this nod side, in time (mjd) order
                grp = order[on_side[order]]
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
        combined using up-the-ramp fitting (see
        :func:`~pypeit.spectrographs.ramp_spectrograph.RampSpectrograph._ramp_fit_image`);
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

        rampfit_file = self.rampfit_path(fil)

        if hdu[0].header.get('RAMPFIT') is None \
                and self._count_reads(hdu) >= self.ramp_min_reads:
            # Multi-read cube: swap in a fresh preprocessed 2D image if one
            # exists in the reduction directory
            if ramp.rampfit_fresh(rampfit_file, fil):
                log.info(f'Loading preprocessed ramp image: {rampfit_file}')
                hdu.close()
                hdu = io.fits_open(rampfit_file)

        head1 = hdu[1].header

        detector_par = self.get_detector_par(det if det is not None else 1, hdu=hdu)

        # get the x and y binning factors...
        binning = head1['CCDSUM']
        xbin, ybin = [int(ibin) for ibin in binning.split(' ')]

        # Need the exposure time
        exptime = self.get_meta_value(hdu, 'exptime')
        gain = detector_par['gain'][0]

        if hdu[0].header.get('RAMPFIT') is not None:
            # Preprocessed 2D count-rate image (e-/s): convert to ADU
            array = hdu[1].data.astype(np.float64) * exptime / gain
            detector_par['ronoise'] = np.atleast_1d(hdu[0].header['RAMPRON'])
        elif self._count_reads(hdu) >= self.ramp_min_reads:
            # Up-the-ramp fitting with jump detection
            rate, sig, eff_ronoise = self._ramp_fit_image(hdu, detector_par)
            detector_par['ronoise'] = np.atleast_1d(eff_ronoise)
            array = rate * exptime / gain
            # Persist the fit so later loads (and other scripts) reuse it.
            # The RampFit directory lives in the reduction directory, which
            # must be writable for the rest of the reduction anyway, so a
            # write failure is left to propagate like any other output.
            ramp.write_rampfit(rampfit_file, rate, hdu, sig, eff_ronoise,
                               self._count_reads(hdu), Path(fil).stat().st_mtime,
                               raw_file=fil)
            log.info(f'Wrote preprocessed ramp image: {rampfit_file}')
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

    def _load_ramp(self, hdu, detector_par):
        """
        Load the non-destructive reads of a raw MMIRS up-the-ramp cube.

        Instrument hook backing the base-class ramp orchestration (see
        :func:`~pypeit.spectrographs.ramp_spectrograph.RampSpectrograph._load_ramp`).
        The reference-pixel-corrected reads (:func:`mmirs_load_ramp`) are
        scaled to electrons by the detector gain, and the covariance is built
        from the per-read time (``GRPTIME``).  The fit itself uses the
        ``fitramp`` algorithm of Brandt (2024,
        https://arxiv.org/abs/2404.01326; reference implementation:
        https://github.com/t-brandt/fitramp) and was inspired by the prototype
        at https://github.com/zhechenghu/mmt-mmirs-up-the-ramp-pypeit.

        Args:
            hdu (`astropy.io.fits.HDUList`_):
                Opened raw MMIRS cube.
            detector_par (:class:`~pypeit.images.detector_container.DetectorContainer`):
                Detector parameters; provides the gain.

        Returns:
            :obj:`tuple`: The reads in electrons (`numpy.ndarray`_, shape
            ``(ngroups, ny, nx)``), the matching
            :class:`~pypeit.ext.fitramp.fitramp.Covar`, and the extension-1
            `astropy.io.fits.Header`_.
        """
        reads, head1 = mmirs_load_ramp(hdu)
        reads *= detector_par['gain'][0]      # ADU -> electrons
        ngroups = reads.shape[0]
        covar = fitramp.Covar([head1['GRPTIME'] * (i + 1)
                               for i in range(ngroups)])
        return reads, covar, head1

    def _count_reads(self, hdu):
        """
        Count the non-destructive reads in a raw MMIRS cube.

        Instrument hook; see
        :func:`~pypeit.spectrographs.ramp_spectrograph.RampSpectrograph._count_reads`.
        """
        return mmirs_count_reads(hdu)


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




def _sexed_ra_hours(s):
    """
    Convert a sexagesimal RA string in *hours* to degrees.

    Parameters
    ----------
    s : :obj:`str`
        Sexagesimal right ascension in hours (e.g. ``'17:22:20.2490'``).

    Returns
    -------
    :obj:`float`
        Right ascension in degrees.
    """
    return Angle(s, unit=units.hourangle).to('deg').value


def _sexed_dec_deg(s):
    """
    Convert a sexagesimal Dec string in *degrees* to degrees.

    Parameters
    ----------
    s : :obj:`str`
        Sexagesimal declination in degrees (e.g. ``'65:56:13.040'``).

    Returns
    -------
    :obj:`float`
        Declination in degrees.
    """
    return Angle(s, unit=units.deg).value


def read_mmirs_maskfile(path):
    """
    Parse an MMT/MMIRS xfitmask ``.msk`` mask-design file.

    The ``.msk`` file (xfitmask output) is plain text with mixed tab/space
    separation and three sections: a mask-level key/value header, a
    ``GuideStars`` block (ignored), and a slit table.  The slit-table columns
    are ``slit ra dec x y target object type height width offset theta bbox
    polygon`` (field indices 0..13), where ``bbox`` (field 12) and ``polygon``
    (field 13) are space-separated numbers within a single tab-delimited field.

    Parameters
    ----------
    path : :obj:`str` or `Path`_
        Path to the ``.msk`` file.

    Returns
    -------
    header : :obj:`dict`
        Mask-level metadata: ``label``, ``ra_deg``, ``dec_deg``, ``pa``,
        ``scale``, ``arc2mm``, ``corners`` (list of 4 floats, mm), and, when
        present in the file, ``grism`` and ``filter``.
    slits : `astropy.table.Table`_
        One row per slit with columns ``slit`` (int), ``ra_deg`` (float),
        ``dec_deg`` (float), ``x_mm`` (float), ``y_mm`` (float), ``target``
        (int), ``object`` (str), ``type`` (str), ``height_mm`` (float),
        ``width_mm`` (float), ``offset_mm`` (float, target displacement from
        the slit center along the slit long axis), ``theta_deg`` (float, slit
        tilt), ``bbox`` (float, shape ``(N, 4)``), ``polygon`` (float, shape
        ``(N, 8)``).

    Raises
    ------
    PypeItError
        If the file has no slit-table header row, or no slit rows can be
        parsed.
    """
    lines = Path(path).read_text().splitlines()

    header = {}
    slit_hdr_idx = None
    for i, line in enumerate(lines):
        fields = re.split(r'\t+', line.rstrip())
        if not fields or fields[0] == '':
            continue
        key = fields[0].strip()
        if key == 'label' and len(fields) > 1:
            header['label'] = fields[1].strip()
        elif key == 'ra' and len(fields) > 1:
            header['ra_deg'] = _sexed_ra_hours(fields[1].strip())
        elif key == 'dec' and len(fields) > 1:
            header['dec_deg'] = _sexed_dec_deg(fields[1].strip())
        elif key in ('pa', 'scale', 'arc2mm') and len(fields) > 1:
            header[key] = float(fields[1].strip())
        elif key in ('grism', 'filter') and len(fields) > 1:
            header[key] = fields[1].strip()
        elif key == 'corners' and len(fields) >= 5:
            header['corners'] = [float(v) for v in fields[1:5]]
        elif key == 'slit' and slit_hdr_idx is None:
            slit_hdr_idx = i

    if slit_hdr_idx is None:
        raise PypeItError(f'No slit table found in mask file {path}')

    rows = []
    for line in lines[slit_hdr_idx + 1:]:
        fields = re.split(r'\t+', line.rstrip())
        if len(fields) < 14 or not re.match(r'^\s*\d+\s*$', fields[0]):
            continue
        try:
            rows.append(dict(
                slit=int(fields[0]),
                ra_deg=_sexed_ra_hours(fields[1].strip()),
                dec_deg=_sexed_dec_deg(fields[2].strip()),
                x_mm=float(fields[3]),
                y_mm=float(fields[4]),
                target=int(fields[5]),
                object=fields[6].strip(),
                type=fields[7].strip(),
                height_mm=float(fields[8]),
                width_mm=float(fields[9]),
                offset_mm=float(fields[10]),
                theta_deg=float(fields[11]),
                bbox=[float(v) for v in fields[12].split()],
                polygon=[float(v) for v in fields[13].split()],
            ))
        except (ValueError, IndexError):
            continue

    if len(rows) == 0:
        raise PypeItError(f'No slit rows parsed from mask file {path}')

    slits = Table()
    slits['slit'] = np.array([r['slit'] for r in rows], dtype=int)
    slits['ra_deg'] = np.array([r['ra_deg'] for r in rows], dtype=float)
    slits['dec_deg'] = np.array([r['dec_deg'] for r in rows], dtype=float)
    slits['x_mm'] = np.array([r['x_mm'] for r in rows], dtype=float)
    slits['y_mm'] = np.array([r['y_mm'] for r in rows], dtype=float)
    slits['target'] = np.array([r['target'] for r in rows], dtype=int)
    slits['object'] = np.array([r['object'] for r in rows], dtype=object)
    slits['type'] = np.array([r['type'] for r in rows], dtype=object)
    slits['height_mm'] = np.array([r['height_mm'] for r in rows], dtype=float)
    slits['width_mm'] = np.array([r['width_mm'] for r in rows], dtype=float)
    slits['offset_mm'] = np.array([r['offset_mm'] for r in rows], dtype=float)
    slits['theta_deg'] = np.array([r['theta_deg'] for r in rows], dtype=float)
    slits['bbox'] = np.array([r['bbox'] for r in rows], dtype=float)
    slits['polygon'] = np.array([r['polygon'] for r in rows], dtype=float)
    return header, slits
