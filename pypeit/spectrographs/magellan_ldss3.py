"""
Module for Magellan/LDSS3 specific methods.

LDSS3-C reads a single CCD through two amplifiers and writes **each amplifier
to its own file**, e.g. ``ccd0042c1.fits`` and ``ccd0042c2.fits``.  PypeIt joins
the two halves internally in :func:`MagellanLDSS3Spectrograph.get_rawimage`, so
no external merging step is required.  Only the amplifier-1 files are ingested
during setup (see :func:`MagellanLDSS3Spectrograph.find_raw_files`) and the
amplifier-2 file is located and read automatically.

.. include:: ../include/links.rst
"""
import re
from pathlib import Path

import numpy as np
from astropy.time import Time

from pypeit import log
from pypeit import PypeItError
from pypeit import telescopes
from pypeit import io
from pypeit.core import framematch
from pypeit.core import parse
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container


class MagellanLDSS3Spectrograph(spectrograph.Spectrograph):
    """
    Class to handle Magellan/LDSS3 data. Works for both MultiSlit and LongSlit modes.
    LDSS3 Technical documentation: https://www.lco.cl/technical-documentation/index-2/
    

    LDS3-C specs: https://www.lco.cl/?epkb_post_type_1=ldss3_c-2
    
    Currently supported grisms: VPH-All, VPH-Blue and VPH-Red.  Other grisms are not supported.
    
    Only 1x1 binning has been tested.
    
    """
    ndet = 1
    name = 'magellan_ldss3'
    telescope = telescopes.MagellanTelescopePar()
    camera = 'LDSS3-C'
    header_name = 'LDSS3-C'
    supported = True
    comment = 'Low Dispersion Survey Spectrograph 3-C; VPH-All, VPH-Blue and VPH-Red'
    url = 'https://www.lco.cl/technical-documentation/ldss-3-user-manual/'

    # Nominal amplifier properties, used when the raw headers are unavailable.
    # Both are superseded by the measured tables below whenever the readout mode
    # is one we have characterised.
    # This is readout mode-dependent
    nominal_gain = np.array([1.65, 1.44])
    nominal_ronoise = np.array([6.4, 6.7])

    # Read noise for LDSS3-C, in e-, for amplifiers 1 and 2 in each readout mode.
    #
    # The ENOISE header card is not reliable. All frames carry
    # the Slow-mode value regardless of the mode actually used.
    #
    # Slow and Turbo are published values and are not verified here.  Fast was
    # measured: 6.37 +- 0.07 and 6.65 +- 0.07 e- from the scatter of bias-pair
    # differences, over 34 pairs at five epochs from 2018 to 2025.  
    # Other modes need to be measured if they are ever used,
    # but the published values are probably good enough for most purposes.
    #
    # The first bias of a sequence should be discarded. Its pedestal is offset and it carries
    # structure the rest do not. 
    #

    ronoise_by_speed = {'slow': (4.4, 5.3),    # from LCO documentation
                        'fast': (6.4, 6.7),    # measured
                        'turbo': (10.0, 10.0)} # from LCO documentation

    # Gain for LDSS3-C, in e-/ADU, for amplifiers 1 and 2 in each readout mode.
    #
    # Slow and Turbo are the values published by LCO and have not been verified recently.
    #
    # Fast is measured, and differs from the published 1.5/1.8.
    # Photon transfer on flat pairs, over eight data sets spanning 2018 to 2025
    # and three grisms, gives 1.648 +- 0.006 and 1.439 +- 0.007; EGAIN reports
    # 1.65/1.47.  An abnormal discountinuity seen at the amplifier boundary is due to a gain mismatch.
    # 
    # Measuring these values accurately is needed to remove any step in the
    # image at the amplifier boundary.  
    # On the VPH-All longslit flat,
    # which has a slit spanning the boundary, the step is
    #
    #     g2/g1 = 0.873 (measured)   ->   +0.8 %
    #     g2/g1 = 0.891 (EGAIN)      ->   +2.8 %
    #     g2/g1 = 1.200 (published)  ->  +32.3 %
    #
    # The value measured here is used as the default, but the EGAIN card is used as a fallback
    # A similar characterization is needed for other readout modes
    # 
    gain_by_speed = {'slow': (0.16, 0.19), # from LCO documentation
                     'fast': (1.65, 1.44), # measured 
                     'turbo': (2.7, 3.1)}  # from LCO documentation

    # Matches an LDSS3 raw filename of the form ``<stem>c<amp>.fits[.gz]``.  The
    # amplifier index is anchored immediately before the extension so that a
    # stem which happens to contain "c1" is not mistaken for the amplifier tag.
    _amp_file_regex = re.compile(r'(?P<stem>.+?)c(?P<amp>\d+)(?P<ext>\.fits(?:\.gz)?)\Z',
                                 re.IGNORECASE)

    @classmethod
    def parse_amp_file(cls, filename):
        """
        Split an LDSS3 raw filename into its exposure stem, amplifier number and
        extension.

        Args:
            filename (:obj:`str`, `Path`_):
                Filename to parse.  Only the base name is inspected.

        Returns:
            :obj:`tuple`: The exposure stem (:obj:`str`), the 1-indexed
            amplifier number (:obj:`int`), and the file extension (:obj:`str`).
            Returns None if the name does not follow the amplifier-split
            convention, which is the case for files that have already been
            merged outside of PypeIt.
        """
        match = cls._amp_file_regex.fullmatch(Path(filename).name)
        if match is None:
            return None
        return match.group('stem'), int(match.group('amp')), match.group('ext')

    @classmethod
    def amp_files(cls, raw_file):
        """
        Find every per-amplifier file belonging to one exposure.

        Args:
            raw_file (:obj:`str`, `Path`_):
                Any one of the amplifier files for the exposure.

        Returns:
            :obj:`list`: `Path`_ objects for the amplifier files, ordered by the
            amplifier index encoded in the filename.  If ``raw_file`` does not
            follow the amplifier-split naming convention it is assumed to be an
            already-merged image and is returned on its own.
        """
        path = Path(raw_file)
        parsed = cls.parse_amp_file(path)
        if parsed is None:
            # Not an amplifier-split name; treat it as a single, merged image
            log.warning(f'File {raw_file} does not follow the LDSS3 amplifier-split naming '
                        'convention; assuming it is an already-merged image.')
            return [path]

        stem, _, _ = parsed
        found = []
        for candidate in path.parent.glob(f'{stem}c*'):
            cand_parsed = cls.parse_amp_file(candidate)
            # Require an exact stem match so that, e.g., "ccd004" does not pick
            # up the amplifier files of "ccd0042"
            if cand_parsed is not None and cand_parsed[0] == stem:
                found.append((cand_parsed[1], candidate))
        if len(found) == 0:
            raise PypeItError(f'No LDSS3 amplifier files found for {raw_file}')
        return [f[1] for f in sorted(found)]

    @classmethod
    def find_raw_files(cls, root, extension=None):
        """
        Find raw observations for this spectrograph in the provided directory.

        Overrides the base class to return only **one file per exposure**.
        LDSS3-C writes each amplifier to its own file; the companion amplifier
        files are found and read automatically by :func:`get_rawimage`.  Without
        this filter, every exposure would be ingested twice and reduced twice.

        Args:
            root (:obj:`str`, `Path`_, :obj:`list`):
                One or more paths to search for files.  See
                :func:`~pypeit.spectrographs.spectrograph.Spectrograph.find_raw_files`.
            extension (:obj:`str`, :obj:`list`, optional):
                One or more file extensions to search on.

        Returns:
            :obj:`list`: `Path`_ objects for the unique raw exposures found.
        """
        files = super().find_raw_files(root, extension=extension)

        # Group the amplifier-split files by exposure and keep the lowest
        # amplifier index in each group.  Grouping (rather than simply keeping
        # amplifier 1) means a directory holding only amplifier-2 files still
        # yields one file per exposure.
        exposures = {}
        keep = []
        for path in files:
            parsed = cls.parse_amp_file(path)
            if parsed is None:
                # Already-merged image; pass it through untouched
                keep.append(path)
                continue
            key = (str(Path(path).parent), parsed[0])
            if key not in exposures or parsed[1] < exposures[key][0]:
                exposures[key] = (parsed[1], path)

        keep += [v[1] for v in exposures.values()]
        ndropped = len(files) - len(keep)
        if ndropped > 0:
            log.info(f'Ignoring duplicate amplifier file(s) for {ndropped} exposure(s). Only c1 '
                     'files are read from the .pypeit and c2 are read automatically.')
        return sorted(keep)

    def get_detector_par(self, det, hdu=None, amp_headers=None):
        """
        Return metadata for the selected detector.

        LDSS3-C has a single CCD read through two amplifiers.

        Args:
            det (:obj:`int`):
                1-indexed detector number.
            hdu (`astropy.io.fits.HDUList`_, optional):
                The open fits file with the raw image of interest.  If not
                provided, frame-dependent parameters are set to a default.
            amp_headers (:obj:`list`, optional):
                `astropy.io.fits.Header`_ objects, one per amplifier actually
                found for this exposure, ordered by amplifier.  When provided,
                the number of amplifiers, the gain and the read noise are taken
                from these headers.  When not provided, nominal two-amplifier
                values are used.

        Returns:
            :class:`~pypeit.images.detector_container.DetectorContainer`:
            Object with the detector metadata.
        """
        # Binning
        binning = '1,1' if hdu is None else self.get_meta_value(self.get_headarr(hdu), 'binning')

        if amp_headers is None:
            gain = self.nominal_gain.copy()
            ronoise = self.nominal_ronoise.copy()
        else:
            # Take the gain from the published table, keyed on the readout mode.
            gain = []
            for i, h in enumerate(amp_headers):
                iamp = int(h.get('OPAMP', i + 1)) - 1
                speed = str(h.get('SPEED', '')).strip().lower()
                table = self.gain_by_speed.get(speed)
                if table is not None and 0 <= iamp < len(table):
                    gain.append(float(table[iamp]))
                else:
                    fallback = float(h.get('EGAIN', self.nominal_gain[min(i, 1)]))
                    log.warning(f'Gain for LDSS3 readout mode {h.get("SPEED")!r} has not '
                                f'been characterised; falling back on the EGAIN card '
                                f'({fallback} e-/ADU). ')
                    gain.append(fallback)
            gain = np.array(gain)
            # Take the read noise from the published table, keyed on the
            # readout mode
            ronoise = []
            for i, h in enumerate(amp_headers):
                iamp = int(h.get('OPAMP', i + 1)) - 1
                speed = str(h.get('SPEED', '')).strip().lower()
                table = self.ronoise_by_speed.get(speed)
                if table is not None and 0 <= iamp < len(table):
                    ronoise.append(float(table[iamp]))
                else:
                    fallback = float(h.get('ENOISE', self.nominal_ronoise[min(i, 1)]))
                    log.warning(f'Unrecognised LDSS3 readout mode {h.get("SPEED")!r}; '
                                f'falling back on the ENOISE card ({fallback} e-), '
                                'which is known to report the Slow-mode value '
                                'regardless of the mode actually used.')
                    ronoise.append(fallback)
            ronoise = np.array(ronoise)

        detector_dict = dict(
                            binning         = binning,
                            det             = 1,
                            dataext         = 0,
                            specaxis        = 0,
                            specflip        = False,
                            spatflip        = False,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.189,    # from the LDSS3-C manual
                            darkcurr        = 25.0,     # from the LDSS3-C manual, in e-/hr/pix
                            # LDSS3-C has a full well of ~205000 e- (10%
                            # non-linear) and stays within 1% of linear below
                            # ~175000 e-. The ADC is
                            # 16-bit, so max ADU is 65535
                            #
                            # The threshold is a single scalar compared in
                            # electrons, while the two amplifiers have different
                            # gains, so it has to be low enough to catch the
                            # amplifier with the *smaller* gain.
                            saturation      = 65535.,
                            nonlinear       = 0.99,
                            mincounts       = -1e10,
                            numamplifiers   = len(gain),
                            gain            = np.atleast_1d(gain),
                            ronoise         = np.atleast_1d(ronoise),
                            )

        # Instantiate
        return detector_container.DetectorContainer(**detector_dict)

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the ``PypeIt``-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=0, card='RA')
        self.meta['dec'] = dict(ext=0, card='DEC')
        self.meta['target'] = dict(ext=0, card='OBJECT')
        self.meta['decker'] = dict(ext=0, card='APERTURE')
        self.meta['binning'] = dict(ext=0, card=None, compound=True)
        # mjd is not in the header, but can be constructed from UT-DATE and UT-TIME
        self.meta['mjd'] = dict(ext=0, card=None, compound=True)
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')

        self.meta['dispname'] = dict(ext=0, card='GRISM')
        self.meta['filter1'] = dict(ext=0, card='FILTER')
        self.meta['idname'] = dict(ext=0, card='EXPTYPE')
        # Distinguishes the two amplifier files of a single exposure
        self.meta['amp'] = dict(ext=0, card='OPAMP')

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
        if meta_key == 'mjd':
            # NOTE: The JD header card is not usable; build the MJD from the UT
            # date and time instead.
            ttime = Time('{:s}T{:s}'.format(headarr[0]['UT-DATE'], headarr[0]['UT-TIME']),
                         format='isot')
            return ttime.mjd
        if meta_key == 'binning':
            # The BINNING card is written in the raw detector frame, which is
            # spatial along NAXIS1 and spectral along NAXIS2, i.e. the opposite
            # of the PypeIt (spectral,spatial) convention.
            binspatial, binspec = parse.parse_binning(headarr[0]['BINNING'])
            return parse.binning2string(binspec, binspatial)
        return None

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.

        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of ``PypeIt`` methods.
        """
        par = super().default_pypeit_par()

        # Wavelengths
        # 1D wavelength solution
        par['calibrations']['wavelengths']['rms_thresh_frac_fwhm'] = 0.5
        par['calibrations']['wavelengths']['sigdetect'] = 6.
        par['calibrations']['wavelengths']['fwhm'] = 5.0
        par['calibrations']['wavelengths']['match_toler'] = 2.5
        par['calibrations']['wavelengths']['n_first'] = 3
        par['calibrations']['wavelengths']['n_final'] = 5
        par['calibrations']['wavelengths']['method'] = 'holy-grail'

        # Tilt and slit parameters
        par['calibrations']['tilts']['tracethresh'] = 20.0
        par['calibrations']['tilts']['spat_order'] = 6
        par['calibrations']['tilts']['spec_order'] = 6

        # edges.  Narrow detector defects otherwise get synced into spurious
        # few-pixel slits; the narrowest real LDSS3 slitlets are ~5 arcsec.
        par['calibrations']['slitedges']['edge_thresh'] = 20.
        par['calibrations']['slitedges']['minimum_slit_length'] = 2.

        # Processing steps
        turn_off = dict(use_biasimage=False, use_darkimage=False)
        par.reset_all_processimages_par(**turn_off)

        # Combine the (usually few) bias and dark frames with a median
        par['calibrations']['biasframe']['process']['combine'] = 'median'
        par['calibrations']['darkframe']['process']['combine'] = 'median'

        # Mask saturated slits; the alignment boxes of a slitmask are routinely
        # saturated in the flats
        par['calibrations']['flatfield']['saturated_slits'] = 'mask'

        # Extraction
        par['reduce']['skysub']['bspline_spacing'] = 0.8
        par['reduce']['extraction']['sn_gauss'] = 4.0
        # Do not perform global sky subtraction for standard stars
        par['reduce']['skysub']['global_sky_std'] = False

        # Flexure.  The default wavelength solution comes from arc lamps taken
        # at a different telescope pointing, so correct the residual spectral
        # flexure against the sky.  Set this to 'skip' when calibrating on sky
        # lines instead (see the magellan_ldss3 documentation).
        par['flexure']['spec_method'] = 'boxcar'

        # cosmic ray rejection parameters for science frames
        par['scienceframe']['process']['sigclip'] = 5.0
        par['scienceframe']['process']['objlim'] = 2.0

        # Set the default exposure time ranges for the frame typing.  The
        # science/standard split is driven entirely by these ranges.
        # The types have to be checked in the .pypeit file anyways
        par['calibrations']['standardframe']['exprng'] = [None, 100]
        par['calibrations']['arcframe']['exprng'] = [None, 60]
        par['calibrations']['tiltframe']['exprng'] = [None, 60]
        par['calibrations']['darkframe']['exprng'] = [1, None]
        par['scienceframe']['exprng'] = [100, None]

        # Sensitivity function parameters
        par['sensfunc']['algorithm'] = 'IR'
        par['sensfunc']['polyorder'] = 7
        par['sensfunc']['IR']['telgridfile'] = 'TellPCA_3000_26000_R15000.fits'

        return par

    def config_specific_par(self, inp, inp_par=None):
        """
        Modify the ``PypeIt`` parameters to hard-wired values used for
        specific instrument configurations.

        Args:
            inp (:obj:`str`, `Path`_, `astropy.io.fits.Header`_, `astropy.table.Table`_):
                File or metadata used to determine the configuration and how to
                adjust the input parameters.
            inp_par (:class:`~pypeit.par.parset.ParSet`, optional):
                Parameter set used for the full run of PypeIt.  If None,
                use :func:`default_pypeit_par`.

        Returns:
            :class:`~pypeit.par.parset.ParSet`: The PypeIt parameter set
            adjusted for configuration specific parameter values.
        """
        # Start with instrument-wide parameters
        par = super().config_specific_par(inp, inp_par=inp_par)

        # Adjust parameters based on settings used
        decker = self.get_meta_value(inp, 'decker')
        grating = self.get_meta_value(inp, 'dispname')
        binning = self.get_meta_value(inp, 'binning')

        # Turn the slit-edge PCA off for the longslit deckers.  ``decker`` is
        # None when this is called on a reduced spec1d/spec2d file.
        if decker is not None and ('center' in decker.lower() or 'long' in decker.lower()):
            par['calibrations']['slitedges']['sync_predict'] = 'nearest'

        # Arc-lamp wavelength solutions.  Sky-line calibration is available as a
        # documented alternative; see doc/spectrographs/magellan_ldss3.rst.
        if grating == 'VPH-Blue':
            par['calibrations']['wavelengths']['lamps'] = ['HeI', 'NeI', 'ArI']
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] \
                    = 'magellan_ldss3_vph_blue_HeINeIArI.fits'
        elif grating == 'VPH-Red':
            par['calibrations']['wavelengths']['lamps'] = ['HeI', 'NeI', 'ArI']
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] \
                    = 'magellan_ldss3_vph_red_HeINeIArI.fits'
        elif grating == 'VPH-All':
            par['calibrations']['wavelengths']['lamps'] = ['HeI', 'NeI', 'ArI']
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] \
                    = 'magellan_ldss3_vph_all_HeINeIArI.fits'
            # NOTE: The bluest arc line detectable with VPH-All is near 3890 A,
            # around pixel 740.  The wavelength solution is therefore
            # extrapolated over the bluest ~740 pixels, where the grism has no
            # throughput and the science frames carry no signal either.
        elif grating is not None:
            log.warning(f'No archived LDSS3 wavelength solution for grism {grating}; '
                        'falling back on holy-grail.')
            par['calibrations']['wavelengths']['lamps'] = ['HeI', 'NeI', 'ArI']
            par['calibrations']['wavelengths']['method'] = 'holy-grail'

        # FWHM, in binned pixels
        if binning is not None:
            par['calibrations']['wavelengths']['fwhm'] = 6.0 / parse.parse_binning(binning)[0]

        return par

    def configuration_keys(self):
        """
        Return the metadata keys that define a unique instrument
        configuration.

        Returns:
            :obj:`list`: List of keywords of data pulled from file headers
            and used to construct the
            :class:`~pypeit.metadata.PypeItMetaData` object.
        """
        # filter1 matters because the order-blocking filter (e.g. BPF-380-590)
        # sets the usable wavelength range for a given grism.
        # observers are recommended to use the same filter for all exposures in a given configuration, 
        # but it is not enforced and this can be overwritten in the .pypeit file if necessary
        return ['dispname', 'decker', 'binning', 'filter1']

    def valid_configuration_values(self):
        """
        Return a fixed set of valid values for metadata used to define a
        unique instrument configuration.

        Frames taken with the grism wheel open are acquisition, through-slit
        or direct images.  PypeIt cannot reduce them, so they are dropped
        during setup rather than being carried through as untyped rows.
        
        Only the three VPH grisms are supported, so the valid values for ``dispname``
        are hard-coded here. 

        Returns:
            :obj:`dict`: Keys are the metadata keywords and values are the
            lists of valid values for that keyword.
        """
        return {'dispname': ['VPH-All', 'VPH-Blue', 'VPH-Red']}

    def config_independent_frames(self):
        """
        Define frame types that are independent of the fully defined
        instrument configuration.

        Bias and dark frames are taken with the grism wheel left wherever it
        happened to be, so their ``dispname`` is meaningless and they must not
        be matched on it.

        Returns:
            :obj:`dict`: Dictionary where the keys are the frame types that
            are configuration independent and the values are the metadata
            keywords that can be used to assign the frames to a configuration
            group.
        """
        return {'bias': 'binning', 'dark': 'binning'}

    def pypeit_file_keys(self):
        """
        Define the list of keys to be output into a standard ``PypeIt`` file.

        Returns:
            :obj:`list`: The list of keywords in the relevant
            :class:`~pypeit.metadata.PypeItMetaData` instance to print to the
            :ref:`pypeit_file`.
        """
        return super().pypeit_file_keys() + ['amp']

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.

        LDSS3's ``EXPTYPE`` card alone is not sufficient to identify frame types.
        The ``OBJECT`` card is the exposure name, which is set by the observer and may not reflect the actual exposure type.
        
        This method implements a set of rules to determine the frame type but it is not guaranteed to be correct.  
        The user should check the frame typing in the .pypeit file and adjust it if necessary based
        on an observing log.
        See :ref:`magellan_ldss3` for the naming conventions this recognises.

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
        # Only ever type the amplifier-1 file of an exposure.  find_raw_files()
        # already filters these out during setup
        primary = np.array([str(a).strip() in ('1', 'None', '') for a in fitstbl['amp']])

        idname = np.array([str(n).strip().lower() for n in fitstbl['idname']])
        target = np.array([str(n).strip().lower() for n in fitstbl['target']])
        grism = np.array([str(n).strip().lower() for n in fitstbl['dispname']])

        # Frames taken with the grism wheel open are acquisition, through-slit
        # or direct images; they are never spectroscopic calibrations or data.
        dispersed = grism != 'open'
        # Explicitly non-spectroscopic object names
        setup_frame = self._has_keyword(target, words=self.setup_words)
        is_arc = self._has_keyword(target, words=self.arc_words)
        is_flat = self._has_keyword(target, words=self.flat_words,
                                    substrings=self.flat_substrings)

        if ftype == 'bias':
            return good_exp & primary & ((idname == 'bias')
                                         | self._has_keyword(target, words=['zero', 'bias']))
        if ftype == 'dark':
            return good_exp & primary & ((idname == 'dark')
                                         | self._has_keyword(target, words=['dark']))
        if ftype in ['arc', 'tilt']:
            # Checked before the flats: an arc frame is sometimes written with
            # EXPTYPE='Flat' and an object name that also names the flat
            return good_exp & primary & dispersed & is_arc
        if ftype in ['pixelflat', 'illumflat', 'trace']:
            return good_exp & primary & dispersed & is_flat & np.logical_not(is_arc)
        if ftype in ['science', 'standard']:
            # The science/standard split is set by the exposure-time ranges in
            # the 'scienceframe' and 'standardframe' parameters.
            is_cal = is_arc | is_flat | self._has_keyword(target,
                                                          words=['zero', 'bias', 'dark'])
            return good_exp & primary & dispersed & np.logical_not(is_cal) \
                    & np.logical_not(setup_frame)

        log.warning('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

    # Object-name keywords used for frame typing.  ``*_words`` are matched
    # against whitespace-, underscore- or hyphen-delimited words, so that short
    # or ambiguous tokens cannot match inside a target designation.
    # ``*_substrings`` are matched anywhere in the name,
    # for distinctive tokens that observers run together with others, such as
    # the 'flat' in 'flatQh'.
    # this is intended to help frame typing, but it is not guaranteed to be correct.  
    # The user should check the frame typing in the .pypeit file and adjust it if necessary based on an observing log.
    arc_words = ['arc', 'arcs', 'henear', 'hene', 'lamp', 'comp']
    flat_words = ['qh', 'ff', 'flatfield']
    flat_substrings = ['flat', 'quartz', 'dome']
    setup_words = ['align', 'thr', 'field', 'focus', 'acq']

    @staticmethod
    def _has_keyword(names, words=None, substrings=None):
        """
        Test each of ``names`` for any of the given keywords.

        Args:
            names (`numpy.ndarray`_):
                Lower-cased strings to test.
            words (:obj:`list`, optional):
                Lower-case words to look for as whitespace-, underscore- or
                hyphen-delimited words.
            substrings (:obj:`list`, optional):
                Lower-case strings to look for anywhere in the name.

        Returns:
            `numpy.ndarray`_: Boolean array, one element per entry in ``names``.
        """
        _words = set() if words is None else set(words)
        _subs = [] if substrings is None else substrings
        return np.array([len(_words & set(re.split(r'[\s_\-]+', n))) > 0
                         or any(s in n for s in _subs) for n in names])

    def get_rawimage(self, raw_file, det):
        """
        Read raw images and generate a few other bits and pieces
        that are key for image processing.

        The two LDSS3 amplifiers live in separate files.  Both are located and
        read here, and joined into a single image laid out as::

            [ overscan_1 | data_1 | flip(data_2) | flip(overscan_2) ]

        along the spatial axis.  If the companion amplifier file is missing, the
        image is built from the single amplifier available and the detector is
        reconfigured to have one amplifier, so that the gain and read noise are
        not applied to a region that holds no data.

        Args:
            raw_file (:obj:`str`, `Path`_):
                File to read.  Any one of the exposure's amplifier files.
            det (:obj:`int`):
                1-indexed detector to read.

        Returns:
            :obj:`tuple`: The detector metadata
            (:class:`~pypeit.images.detector_container.DetectorContainer`), the
            raw image (`numpy.ndarray`_), the opened file
            (`astropy.io.fits.HDUList`_), the exposure time (:obj:`float`), and
            the data-section and overscan-section amplifier images
            (`numpy.ndarray`_).
        """
        self._check_extensions(raw_file)
        amp_files = self.amp_files(raw_file)

        # Read each amplifier and sort by the OPAMP header card, which is
        # authoritative; the filename index is only a fallback.
        amps = []
        for ifile in amp_files:
            data, overscan, header, nx, nxb = ldss3_read_amp(ifile)
            parsed = self.parse_amp_file(ifile)
            opamp = header.get('OPAMP', parsed[1] if parsed is not None else len(amps) + 1)
            amps.append((int(opamp), data, overscan, header, nx, nxb))
        amps.sort(key=lambda a: a[0])

        namp = len(amps)
        if namp == 1:
            log.warning(f'Only one amplifier file found for {Path(raw_file).name}; '
                        'proceeding with a single amplifier.')
        elif namp != 2:
            raise PypeItError(f'Expected 1 or 2 LDSS3 amplifier files, found {namp} for '
                              f'{raw_file}.')

        # Warn about binning we have not been able to validate
        binning = self.get_meta_value([amps[0][3]], 'binning')
        if binning != '1,1':
            log.warning(f'LDSS3 support has only been verified for 1x1 binning; found '
                        f'{binning}. Check the data and overscan sections carefully.')

        # The image is opened again here so that the returned HDUList is the one
        # PypeIt associates with the file it was given.
        hdu = io.fits_open(amp_files[0])

        detector_par = self.get_detector_par(det if det is not None else 1, hdu=hdu,
                                             amp_headers=[a[3] for a in amps])

        # Assemble.  Work in the raw detector frame, which is spatial along the
        # first axis and spectral along the second, then transpose at the end to
        # give the PypeIt (nspec,nspat) convention.
        nspec = amps[0][1].shape[1]
        nspat = sum(a[4] for a in amps)
        array = np.zeros((nspat, nspec), dtype=float)
        rawdatasec_img = np.zeros_like(array, dtype=int)
        oscansec_img = np.zeros_like(array, dtype=int)

        # Amplifier 1: overscan at the low-spatial edge, then the data
        _, data1, overscan1, _, _, nxb1 = amps[0]
        ndata1 = data1.shape[0]
        array[:nxb1, :] = overscan1
        oscansec_img[:nxb1, :] = 1
        array[nxb1:nxb1+ndata1, :] = data1
        rawdatasec_img[nxb1:nxb1+ndata1, :] = 1

        # Amplifier 2: flipped about the spatial axis so that it abuts amplifier
        # 1 at the physical centre of the CCD, then its overscan at the high edge
        if namp > 1:
            _, data2, overscan2, _, _, nxb2 = amps[1]
            ndata2 = data2.shape[0]
            start = nxb1 + ndata1
            array[start:start+ndata2, :] = np.flipud(data2)
            rawdatasec_img[start:start+ndata2, :] = 2
            array[start+ndata2:start+ndata2+nxb2, :] = np.flipud(overscan2)
            oscansec_img[start+ndata2:start+ndata2+nxb2, :] = 2

        exptime = self.get_meta_value([amps[0][3]], 'exptime')

        return detector_par, array.T, hdu, exptime, rawdatasec_img.T, oscansec_img.T

    def bpm(self, filename, det, shape=None, msbias=None):
        """
        Generate a default bad-pixel mask.

        Args:
            filename (:obj:`str` or None):
                An example file to use to get the image shape.
            det (:obj:`int`):
                1-indexed detector number to use when getting the image
                shape from the example file.
            shape (tuple, optional):
                Processed image shape.  Required if filename is None.
                Ignored if filename is not None.
            msbias (`numpy.ndarray`_, optional):
                Processed bias frame used to identify bad pixels.

        Returns:
            `numpy.ndarray`_: An integer array with a masked value set
            to 1 and an unmasked value set to 0.
        """
        bpm_img = super().bpm(filename, det, shape=shape, msbias=msbias)

        binspat = parse.parse_binning(self.get_meta_value(filename, 'binning'))[1] \
                if filename is not None else 1

        # Bad columns, as (first, last) inclusive spatial pixel ranges of the
        # assembled and trimmed unbinned frame.  Measured from the column
        # median of combined bias frames, and confirmed on two epochs four years
        # apart (2018 and 2022).  The first and last entries are the outer edges
        # of amplifiers 1 and 2.
        for c1, c2 in [(0, 11), (443, 443), (608, 608), (1413, 1413), (1492, 1494),
                       (1549, 1551), (1602, 1606), (1635, 1640), (1688, 1690),
                       (1695, 1695), (1999, 2000), (2034, 2047)]:
            bpm_img[:, c1 // binspat:c2 // binspat + 1] = 1

        return bpm_img


def ldss3_read_amp(fil):
    """
    Read a single amplifier file of LDSS3 data.

    Args:
        fil (:obj:`str`, `Path`_):
            Filename.

    Returns:
        :obj:`tuple`: The data section and overscan section of the amplifier,
        each a `numpy.ndarray`_ ordered as (spatial,spectral); the primary
        `astropy.io.fits.Header`_; and the total and overscan widths along the
        spatial axis (:obj:`int`).
    """
    log.info(f'Reading LDSS3 amplifier file: {fil}')
    hdu = io.fits_open(fil)
    head = hdu[0].header

    # NOTE: DATASEC and BIASSEC are written in binned pixels -- they sum to
    # NAXIS1 -- so no binning correction is applied here.
    x1, x2, _, _ = np.array(parse.load_sections(head['DATASEC'], fmt_iraf=False)).flatten()
    b1, b2, _, _ = np.array(parse.load_sections(head['BIASSEC'], fmt_iraf=False)).flatten()
    _, _, y1, y2 = np.array(parse.load_sections(head['DATASEC'], fmt_iraf=False)).flatten()

    nxb = b2 - b1 + 1
    nx = (x2 - x1 + 1) + nxb
    ny = y2 - y1 + 1

    # Transpose to (spatial,spectral) and drop the unused overscan rows
    array = hdu[0].data.T[:, :ny] * 1.0
    data = array[:nx-nxb, :]
    overscan = array[nx-nxb:nx, :]

    hdu.close()
    return data, overscan, head, nx, nxb
