"""
Module for Magellan/LDSS3 specific methods.

LDSS3-C reads its single CCD through two amplifiers and writes **each amplifier
to its own file**, e.g. ``ccd0042c1.fits`` and ``ccd0042c2.fits``.  PypeIt joins
the two halves internally in :func:`MagellanLDSS3Spectrograph.get_rawimage`, so
no external merging step is required.  Only the amplifier-1 files are ingested
during setup (see :func:`MagellanLDSS3Spectrograph.find_raw_files`); the
amplifier-2 file is located and read automatically.

.. include:: ../include/links.rst
"""
import re
from pathlib import Path

import numpy as np
from astropy.time import Time
from astropy.io import fits

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
    Child to handle Magellan/LDSS3 specific code
    """
    ndet = 1
    name = 'magellan_ldss3'
    telescope = telescopes.MagellanTelescopePar()
    camera = 'LDSS3-C'
    header_name = 'LDSS3-C'
    supported = True
    comment = 'Low Dispersion Survey Spectrograph 3-C; VPH-All, VPH-Blue and VPH-Red'
    url = 'https://www.lco.cl/technical-documentation/ldss-3-user-manual/'

    # Nominal amplifier properties, used when the raw headers are unavailable
    # (e.g., when building the documentation).  The true values are read from
    # the EGAIN/ENOISE cards of each amplifier file; they depend on the readout
    # speed (the SPEED header card).
    nominal_gain = np.array([1.65, 1.47])
    nominal_ronoise = np.array([4.67, 5.06])

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
            log.info(f'Ignoring {ndropped} companion amplifier file(s); they are read '
                     'automatically with their amplifier-1 counterpart.')
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
            # EGAIN/ENOISE track the readout speed (the SPEED card), so prefer
            # the header values over the nominal ones.
            gain = np.array([float(h.get('EGAIN', self.nominal_gain[i]))
                             for i, h in enumerate(amp_headers)])
            ronoise = np.array([float(h.get('ENOISE', self.nominal_ronoise[i]))
                                for i, h in enumerate(amp_headers)])

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
                            platescale      = 0.189,
                            darkcurr        = 25.0,
                            saturation      = 205000.,
                            nonlinear       = 0.85,
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

        # edges
        par['calibrations']['slitedges']['edge_thresh'] = 20.

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
            # The archived VPH-All solution extends well blueward of the real
            # throughput of the grism; restrict it when it is read.
            par['calibrations']['wavelengths']['wvrng_arxiv'] = [4000., 10500.]
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
        return ['dispname', 'decker', 'binning', 'filter1']

    def valid_configuration_values(self):
        """
        Return a fixed set of valid values for metadata used to define a
        unique instrument configuration.

        Frames taken with the grism wheel open are acquisition, through-slit
        or direct images.  PypeIt cannot reduce them, so they are dropped
        during setup rather than being carried through as untyped rows.

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

        LDSS3's ``EXPTYPE`` card alone is not sufficient: arc and flat frames
        both appear as ``EXPTYPE='Flat'`` *and* ``EXPTYPE='Object'`` depending
        on how the observation was taken.  The ``OBJECT`` card carries the
        observer's intent, so the two are combined.  See
        :ref:`magellan_ldss3` for the naming conventions this recognises.

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
        # already filters these out during setup; this guards hand-written
        # PypeIt files, where listing both halves would reduce every exposure
        # twice.
        primary = np.array([str(a).strip() in ('1', 'None', '') for a in fitstbl['amp']])

        idname = np.array([str(n).strip().lower() for n in fitstbl['idname']])
        target = np.array([str(n).strip().lower() for n in fitstbl['target']])
        grism = np.array([str(n).strip().lower() for n in fitstbl['dispname']])

        # Frames taken with the grism wheel open are acquisition, through-slit
        # or direct images; they are never spectroscopic calibrations or data.
        dispersed = grism != 'open'
        # Explicitly non-spectroscopic object names
        setup_frame = self._has_keyword(target, ['align', 'thr', 'field', 'focus', 'acq'])

        if ftype == 'bias':
            return good_exp & primary & ((idname == 'bias')
                                         | self._has_keyword(target, ['zero', 'bias']))
        if ftype == 'dark':
            return good_exp & primary & ((idname == 'dark')
                                         | self._has_keyword(target, ['dark']))
        if ftype in ['arc', 'tilt']:
            # Checked before the flats: an arc frame is sometimes written with
            # EXPTYPE='Flat'
            return good_exp & primary & dispersed \
                    & self._has_keyword(target, ['arc', 'henear', 'hene', 'lamp', 'comp'])
        if ftype in ['pixelflat', 'illumflat', 'trace']:
            is_flat = self._has_keyword(target, ['flat', 'qh', 'quartz', 'dome'])
            is_arc = self._has_keyword(target, ['arc', 'henear', 'hene', 'lamp', 'comp'])
            return good_exp & primary & dispersed & is_flat & np.logical_not(is_arc)
        if ftype in ['science', 'standard']:
            # The science/standard split is set by the exposure-time ranges in
            # the 'scienceframe' and 'standardframe' parameters.
            is_cal = self._has_keyword(target, ['arc', 'henear', 'hene', 'lamp', 'comp',
                                                'flat', 'qh', 'quartz', 'dome',
                                                'zero', 'bias', 'dark'])
            return good_exp & primary & dispersed & np.logical_not(is_cal) \
                    & np.logical_not(setup_frame)

        log.warning('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

    @staticmethod
    def _has_keyword(names, keywords):
        """
        Test whether any of ``keywords`` appears as a whitespace-delimited word
        in each of ``names``.

        Matching on words rather than substrings avoids false positives such as
        the ``arc`` in ``search`` or the ``qh`` in a target designation.

        Args:
            names (`numpy.ndarray`_):
                Lower-cased strings to test.
            keywords (:obj:`list`):
                Lower-case words to look for.

        Returns:
            `numpy.ndarray`_: Boolean array, one element per entry in ``names``.
        """
        _keywords = set(keywords)
        return np.array([len(_keywords & set(re.split(r'[\s_\-]+', n))) > 0 for n in names])

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

        # Bad columns on amplifier 1, in spatial pixels of the assembled and
        # trimmed frame.  These are for unbinned data.
        # TODO: The amplifier-2 bad columns still need to be re-derived from a
        # flat field.  Amplifier 2 is flipped during assembly, so the column
        # indices measured on the raw c2 frame do not carry over directly.
        binspec, binspat = parse.parse_binning(self.get_meta_value(filename, 'binning')) \
                if filename is not None else (1, 1)
        for col in [443]:
            bpm_img[:, col // binspat] = 1

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
