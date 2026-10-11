"""
Module for Magellan/LDSS3 specific methods.

.. include:: ../include/links.rst
"""
import re
import copy
from pathlib import Path
from dataclasses import dataclass

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


@dataclass(frozen=True)
class LDSS3AmplifierProperties:
    """Gain and read noise of one LDSS3 amplifier in one readout mode."""
    speed: str          # SPEED header value: 'Slow', 'Fast' or 'Turbo'
    readout_gain: str   # GAIN header value: 'Low' or 'High'
    opamp: int          # OPAMP header value (1 = c1 file, 2 = c2 file)
    gain: float         # e-/ADU
    ronoise: float      # e-
    det: int = 1         # PypeIt detector number (LDSS3-C has a single detector)
    source: str = ''    # where the numbers come from
    # e-; the response departs from linear by >1% above this level.
    # LDSS3 User Manual, Section 6.1.
    full_well_non_linearity: float = 175000.
    saturation_level: float = 65535.  # ADU, 16-bit ADC

    @property
    def nonlinear_relative_to_saturation(self):
        """
        Return the non-linearity fraction relative to the ADC saturation level for this amplifier.
        This ensures that pixels above the non-linearity threshold are masked even if they are
        below the saturation level.

        Not used for now, but could be used to set the nonlinear parameter in the detector container
        if the non-linearity threshold is characterized in the future.
        """
        return self.full_well_non_linearity / self.gain / self.saturation_level

    def matches(self, speed, readout_gain):
        """
        Check whether this entry applies to a readout mode.

        Args:
            speed (:obj:`str`):
                Value of the ``SPEED`` header card.
            readout_gain (:obj:`str`):
                Value of the ``GAIN`` header card.

        Returns:
            :obj:`bool`: True if both values match this entry, ignoring case
            and surrounding whitespace.
        """
        return (self.speed.lower() == str(speed).strip().lower()
                and self.readout_gain.lower() == str(readout_gain).strip().lower())


class MagellanLDSS3Spectrograph(spectrograph.Spectrograph):
    """
    Class to handle Magellan/LDSS3 data. Supports both multislit and longslit modes.

    LDSS3 has one detector read by two amplifiers, written to separate files.
    Both amplifiers can be reduced together, or either one alone:

    - Two-amplifier reduction: the c1 and c2 files of each exposure are in the
      same raw data directory.
    - Single-amplifier reduction: the raw data directory holds only c1 or only
      c2 files.  Keep the c1 and c2 files in separate directories, because a
      companion file found next to the listed file will always be joined.

    """
    ndet = 1
    name = 'magellan_ldss3'
    telescope = telescopes.MagellanTelescopePar()
    camera = 'LDSS3-C'
    header_name = 'LDSS3-C'
    supported = True
    comment = 'VPH-All, VPH-Blue and VPH-Red grisms; see :doc:`magellan_ldss3`'
    url = 'https://www.lco.cl/technical-documentation/ldss-3-user-manual'

    # Gain (e-/ADU) and read noise (e-) per readout mode, Low gain setting only.
    # Other settings fall back on the EGAIN and ENOISE header cards; see
    # get_detector_par().
    LDSS3_AMPLIFIERS_PROPERTIES = [
        # Slow/Low: LDSS3 User Manual, Section 6.1.  At these gains the 16-bit
        # ADC saturates near 1e4 e-, well below the full well.
        LDSS3AmplifierProperties(speed='Slow', readout_gain='Low', opamp=1, gain=0.15, ronoise=4.4),
        LDSS3AmplifierProperties(speed='Slow', readout_gain='Low', opamp=2, gain=0.13, ronoise=5.3),

        # Fast/Low: measured.  Gain from photon transfer on flat pairs (eight
        # data sets, 2018-2025, three grisms): 1.648 +/- 0.006 and 1.439 +/- 0.007.
        # Read noise from bias-pair differences: 6.37 +/- 0.07 and 6.65 +/- 0.07.
        # The User Manual gives 1.67/1.45 e-/ADU and 7.0/7.2 e-.  The gain ratio
        # sets the step at the amplifier boundary.
        LDSS3AmplifierProperties(speed='Fast',  readout_gain='Low', opamp=1, gain=1.65, ronoise=6.4,
                                 source='Measured'),
        LDSS3AmplifierProperties(speed='Fast',  readout_gain='Low', opamp=2, gain=1.44, ronoise=6.7,
                                 source='Measured'),

        # Turbo/Low: LDSS3 User Manual, Section 6.1.  The read noise is
        # approximate ("not measured lately") and Turbo is not meant for science.
        LDSS3AmplifierProperties(speed='Turbo', readout_gain='Low', opamp=1, gain=3.04, ronoise=10.0),
        LDSS3AmplifierProperties(speed='Turbo', readout_gain='Low', opamp=2, gain=2.63, ronoise=10.0)
        ]

    # Expected shape of the 2-amplifier combined image, in unbinned pixels
    n_spat_unbinned = 2048
    n_spec_unbinned = 4096

    # Matches an LDSS3 raw filename of the form ``ccd<number>c<amp>.fits[.gz]``.
    # For an example filename ``ccd0042c1.fits``, the stem is ``ccd0042``
    # (the ``ccd`` prefix is kept in the stem), the amplifier is ``1``, and the
    # extension is ``.fits``.
    _amp_file_regex = re.compile(r'(?P<stem>ccd\d+)c(?P<amp>\d+)(?P<ext>\.fits(?:\.gz)?)\Z',
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
            Returns None if the name does not follow the LDSS3 file naming convention.
        """
        match = cls._amp_file_regex.fullmatch(Path(filename).name)
        if match is None:
            return None
        return match.group('stem'), int(match.group('amp')), match.group('ext')

    @classmethod
    def amp_files(cls, raw_file):
        """
        Find both amplifier files associated with one exposure.

        The two amplifier files are assumed to be in the same directory
        and have the same extension, differing only in the amplifier index.

        Args:
            raw_file (:obj:`str`, `Path`_):
                Any one of the amplifier files for the exposure.

        Returns:
            :obj:`list`: `Path`_ objects for the amplifier files, sorted by amplifier number.
            If raw_file does not follow the LDSS3 file naming convention, or if the companion
            amplifier file is missing, the list will contain None for the missing amplifier(s).
        """

        path = Path(raw_file)
        parsed = cls.parse_amp_file(path)
        amps_files = [None, None]

        if parsed is None:
            # Filename was not matched
            log.info(f'{raw_file} does not follow the LDSS3 file naming convention '
                     '(ccd<number>c<amp>.fits[.gz]).')
            return amps_files

        stem, amp, file_ext = parsed

        opamp = fits.getheader(path, ext=0)['OPAMP']
        # Check that the OPAMP header card matches the amplifier number in the filename.
        if opamp != amp:
            # Skip this exposure if not.
            log.warning(f'File {raw_file} has OPAMP={opamp} but filename indicates {amp}; skipping.')
            return amps_files

        # Save the file. The list is 0-indexed, but the amplifier numbers are 1-indexed.
        amps_files[amp-1] = path

        # Companion amplifier (1-indexed, 1 or 2).
        other_amp = 3-amp

        # The companion amplifier file is expected to be in the same directory.
        other_amp_path = path.parent / f'{stem}c{other_amp}{file_ext}'

        # Check that the companion amplifier file exists.
        if not other_amp_path.exists():
            log.info(f'File {other_amp_path} not found for {raw_file}')
            return amps_files

        # Check that the OPAMP header card matches too.
        other_opamp = fits.getheader(other_amp_path, ext=0)['OPAMP']
        if other_opamp != other_amp:
            # Skip this exposure if not.
            log.warning(f'File {other_amp_path} has OPAMP={other_opamp} but the filename '
                        f'indicates {other_amp}; skipping.')
            return amps_files

        # Save the companion amplifier file
        amps_files[other_amp-1] = other_amp_path

        return amps_files

    @classmethod
    def find_raw_files(cls, root, extension=None):
        """
        Find raw observations for this spectrograph in the provided directory.

        Overrides the base class to return one file per exposure.  All LDSS3
        raw files must be in the same directory.  The amplifier files are
        associated by
        :func:`~pypeit.spectrographs.magellan_ldss3.MagellanLDSS3Spectrograph.amp_files`.

        Only one group of exposures is returned, chosen in this order:

        1. Exposures with both amplifier files; the ``c1`` file is returned.
        2. Otherwise, exposures with only a ``c1`` file.
        3. Otherwise, exposures with only a ``c2`` file.

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
        parents = {Path(f).parent for f in files}
        if len(parents) > 1:
            raise PypeItError(
                f'LDSS3 raw files were found in {len(parents)} directories: '
                f'{", ".join(sorted(str(p) for p in parents))}. All amplifier files for a '
                f'given exposure must be in the same directory; a symlink is enough.')

        # One entry per exposure.
        exposures = {}  # stem -> [c1, c2]
        for file in files:
            path = Path(file)
            parsed = cls.parse_amp_file(path)
            if parsed is None:
                continue
            stem, amp, file_ext = parsed
            if stem not in exposures:
                exposures[stem] = cls.amp_files(path)

        # (c1 present, c2 present), description, index of file to keep, detector mode
        _amp_priority = (
            ((True, True),  'with both amplifiers',  0, 'TWO-AMPLIFIER (whole detector)'),
            ((True, False), 'with amplifier 1 only', 0, 'SINGLE-AMPLIFIER (half detector)'),
            ((False, True), 'with amplifier 2 only', 1, 'SINGLE-AMPLIFIER (half detector)'),
        )

        # Reduce only the highest-priority group that exists, so every frame in a
        # setup can be associated with the same calibrations.
        for present, label, i, mode in _amp_priority:
            usable = sorted(v[i] for v in exposures.values()
                            if tuple(map(bool, v)) == present)
            if not usable:
                continue
            if nskip := len(exposures) - len(usable):
                log.warning(f'Found {len(usable)} exposure(s) {label}; the '
                            f'remaining {nskip} of {len(exposures)} will be skipped.')
            msg = f'Found {len(usable)} exposure(s) {label} as {mode} images.'
            if present == (True, True):
                log.info(msg)
            else:
                log.info(f'No complete c1/c2 pairs found.  {msg}')
            return usable

        log.warning('No usable LDSS3 exposures found.')
        return []

    def get_detector_par(self, det, hdu=None, amp_headers=None):
        """
        Return metadata for the selected detector.

        LDSS3-C has a single CCD read by two amplifiers.

        The gain and read noise depend on the readout mode,
        which is the same for both amplifiers.

        This function accepts either a single ``hdu`` or a list of amplifier
        headers (``amp_headers``). If both are provided, ``amp_headers`` is used.

        If two amplifier headers are provided, assumes a 2-amplifier reduction
        and returns the gain and read noise for each amplifier.

        If one amplifier header (or hdu) is provided, assumes a 1-amplifier reduction
        and returns the gain and read noise for that amplifier.

        Args:
            det (:obj:`int`):
                1-indexed detector number. Only det=1 is supported for LDSS3-C.
            hdu (`astropy.io.fits.HDUList`_, optional):
                The opened FITS file. Can be any of the amplifier files for the exposure.
            amp_headers (:obj:`list`, optional):
                `astropy.io.fits.Header`_ objects, one per amplifier.

        Returns:
            :class:`~pypeit.images.detector_container.DetectorContainer`:
            Object with the detector metadata.
        """
        if det != 1:
            # Only one detector is supported for LDSS3-C.
            raise PypeItError(f'LDSS3 has only one detector, but det={det} was requested.')

        if amp_headers is not None:
            if len(amp_headers) not in (1, 2):
                raise PypeItError('LDSS3 amplifier headers must be provided for 1 or 2 amplifiers, '
                                  f'but {len(amp_headers)} were given.')
            headers = list(amp_headers)
        elif hdu is not None:
            headers = [hdu[0].header]
        else:
            log.warning('No header or hdu provided; assuming Fast speed, Low gain, and 1x1 binning.')
            headers = []

        # The readout speed, gain and binning are the same for both amplifiers
        if headers:
            readout_speed, readout_gain = self.get_meta_value(headers, 'readmode').split(',')
            binning = self.get_meta_value(headers, 'binning')
        else:
            # Assume Fast speed and Low gain if a header is not available.
            readout_speed, readout_gain, binning = 'fast', 'low', '1,1'

        # The OPAMP header card is 1-indexed
        opamps = sorted(int(h.get('OPAMP', i + 1)) for i, h in enumerate(headers)) \
                    if headers else [1, 2]
        # Amplifier properties.
        amps = sorted((a for a in self.LDSS3_AMPLIFIERS_PROPERTIES
                       if a.matches(readout_speed, readout_gain) and a.opamp in opamps),
                      key=lambda a: a.opamp)
        # In case no amplifier properties are found.
        amps = amps if amps else None

        # The non-linearity threshold, defined relative to the saturation level.
        # NOTE: The LDSS3 non-linearity is expressed in terms
        # of e- in 6.1, and in ADUs, in 6.2.
        # The e- values imply non-linearity is not reached below the saturation level
        # for the Slow and Fast readout modes,
        # but 6.2 says that linearity is better than 1% below 40,000 ADUs
        # (0.61 of the ADC saturation).
        # TODO: This value will be adopted to avoid masking too many pixels,
        # but a good characterization of the non-linearity threshold is needed.
        nonlinear = 0.95

        if amps is None:
            log.warning(f'LDSS3 readout mode {readout_speed}/{readout_gain} is not characterized; '
                        'taking the gain and read noise from the EGAIN and ENOISE header keys.')
            gain = [h['EGAIN'] for h in headers]
            ronoise = [h['ENOISE'] for h in headers]
            # The same for every amplifier and readout mode.
            saturation = self.LDSS3_AMPLIFIERS_PROPERTIES[0].saturation_level
        elif len(amps) == len(opamps):
            # All found amplifiers are characterized
            gain = [amp.gain for amp in amps]
            ronoise = [amp.ronoise for amp in amps]
            saturation = amps[0].saturation_level
        else:
            raise PypeItError(f'LDSS3 readout mode {readout_speed}/{readout_gain} is tabulated for '
                              f'{len(amps)} amplifier(s); '
                              f'expected {len(opamps)}, with OPAMP {", ".join(str(o) for o in opamps)}.')

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
                            platescale      = 0.189,    # From the LDSS3-C User Manual
                            darkcurr        = 25.0,     # From the LDSS3-C User Manual, in e-/hr/pix
                            saturation      = saturation, # 16-bit ADC
                            nonlinear       = nonlinear,
                            mincounts       = -1e10,
                            numamplifiers   = len(gain), # 2 if two-amplifier reduction, 1 if single-amplifier reduction
                            gain            = np.atleast_1d(gain),
                            ronoise         = np.atleast_1d(ronoise),
                            )

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
        self.meta['instrument'] = dict(ext=0, card='INSTRUME')
        self.meta['binning'] = dict(ext=0, card=None, compound=True)
        # mjd is not in the header, but can be constructed from UT-DATE and UT-TIME
        self.meta['mjd'] = dict(ext=0, card=None, compound=True)
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')

        self.meta['dispname'] = dict(ext=0, card='GRISM')
        self.meta['filter1'] = dict(ext=0, card='FILTER')
        self.meta['idname'] = dict(ext=0, card='EXPTYPE')
        self.meta['readmode'] = dict(ext=0, card=None, compound=True)

    def compound_meta(self, headarr, meta_key):
        """
        Methods to generate metadata requiring interpretation of the header
        data, instead of simply reading the value of a header card.

        If ``headarr`` is a list of headers, the first header is used to construct the metadata value.

        Args:
            headarr (:obj:`list`):
                List of `astropy.io.fits.Header`_ objects.
            meta_key (:obj:`str`):
                Metadata keyword to construct.

        Returns:
            :obj:`object`: Metadata value read from the header(s).
        """
        if meta_key == 'mjd':
            # Build the MJD from the UT date and time.
            ttime = Time('{:s}T{:s}'.format(headarr[0]['UT-DATE'], headarr[0]['UT-TIME']),
                         format='isot')
            return ttime.mjd

        if meta_key == 'binning':
            # The BINNING card is written in the raw detector frame, which is
            # spatial along NAXIS1 and spectral along NAXIS2.
            binspatial, binspec = parse.parse_binning(headarr[0]['BINNING'])
            return parse.binning2string(binspec, binspatial)

        if meta_key == 'readmode':
            speed, gain = str(headarr[0]['SPEED']).strip().lower(), str(headarr[0]['GAIN']).strip().lower()
            return f'{speed},{gain}'

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

        # Slit edges.  Narrow detector defects otherwise get synced into spurious
        # few-pixel slits; the narrowest real LDSS3 slitlets are ~5 arcsec.
        par['calibrations']['slitedges']['edge_thresh'] = 20.
        par['calibrations']['slitedges']['minimum_slit_length'] = 2.

        # Processing steps
        # TODO: is there a reason why use_biasimage=False?
        turn_off = dict(use_darkimage=False, use_biasimage=False)
        par.reset_all_processimages_par(**turn_off)

        # Combine the (usually few) bias and dark frames with a median
        par['calibrations']['biasframe']['process']['combine'] = 'median'
        par['calibrations']['darkframe']['process']['combine'] = 'median'

        # Mask saturated slits; the alignment boxes of a slitmask are routinely
        # saturated in the flats.
        par['calibrations']['flatfield']['saturated_slits'] = 'mask'

        # Extraction
        par['reduce']['skysub']['bspline_spacing'] = 0.8
        par['reduce']['extraction']['sn_gauss'] = 4.0
        # Do not perform global sky subtraction for standard stars
        par['reduce']['skysub']['global_sky_std'] = False

        # Flexure.  The default wavelength solution comes from arc lamps taken
        # at a different telescope pointing, so correct the residual spectral
        # flexure against the sky.  May set this to 'skip' when calibrating on sky
        # lines instead (see the magellan_ldss3 documentation).
        par['flexure']['spec_method'] = 'boxcar'

        # Cosmic ray rejection parameters for science frames
        par['scienceframe']['process']['sigclip'] = 5.0
        par['scienceframe']['process']['objlim'] = 2.0

        # Set the default exposure time ranges for the frame typing.  The
        # science/standard split is driven entirely by these ranges.
        # The types have to be checked in the .pypeit file against an observing log.
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
                use :func:`~pypeit.spectrographs.magellan_ldss3.MagellanLDSS3Spectrograph.default_pypeit_par`.

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

        par['calibrations']['wavelengths']['lamps'] = ['HeI', 'NeI', 'ArI']
        par['calibrations']['wavelengths']['method'] = 'full_template'

        match grating:
            case 'VPH-Blue':
                par['calibrations']['wavelengths']['reid_arxiv'] = 'magellan_ldss3_vph_blue_HeINeIArI.fits'
            case 'VPH-Red':
                par['calibrations']['wavelengths']['reid_arxiv'] = 'magellan_ldss3_vph_red_HeINeIArI.fits'
            case 'VPH-All':
                par['calibrations']['wavelengths']['reid_arxiv'] = 'magellan_ldss3_vph_all_HeINeIArI.fits'
            case _:
                par['calibrations']['wavelengths']['method'] = 'holy-grail'
                log.warning(f'No archived LDSS3 wavelength solution for grism {grating}; '
                            'falling back on holy-grail.')
        # NOTE: A wavelength solution using sky lines on deep observations
        # is also available for VPH-Red, with the reid_arxiv=magellan_ldss3_vph_red_sky.fits
        # and lamps=OH_LDSS3_vac.

        # FWHM, in binned pixels
        if binning is not None:
            par['calibrations']['wavelengths']['fwhm'] = 5.0 / parse.parse_binning(binning)[0]

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
        # Observers are recommended to use the same filter for all exposures in a given configuration,
        # but it is not enforced and this can be overridden in the .pypeit file if necessary.
        return ['dispname', 'decker', 'binning', 'filter1']

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
        # Do not propagate ENOISE, EGAIN because values are set from the table.
        return ['GRISM', 'APERTURE', 'BINNING', 'FILTER', 'SPEED', 'GAIN']

    def modify_config(self, row, cfg):
        """
        Modify the configuration requirement for one frame.

        Arcs and flats are often taken without the order-blocking filter used
        for the science frames.  For any frame that is not a science or standard
        frame, the filter is ignored when matching it to a configuration.

        Args:
            row (`astropy.table.Row`_):
                The table row with the metadata for one frame.
            cfg (:obj:`dict`):
                Dictionary with metadata associated to a specific configuration.

        Returns:
            :obj:`dict`: The configuration to match this frame against.
        """
        ftype = str(row['frametype'])
        if 'science' in ftype or 'standard' in ftype:
            return cfg
        _cfg = copy.deepcopy(cfg)
        # Match on this frame's own filter, so the filter never blocks a match.
        _cfg['filter1'] = row['filter1']
        return _cfg

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.

        LDSS3's ``EXPTYPE`` card alone is not sufficient to identify frame types.

        This method implements a set of rules to determine the frame type, but it is not
        guaranteed to be correct.

        NOTE: The user should always double-check the frame typing in the .pypeit
        file and adjust it if necessary against an observing log.

        Args:
            ftype (:obj:`str`):
                Type of frame to check. Must be a valid frame type; see
                frame-type :ref:`frame_type_defs`.
            fitstbl (`astropy.table.Table`_):
                The table with the metadata for one or more frames to check.
            exprng (:obj:`list`, optional):
                Range in the allowed exposure time for a frame of type
                ``ftype``. See
                :func:`~pypeit.core.framematch.check_frame_exptime`.

        Returns:
            `numpy.ndarray`_: Boolean array with the flags selecting the
            exposures in ``fitstbl`` that are ``ftype`` type frames.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        idname = np.array([str(n).strip().lower() for n in fitstbl['idname']])
        target = np.array([str(n).strip().lower() for n in fitstbl['target']])
        grism = np.array([str(n).strip().lower() for n in fitstbl['dispname']])

        # Open grism means no disperser is in the beam.
        dispersed = grism != 'open'

        setup_frame = self._has_keyword(target, self.setup_keywords)

        is_arc = self._has_keyword(target, self.arc_keywords)
        is_flat = self._has_keyword(target, self.flat_keywords)
        # Bias and dark frames are identified from EXPTYPE alone.
        is_bias = idname == 'bias'
        is_dark = idname == 'dark'

        if ftype == 'bias':
            return good_exp & is_bias
        if ftype == 'dark':
            return good_exp & is_dark
        if ftype in ['arc', 'tilt']:
            # Checked before the flats because arc frames are sometimes written with
            # EXPTYPE='Flat'.
            return good_exp & dispersed & is_arc
        if ftype in ['pixelflat', 'illumflat', 'trace']:
            return good_exp & dispersed & is_flat & np.logical_not(is_arc)
        if ftype in ['science', 'standard']:
            # The science/standard split is set by the exposure-time ranges in
            # the 'scienceframe' and 'standardframe' parameters.
            is_cal = is_arc | is_flat | is_bias | is_dark
            return (good_exp & dispersed & np.logical_not(is_cal)
                    & np.logical_not(setup_frame))

        log.warning(f'Cannot determine if frames are of type {ftype}.')
        return np.zeros(len(fitstbl), dtype=bool)

    # OBJECT keywords used for frame typing.  They are matched as substrings,
    # ignoring case, so 'arc' also covers 'arcs' and 'wave' covers 'wavecal'.
    arc_keywords = ['arc', 'hene', 'hear', 'wave']
    # 'qh' and 'ql' also cover the combined 'qlqh' and 'qhql' lamp names.  A
    # science or standard target containing any of these is typed as a flat.
    flat_keywords = ['flat', 'quartz', 'dome', 'qh', 'ql']
    # Frames whose OBJECT contains one of these are never typed as science or
    # standard.  The setup keywords are used to identify setup frames that might
    # be taken with a disperser in the beam.
    setup_keywords = ['align', 'focus', 'acq']

    @staticmethod
    def _has_keyword(names, keywords):
        """
        Test each of ``names`` for any of the given keywords.

        Args:
            names (`numpy.ndarray`_):
                Lower-cased strings to test.
            keywords (:obj:`list`):
                Lower-case strings to look for anywhere in the name.

        Returns:
            `numpy.ndarray`_: Boolean array, one element per entry in ``names``.
        """
        return np.array([any(k in n for k in keywords) for n in names])

    def get_rawimage(self, raw_file, det):
        """
        Take a raw LDSS3 amplifier file and return the combined image for the detector.

        The two amplifier files are located by
        :func:`~pypeit.spectrographs.magellan_ldss3.MagellanLDSS3Spectrograph.amp_files`.
        Amplifier 2 is flipped
        about the spatial axis so that it abuts amplifier 1 at the center of the CCD, and
        the two are joined along the spatial axis with the overscan regions at the outer
        edges and the data sections in the middle::

            [ overscan_1 | data_1 | flip(data_2) | flip(overscan_2) ]

        If one amplifier is missing, its data and overscan sections are left as
        zeros and are trimmed away during processing.  If only amplifier 2 is
        present, it is labeled as amplifier 1 in the data-section and
        overscan-section images, because PypeIt numbers amplifiers from 1.

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
        amp_files = [f for f in self.amp_files(raw_file) if f is not None]
        if not amp_files:
            raise PypeItError(f'No amplifier file found for {raw_file}.')

        # amps is a dict indexed by the amplifier number (which are 1-indexed).
        amps = {}
        for amp_file in amp_files:
            data, overscan, header, hdu = ldss3_read_amp(amp_file)
            opamp = int(header['OPAMP'])
            amps[opamp] = {'data': data, 'overscan': overscan,
                                'header': header, 'hdu': hdu}
        # Use amplifier 1 as the reference for the detector parameters,
        # or amplifier 2 if amplifier 1 is missing.
        opamp = min(amps.keys())

        # Warn about binning we have not been able to validate.
        binning = self.get_meta_value([amps[opamp]['header']], 'binning')
        if binning != '1,1':
            log.warning('LDSS3 support has only been verified for 1x1 binning; found '
                        f'{binning}. Check the data and overscan sections carefully.')

        # Assemble.  Work in the raw detector frame, which is spatial along the
        # first axis and spectral along the second, then transpose at the end to
        # give the PypeIt (nspec, nspat) convention.
        ndata, nspec = amps[opamp]['data'].shape
        nx_oscan = amps[opamp]['overscan'].shape[0]
        array = np.zeros((2 * (ndata + nx_oscan), nspec), dtype=float)
        rawdatasec_img = np.zeros_like(array, dtype=int)
        oscansec_img = np.zeros_like(array, dtype=int)

        # Fill in the data and overscan sections for each amplifier.  The
        # amplifier 2 data and overscan are flipped about the spatial axis.
        if 1 in amps:
            array[:nx_oscan, :] = amps[1]['overscan']
            oscansec_img[:nx_oscan, :] = 1
            array[nx_oscan:nx_oscan+ndata, :] = amps[1]['data']
            rawdatasec_img[nx_oscan:nx_oscan+ndata, :] = 1

        if 2 in amps:
            # If only amplifier 2 is present, set the amplifier id to 1 in the data-section and overscan-section images.
            # This is needed because PypeIt expects the amplifier id to be 1-indexed.
            # Cannot set amp_id=2 if amp_id=1 does not exist.
            amp2_id = 2 if 1 in amps else 1

            start = nx_oscan + ndata
            array[start:start+ndata, :] = np.flipud(amps[2]['data'])
            rawdatasec_img[start:start+ndata, :] = amp2_id
            array[start+ndata:start+ndata+nx_oscan, :] = np.flipud(amps[2]['overscan'])
            oscansec_img[start+ndata:start+ndata+nx_oscan, :] = amp2_id

        detector_par = self.get_detector_par(det,
                                             amp_headers=[amps[i]['header'] for i in amps])

        exptime = self.get_meta_value([amps[opamp]['header']], 'exptime')
        hdu = amps[opamp]['hdu']  # return the first amplifier's HDU for the header

        # Close the other amplifier HDUs to avoid leaving open file handles.
        for i in amps:
            if i != opamp:
                amps[i]['hdu'].close()

        return detector_par, array.T, hdu, exptime, rawdatasec_img.T, oscansec_img.T

    def bpm(self, filename, det, shape=None, msbias=None):
        """
        Generate a default bad-pixel mask.

        NOTE: Defined for the two-amplifier joined image. When only
        one amplifier file is present, the bad-pixel mask is trimmed to the
        data section of that amplifier.

        Args:
            filename (:obj:`str` or None):
                Either amplifier file, used to get the spatial binning and, for a
                single-amplifier image, identify the amplifier number.
            det (:obj:`int`):
                1-indexed detector number to use when getting the image
                shape from the input file.
            shape (:obj:`tuple`, optional):
                Processed image shape, (nspec, nspat).  Required if filename is
                None.
            msbias (`numpy.ndarray`_, optional):
                Processed bias frame used to identify bad pixels. Not used.

        Returns:
            `numpy.ndarray`_: An integer array with a masked value set
            to 1 and an unmasked value set to 0.
        """
        bpm_img = super().bpm(filename, det, shape=shape, msbias=msbias)
        log.info("Using hard-coded BPM")

        if filename:
            # Get the binning from the header of the example file
            binning = parse.parse_binning(self.get_meta_value(filename, 'binning'))
            _, binspat = binning
        else:
            # Guess the binning from the shape of bpm_img
            bpm_shape = bpm_img.shape
            if self.n_spat_unbinned % bpm_shape[1] != 0:
                raise PypeItError(f'Cannot determine the spatial binning: bpm shape {bpm_shape} '
                                  'does not divide the unbinned width '
                                  f'{self.n_spat_unbinned}.')
            binspat = self.n_spat_unbinned // bpm_shape[1]

        # Bad columns, measured from the column
        # median of combined bias frames.
        # The first and last entries are the outer edges
        # of amplifiers 1 and 2.
        nspat_unbinned = bpm_img.shape[1] * binspat
        offset = 0

        if nspat_unbinned < self.n_spat_unbinned:
            # This means that only one amplifier is present, so we need to determine which one it is.
            parsed = None if filename is None else self.parse_amp_file(filename)
            if parsed is None:
                raise PypeItError('Cannot identify the amplifier of a single-amplifier LDSS3 '
                                  'image without an LDSS3 amplifier file.')
            amp = parsed[1]
            # If only amplifier 2 is present, the bad columns are offset by the width of amplifier 1.
            offset = (amp - 1) * nspat_unbinned

        # These coordinates are defined in the unbinned, two-amplifier frame,
        # so we need to adjust for the offset of the single amplifier if only one is present
        # because Pypeit trims the bpm to the data section of the amplifier.
        # The offset is 0 for amplifier 1 and nspat_unbinned for amplifier 2.
        for c1, c2 in [(0, 11), (443, 443), (608, 608), (1413, 1413), (1492, 1494),
                       (1549, 1551), (1602, 1606), (1635, 1640), (1688, 1690),
                       (1695, 1695), (1999, 2000), (2034, 2047)]:
            lo = max((c1 - offset) // binspat, 0)
            hi = min((c2 - offset) // binspat + 1, bpm_img.shape[1])
            if lo < hi:
                bpm_img[:, lo:hi] = 1

        return bpm_img

def ldss3_read_amp(filename):
    """
    Read a single amplifier file of LDSS3 data.

    Args:
        filename (:obj:`str`, `Path`_):
            Filename.

    Returns:
        :obj:`tuple`: The data section and overscan section of the amplifier,
        each a `numpy.ndarray`_ ordered as (spatial, spectral); the primary
        `astropy.io.fits.Header`_; and the opened `astropy.io.fits.HDUList`_.
    """
    log.debug(f'Reading LDSS3 amplifier file: {filename}')
    hdu = io.fits_open(filename)
    array, header = hdu[0].data, hdu[0].header

    array = array.T * 1.0
    data = array[parse.sec2slice(header['DATASEC'], one_indexed=True, include_end=True)]
    overscan = array[parse.sec2slice(header['BIASSEC'], one_indexed=True, include_end=True)]

    return data, overscan, header, hdu
