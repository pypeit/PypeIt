"""
Module for LDT/DeVeny specific methods.

.. include:: ../include/links.rst
"""
import numpy as np

from astropy.time import Time

from pypeit import io
from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.core import parse
from pypeit.images import detector_container
from pypeit.spectrographs import spectrograph


class LDTDeVenySpectrograph(spectrograph.Spectrograph):
    """
    Child to handle LDT/DeVeny specific code
    """
    ndet = 1
    name = 'ldt_deveny'
    telescope = telescopes.LDTTelescopePar()
    camera = 'DeVeny'
    url = 'https://lowell.edu/research/telescopes-and-facilities/ldt/deveny-optical-spectrograph/'
    header_name = 'Deveny'
    comment = 'LDT DeVeny Optical Spectrograph'
    supported = True

    # Parameters equal to the PypeIt defaults, shown here for completeness
    # pypeline = 'MultiSlit'

    def get_detector_par(self, det, hdu=None):
        """
        Return metadata for the selected detector.

        .. warning::

            Many of the necessary detector parameters are read from the file
            header, meaning the ``hdu`` argument is effectively **required** for
            LTD/DeVeny.  The optional use of ``hdu`` is only viable for
            automatically generated documentation.

        Args:
            det (:obj:`int`):
                1-indexed detector number.
            hdu (`astropy.io.fits.HDUList`_, optional):
                The open fits file with the raw image of interest.

        Returns:
            :class:`~pypeit.images.detector_container.DetectorContainer`:
            Object with the detector metadata.
        """
        if hdu is None:
            binning = '1,1'                 # Most common use mode
            gain = np.atleast_1d(1.52)      # Hardcoded in the header
            ronoise = np.atleast_1d(4.9)    # Hardcoded in the header
            datasec = np.atleast_1d('[5:512,53:2095]')   # For 1x1 binning
            oscansec = np.atleast_1d('[5:512,5:48]')     # For 1x1 binning
        else:
            binning = self.get_meta_value(self.get_headarr(hdu), 'binning')
            gain = np.atleast_1d(hdu[0].header['GAIN'])
            ronoise = np.atleast_1d(hdu[0].header['RDNOISE'])
            datasec = self.rotate_trimsections(hdu[0].header['TRIMSEC'],
                                               hdu[0].header['NAXIS1'])
            oscansec = self.rotate_trimsections(hdu[0].header['BIASSEC'],
                                               hdu[0].header['NAXIS1'])

        # Detector
        detector_dict = dict(
            binning         = binning,
            det             = 1,
            dataext         = 0,
            specaxis        = 1,        # Native spectrum is along the x-axis
            specflip        = True,     # DeVeny CCD has blue at the right
            spatflip        = False,
            platescale      = 0.34,     # Arcsec / pixel
            darkcurr        = 4.5,      # Electrons per hour
            saturation      = 65535.,   # 16-bit ADC
            nonlinear       = 0.97,     # Linear to ~97% of saturation
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = gain,     # See above
            ronoise         = ronoise,  # See above
            # Data & Overscan Sections -- Edge tracing can handle slit edges
            datasec         = datasec,  # See above
            oscansec        = oscansec  # See above
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
        self.meta['target'] = dict(ext=0, card='OBJNAME')
        self.meta['dispname'] = dict(card=None, compound=True)
        self.meta['decker'] = dict(card=None, compound=True)
        self.meta['binning'] = dict(card=None, compound=True)
        self.meta['mjd'] = dict(card=None, compound=True)
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['instrument'] = dict(ext=0, card='INSTRUME')

        # Extras for config and frametyping
        self.meta['idname'] = dict(ext=0, card='IMAGETYP')
        self.meta['dispangle'] = dict(ext=0, card='GRANGLE', rtol=1e-3)
        self.meta['cenwave'] = dict(card=None, compound=True, rtol=2.0)
        self.meta['filter1'] = dict(card=None, compound=True)
        self.meta['slitwid'] = dict(ext=0, card='SLITASEC')
        self.meta['lampstat01'] = dict(card=None, compound=True)

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
        if meta_key == 'binning':
            # Binning in lois headers is space-separated (like Gemini).
            binspec, binspatial = parse.parse_binning(headarr[0]['CCDSUM'])
            return parse.binning2string(binspec, binspatial)

        if meta_key == 'mjd':
            # Use astropy to convert 'DATE-OBS' into a mjd.
            ttime = Time(headarr[0]['DATE-OBS'], format='isot')
            return ttime.mjd

        if meta_key == 'lampstat01':
            # The spectral comparison lamps turned on are listed in `LAMPCAL`, but
            #  if no lamps are on, then this string is blank.  Return either the
            #  populated `LAMPCAL` string, or 'off' to ensure a positive entry for
            #  `lampstat01`.
            lampcal = headarr[0]['LAMPCAL'].strip()
            return 'off' if lampcal == '' else lampcal

        if meta_key == 'dispname':
            # Convert older FITS keyword GRATING (gpmm/blaze) into the newer
            #  Grating ID names (DVx) for easier identification of disperser.
            gratings = {"150/5000":"DV1", "300/4000":"DV2", "300/6750":"DV3",
                        "400/8500":"DV4", "500/5500":"DV5", "600/4900":"DV6",
                        "600/6750":"DV7", "831/8000":"DV8", "1200/5000":"DV9",
                        "2160/5000":"DV10", "UNKNOWN":"DVxx"}
            if headarr[0]['GRATING'] not in gratings:
                raise ValueError(f"Grating value {headarr[0]['GRATING']} not recognized.")
            return f"{gratings[headarr[0]['GRATING']]} ({headarr[0]['GRATING']})"

        if meta_key == 'decker':
            # Provide a stub for future inclusion of a decker on LDT/DeVeny.
            return headarr[0]['DECKER'] if 'DECKER' in headarr[0].keys() else 'None'

        if meta_key == 'filter1':
            # Remove the parenthetical knob position to leave just the filter name
            return headarr[0]['FILTREAR'].split()[0].upper()

        if meta_key == 'cenwave':
            # The central wavelength is more descriptive of the grating angle position
            #  than just the angle value.  Use the DeVeny grating angle formula to
            #  return the central wavelength of the configuration.

            # Extract lines/mm, catch 'UNKNOWN' grating
            lpmm = (
                np.inf
                if headarr[0]["GRATING"] == "UNKNOWN"
                else float(headarr[0]["GRATING"].split("/")[0])
            )

            # DeVeny Fixed Optical Angles in radians
            CAMCOL = np.deg2rad(55.00)  # Camera-to-Collimator Angle
            COLL = np.deg2rad(10.00)    # Collimator-to-Grating Angle
            # Grating angle in radians
            theta = np.deg2rad(float(headarr[0]['GRANGLE']))
            # Wavelength in Angstroms
            wavelen = (
                (np.sin(COLL + theta) + np.sin(COLL + theta - CAMCOL))  # Angles
                * 1.0e7                                                 # Angstroms/mm
                / lpmm                                                  # lines / mm
            )
            # Round the wavelength to the nearest 5A
            return np.around(wavelen / 5, decimals=0) * 5

        msgs.error(f"Not ready for compound meta {meta_key} for LDT/DeVeny")

    def configuration_keys(self):
        """
        Return the metadata keys that define a unique instrument
        configuration.

        This list is used by :class:`~pypeit.metadata.PypeItMetaData` to
        identify the unique configurations among the list of frames read
        for a given reduction.

        Returns:
            :obj:`list`: List of keywords of data pulled from file headers
            and used to constuct the :class:`~pypeit.metadata.PypeItMetaData`
            object.
        """
        return ['dispname', 'cenwave', 'filter1', 'binning']

    def raw_header_cards(self):
        """
        Return additional raw header cards to be propagated in
        downstream output files for configuration identification.

        The list of raw data FITS keywords should be those used to populate
        the :meth:`~pypeit.spectrograph.Spectrograph.configuration_keys`
        or are used in :meth:`~pypeit.spectrograph.Spectrograph.config_specific_par`
        for a particular spectrograph, if different from the name of the
        PypeIt metadata keyword.

        This list is used by :meth:`~pypeit.spectrograph.Spectrograph.subheader_for_spec`
        to include additional FITS keywords in downstream output files.

        Returns:
            :obj:`list`: List of keywords from the raw data files that should
            be propagated in output files.
        """
        return ['GRATING', 'GRANGLE', 'FILTREAR', 'CCDSUM']

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.

        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of ``PypeIt`` methods.
        """
        par = super().default_pypeit_par()

        # Turn off illumflat unless/until we can deal properly with flexure in
        #   the spatial direction.  All other defaults OK (as of v1.7.0)
        set_use = dict(use_illumflat=False)
        par.reset_all_processimages_par(**set_use)

        # Make a bad pixel mask
        par['calibrations']['bpm_usebias'] = True

        # Wavelength Calibration Parameters
        # Do not sigmaclip the arc frames for better MasterArc and better wavecalib
        par['calibrations']['arcframe']['process']['clip'] = False
        # Do not sigmaclip the tilt frames
        par['calibrations']['tiltframe']['process']['clip'] = False
        # Arc lamps list from header -- instead of defining the full list here
        par['calibrations']['wavelengths']['lamps'] = ['use_header']
        #par['calibrations']['wavelengths']['lamps'] = ['NeI', 'ArI', 'CdI', 'HgI']
        # Set this as default... but use `holy-grail` for DV4, DV8
        par['calibrations']['wavelengths']['method'] = 'full_template'
        # The DeVeny arc line FWHM varies based on slitwidth used
        par['calibrations']['wavelengths']['fwhm_fromlines'] = True
        par['calibrations']['wavelengths']['nsnippet'] = 1  # Default: 2
        # Because of the wide wavelength range, solution more non-linear; user higher orders
        par['calibrations']['wavelengths']['n_first'] = 3  # Default: 2
        par['calibrations']['wavelengths']['n_final'] = 5  # Default: 4

        # Slit-edge settings for long-slit data (DeVeny's slit is > 90" long)
        par['calibrations']['slitedges']['bound_detector'] = True
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'
        par['calibrations']['slitedges']['minimum_slit_length'] = 90.

        # For the tilts, our lines are not as well-behaved as others',
        #   possibly due to the Wynne type E camera.
        par['calibrations']['tilts']['spat_order'] = 4  # Default: 3
        par['calibrations']['tilts']['spec_order'] = 5  # Default: 4

        # Cosmic ray rejection parameters for science frames
        par['scienceframe']['process']['sigclip'] = 5.0  # Default: 4.5
        par['scienceframe']['process']['objlim'] = 2.0   # Default: 3.0

        # Reduction and Extraction Parameters -- Look for fainter objects
        par['reduce']['findobj']['snr_thresh'] = 5.0   # Default: 10.0

        # Flexure Correction Parameters
        par['flexure']['spec_method'] = 'boxcar'  # Default: 'skip'

        # Sensitivity Function Parameters
        par['sensfunc']['polyorder'] = 7  # Default: 5

        return par

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
        if ftype == 'bias':
            return fitstbl['idname'] == 'BIAS'
        if ftype in ['arc', 'tilt']:
            # FOCUS frames should have frametype None, BIAS is bias regardless of lamp status
            return (
                good_exp
                & (fitstbl['lampstat01'] != 'off')
                & (fitstbl['idname'] != 'FOCUS')
                & (fitstbl['idname'] != 'BIAS')
            )
        if ftype in ['trace', 'pixelflat']:
            return (
                good_exp
                & (fitstbl['idname'] == 'DOME FLAT')
                & (fitstbl['lampstat01'] == 'off')
            )
        if ftype in ['illumflat','sky']:
            return (
                good_exp
                & (fitstbl['idname'] == 'SKY FLAT')
                & (fitstbl['lampstat01'] == 'off')
            )
        if ftype == 'science':
            # Both OBJECT and STANDARD frames should be processed as science frames
            return (
                good_exp
                & ((fitstbl['idname'] == 'OBJECT') | (fitstbl['idname'] == 'STANDARD'))
                & (fitstbl['lampstat01'] == 'off')
            )
        if ftype == 'standard':
            return (
                good_exp
                & (fitstbl['idname'] == 'STANDARD')
                & (fitstbl['lampstat01'] == 'off')
            )
        if ftype == 'dark':
            return (
                good_exp
                & (fitstbl['idname'] == 'DARK')
                & (fitstbl['lampstat01'] == 'off')
            )
        if ftype in ['pinhole','align']:
            # Don't types pinhole or align frames
            return np.zeros(len(fitstbl), dtype=bool)
        msgs.warn(f"Cannot determine if frames are of type {ftype}")
        return np.zeros(len(fitstbl), dtype=bool)

    def pypeit_file_keys(self):
        """
        Define the list of keys to be output into a standard ``PypeIt`` file.

        Returns:
            :obj:`list`: The list of keywords in the relevant
            :class:`~pypeit.metadata.PypeItMetaData` instance to print to the
            :ref:`pypeit_file`.
        """
        return super().pypeit_file_keys() + ['dispangle','slitwid','lampstat01']

    def get_lamps(self, fitstbl):
        """
        Extract the list of arc lamps used from header

        .. note::

            There are some faint Cd and Hg lines in the DV9 spectra that are
            helpful for nailing down the wavelength calibration for that
            grating, but these lines are too faint / close to other lines for
            use with other gratings.  Therefore, use a grating-specific line
            list for Cd and Hg with DV9, but use the usual ion lists for
            everything else.

        Args:
            fitstbl (`astropy.table.Table`_):
                The table with the metadata for one or more arc frames.
        Returns:
            lamps (:obj:`list`) : List used arc lamps
        """
        grating = fitstbl['dispname'][0].split()[0]              # Get the DVn specifier
        return [
            f"{lamp.strip()}_DeVeny1200"                         # Instrument-specific list
            if grating == "DV9" and lamp.strip() in ["Cd", "Hg"] # Under these conditions
            else f"{lamp.strip()}I"                              # Otherwise, the usuals
            for lamp in np.unique(
                np.concatenate([lname.split(",") for lname in fitstbl["lampstat01"]])
            )
        ]

    def config_specific_par(self, scifile, inp_par=None):
        """
        Modify the ``PypeIt`` parameters to hard-wired values used for
        specific instrument configurations.

        Args:
            scifile (:obj:`str`):
                File to use when determining the configuration and how
                to adjust the input parameters.
            inp_par (:class:`~pypeit.par.parset.ParSet`, optional):
                Parameter set used for the full run of PypeIt.  If None,
                use :func:`default_pypeit_par`.

        Returns:
            :class:`~pypeit.par.parset.ParSet`: The PypeIt parameter set
            adjusted for configuration specific parameter values.
        """
        # Start with instrument wide
        par = super().config_specific_par(scifile, inp_par=inp_par)

        # Set parameters based on grating used:
        grating = self.get_meta_value(scifile, 'dispname')

        if grating == 'DV1 (150/5000)':
            # Use this `reid_arxiv` with the `full-template` method:
            par['calibrations']['wavelengths']['reid_arxiv'] = 'ldt_deveny_150_HgCdAr.fits'
            # Because of the wide wavelength range, split DV1 arcs in half for reidentification
            par['calibrations']['wavelengths']['nsnippet'] = 2

        elif grating in ['DV2 (300/4000)', 'DV3 (300/6750)']:
            # Use this `reid_arxiv` with the `full-template` method:
            par['calibrations']['wavelengths']['reid_arxiv'] = 'ldt_deveny_300_HgCdAr.fits'

        elif grating == 'DV4 (400/8000)':
            # We don't have a good `reid_arxiv`` for this grating yet; use `holy-grail`
            #  and it's associated tweaks in parameters
            par['calibrations']['wavelengths']['method'] = 'holy-grail'
            par['calibrations']['wavelengths']['sigdetect'] = 10.0  # Default: 5.0
            par['calibrations']['wavelengths']['rms_threshold'] = 0.5  # Default: 0.15

        elif grating == 'DV5 (500/5500)':
            # Use this `reid_arxiv` with the `full-template` method:
            par['calibrations']['wavelengths']['reid_arxiv'] = 'ldt_deveny_500_HgCdAr.fits'
            par['calibrations']['wavelengths']['n_first'] = 2  # Default: 3

        elif grating in ['DV6 (600/4900)', 'DV7 (600/6750)']:
            # Use this `reid_arxiv` with the `full-template` method:
            par['calibrations']['wavelengths']['reid_arxiv'] = 'ldt_deveny_600_HgCdAr.fits'
            # Given the narrow range of wavelengths, use a lower initial order of the fit
            par['calibrations']['wavelengths']['n_first'] = 2  # Default: 3

        elif grating == 'DV8 (831/8000)':
            # We don't have a good `reid_arxiv`` for this grating yet; use `holy-grail`
            #  and it's associated tweaks in parameters
            par['calibrations']['wavelengths']['method'] = 'holy-grail'
            par['calibrations']['wavelengths']['sigdetect'] = 10.0  # Default: 5.0
            par['calibrations']['wavelengths']['rms_threshold'] = 0.5  # Default: 0.15
            # Given the narrow range of wavelengths, use a lower final order of the fit
            par['calibrations']['wavelengths']['n_first'] = 2  # Default: 3
            par['calibrations']['wavelengths']['n_final'] = 4  # Default: 5

        elif grating == 'DV9 (1200/5000)':
            # Use this `reid_arxiv` with the `full-template` method:
            par['calibrations']['wavelengths']['reid_arxiv'] = 'ldt_deveny_1200_HgCdAr.fits'
            # Given the narrow range of wavelengths, use a lower final order of the fit
            par['calibrations']['wavelengths']['n_first'] = 2  # Default: 3
            par['calibrations']['wavelengths']['n_final'] = 4  # Default: 5

        elif grating == 'DV10 (2160/5000)':
            # Presently unsupported; no parameter changes
            pass

        else:
            pass

        return par

    def get_rawimage(self, raw_file, det):
        """
        Read raw images and generate a few other bits and pieces
        that are key for image processing.

        For LDT/DeVeny, the LOIS control system automatically adjusts the
        DATASEC and OSCANSEC regions if the CCD is used in a binning other
        than 1x1.  The get_rawimage() method in the base class assumes these
        sections are fixed and adjusts them based on the binning -- incorrect
        for this instrument.

        This method is a stripped-down version of the base class method and
        additionally does NOT send the binning to parse.sec2slice().

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
            Exposure time *in seconds*.
        rawdatasec_img : `numpy.ndarray`_
            Data (Science) section of the detector as provided by setting the
            (1-indexed) number of the amplifier used to read each detector
            pixel. Pixels unassociated with any amplifier are set to 0.
        oscansec_img : `numpy.ndarray`_
            Overscan section of the detector as provided by setting the
            (1-indexed) number of the amplifier used to read each detector
            pixel. Pixels unassociated with any amplifier are set to 0.
        """
        # Open
        hdu = io.fits_open(raw_file)

        # Grab the DetectorContainer and extract the raw image
        detector = self.get_detector_par(det, hdu=hdu)
        raw_img = hdu[detector['dataext']].data.astype(float)

        # Exposure time (used by RawImage) from the header
        headarr = self.get_headarr(hdu)
        exptime = self.get_meta_value(headarr, 'exptime')

        for section in ['datasec', 'oscansec']:
            # Get the data section from Detector
            image_sections = detector[section]

            # Initialize the image (0 means no amplifier)
            pix_img = np.zeros(raw_img.shape, dtype=int)
            for i in range(detector['numamplifiers']):

                if image_sections is not None:
                    # Convert the (FITS) data section from a string to a slice
                    # DO NOT send the binning (default: None)
                    datasec = parse.sec2slice(image_sections[i], one_indexed=True,
                                              include_end=True, require_dim=2)
                    # Assign the amplifier
                    pix_img[datasec] = i+1

            # Finish
            if section == 'datasec':
                rawdatasec_img = pix_img.copy()
            else:
                oscansec_img = pix_img.copy()

        # Return
        return detector, raw_img, hdu, exptime, rawdatasec_img, oscansec_img

    def rotate_trimsections(self, section_string, nspecpix):
        """
        In order to orient LDT/DeVeny images into the PypeIt-standard
        configuration, frames are essentially rotated 90º clockwise.  As such,
        x' = y and y' = -x.

        The TRIMSEC / BIASSEC FITS keywords in LDT/DeVeny data specify the
        proper regions to be trimmed for the data and overscan arrays,
        respectively, in the native orientation.  This method performs the
        rotation and returns the section in the Numpy image section required
        by the PypeIt processing routines.

        The LDT/DeVeny FITS header lists the sections as '[SPEC_SEC,SPAT_SEC]'.

        Args:
            section_string (:obj:`str`):
                The FITS keyword string to be parsed / translated
            nspecpix (:obj:`int`):
                The total number of pixels in the spectral direction
        Returns:
            section (:obj:`numpy.ndarray`): Numpy image section needed by PypeIt
        """
        # Split out the input section into spectral and spatial pieces
        spec_sec, spat_sec = section_string.strip('[]').split(',')

        # The spatial section is unchanged, but the spectral section flips
        #  Add 1 because the pixels are 1-indexed (FITS standard)
        # TODO: Why does this need to be cast specifically as int32?
        y2p, y1p = nspecpix - np.int32(spec_sec.split(':')) + 1

        # Return the PypeIt-standard Numpy array
        return np.atleast_1d(f"[{spat_sec},{y1p}:{y2p}]")
