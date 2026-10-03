"""
Module for Subaru MOIRCS

.. include:: ../include/links.rst
"""
import re
from pathlib import Path

import numpy as np
from astropy.io import fits
from IPython import embed

from pypeit import log, PypeItError
from pypeit import dataPaths
from pypeit import telescopes
from pypeit.core import framematch, parse
from pypeit.images import detector_container
from pypeit.spectrographs import spectrograph


class SubaruMOIRCSSpectrograph(spectrograph.Spectrograph):
    """
    Child of Spectrograph to handle Subaru/MOIRCS specific code.

    MOIRCS has two optically independent channels, each imaged onto its own
    Hawaii-2RG detector and written to its own FITS file.  For each exposure,
    the chip-1 file (odd frame number) is the file listed in the PypeIt file,
    and the chip-2 file (the next frame number, same ``EXP-ID``) is located and
    read automatically when ``det = 2`` is requested; see
    :func:`companion_file`.
    """

    ndet = 2
    telescope = telescopes.SubaruTelescopePar()
    url = "https://www.naoj.org/Instruments/MOIRCS/index.html"

    name = "subaru_moircs"
    camera = "MOIRCS"
    header_name = "MOIRCS"
    supported = True
    comment = "Without using the mask definition files"

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.

        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        # Wavelengths: OH sky lines in the science frames.  Grism-specific
        # choices (e.g. templates) are set in config_specific_par.
        par["calibrations"]["wavelengths"]["lamps"] = ["OH_NIRES"]
        par["calibrations"]["wavelengths"]["method"] = "holy-grail"
        par["calibrations"]["wavelengths"]["rms_thresh_frac_fwhm"] = 0.11
        par["calibrations"]["wavelengths"]["sigdetect"] = 5.0
        par["calibrations"]["wavelengths"]["fwhm"] = 5.0
        par["calibrations"]["wavelengths"]["n_final"] = 4

        # Slit edges.  With these values, edge tracing recovers every slit
        # in the HK500 dev-suite mask (17 + 3 boxes on chip 1, 15 + 4 on
        # chip 2).  PCA is not used because the slits have very different
        # spectral extents.
        par["calibrations"]["slitedges"]["edge_thresh"] = 50.0
        par["calibrations"]["slitedges"]["sync_predict"] = "nearest"
        par["calibrations"]["slitedges"]["fit_order"] = 3
        par["calibrations"]["slitedges"]["max_shift_adj"] = 0.5
        # Alignment-star boxes are ~4.4 arcsec; flag them as boxes so they
        # are not reduced as science slits
        par['calibrations']['slitedges']['minimum_slit_length_sci'] = 5.
        # Remove slits that are too short
        par['calibrations']['slitedges']['minimum_slit_length'] = 3.

        # Tilts, from the OH lines
        par["calibrations"]["tilts"]["tracethresh"] = 25.0
        par["calibrations"]["tilts"]["spat_order"] = 3
        par["calibrations"]["tilts"]["spec_order"] = 4

        # No bias, overscan, or dark frames
        turn_off = dict(use_biasimage=False, use_overscan=False,
                        use_darkimage=False)
        par.reset_all_processimages_par(**turn_off)

        # Science-frame processing
        par["scienceframe"]["process"]["sigclip"] = 20.0
        par["scienceframe"]["process"]["satpix"] = "nothing"

        # Object finding: MOS targets are typically faint in the A-B
        # images, which are missed with the default threshold (10)
        par["reduce"]["findobj"]["snr_thresh"] = 5.0

        # Sky subtraction and extraction
        par["reduce"]["skysub"]["bspline_spacing"] = 0.8
        par["reduce"]["extraction"]["sn_gauss"] = 4.0

        # The wavelength solution comes from the science frames themselves,
        # so no spectral flexure correction is needed (as for other NIR
        # spectrographs calibrated on OH lines)
        par["flexure"]["spec_method"] = "skip"

        # Set the default exposure time ranges for the frame typing
        par["calibrations"]["standardframe"]["exprng"] = [None, 20]
        par["calibrations"]["arcframe"]["exprng"] = [1, None]
        par["calibrations"]["darkframe"]["exprng"] = [1, None]
        par["scienceframe"]["exprng"] = [20, None]

        # Sensitivity function parameters (not yet tested on MOIRCS data)
        par["sensfunc"]["extrap_blu"] = 0.0
        par["sensfunc"]["extrap_red"] = 0.0
        par["fluxcalib"]["extrap_sens"] = True
        par["sensfunc"]["algorithm"] = "IR"
        par["sensfunc"]["polyorder"] = 13
        par["sensfunc"]["IR"]["maxiter"] = 2
        par["sensfunc"]["IR"]["telgridfile"] \
            = "TelFit_MaunaKea_3100_26100_R20000.fits"

        return par

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)
        self.meta["ra"] = dict(
            ext=0, card="RA", required_ftypes=["science", "standard"]
        )  # Need to convert to : separated
        self.meta["dec"] = dict(
            ext=0, card="DEC", required_ftypes=["science", "standard"]
        )
        self.meta["target"] = dict(ext=0, card="OBJECT")
        self.meta["binning"] = dict(ext=0, card=None, compound=True)

        self.meta["mjd"] = dict(ext=0, card="MJD")
        self.meta["exptime"] = dict(ext=0, card="EXPTIME")
        self.meta["airmass"] = dict(ext=0, card="AIRMASS")
        #
        self.meta["decker"] = dict(ext=0, card="SLIT")

        # Extras for config and frametyping
        self.meta["dispname"] = dict(
            ext=0, card="DISPERSR", required_ftypes=["science", "standard"]
        )
        # DATA-TYP, refined using OBJECT for the dome-flat variants
        self.meta["idname"] = dict(ext=0, card=None, compound=True)
        self.meta["lampstat01"] = dict(ext=0, card=None, compound=True)
        # Chip ID; used only to keep chip-2 files out of the metadata table
        self.meta["detector"] = dict(ext=0, card="DET-ID")
        self.meta["instrument"] = dict(ext=0, card="INSTRUME")
        self.meta["frameno"] = dict(ext=0, card="FRAMEID")

        # Dithering
        self.meta["dithpat"] = dict(ext=0, card=None, compound=True)
        self.meta["dithpos"] = dict(ext=0, card=None, compound=True)
        self.meta["dithoff"] = dict(ext=0, card=None, compound=True)

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
        hdr = headarr[0]
        if meta_key == "binning":
            # BIN-FCT1 is along x, which is the dispersion axis
            # (DISPAXIS = 1, specaxis = 1); BIN-FCT2 is along the slit.
            binspec = hdr.get("BIN-FCT1", 1)
            binspatial = hdr.get("BIN-FCT2", 1)
            return parse.binning2string(binspec, binspatial)

        if meta_key == "idname":
            # Lamp-on flats, lamp-off flats and the mask image all have
            # DATA-TYP = DOMEFLAT; only OBJECT distinguishes them.
            datatyp = str(hdr.get("DATA-TYP", "")).strip()
            obj = str(hdr.get("OBJECT", "")).strip().upper()
            if datatyp == "DOMEFLAT" and obj in ["DOMEFLAT_OFF",
                                                 "MASKIMAGE"]:
                return obj
            # Internal-lamp frames have DATA-TYP = INSTFLAT; OBJECT gives
            # the lamp (e.g. TH-AR for the ThAr arcs)
            if datatyp == "INSTFLAT" and obj != "":
                return obj
            return datatyp

        if meta_key == "lampstat01":
            # Dome-flat lamp status; only lamp-on dome flats are 'on'
            datatyp = str(hdr.get("DATA-TYP", "")).strip()
            obj = str(hdr.get("OBJECT", "")).strip().upper()
            return "on" if datatyp == "DOMEFLAT" and obj == "DOMEFLAT" \
                else "off"

        if meta_key in ["dithpat", "dithpos", "dithoff"]:
            return self._parse_dither(hdr, meta_key)

        raise PypeItError(f"Not ready for compound meta {meta_key}")

    @staticmethod
    def _parse_dither(hdr, meta_key):
        """
        Interpret the MOIRCS dither header cards.

        The MOIRCS dither cards are ``K_DITPAT`` (pattern name), ``K_DITCNT``
        (1-indexed position within the pattern) and ``K_DITWID`` (dither
        length in arcsec).  Only the two-position ``LINE2`` pattern has been
        seen so far: position 1 is taken to be A and position 2 to be B, with
        offsets of +/- ``K_DITWID``/2 along the slit.  The sign convention of
        the offset has not been verified.  For any other pattern, the
        position is reported as ``P<N>`` (so it is not paired automatically
        by :func:`get_comb_group`) and the offset is 0.

        Generated by JXP and Claude.

        Args:
            hdr (`astropy.io.fits.Header`_):
                Primary header of the chip-1 file.
            meta_key (:obj:`str`):
                One of ``dithpat``, ``dithpos``, or ``dithoff``.

        Returns:
            :obj:`str` or :obj:`float`: The dither pattern name, the dither
            position (``A``, ``B``, ``P<N>``, or ``none``) or the dither
            offset in arcsec.
        """
        pattern = str(hdr.get("K_DITPAT", "NONE")).strip()
        count = hdr.get("K_DITCNT", 0)
        width = hdr.get("K_DITWID", 0.0)
        # Not dithered (e.g. calibrations)
        no_dither = pattern.upper() in ["NONE", ""] or count in [None, 0]

        if meta_key == "dithpat":
            return "none" if no_dither else pattern
        if meta_key == "dithpos":
            if no_dither:
                return "none"
            if pattern == "LINE2" and int(count) in [1, 2]:
                return "A" if int(count) == 1 else "B"
            return f"P{int(count)}"
        # dithoff
        if no_dither or pattern != "LINE2" or width is None:
            return 0.0
        return 0.5 * float(width) * (1. if int(count) == 1 else -1.)

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
        good_exp = framematch.check_frame_exptime(fitstbl["exptime"], exprng)
        if ftype in ["science", "standard"]:
            # Science and standards are separated by exposure time only
            # (see default_pypeit_par).  Standards may have DATA-TYP =
            # STANDARD_STAR (e.g. VB_K data) or OBJECT.
            return good_exp & ((fitstbl["idname"] == "OBJECT")
                               | (fitstbl["idname"] == "STANDARD_STAR"))
        if ftype in ["arc", "tilt"]:
            # Wavelengths and tilts come from the OH lines in the science
            # frames
            return good_exp & (fitstbl["idname"] == "OBJECT")
        if ftype == "bias":
            return good_exp & (fitstbl["idname"] == "BIAS")
        if ftype in ["pixelflat", "trace", "illumflat"]:
            # Lamp-on dome flats
            return good_exp & (fitstbl["idname"] == "DOMEFLAT") \
                & (fitstbl["lampstat01"] == "on")
        if ftype == "lampoffflats":
            return good_exp & (fitstbl["idname"] == "DOMEFLAT_OFF")
        # NOTE: Mask images (idname = MASKIMAGE) are deliberately not typed.
        # NOTE: ThAr arcs (idname = TH-AR) are deliberately not typed as
        # 'arc' either: the OH lines in the science frames are the default
        # wavelength calibration, and PypeIt would combine all 'arc' frames
        # in a calibration group.  The arcs share the science setup, and
        # pypeit_setup writes their (untyped) rows commented out.  To use
        # them, uncomment the rows, set their frametype to 'arc', remove
        # 'arc' from the science frames (keep 'tilt' there: the OH lines
        # trace the tilts better than the sparse K-band ThAr lines), and
        # set the wavelength-calibration lamps accordingly.

        log.debug(f"Cannot determine if frames are of type {ftype}.")
        return np.zeros(len(fitstbl), dtype=bool)

    def valid_configuration_values(self):
        """
        Return a fixed set of valid values for any/all of the configuration
        keys.

        For MOIRCS, this is used to remove the chip-2 files from the
        metadata table built by ``pypeit_setup``.  Each chip-2 file is read
        through its chip-1 companion; see :func:`get_rawimage`.

        Returns:
            :obj:`dict`: A dictionary with any/all of the configuration keys
            and their associated discrete set of valid values.
        """
        return {"detector": ["1"]}

    @staticmethod
    def _chip1_header(raw_file):
        """
        Read the primary header of a file and check it is a chip-1 file.

        Generated by JXP and Claude.

        Args:
            raw_file (:obj:`str`, `Path`_):
                File to check.

        Returns:
            `astropy.io.fits.Header`_: The primary header.

        Raises:
            :class:`~pypeit.PypeItError`: Raised if ``DET-ID`` is not 1.
        """
        hdr = fits.getheader(raw_file, 0)
        if int(hdr.get("DET-ID", -1)) != 1:
            raise PypeItError(
                f"{Path(raw_file).name} is not a MOIRCS chip-1 file "
                f"(DET-ID = {hdr.get('DET-ID')}).  List only the chip-1 "
                "files in the PypeIt file; chip 2 is read automatically.")
        return hdr

    @staticmethod
    def companion_file(raw_file):
        """
        Find and validate the chip-2 file associated with a chip-1 file.

        The chip-2 file has the next frame number (e.g.,
        ``MCSP00237323.fits`` -> ``MCSP00237324.fits``) in the same
        directory, with the same ``EXP-ID`` and ``DET-ID = 2``.

        Generated by JXP and Claude.

        Args:
            raw_file (:obj:`str`, `Path`_):
                Path to the chip-1 file.

        Returns:
            `Path`_: Path to the chip-2 file.

        Raises:
            :class:`~pypeit.PypeItError`: Raised if ``raw_file`` is not a
            chip-1 file, or if the chip-2 file is missing or does not match.
        """
        raw_file = Path(raw_file)
        hdr1 = SubaruMOIRCSSpectrograph._chip1_header(raw_file)

        # Split the name into prefix, frame number, and extension(s)
        root, ext = raw_file.name.split(".", 1)
        match = re.fullmatch(r"(\D*)(\d+)", root)
        if match is None:
            raise PypeItError(
                f"Cannot parse a frame number from {raw_file.name}.")
        prefix, number = match.groups()
        name2 = f"{prefix}{int(number) + 1:0{len(number)}d}.{ext}"
        file2 = raw_file.with_name(name2)
        if not file2.exists():
            raise PypeItError(
                f"Missing the chip-2 file {name2} for {raw_file.name}; it "
                f"must be in the same directory ({raw_file.parent}).")

        # Validate the pairing
        hdr2 = fits.getheader(file2, 0)
        if int(hdr2.get("DET-ID", -1)) != 2 \
                or hdr2.get("EXP-ID") != hdr1.get("EXP-ID"):
            raise PypeItError(
                f"{name2} is not the chip-2 companion of {raw_file.name}: "
                f"DET-ID = {hdr2.get('DET-ID')}, EXP-ID = "
                f"{hdr2.get('EXP-ID')} (expected 2 and "
                f"{hdr1.get('EXP-ID')}).")
        return file2

    def get_rawimage(self, raw_file, det, **kwargs):
        """
        Read a raw MOIRCS image and return the data and relevant metadata.

        ``raw_file`` is always the chip-1 file.  For ``det = 2`` the chip-2
        file is located with :func:`companion_file` and read instead.
        Everything else follows
        :func:`~pypeit.spectrographs.spectrograph.Spectrograph.get_rawimage`.

        Generated by JXP and Claude.

        Args:
            raw_file (:obj:`str`, `Path`_):
                The chip-1 file of the exposure.
            det (:obj:`int`):
                1-indexed detector to read (1 or 2).
            **kwargs:
                Passed to the base-class method.

        Returns:
            tuple: See
            :func:`~pypeit.spectrographs.spectrograph.Spectrograph.get_rawimage`.
            The returned ``hdu`` is the file actually read.
        """
        # The listed file must be a chip-1 file; for det=2, read its
        # companion instead
        if det == 2:
            _raw_file = self.companion_file(raw_file)
        else:
            self._chip1_header(raw_file)
            _raw_file = raw_file
        return super().get_rawimage(_raw_file, det, **kwargs)

    def get_detector_par(self, det, hdu=None):
        """
        Return metadata for the selected detector.

        The chip is selected by ``det``.  When called from
        :func:`get_rawimage`, ``hdu`` is the file actually read, i.e. the
        chip-2 companion file for ``det = 2``, so frame-dependent values
        (binning, read noise) come from the right chip.

        The read noise depends on the number of Fowler samples
        (``DET-NSMP``; e.g. 10 for science frames, 1 for flats and
        standards): ``17.5/sqrt(DET-NSMP)`` e-.  ``DET-NSMP = 10`` is
        assumed when ``hdu`` is None or the card is missing.

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
        # Binning
        binning = "1,1" if hdu is None \
            else self.get_meta_value(self.get_headarr(hdu), "binning")

        # Read noise for the number of Fowler samples of this frame.
        # TODO - The single-read noise (17.5 e-) is to be confirmed by the
        # instrument scientist.
        nsmp = 10 if hdu is None else hdu[0].header.get("DET-NSMP", 10)
        ronoise = 17.5 / np.sqrt(max(int(nsmp), 1))

        # TODO - Confirm gain, read noise, dark current and saturation with
        # the instrument scientist

        # CHIP1
        detector_dict1 = dict(
            binning=binning,
            det=1,
            dataext=0,
            specaxis=1,
            specflip=True,
            spatflip=False,
            platescale=0.116,
            darkcurr=18.0,  # e-/pixel/hour (<0.005 e-/pixel/sec)
            saturation=3.3e4,  # from MOIRCS website
            nonlinear=1.0,
            mincounts=-1e10,
            numamplifiers=1,
            gain=np.atleast_1d(2.07),  # e-/ADU for chip 1
            ronoise=np.atleast_1d(ronoise),  # 17.5/sqrt(DET-NSMP)
            # Effective area (EFP-MIN/EFP-RNG); excludes reference pixels
            datasec=np.atleast_1d("[5:2044,5:2044]"),
        )

        # CHIP2
        detector_dict2 = dict(
            binning=binning,
            det=2,
            dataext=0,
            specaxis=1,
            specflip=False,
            spatflip=False,
            platescale=0.116,
            darkcurr=18.0,  # e-/pixel/hour (<0.005 e-/pixel/sec)
            saturation=3.3e4,  # from MOIRCS website
            nonlinear=1.0,
            mincounts=-1e10,
            numamplifiers=1,
            gain=np.atleast_1d(1.99),  # e-/ADU for chip 2
            ronoise=np.atleast_1d(ronoise),  # 17.5/sqrt(DET-NSMP)
            # Effective area (EFP-MIN/EFP-RNG); excludes reference pixels
            datasec=np.atleast_1d("[5:2044,5:2044]"),
        )
        # Finish
        if det == 1:
            return detector_container.DetectorContainer(**detector_dict1)
        if det == 2:
            return detector_container.DetectorContainer(**detector_dict2)
        raise PypeItError(f"Unknown MOIRCS detector: {det=}!")

    def bpm(self, filename, det, shape=None, msbias=None):
        """
        Generate the default bad-pixel mask.

        The static masks are derived from the NAOJ MOIRCS bad-pixel masks
        (``mcsbadpix_oct2016``), with the imaging beam-splitter shadow
        removed.  They are stored in the raw orientation and are trimmed and
        re-oriented here exactly as the raw images are.

        Generated by JXP and Claude.

        Args:
            filename (:obj:`str` or None):
                An example file to use to get the image shape.
            det (:obj:`int`):
                1-indexed detector number.
            shape (:obj:`tuple`, optional):
                Processed image shape.  Required if ``filename`` is None;
                ignored otherwise.
            msbias (`numpy.ndarray`_, optional):
                Processed bias frame used to identify bad pixels.

        Returns:
            `numpy.ndarray`_: An integer array with a masked value set to 1
            and an unmasked value set to 0.
        """
        # Empty BPM with the processed shape (and bias-based pixels, if any)
        bpm_img = super().bpm(filename, det, shape=shape, msbias=msbias)

        # Static mask in the raw orientation
        bpm_file = dataPaths.static_calibs.get_file_path(
            f'subaru_moircs/bpm_moircs_det{det}.fits.gz')
        static_raw = fits.getdata(bpm_file).astype(bool)

        # Trim to the data section and re-orient, as for the raw images
        detpar = self.get_detector_par(det)
        datasec = parse.sec2slice(detpar['datasec'][0], one_indexed=True,
                                  include_end=True, require_dim=2)
        static = self.orient_image(detpar, static_raw[datasec])
        if static.shape != bpm_img.shape:
            # E.g. binned data; the static mask is for unbinned frames
            log.warning(f'Static BPM shape {static.shape} does not match '
                        f'the image shape {bpm_img.shape}; not applied.')
            return bpm_img
        bpm_img[static] = 1
        return bpm_img

    def config_specific_par(self, inp, inp_par=None):
        """
        Modify the PypeIt parameters to hard-wired values used for
        specific instrument configurations.

        Args:
            inp (:obj:`str`, :obj:`list`, `Path`_, `astropy.io.fits.Header`_, `astropy.table.Table`_):
                Input filename, an `astropy.io.fits.Header`_ object, or a
                list of `astropy.io.fits.Header`_ objects.  Or a row from
                the metadata table.
            inp_par (:class:`~pypeit.par.parset.ParSet`, optional):
                Parameter set used for the full run of PypeIt.  If None,
                use :func:`default_pypeit_par`.

        Returns:
            :class:`~pypeit.par.parset.ParSet`: The PypeIt parameter set
            adjusted for configuration specific parameter values.
        """
        # Start with instrument wide
        par = super().config_specific_par(inp, inp_par=inp_par)

        # Grism-specific wavelength calibration.  Other grisms fall back to
        # holy-grail on the OH lines until templates are built for them.
        dispname = self.get_meta_value(inp, 'dispname')
        if dispname == 'HK500':
            # Archive of 19 holy-grail OH solutions (both detectors) from
            # the dev-suite HK500 mask.  Reidentifying against many slits
            # copes with the different spectral coverage of each slit
            # better than a single full_template (tested; see the dev
            # suite pypeitdev/subaru_moircs logs).
            par['calibrations']['wavelengths']['method'] = 'reidentify'
            par['calibrations']['wavelengths']['reid_arxiv'] \
                = 'subaru_moircs_HK500.fits'
        elif dispname == 'VB_K':
            # Holy-grail on the OH lines solves every science slit of the
            # VB_K test mask.  OH_MOSFIRE_K gives the same solutions as
            # OH_NIRES (median |dlambda| < 0.15 A) with a lower rms.  The
            # measured line FWHM is ~6 px, set by the slit width (~0.7
            # arcsec, R ~ 1900); R ~ 2600 (4.3 px) is for narrower slits.
            # See the dev suite pypeitdev/subaru_moircs_vbk logs.
            par['calibrations']['wavelengths']['lamps'] = ['OH_MOSFIRE_K']
            par['calibrations']['wavelengths']['fwhm'] = 6.0

        return par

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
        return ["dispname", "decker", "binning"]

    def raw_header_cards(self):
        """
        Return additional raw header cards to be propagated in
        downstream output files for configuration identification.

        Returns:
            :obj:`list`: List of keywords from the raw data files that should
            be propagated in output files.
        """
        return ["DISPERSR", "SLIT", "BIN-FCT1", "BIN-FCT2"]

    def pypeit_file_keys(self):
        """
        Define the list of keys to be output into a standard PypeIt file.

        Returns:
            :obj:`list`: The list of keywords in the relevant
            :class:`~pypeit.metadata.PypeItMetaData` instance to print to the
            :ref:`pypeit_file`.
        """
        return super().pypeit_file_keys() \
            + ["lampstat01", "dithpat", "dithpos", "dithoff", "frameno"]

    def get_comb_group(self, fitstbl):
        """
        Automatically assign combination groups and background images by
        parsing the dither positions.

        Within each setup, every science (or standard) frame at dither
        position A is paired with the B frame closest in time, and vice
        versa, following the fall-back logic in
        :func:`~pypeit.spectrographs.keck_mosfire.KeckMOSFIRESpectrograph.get_comb_group`.
        Each frame keeps its own ``comb_id``; only ``bkg_id`` is set.  This
        covers AB, BA, ABBA and longer AB sequences.  Frames at other
        positions are left untouched.

        Generated by JXP and Claude.

        Args:
            fitstbl (`astropy.table.Table`_):
                The table with the metadata for all the frames.

        Returns:
            `astropy.table.Table`_: The modified table.
        """
        for ftype in ["science", "standard"]:
            is_type = np.array([ftype in _ft for _ft in fitstbl["frametype"]])
            for setup in np.unique(fitstbl["setup"][is_type]):
                in_cfg = is_type & np.array(
                    [setup in _s for _s in fitstbl["setup"]])
                for dpat in np.unique(fitstbl["dithpat"][in_cfg]):
                    if dpat == "none":
                        continue
                    in_pat = in_cfg & (fitstbl["dithpat"] == dpat)
                    is_a = np.where(in_pat & (fitstbl["dithpos"] == "A"))[0]
                    is_b = np.where(in_pat & (fitstbl["dithpos"] == "B"))[0]
                    if is_a.size == 0 or is_b.size == 0:
                        continue
                    # Pair each frame with the closest (in time) frame at
                    # the other position
                    for this, other in [(is_a, is_b), (is_b, is_a)]:
                        for i in this:
                            dt = np.absolute(fitstbl["mjd"][other]
                                             - fitstbl["mjd"][i])
                            j = other[np.argmin(dt)]
                            fitstbl["bkg_id"][i] = fitstbl["comb_id"][j]
        return fitstbl
