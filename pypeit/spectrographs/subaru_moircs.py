"""
Module for Subaru MOIRCS

.. include:: ../include/links.rst
"""
import numpy as np
from astropy import units
from astropy.coordinates import SkyCoord
from astropy.io import fits
from IPython import embed

from pypeit import msgs, telescopes
from pypeit.core import framematch, meta, parse
from pypeit.images import detector_container
from pypeit.spectrographs import spectrograph


class SubaruMOIRCSSpectrograph(spectrograph.Spectrograph):
    """
    Child of Spectrograph to handle Subaru/MOIRCS specific code
    """

    ndet = 1  # Because each detector is written to a separate FITS file
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

        # NOTE: These parameters are copied from Keck/MOSFIRE
        # TODO: review each of these parameters and update for MOIRCS

        # Wavelengths
        # 1D wavelength solution
        # 0.20  # Might be grating dependent..
        par["calibrations"]["wavelengths"][
            "rms_thresh_frac_fwhm"
        ] = 0.11  # might be grism dependent
        par["calibrations"]["wavelengths"]["sigdetect"] = 5.0
        par["calibrations"]["wavelengths"]["fwhm"] = 5.0
        par["calibrations"]["wavelengths"]["n_final"] = 4
        par["calibrations"]["wavelengths"]["lamps"] = ["OH_NIRES"]
        par["calibrations"]["wavelengths"]["method"] = "holy-grail"
        # Reidentification parameters
        # par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_nires.fits'
        par["calibrations"]["slitedges"]["edge_thresh"] = 50.0
        par["calibrations"]["slitedges"]["sync_predict"] = "nearest"

        # Flats
        turn_off = dict(use_biasimage=False, use_overscan=False, use_darkimage=False)
        par.reset_all_processimages_par(**turn_off)

        # Extraction
        par["reduce"]["skysub"]["bspline_spacing"] = 0.8
        par["reduce"]["extraction"]["sn_gauss"] = 4.0

        # Flexure
        par["flexure"]["spec_method"] = "boxcar"

        # Adjustments to slit and tilts for NIR
        par["calibrations"]["slitedges"]["edge_thresh"] = 50.0
        par["calibrations"]["slitedges"]["fit_order"] = 3
        par["calibrations"]["slitedges"]["max_shift_adj"] = 0.5

        # lris alignment boxes are typically 4 arcsec
        par['calibrations']['slitedges']['minimum_slit_length_sci'] = 5.
        # Remove slits that are too short
        par['calibrations']['slitedges']['minimum_slit_length'] = 3.

        # Tilt parameters
        par["calibrations"]["tilts"]["tracethresh"] = 25.0
        par["calibrations"]["tilts"]["spat_order"] = 3
        par["calibrations"]["tilts"]["spec_order"] = 4

        par["scienceframe"]["process"]["sigclip"] = 20.0
        par["scienceframe"]["process"]["satpix"] = "nothing"

        # Set the default exposure time ranges for the frame typing
        par["calibrations"]["standardframe"]["exprng"] = [None, 20]
        par["calibrations"]["arcframe"]["exprng"] = [1, None]
        par["calibrations"]["darkframe"]["exprng"] = [1, None]
        par["scienceframe"]["exprng"] = [20, None]

        # Sensitivity function parameters
        par["sensfunc"][
            "extrap_blu"
        ] = 0.0  # Y-band contaminated by higher order so don't extrap much
        par["sensfunc"]["extrap_red"] = 0.0
        par["fluxcalib"]["extrap_sens"] = True
        par["sensfunc"]["extrap_red"] = 0.0
        par["sensfunc"]["algorithm"] = "IR"
        par["sensfunc"]["polyorder"] = 13
        par["sensfunc"]["IR"]["maxiter"] = 2
        par["sensfunc"]["IR"]["telgridfile"] = "TelFit_MaunaKea_3100_26100_R20000.fits"

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
        self.meta["binning"] = dict(ext=0, card=None, default="1,1")

        self.meta["mjd"] = dict(ext=0, card="MJD")
        self.meta["exptime"] = dict(ext=0, card="EXPTIME")
        self.meta["airmass"] = dict(ext=0, card="AIRMASS")
        #
        self.meta["decker"] = dict(ext=0, card="SLIT")

        # Extras for config and frametyping
        self.meta["dispname"] = dict(
            ext=0, card="DISPERSR", required_ftypes=["science", "standard"]
        )
        # mosfire doe not have this and it's not needed
        # # TODO - FIX THIS!!
        # self.meta["dispangle"] = dict(
        #     ext=0, card="BZERO", rtol=2.0
        # )  # , required_ftypes=['science', 'standard'])
        self.meta["idname"] = dict(ext=0, card="DATA-TYP")
        self.meta["detector"] = dict(ext=0, card="DET-ID")
        self.meta["instrument"] = dict(ext=0, card="INSTRUME")

        # Dithering
        self.meta['dithpat'] = dict(ext=0, card='K_DITPAT')
        self.meta['dithpos'] = dict(card=None, compound=True)
        self.meta['dithoff'] = dict(card=None, compound=True)

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
        # if meta_key == "binning":
        #     binspatial = headarr[0]["BIN-FCT1"]  # X
        #     binspec = headarr[0]["BIN-FCT2"]  # Y
        #     # TODO -- CHECK THE FOLLOWING
        #     binning = parse.binning2string(binspec, binspatial)
        #     return binning
        if meta_key == "dithpos":
            pos_key = headarr[0].get("K_DITCNT")
            if pos_key is not None and pos_key in [1.,2.]:
                return 'A' if pos_key == 1. else 'B'
            else:
                return pos_key

        if meta_key == "dithoff":
            ditwit = headarr[0].get("K_DITWID")
            ditpos = self.get_meta_value(headarr[0], 'dithpos')
            if ditpos in ['A', 'B']:
                return ditwit/2 if ditpos == 'A' else -ditwit/2
            else:
                return ditwit/2
        else:
            msgs.error("Not ready for this compound meta")


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
        if ftype == "science":
            return good_exp & (fitstbl["idname"] == "OBJECT")
        if ftype == "standard":
            return good_exp & (fitstbl["idname"] == "OBJECT")
        if ftype == "bias":
            return good_exp & (fitstbl["idname"] == "BIAS")
        if ftype == "lampoffflats":
            return good_exp & (fitstbl["idname"] == "DOMEFLAT") & (fitstbl["target"] == "DOMEFLAT_OFF")
        if ftype in ["pixelflat", "trace", "illumflat"]:
            # Flats and trace frames are typed together
            # TODO -- Are there internal flats?
            return good_exp & (fitstbl["idname"] == "DOMEFLAT") & (fitstbl["target"] == "DOMEFLAT")
        if ftype in ["arc", "tilt"]:
            # can be an arc image
            is_arc = (fitstbl["idname"] == "COMPARISON")
            # or a science image
            is_obj = (fitstbl["idname"] == "OBJECT")
            return good_exp & (is_arc | is_obj)

        msgs.warn("Cannot determine if frames are of type {0}.".format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

    def get_detector_par(self, det, hdu=None):
        """
        Return metadata for the selected detector.

        Args:
            det (:obj:`int`):
                1-indexed detector number.  MORICS writes each of the two detectors
                to separate files.  Only the primary HDU will be used.
            hdu (`astropy.io.fits.HDUList`_, optional):
                The open fits file with the raw image of interest.  If not
                provided, frame-dependent parameters are set to a default.

        Returns:
            :class:`~pypeit.images.detector_container.DetectorContainer`:
            Object with the detector metadata.
        """

        if hdu is None:
            chip = "1" if det == 1 else "2"
        else:
            # Binning
            # TODO: Could this be detector dependent??
            binning = self.get_meta_value(self.get_headarr(hdu), "binning")
            chip = self.get_meta_value(self.get_headarr(hdu), "detector")

        # TODO - UPDATE ALL OF THIS

        # CHIP1
        detector_dict1 = dict(
            binning="1,1",
            det=1,  # because each detector is written to a separate FITS file
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
            ronoise=np.atleast_1d(5.534),  # for NDR=10 (17.5/sqrt(10))
            datasec=np.atleast_1d("[5:2044,5:2044]"),  # copied from the MOSFIRE config
        )

        # CHIP2
        detector_dict2 = dict(
            binning="1,1",
            det=1,  # because each detector is written to a separate FITS file
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
            ronoise=np.atleast_1d(5.534),  # for NDR=10 (17.5/sqrt(10))
            datasec=np.atleast_1d("[5:2044,5:2044]"),  # copied from the MOSFIRE config
        )
        # Finish
        if chip == "1":
            return detector_container.DetectorContainer(**detector_dict1)
        elif chip == "2":
            return detector_container.DetectorContainer(**detector_dict2)
        else:
            msgs.error(f"Unknown detector chip: {chip=}!")

    def config_specific_par(self, scifile, inp_par=None):
        """
        Modify the PypeIt parameters to hard-wired values used for
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

        return par

    def config_independent_frames(self):
        """
        Define frame types that are independent of the fully defined
        instrument configuration.

        This method returns a dictionary where the keys of the dictionary are
        the list of configuration-independent frame types. The value of each
        dictionary element can be set to one or more metadata keys that can
        be used to assign each frame type to a given configuration group. See
        :func:`~pypeit.metadata.PypeItMetaData.set_configurations` and how it
        interprets the dictionary values, which can be None.

        Returns:
            :obj:`dict`: Dictionary where the keys are the frame types that
            are configuration-independent and the values are the metadata
            keywords that can be used to assign the frames to a configuration
            group.
        """
        return {"bias": "detector", "dark": "detector"}

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
        # return ['dispname', 'dispangle', 'decker', 'detector']
        # TODO -- Consider dispangle
        return ["dispname", "decker", "detector"]

    def pypeit_file_keys(self):
        """
        Define the list of keys to be output into a standard PypeIt file.

        Returns:
            :obj:`list`: The list of keywords in the relevant
            :class:`~pypeit.metadata.PypeItMetaData` instance to print to the
            :ref:`pypeit_file`.
        """

        pypeit_keys = super().pypeit_file_keys()
        return pypeit_keys + ['dithpat', 'dithpos', 'dithoff']

    # TODO -- Convert this into get_comb_group()
    def parse_dither_pattern(self, file_list, ext=None):
        """
        Parse headers from a file list to determine the dither pattern.

        Parameters
        ----------
        file_list (list of strings):
            List of files for which dither pattern is desired
        ext (int, optional):
            Extension containing the relevant header for these files. Default=None. If None, code uses
            self.primary_hdrext

        Returns
        -------
        dither_pattern, dither_id, offset_arcsec

        dither_pattern (str `numpy.ndarray`_):
            Array of dither pattern names
        dither_id (str `numpy.ndarray`_):
            Array of dither pattern IDs
        offset_arc (float `numpy.ndarray`_):
            Array of dither pattern offsets
        """
        nfiles = len(file_list)
        offset_arcsec = np.zeros(nfiles)
        dither_pattern = []
        dither_id = []
        for ifile, file in enumerate(file_list):
            hdr = fits.getheader(file, self.primary_hdrext if ext is None else ext)
            dither_pattern.append(self.get_meta_value(hdr, 'dithpat'))
            dither_id.append(self.get_meta_value(hdr, 'dithpos'))
            offset_arcsec[ifile] = self.get_meta_value(hdr, 'dithoff')
        return np.array(dither_pattern), np.array(dither_id), np.array(offset_arcsec)


