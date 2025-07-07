"""
Module for LDT/RIMAS specific methods.

The Rapid infrared IMAger Spectrometer (RIMAS) was built at the NASA Goddard
Space Flight Center in partnership with the Astronomy Department of the
University of Maryland at College Park.  RIMAS was designed for Target of
Opportunity (ToO) observations, and thus once commissioned RIMAS is expected to
remain on Port A of the LDT Instrument Cube for some time to come in order to
support ToO observations.  RIMAS arrived at LDT in May 2025, and is undergoing
commissioning.

RIMAS is a near-infrared imager spectrometer (R≈25, R≈250, and R≈4000) designed
to observe the afterglows of high-redshift GRBs in the Y, J, H and K bands.  It
is a fully cryogenic instrument operating along two optical arms: one covering
H and K bands (HK) and the other covering Y and J bands (YJ). RIMAS is designed
for photometry, low resolution spectroscopy (R ≈ 25 and R ≈ 250), and moderate
resolution spectroscopy (R ≈ 4000). It accomplishes this using a dichroic
mirror, H, K, Y and J broadband filters, two VPHs, ZnSe grisms with cross-
dispersers, and H4RG-10 detectors.  RIMAS is a flexible tool with both imaging
and spectral modes available.

The design for the science optics is broken into three lens assemblies, one
collimator and a camera for each of two optical arms.  Each assembly is
composed of five elements to minimize image aberrations.  NIR light from LDT is
directed to RIMAS by positioning a dichroic in the instrument cube (Figure 1).
The beam is then passed through RIMAS's dewar window, after which all optics
are cooled to ~70 K.  Once the beam has been collimated a dichroic mirror is
used to divide the wavelength coverage into two optical arms, “YJ” (0.9 -
1.4 μm) and “HK” (1.4 - 2.4 μm).  Before being refocused by the cameras, the
beams are either filtered or dispersed by additional optics located on wheels.  
A Teledyne H4RG-10 HgCdTe detector (4096 x 4096 pixels) is positioned at each
arm's focal plane.

.. include:: ../include/links.rst
"""

import astropy.table
import astropy.time
import numpy as np

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.core import parse
from pypeit.images import detector_container
from pypeit.spectrographs import spectrograph


class LDTRIMASSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle LDT/RIMAS specific code

    This class contains the common methods for all 4 operating modes of RIMAS
    data reduction: 2 arms x [longslit, echelle].
    """

    telescope = telescopes.LDTTelescopePar()
    url = "https://lowell.edu/research/telescopes-and-facilities/ldt/rimas/"
    header_name = "RIMAS"
    supported = True

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}

        # Required (core)
        self.meta["ra"] = dict(ext=0, card="RA")
        self.meta["dec"] = dict(ext=0, card="DEC")
        self.meta["target"] = dict(ext=0, card="OBJNAME")
        self.meta["dispname"] = dict(card=None, compound=True)
        self.meta["decker"] = dict(ext=0, card="FILTER4")  # SLIT filter wheel
        self.meta["binning"] = dict(card=None, compound=True)
        self.meta["mjd"] = dict(card=None, compound=True)
        self.meta["airmass"] = dict(ext=0, card="SECZ")
        self.meta["exptime"] = dict(ext=0, card="EXPTIMEE")
        self.meta["instrument"] = dict(ext=0, card="INSTRUME")

        # Extras for config and frametyping
        # NOTE: `rtol` is _relative_ tolerance (e.g. 1 part in 1,000)
        self.meta["arm"] = dict(card=None, compound=True)
        self.meta["idname"] = dict(ext=0, card="OBJTYPE")
        self.meta["filter1"] = dict(ext=0, card="FILTER3")  # AUX filter wheel
        self.meta["lampstat01"] = dict(card=None, compound=True)
        self.meta["slitwid"] = dict(card=None, compound=True)

    def compound_meta(self, headarr: list, meta_key: str) -> object:
        """
        Methods to generate metadata requiring interpretation of the header
        data, instead of simply reading the value of a header card.

        Args:
            headarr (:obj:`list`):
                List of `astropy.io.fits.Header`_ objects.
            meta_key (:obj:`str`):
                Metadata keyword to construct.

        Returns:
            :obj:`object`: Metadata value read from the header(s).
        """
        if meta_key == "binning":
            # Binning in RIMAS headers given as separate values
            return parse.binning2string(headarr[0]["BINX"], headarr[0]["BINY"])

        if meta_key == "mjd":
            # Use custom scrubber + AstroPy to convert 'DATE-OBS' into a mjd.
            ttime = self.scrub_isot_dateobs(headarr[0]["DATE-BEG"])
            return ttime.mjd

        if meta_key == "lampstat01":
            # The spectral comparison lamps turned on are listed in `LAMPCAL`, but
            #  if no lamps are on, then this string is blank.  Return either the
            #  populated `LAMPCAL` string, or 'off' to ensure a positive entry for
            #  `lampstat01`.
            lampcal = ""
            return "off" if lampcal == "" else lampcal

        if meta_key == "arm":
            # Return the ARM name from the filename
            # TODO: Have the RIMAS folks add a FITS keyword indicating the arm
            return "YJ" if "YJ" in headarr[0]["ASDFNAME"] else "HK"

        if meta_key == "dispname":
            # Return FILTER1 (Filter Wheel Name FW YJ) for YJ frames and
            #   FILTER2 (Filter Wheel Name FW HK) for HK frames
            # TODO: Have the RIMAS folks add a FITS keyword indicating the arm
            return (
                headarr[0]["FILTER1"]
                if "YJ" in headarr[0]["ASDFNAME"]
                else headarr[0]["FILTER2"]
            )

        if meta_key == "slitwid":
            # Convert the decker into a slitwidth in arcseconds
            match headarr[0]["FILTER4"].strip():
                case "250 um":
                    return 2.5
                case "130 um":
                    return 1.2
                case "80 um":
                    return 0.6
                case _:
                    return 0

        if meta_key == "target":
            # Revert to TCS's SCITARG if target not set in LOUI for OBJECT frames
            return (
                headarr[0]["SCITARG"].strip()
                if (
                    headarr[0]["IMAGETYP"].strip() == "OBJECT"
                    and headarr[0]["OBJNAME"].strip() in ["UNKNOWN", ""]
                )
                else headarr[0]["OBJNAME"].strip()
            )

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
        return ["arm", "dispname", "filter1"]

    def raw_header_cards(self):
        """
        Return additional raw header cards to be propagated in
        downstream output files for configuration identification.

        The list of raw data FITS keywords should be those used to populate
        the :meth:`~pypeit.spectrographs.spectrograph.Spectrograph.configuration_keys`
        or are used in :meth:`~pypeit.spectrographs.spectrograph.Spectrograph.config_specific_par`
        for a particular spectrograph, if different from the name of the
        PypeIt metadata keyword.

        This list is used by
        :meth:`~pypeit.spectrographs.spectrograph.Spectrograph.subheader_for_spec`
        to include additional FITS keywords in downstream output files.

        Returns:
            :obj:`list`: List of keywords from the raw data files that should
            be propagated in output files.
        """
        return ["ASDFNAME", "FILTER1", "FILTER2", "FILTER3", "FILTER4"]

    def pypeit_file_keys(self):
        """
        Define the list of keys to be output into a standard PypeIt file.

        Returns:
            :obj:`list` : The list of keywords in the relevant
            :class:`~pypeit.metadata.PypeItMetaData` instance to print to the
            :ref:`pypeit_file`.
        """
        return super().pypeit_file_keys() + ["slitwid", "lampstat01", "dither"]

    def check_frame_type(self, ftype: str, fitstbl: astropy.table.Table, exprng=None):
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
                :func:`~pypeit.core.framematch.check_frame_exptime`.

        Returns:
            `numpy.ndarray`_: Boolean array with the flags selecting the
            exposures in ``fitstbl`` that are ``ftype`` type frames.
        """
        good_exp = framematch.check_frame_exptime(fitstbl["exptime"], exprng)
        if ftype == "bias":
            return fitstbl["idname"] == "BIAS"
        if ftype in ["arc", "tilt"]:
            # FOCUS frames should have frametype None, BIAS is bias regardless of lamp status
            return (
                good_exp
                & (fitstbl["lampstat01"] != "off")
                & (fitstbl["idname"] != "FOCUS")
                & (fitstbl["idname"] != "BIAS")
            )
        if ftype in ["trace", "pixelflat"]:
            return (
                good_exp
                & (fitstbl["idname"] == "DOME FLAT")
                & (fitstbl["lampstat01"] == "off")
            )
        if ftype == "illumflat":
            return (
                good_exp
                & (fitstbl["idname"] == "SKY FLAT")
                & (fitstbl["lampstat01"] == "off")
            )
        if ftype == "science":
            return (
                good_exp
                & (fitstbl["idname"] == "OBJECT")
                & (fitstbl["lampstat01"] == "off")
            )
        if ftype == "standard":
            return (
                good_exp
                & (fitstbl["idname"] == "STANDARD")
                & (fitstbl["lampstat01"] == "off")
            )
        if ftype == "dark":
            return (
                good_exp
                & (fitstbl["idname"] == "DARK")
                & (fitstbl["lampstat01"] == "off")
            )
        if ftype in ["pinhole", "align", "sky", "lampoffflats", "scattlight"]:
            # DeVeny doesn't have any of these types of frames
            return np.zeros(len(fitstbl), dtype=bool)
        msgs.warn(f"Cannot determine if frames are of type {ftype}")
        return np.zeros(len(fitstbl), dtype=bool)

    @staticmethod
    def scrub_isot_dateobs(dt_str: str):
        """Scrub the input ``DATE-OBS`` for ingestion by AstroPy Time

        The main issue this method addresses is that sometimes the LOIS
        software at LDT has roundoff abnormalities in the time string written
        to the ``DATE-OBS`` header keyword.  For example, in one header
        ``2020-01-30T13:17:010.0`` was written, where the seconds has 3 digits
        -- presumably the seconds field was constructed with a leading zero
        because ``sec`` was < 10, but when rounded for printing
        yielded "10.00", producing a complete seconds field of ``010.00``.

        This abnormality, along with a seconds field equaling ``60.00``, causes
        AstroPy's Time parser to freak out with a ``ValueError``.  This
        method attempts to return the `astropy.time.Time`_ object directly, but
        then scrubs any values that cause a ``ValueError``.

        The scrubbing consists of deconstructing the string into its components,
        then carefully reconstructing it into proper ISO 8601 format.  Also,
        some recursive edge-case catching is done, but at some point you just
        have to give up and go buy a lottery ticket.

        If you have a truly bizarre ``DATE-OBS`` string, simply edit that keyword
        in the FITS header and then re-run PypeIt.

        Parameters
        ----------
        dt_str : :obj:`str`
            Input datetime string from the ``DATE-OBS`` header keyword

        Returns
        -------
        `astropy.time.Time`_
            The AstroPy Time object corresponding to the ``DATE-OBS`` input string
        """
        # Clean all leading / trailing whitespace
        dt_str = dt_str.strip()

        # Attempt to directly return the AstroPy Time object
        try:
            return astropy.time.Time(dt_str, format="isot")
        except ValueError:
            # Split out all pieces of the datetime, and recompile
            date, time = dt_str.split("T")
            yea, mon, day = date.split("-")
            hou, mnt, sec = time.split(":")
            # Check if the seconds is exactly equal to 60... increment minute
            if sec == "60.00":
                sec = "00.00"
                if mnt != "59":
                    mnt = int(mnt) + 1
                else:
                    mnt = "00"
                    if hou != "23":
                        hou = int(hou) + 1
                    else:
                        hou = "00"
                        # If the edge cases go past here, go buy a lottery ticket!
                        day = int(day) + 1
            # Reconstitute the DATE-OBS string, and return the Time() object
            date = f"{int(yea):04d}-{int(mon):02d}-{int(day):02d}"
            time = f"{int(hou):02d}:{int(mnt):02d}:{float(sec):09.6f}"
            return astropy.time.Time(f"{date}T{time}", format="isot")


class LDTRIMASYJSpectrograph(LDTRIMASSpectrograph):
    """
    Child to handle common aspects of the LDT/RIMAS YJ Arm
    """

    ndet = 1
    camera = "RIMAS_YJ"

    def get_detector_par(self, _, hdu=None):
        """
        Return metadata for the LDT/RIMAS YJ detector.

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
            binning = "1,1"  # Most common use mode
            gain = np.atleast_1d(1.52)  # Hardcoded in the header
            ronoise = np.atleast_1d(4.9)  # Hardcoded in the header
            datasec = np.atleast_1d("[5:512,53:2095]")  # For 1x1 binning
            oscansec = np.atleast_1d("[5:512,5:48]")  # For 1x1 binning
        else:
            binning = self.get_meta_value(self.get_headarr(hdu), "binning")
            gain = np.atleast_1d(hdu[0].header["GAIN"])
            ronoise = np.atleast_1d(hdu[0].header["RDNOISE"])
            datasec = hdu[0].header["TRIMSEC"]
            oscansec = hdu[0].header["BIASSEC"]

        # Detector
        detector_dict = dict(
            binning=binning,
            det=1,  # The YJ channel has but one detector
            dataext=0,  # Image is in extension 0
            specaxis=1,  # Native spectrum is along the x-axis
            specflip=True,  # RIMAS IR FPAs have blue at the right
            spatflip=False,
            platescale=0.19,  # Arcsec / pixel
            darkcurr=0,  # e-/pixel/hour
            saturation=0,  # 16-bit ADC
            nonlinear=0,  # Linear to ~97% of saturation
            mincounts=-1e10,
            numamplifiers=1,
            gain=gain,  # See above
            ronoise=ronoise,  # See above
            # Data & Overscan Sections -- Edge tracing can handle slit edges
            datasec=datasec,  # See above
            oscansec=oscansec,  # See above
        )
        return detector_container.DetectorContainer(**detector_dict)


class LDTRIMASHKSpectrograph(LDTRIMASSpectrograph):
    """
    Child to handle common aspects of the LDT/RIMAS HK Arm
    """

    ndet = 1
    camera = "RIMAS_HK"

    def get_detector_par(self, _, hdu=None):
        """
        Return metadata for the LDT/RIMAS YJ detector.

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
            dataext = 0  # Raw data
            binning = "1,1"  # Most common use mode
            gain = np.atleast_1d(1.52)  # Hardcoded in the header
            ronoise = np.atleast_1d(4.9)  # Hardcoded in the header
            datasec = np.atleast_1d("[5:512,53:2095]")  # For 1x1 binning
            oscansec = np.atleast_1d("[5:512,5:48]")  # For 1x1 binning
        else:
            # If file is post-processed, data extension is specified.  Raw is 0.
            dataext = hdu[0].header.get("POST_EXT", 0)
            binning = self.get_meta_value(self.get_headarr(hdu), "binning")
            gain = np.atleast_1d(hdu[0].header["GAIN"])
            ronoise = np.atleast_1d(hdu[0].header["RDNOISE"])
            datasec = hdu[0].header["TRIMSEC"]
            oscansec = hdu[0].header["BIASSEC"]

        # Detector
        detector_dict = dict(
            binning=binning,
            det=1,  # DeVeny has but one detector
            dataext=dataext,
            specaxis=1,  # Native spectrum is along the x-axis
            specflip=True,  # DeVeny CCD has blue at the right
            spatflip=False,
            platescale=0.34,  # Arcsec / pixel
            darkcurr=4.5,  # e-/pixel/hour
            saturation=65535.0,  # 16-bit ADC
            nonlinear=0.97,  # Linear to ~97% of saturation
            mincounts=-1e10,
            numamplifiers=1,
            gain=gain,  # See above
            ronoise=ronoise,  # See above
            # Data & Overscan Sections -- Edge tracing can handle slit edges
            datasec=datasec,  # See above
            oscansec=oscansec,  # See above
        )
        return detector_container.DetectorContainer(**detector_dict)


# Actual Operational Modes ===================================================#
class LDTRIMASLowYJSpectrograph(LDTRIMASYJSpectrograph):
    """
    Child to handle LDT/RIMAS YJ Arm, lowres-specific code
    """

    name = "ldt_rimas_yj_low"
    comment = "LDT Rapid infrared IMAger Spectrometer, YJ Arm Low-Res Gratings"
    pypeline = "MultiSlit"

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.

        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        # No bias or overscan for IRFPAs
        set_procpars = {
            "use_biasimage": False,
            "use_overscan": False,
            "use_darkimage": True,
            "use_illumflat": False,
        }
        par.reset_all_processimages_par(**set_procpars)

        # Do not use Dark image for Dark frame procesing
        par["calibrations"]["darkframe"]["process"]["use_darkimage"] = False

        # Slit-edge settings for NIHTS' slitlets
        par["calibrations"]["slitedges"]["edge_thresh"] = 15.0  # Default: 20.0
        par["calibrations"]["slitedges"]["fit_order"] = 2  # Default: 5
        par["calibrations"]["slitedges"]["max_nudge"] = 5  # Default: None
        par["calibrations"]["slitedges"]["minimum_slit_length"] = 6.0  # Default: None
        par["calibrations"]["slitedges"]["smash_range"] = [0.2, 0.5]  # Default: None
        par["calibrations"]["slitedges"]["sync_predict"] = "nearest"  # Default: 'pca'
        par["calibrations"]["slitedges"]["trace_thresh"] = 50  # Default: None
        par["calibrations"]["slitedges"]["trim_spec"] = [0, 50]  # Default: None

        # Only use LONG arc frames
        par["calibrations"]["arcframe"]["exprng"] = [30, None]
        # For processing the arc frame, these settings allow for the combination of
        #   of frames from different lamps into a comprehensible Master
        par["calibrations"]["arcframe"]["process"]["clip"] = False
        par["calibrations"]["arcframe"]["process"]["combine"] = "mean"
        # par['calibrations']['arcframe']['process']['subtract_continuum'] = True
        par["calibrations"]["tiltframe"]["process"]["clip"] = False
        par["calibrations"]["tiltframe"]["process"]["combine"] = "mean"
        # par['calibrations']['tiltframe']['process']['subtract_continuum'] = True

        # Wavelength Calibration Parameters
        # Arc lamps list from header -- instead of defining the full list here
        par["calibrations"]["wavelengths"]["lamps"] = ["XeI"]
        # Set this as default... but use `holy-grail` for DV4, DV8
        par["calibrations"]["wavelengths"][
            "method"
        ] = "holy-grail"  #'full_template'  # Default: 'holy-grail'
        # Reidentification parameters
        par["calibrations"]["wavelengths"]["reid_arxiv"] = "ldt_nihts.fits"
        # The DeVeny arc line FWHM varies based on slitwidth used
        par["calibrations"]["wavelengths"]["fwhm_fromlines"] = True  # Default: True
        par["calibrations"]["wavelengths"]["nsnippet"] = 1  # Default: 2

        # # For the tilts, our lines are not as well-behaved as others',
        # #   possibly due to the Wynne version E camera.
        # par["calibrations"]["tilts"]["spat_order"] = 4  # Default: 3
        # par["calibrations"]["tilts"]["spec_order"] = 5  # Default: 4

        # Flat-field parameter modification
        par["calibrations"]["flatfield"]["pixelflat_min_wave"] = 3000.0  # Default: None
        par["calibrations"]["flatfield"]["slit_illum_finecorr"] = False  # Default: True
        par["calibrations"]["flatfield"]["spec_samp_fine"] = 30  # Default: 1.2
        par["calibrations"]["flatfield"]["tweak_slits"] = False  # Default: True

        # Cosmic ray rejection parameters for science frames
        par["scienceframe"]["process"]["sigclip"] = 5.0  # Default: 4.5
        par["scienceframe"]["process"]["objlim"] = 2.0  # Default: 3.0

        # Object Finding, Extraction, and Sky Subtraction Parameters
        assumed_seeing = 1.5  # arcsec
        par["reduce"]["findobj"]["trace_npoly"] = 3  # Default: 5
        par["reduce"]["findobj"]["snr_thresh"] = 50.0  # Default: 10.0
        par["reduce"]["findobj"]["maxnumber_std"] = 1  # Default: 5
        par["reduce"]["findobj"]["maxnumber_sci"] = 5  # Default: 10
        par["reduce"]["findobj"]["find_fwhm"] = np.round(
            assumed_seeing / 0.34, 1
        )  # Default: 5.0 pix
        par["reduce"]["findobj"]["find_trim_edge"] = [0, 0]  # Default: [5, 5]
        # Boxcar width = ±3σ of Gaussian profile = >99% enclosed flux; radius = 1.28 * seeing
        par["reduce"]["extraction"]["boxcar_radius"] = np.round(
            assumed_seeing * 1.28, 1
        )  # Default: 1.5"
        par["reduce"]["extraction"]["use_2dmodel_mask"] = False  # Default: True
        par["reduce"]["skysub"]["sky_sigrej"] = 4.0  # Default: 3.0

        # Flexure Correction Parameters
        par["flexure"]["spec_method"] = "boxcar"  # Default: 'skip'
        par["flexure"]["spec_maxshift"] = 30  # Default: 20

        # Sensitivity Function Parameters
        par["sensfunc"]["UVIS"]["nresln"] = 15  # Default: 20
        par["sensfunc"]["UVIS"]["polycorrect"] = False  # Default: True

        return par

    def config_specific_par(self, scifile, inp_par=None):
        """
        Modify the PypeIt parameters to hard-wired values used for
        specific instrument configurations.

        In this case, choose between the two low-resolution longslit modes
        available for the YJ arm of RIMAS.

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
        # Start with instrument-wide parameters
        par = super().config_specific_par(scifile, inp_par=inp_par)

        # Adjust parameters based on DeVeny grating used
        grating = self.get_meta_value(scifile, "dispname")

        if grating == "DV1 (150/5000)":
            # Use this `reid_arxiv` with the `full-template` method:
            par["calibrations"]["wavelengths"][
                "reid_arxiv"
            ] = "ldt_deveny_150_HgCdAr.fits"
            # Because of the wide wavelength range, split DV1 arcs in half for reidentification
            par["calibrations"]["wavelengths"]["nsnippet"] = 2
            # Higher order wavelength fits because of larger span
            par["calibrations"]["wavelengths"]["n_first"] = 3  # Default: 2
            par["calibrations"]["wavelengths"]["n_final"] = 5  # Default: 4
            # The approximate resolution of this grating
            par["sensfunc"]["UVIS"]["resolution"] = 400

        elif grating == "DV2 (300/4000)":
            # Use this `reid_arxiv` with the `full-template` method:
            par["calibrations"]["wavelengths"][
                "reid_arxiv"
            ] = "ldt_deveny_300_HgCdAr.fits"
            # Higher order wavelength fits because of larger span
            par["calibrations"]["wavelengths"]["n_first"] = 3  # Default: 2
            par["calibrations"]["wavelengths"]["n_final"] = 5  # Default: 4
            # The approximate resolution of this grating
            par["sensfunc"]["UVIS"]["resolution"] = 800

        else:
            pass

        # Adjust parameters based on CCD binning
        binspec, binspat = parse.parse_binning(self.get_meta_value(scifile, "binning"))
        par["reduce"]["findobj"][
            "find_fwhm"
        ] /= binspat  # Specified in pixels and not arcsec
        par["flexure"]["spec_maxshift"] //= binspec  # Must be an integer
        par["sensfunc"]["UVIS"]["resolution"] /= binspec

        # SlitEdges Exclusion Regions (30 pixels at each edge) -- adjust based on binning
        excl_l, excl_r, last = np.array([30, 485, 515], dtype=int) // binspat
        par["calibrations"]["slitedges"][
            "exclude_regions"
        ] = f"1:0:{excl_l},1:{excl_r}:{last}"

        return par


class LDTRIMASLowHKSpectrograph(LDTRIMASHKSpectrograph):
    """
    Child to handle LDT/RIMAS HK Arm, lowres-specific code
    """

    name = "ldt_rimas_hk_low"
    comment = "LDT Rapid infrared IMAger Spectrometer, HK Arm Low-Res Gratings"
    pypeline = "MultiSlit"

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.

        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        # Turn off illumflat unless/until we can deal properly with flexure in
        #   the spatial direction.  All other defaults OK (as of v1.7.0)
        #   Also, use an order=1 chebyshev polynomial for fitting the overscan
        #   rather a SavGol filter -- more appropriate for this CCD.
        set_procpars = {
            "use_illumflat": False,
            "overscan_method": "chebyshev",
            "overscan_par": 1,
        }
        par.reset_all_processimages_par(**set_procpars)

        # For processing the arc frame, these settings allow for the combination of
        #   of frames from different lamps into a comprehensible Master
        par["calibrations"]["arcframe"]["process"]["clip"] = False
        par["calibrations"]["arcframe"]["process"]["combine"] = "mean"
        # par['calibrations']['arcframe']['process']['subtract_continuum'] = True
        par["calibrations"]["tiltframe"]["process"]["clip"] = False
        par["calibrations"]["tiltframe"]["process"]["combine"] = "mean"
        # par['calibrations']['tiltframe']['process']['subtract_continuum'] = True

        # Make a bad pixel mask
        par["calibrations"]["bpm_usebias"] = True

        # Wavelength Calibration Parameters
        # Arc lamps list from header -- instead of defining the full list here
        par["calibrations"]["wavelengths"]["lamps"] = ["use_header"]
        # Set this as default... but use `holy-grail` for DV4, DV8
        par["calibrations"]["wavelengths"][
            "method"
        ] = "full_template"  # Default: 'holy-grail'
        # The DeVeny arc line FWHM varies based on slitwidth used
        par["calibrations"]["wavelengths"]["fwhm"] = 3.0  # Default: 4.0
        par["calibrations"]["wavelengths"]["nsnippet"] = 1  # Default: 2

        # Slit-edge settings for long-slit data (DeVeny's slit is > 90" long)
        par["calibrations"]["slitedges"]["bound_detector"] = True  # Defualt: False
        par["calibrations"]["slitedges"]["sync_predict"] = "nearest"  # Default: 'pca'
        par["calibrations"]["slitedges"]["minimum_slit_length"] = 170.0  # Default: None
        par["calibrations"]["slitedges"]["max_nudge"] = 5  # Default: None

        # Flat-field parameter modification
        par["calibrations"]["flatfield"]["pixelflat_min_wave"] = 3000.0  # Default: None
        par["calibrations"]["flatfield"]["slit_illum_finecorr"] = False  # Default: True
        par["calibrations"]["flatfield"]["spec_samp_fine"] = 30  # Default: 1.2
        par["calibrations"]["flatfield"]["tweak_slits"] = False  # Default: True

        # For the tilts, our lines are not as well-behaved as others',
        #   possibly due to the Wynne version E camera.
        par["calibrations"]["tilts"]["spat_order"] = 4  # Default: 3
        par["calibrations"]["tilts"]["spec_order"] = 5  # Default: 4

        # Cosmic ray rejection parameters for science frames
        par["scienceframe"]["process"]["sigclip"] = 5.0  # Default: 4.5
        par["scienceframe"]["process"]["objlim"] = 2.0  # Default: 3.0

        # Object Finding, Extraction, and Sky Subtraction Parameters
        assumed_seeing = 1.5  # arcsec
        par["reduce"]["findobj"]["trace_npoly"] = 3  # Default: 5
        par["reduce"]["findobj"]["snr_thresh"] = 50.0  # Default: 10.0
        par["reduce"]["findobj"]["maxnumber_std"] = 1  # Default: 5
        par["reduce"]["findobj"]["maxnumber_sci"] = 5  # Default: 10
        par["reduce"]["findobj"]["find_fwhm"] = np.round(
            assumed_seeing / 0.34, 1
        )  # Default: 5.0 pix
        par["reduce"]["findobj"]["find_trim_edge"] = [0, 0]  # Default: [5, 5]
        # Boxcar width = ±3σ of Gaussian profile = >99% enclosed flux; radius = 1.28 * seeing
        par["reduce"]["extraction"]["boxcar_radius"] = np.round(
            assumed_seeing * 1.28, 1
        )  # Default: 1.5"
        par["reduce"]["extraction"]["use_2dmodel_mask"] = False  # Default: True
        par["reduce"]["skysub"]["sky_sigrej"] = 4.0  # Default: 3.0

        # Flexure Correction Parameters
        par["flexure"]["spec_method"] = "boxcar"  # Default: 'skip'
        par["flexure"]["spec_maxshift"] = 30  # Default: 20

        # Sensitivity Function Parameters
        par["sensfunc"]["UVIS"]["nresln"] = 15  # Default: 20
        par["sensfunc"]["UVIS"]["polycorrect"] = False  # Default: True

        return par

    def config_specific_par(self, scifile, inp_par=None):
        """
        Modify the PypeIt parameters to hard-wired values used for
        specific instrument configurations.

        In this case, choose between the two low-resolution longslit modes
        available for the YJ arm of RIMAS.

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
        # Start with instrument-wide parameters
        par = super().config_specific_par(scifile, inp_par=inp_par)

        # Adjust parameters based on DeVeny grating used
        grating = self.get_meta_value(scifile, "dispname")

        if grating == "DV1 (150/5000)":
            # Use this `reid_arxiv` with the `full-template` method:
            par["calibrations"]["wavelengths"][
                "reid_arxiv"
            ] = "ldt_deveny_150_HgCdAr.fits"
            # Because of the wide wavelength range, split DV1 arcs in half for reidentification
            par["calibrations"]["wavelengths"]["nsnippet"] = 2
            # Higher order wavelength fits because of larger span
            par["calibrations"]["wavelengths"]["n_first"] = 3  # Default: 2
            par["calibrations"]["wavelengths"]["n_final"] = 5  # Default: 4
            # The approximate resolution of this grating
            par["sensfunc"]["UVIS"]["resolution"] = 400

        elif grating == "DV2 (300/4000)":
            # Use this `reid_arxiv` with the `full-template` method:
            par["calibrations"]["wavelengths"][
                "reid_arxiv"
            ] = "ldt_deveny_300_HgCdAr.fits"
            # Higher order wavelength fits because of larger span
            par["calibrations"]["wavelengths"]["n_first"] = 3  # Default: 2
            par["calibrations"]["wavelengths"]["n_final"] = 5  # Default: 4
            # The approximate resolution of this grating
            par["sensfunc"]["UVIS"]["resolution"] = 800

        else:
            pass

        # Adjust parameters based on CCD binning
        binspec, binspat = parse.parse_binning(self.get_meta_value(scifile, "binning"))
        par["reduce"]["findobj"][
            "find_fwhm"
        ] /= binspat  # Specified in pixels and not arcsec
        par["flexure"]["spec_maxshift"] //= binspec  # Must be an integer
        par["sensfunc"]["UVIS"]["resolution"] /= binspec

        # SlitEdges Exclusion Regions (30 pixels at each edge) -- adjust based on binning
        excl_l, excl_r, last = np.array([30, 485, 515], dtype=int) // binspat
        par["calibrations"]["slitedges"][
            "exclude_regions"
        ] = f"1:0:{excl_l},1:{excl_r}:{last}"

        return par


class LDTRIMASEchelleYJSpectrograph(LDTRIMASYJSpectrograph):
    """
    Child to handle LDT/RIMAS YJ Arm, echelle-specific code
    """

    name = "ldt_rimas_yj_med"
    comment = "LDT Rapid infrared IMAger Spectrometer, YJ Arm Med-Res Echelle Grism"
    pypeline = "Echelle"
    ech_fixed_format = True

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.

        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        # Turn off illumflat unless/until we can deal properly with flexure in
        #   the spatial direction.  All other defaults OK (as of v1.7.0)
        #   Also, use an order=1 chebyshev polynomial for fitting the overscan
        #   rather a SavGol filter -- more appropriate for this CCD.
        set_procpars = {
            "use_illumflat": False,
            "overscan_method": "chebyshev",
            "overscan_par": 1,
        }
        par.reset_all_processimages_par(**set_procpars)

        # For processing the arc frame, these settings allow for the combination of
        #   of frames from different lamps into a comprehensible Master
        par["calibrations"]["arcframe"]["process"]["clip"] = False
        par["calibrations"]["arcframe"]["process"]["combine"] = "mean"
        # par['calibrations']['arcframe']['process']['subtract_continuum'] = True
        par["calibrations"]["tiltframe"]["process"]["clip"] = False
        par["calibrations"]["tiltframe"]["process"]["combine"] = "mean"
        # par['calibrations']['tiltframe']['process']['subtract_continuum'] = True

        # Make a bad pixel mask
        par["calibrations"]["bpm_usebias"] = True

        # Wavelength Calibration Parameters
        # Arc lamps list from header -- instead of defining the full list here
        par["calibrations"]["wavelengths"]["lamps"] = ["use_header"]
        # Set this as default... but use `holy-grail` for DV4, DV8
        par["calibrations"]["wavelengths"][
            "method"
        ] = "full_template"  # Default: 'holy-grail'
        # The DeVeny arc line FWHM varies based on slitwidth used
        par["calibrations"]["wavelengths"]["fwhm"] = 3.0  # Default: 4.0
        par["calibrations"]["wavelengths"]["nsnippet"] = 1  # Default: 2

        # Slit-edge settings for long-slit data (DeVeny's slit is > 90" long)
        par["calibrations"]["slitedges"]["bound_detector"] = True  # Defualt: False
        par["calibrations"]["slitedges"]["sync_predict"] = "nearest"  # Default: 'pca'
        par["calibrations"]["slitedges"]["minimum_slit_length"] = 170.0  # Default: None
        par["calibrations"]["slitedges"]["max_nudge"] = 5  # Default: None

        # Flat-field parameter modification
        par["calibrations"]["flatfield"]["pixelflat_min_wave"] = 3000.0  # Default: None
        par["calibrations"]["flatfield"]["slit_illum_finecorr"] = False  # Default: True
        par["calibrations"]["flatfield"]["spec_samp_fine"] = 30  # Default: 1.2
        par["calibrations"]["flatfield"]["tweak_slits"] = False  # Default: True

        # For the tilts, our lines are not as well-behaved as others',
        #   possibly due to the Wynne version E camera.
        par["calibrations"]["tilts"]["spat_order"] = 4  # Default: 3
        par["calibrations"]["tilts"]["spec_order"] = 5  # Default: 4

        # Cosmic ray rejection parameters for science frames
        par["scienceframe"]["process"]["sigclip"] = 5.0  # Default: 4.5
        par["scienceframe"]["process"]["objlim"] = 2.0  # Default: 3.0

        # Object Finding, Extraction, and Sky Subtraction Parameters
        assumed_seeing = 1.5  # arcsec
        par["reduce"]["findobj"]["trace_npoly"] = 3  # Default: 5
        par["reduce"]["findobj"]["snr_thresh"] = 50.0  # Default: 10.0
        par["reduce"]["findobj"]["maxnumber_std"] = 1  # Default: 5
        par["reduce"]["findobj"]["maxnumber_sci"] = 5  # Default: 10
        par["reduce"]["findobj"]["find_fwhm"] = np.round(
            assumed_seeing / 0.34, 1
        )  # Default: 5.0 pix
        par["reduce"]["findobj"]["find_trim_edge"] = [0, 0]  # Default: [5, 5]
        # Boxcar width = ±3σ of Gaussian profile = >99% enclosed flux; radius = 1.28 * seeing
        par["reduce"]["extraction"]["boxcar_radius"] = np.round(
            assumed_seeing * 1.28, 1
        )  # Default: 1.5"
        par["reduce"]["extraction"]["use_2dmodel_mask"] = False  # Default: True
        par["reduce"]["skysub"]["sky_sigrej"] = 4.0  # Default: 3.0

        # Flexure Correction Parameters
        par["flexure"]["spec_method"] = "boxcar"  # Default: 'skip'
        par["flexure"]["spec_maxshift"] = 30  # Default: 20

        # Sensitivity Function Parameters
        par["sensfunc"]["UVIS"]["nresln"] = 15  # Default: 20
        par["sensfunc"]["UVIS"]["polycorrect"] = False  # Default: True

        return par


class LDTRIMASEchelleHKSpectrograph(LDTRIMASHKSpectrograph):
    """
    Child to handle LDT/RIMAS HK Arm, echelle-specific code
    """

    name = "ldt_rimas_hk_med"
    comment = "LDT Rapid infrared IMAger Spectrometer, HK Arm Med-Res Echelle Grism"
    pypeline = "Echelle"
    ech_fixed_format = True

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.

        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        # Turn off illumflat unless/until we can deal properly with flexure in
        #   the spatial direction.  All other defaults OK (as of v1.7.0)
        #   Also, use an order=1 chebyshev polynomial for fitting the overscan
        #   rather a SavGol filter -- more appropriate for this CCD.
        set_procpars = {
            "use_illumflat": False,
            "overscan_method": "chebyshev",
            "overscan_par": 1,
        }
        par.reset_all_processimages_par(**set_procpars)

        # For processing the arc frame, these settings allow for the combination of
        #   of frames from different lamps into a comprehensible Master
        par["calibrations"]["arcframe"]["process"]["clip"] = False
        par["calibrations"]["arcframe"]["process"]["combine"] = "mean"
        # par['calibrations']['arcframe']['process']['subtract_continuum'] = True
        par["calibrations"]["tiltframe"]["process"]["clip"] = False
        par["calibrations"]["tiltframe"]["process"]["combine"] = "mean"
        # par['calibrations']['tiltframe']['process']['subtract_continuum'] = True

        # Make a bad pixel mask
        par["calibrations"]["bpm_usebias"] = True

        # Wavelength Calibration Parameters
        # Arc lamps list from header -- instead of defining the full list here
        par["calibrations"]["wavelengths"]["lamps"] = ["use_header"]
        # Set this as default... but use `holy-grail` for DV4, DV8
        par["calibrations"]["wavelengths"][
            "method"
        ] = "full_template"  # Default: 'holy-grail'
        # The DeVeny arc line FWHM varies based on slitwidth used
        par["calibrations"]["wavelengths"]["fwhm"] = 3.0  # Default: 4.0
        par["calibrations"]["wavelengths"]["nsnippet"] = 1  # Default: 2

        # Slit-edge settings for long-slit data (DeVeny's slit is > 90" long)
        par["calibrations"]["slitedges"]["bound_detector"] = True  # Defualt: False
        par["calibrations"]["slitedges"]["sync_predict"] = "nearest"  # Default: 'pca'
        par["calibrations"]["slitedges"]["minimum_slit_length"] = 170.0  # Default: None
        par["calibrations"]["slitedges"]["max_nudge"] = 5  # Default: None

        # Flat-field parameter modification
        par["calibrations"]["flatfield"]["pixelflat_min_wave"] = 3000.0  # Default: None
        par["calibrations"]["flatfield"]["slit_illum_finecorr"] = False  # Default: True
        par["calibrations"]["flatfield"]["spec_samp_fine"] = 30  # Default: 1.2
        par["calibrations"]["flatfield"]["tweak_slits"] = False  # Default: True

        # For the tilts, our lines are not as well-behaved as others',
        #   possibly due to the Wynne version E camera.
        par["calibrations"]["tilts"]["spat_order"] = 4  # Default: 3
        par["calibrations"]["tilts"]["spec_order"] = 5  # Default: 4

        # Cosmic ray rejection parameters for science frames
        par["scienceframe"]["process"]["sigclip"] = 5.0  # Default: 4.5
        par["scienceframe"]["process"]["objlim"] = 2.0  # Default: 3.0

        # Object Finding, Extraction, and Sky Subtraction Parameters
        assumed_seeing = 1.5  # arcsec
        par["reduce"]["findobj"]["trace_npoly"] = 3  # Default: 5
        par["reduce"]["findobj"]["snr_thresh"] = 50.0  # Default: 10.0
        par["reduce"]["findobj"]["maxnumber_std"] = 1  # Default: 5
        par["reduce"]["findobj"]["maxnumber_sci"] = 5  # Default: 10
        par["reduce"]["findobj"]["find_fwhm"] = np.round(
            assumed_seeing / 0.34, 1
        )  # Default: 5.0 pix
        par["reduce"]["findobj"]["find_trim_edge"] = [0, 0]  # Default: [5, 5]
        # Boxcar width = ±3σ of Gaussian profile = >99% enclosed flux; radius = 1.28 * seeing
        par["reduce"]["extraction"]["boxcar_radius"] = np.round(
            assumed_seeing * 1.28, 1
        )  # Default: 1.5"
        par["reduce"]["extraction"]["use_2dmodel_mask"] = False  # Default: True
        par["reduce"]["skysub"]["sky_sigrej"] = 4.0  # Default: 3.0

        # Flexure Correction Parameters
        par["flexure"]["spec_method"] = "boxcar"  # Default: 'skip'
        par["flexure"]["spec_maxshift"] = 30  # Default: 20

        # Sensitivity Function Parameters
        par["sensfunc"]["UVIS"]["nresln"] = 15  # Default: 20
        par["sensfunc"]["UVIS"]["polycorrect"] = False  # Default: True

        return par
