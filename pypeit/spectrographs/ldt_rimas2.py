# pylint: disable=use-dict-literal
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
mirror, Y, J, H and K broadband filters, two VPHs, ZnSe grisms with cross-
dispersers, and Teledyne H4RG-10 detectors.  RIMAS is a flexible tool with both
imaging and spectral modes available.

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

import dataclasses
import pathlib

import astropy.io.fits
import astropy.table
import astropy.time
import numpy as np

from pypeit import log
from pypeit import PypeItError
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.core import parse
from pypeit.images import detector_container
from pypeit.par import parset
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph

# Module-Level Constants
# TODO: Update these with commissioning data after RIMAS returns in 2026B!
RIMAS_FLAT_EXPTIME_RANGES = {
    ("Vph30", "0.6''", "YJ"): [5.0, 15.0],
    ("Vph30", "0.6''", "HK"): [None, 5.0],
    ("Vph30", "1.0''", "YJ"): [5.0, 15.0],
    ("Vph30", "1.0''", "HK"): [5.0, 15.0],
    ("Vph30", "2.0''", "YJ"): [None, 10.0],
    ("Vph30", "2.0''", "HK"): [None, 10.0],
    ("Vph30", "1.2'' long", "YJ"): [15.0, 45.0],
    ("Vph30", "1.2'' long", "HK"): [None, 10.0],
    ("Vph300", "0.6''", "YJ"): [45.0, None],
    ("Vph300", "0.6''", "HK"): [10.0, 45.0],
    ("Vph300", "1.0''", "YJ"): [30.0, None],
    ("Vph300", "1.0''", "HK"): [10.0, 30.0],
    ("Vph300", "2.0''", "YJ"): [15.0, 45.0],
    ("Vph300", "2.0''", "HK"): [None, 15.0],
    ("Vph300", "1.2'' long", "YJ"): [30.0, None],
    ("Vph300", "1.2'' long", "HK"): [None, 30.0],
    ("Grism", "0.6''", "YJ"): [90.0, None],
    ("Grism", "0.6''", "HK"): [30.0, 90.0],
    ("Grism", "1.0''", "YJ"): [90.0, None],
    ("Grism", "1.0''", "HK"): [15.0, 60.0],
    ("Grism", "2.0''", "YJ"): [90.0, None],
    ("Grism", "2.0''", "HK"): [10.0, 30.0],
}

RIMAS_DITHER_CORRECTION_FILE = "rimas_dither_corrections.ecsv"


@dataclasses.dataclass
class EchelleProps:
    """Fixed-format echelle geometry for one RIMAS configuration.

    Attributes
    ----------
    norders : int
        Number of echelle orders.
    orders : `numpy.ndarray`_
        Order numbers.
    order_spat_pos : `numpy.ndarray`_
        Fractional spatial positions of each order.
    spec_min : `numpy.ndarray`_
        Minimum spectral pixel for each order.
    spec_max : `numpy.ndarray`_
        Maximum spectral pixel for each order.
    platescale : float
        Spatial platescale in arcsec/pixel.
    dloglam : float
        Base-10 logarithmic wavelength step.
    loglam_minmax : tuple
        Minimum and maximum base-10 logarithmic wavelength.
    """

    norders: int
    orders: np.ndarray
    order_spat_pos: np.ndarray
    spec_min: np.ndarray
    spec_max: np.ndarray
    platescale: float
    dloglam: float
    loglam_minmax: tuple[float, float]


# Base LDT/RIMAS Class =======================================================#
class LDTRIMASSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle LDT/RIMAS specific code

    This class contains the common methods for the operating modes of RIMAS
    data reduction: single-order, echelle.
    """

    telescope = telescopes.LDTTelescopePar()
    url = "https://lowell.edu/research/telescopes-and-facilities/ldt/rimas/"
    header_name = "RIMAS"
    supported = True
    ndet = 2
    camera = "RIMAS"
    allowed_extensions = [".fits"]
    add_bg_columns = True

    def __init__(self):
        super().__init__()
        self.arm = None
        self.decker = None

    def get_detector_par(
        self, det: int, hdu: astropy.io.fits.HDUList | None = None
    ) -> detector_container.DetectorContainer:
        """
        Return metadata for the selected detector.

        .. warning::

            Many of the necessary detector parameters are read from the file
            header, meaning the ``hdu`` argument is effectively **required** for
            LDT/RIMAS.  The optional use of ``hdu`` is only viable for
            automatically generated documentation.

        Parameters
        ----------
        det : int
            1-indexed detector number.
        hdu : `astropy.io.fits.HDUList`_, optional
            The open FITS file with the raw image of interest.

        Returns
        -------
        :class:`~pypeit.images.detector_container.DetectorContainer`
            Object with the detector metadata.
        """
        if hdu is None:
            binning = "1,1"  # Most common use mode
            gain = np.atleast_1d(1.8)  # Hardcoded in the header
            ronoise = np.atleast_1d(4.9)  # Hardcoded in the header
            datasec = np.atleast_1d("[:,:]")  # For 1x1 binning
            disp = np.atleast_1d("Grism")

        else:
            binning = self.get_meta_value(self.get_headarr(hdu), "binning")
            gain = np.atleast_1d(hdu[0].header["GAIN0"])
            ronoise = np.atleast_1d(4.9)
            datasec = np.atleast_1d("[:,:]")
            disp = np.atleast_1d(
                hdu[0].header[
                    "FILTER1" if hdu[0].header["CAMNAME"] == "YJ" else "FILTER2"
                ]
            )

        # TODO:
        # Because of the VPH300 YJ issue, the detector SPECFLIP needs to be
        #   cased out based on arm, grating, etc.

        # Detector 1
        detector_dict1 = dict(
            binning=binning,
            det=1,  # The YJ channel is DETECTOR #1 (FITS HEADER CHIP = 1)
            dataext=0,
            specaxis=1,  # Native spectrum is along the x-axis
            specflip=self.get_specflip(1, disp),  # RIMAS IR FPAs have blue at the right
            spatflip=False,
            platescale=0.19,  # Arcsec / pixel
            darkcurr=0,  # e-/pixel/hour
            saturation=65535,  # 16-bit ADC
            nonlinear=0.97,  # Linear to ~97% of saturation
            mincounts=-1e10,
            numamplifiers=1,
            gain=gain,  # See above
            ronoise=ronoise,  # See above
            # Data & Overscan Sections -- Edge tracing can handle slit edges
            datasec=datasec,  # See above
        )
        # Detector 2
        detector_dict2 = detector_dict1.copy()
        detector_dict2.update(
            dict(
                det=2,  # The HK channel is DETECTOR #2 (FITS HEADER CHIP = 2)
                specflip=self.get_specflip(2, disp),
                darkcurr=0.0,  # e-/pixel/hour
            )
        )

        detectors = [detector_dict1, detector_dict2]
        # Return
        return detector_container.DetectorContainer(**detectors[det - 1])

    @staticmethod
    def get_specflip(det: int, disp: str | np.ndarray) -> bool:
        """Get the spectral flip based on grating and arm.

        Parameters
        ----------
        det : int
            Detector number (YJ = 1, HK = 2).
        disp : str or `numpy.ndarray`_
            Disperser for this arm.

        Returns
        -------
        bool
            Whether the spectral direction is flipped from red to blue.
        """
        # The YJ Vph300 grating is installed such that BLUE -> RED
        if (det == 1 and disp == "Vph300") or (det == 2 and disp == "Vph30"):
            return False

        # Most everything in RIMAS is RED -> BLUE
        return True

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
        self.meta["decker"] = dict(card=None, compound=True)  # SLIT filter wheel
        self.meta["binning"] = dict(card=None, compound=True)
        self.meta["rawshape"] = dict(card=None, compound=True)
        self.meta["mjd"] = dict(card=None, compound=True)
        self.meta["airmass"] = dict(ext=0, card="AIRMASS")
        self.meta["exptime"] = dict(card=None, compound=True)
        self.meta["instrument"] = dict(ext=0, card="INSTRUME")

        # Extras for config and frametyping
        # NOTE: `rtol` is _relative_ tolerance (e.g. 1 part in 1,000)
        self.meta["arm"] = dict(ext=0, card="CAMNAME")
        self.meta["idname"] = dict(card=None, compound=True)
        self.meta["filter1"] = dict(ext=0, card="FILTER3")  # AUX filter wheel
        self.meta["lampstat01"] = dict(card=None, compound=True)
        self.meta["slitwid"] = dict(card=None, compound=True)
        self.meta["dithpat"] = dict(ext=0, card="DITHTYP")
        self.meta["dithpos"] = dict(card=None, compound=True)
        self.meta["dithoff"] = dict(ext=0, card="DITHRAD")

    def compound_meta(
        self, headarr: list[astropy.io.fits.Header], meta_key: str
    ) -> object:
        """
        Methods to generate metadata requiring interpretation of the header
        data, instead of simply reading the value of a header card.

        Parameters
        ----------
        headarr : list
            List of `astropy.io.fits.Header`_ objects.
        meta_key : str
            Metadata keyword to construct.

        Returns
        -------
        object
            Metadata value read from the header(s).
        """
        if meta_key == "dispname":
            # Return FILTER1 (Filter Wheel Name FW YJ) for YJ frames and
            #   FILTER2 (Filter Wheel Name FW HK) for HK frames
            return headarr[0]["FILTER1" if headarr[0]["CAMNAME"] == "YJ" else "FILTER2"]

        if meta_key == "decker":
            # Older RIMAS headers used alternate names for the 1.2" long slit.
            decker = headarr[0]["FILTER4"].strip()
            if decker in ["long", "1.2'''' long"]:
                return "1.2'' long"
            return decker

        if meta_key == "idname":
            # Force uppercase to match other LDT instruments
            return headarr[0]["OBJTYPE"].upper()

        if meta_key == "binning":
            # Binning in RIMAS headers given as separate values
            if "BINX" in headarr[0]:
                return parse.binning2string(headarr[0]["BINX"], headarr[0]["BINY"])
            if "XBINNING" in headarr[0]:
                return parse.binning2string(
                    headarr[0]["XBINNING"], headarr[0]["YBINNING"]
                )
            if (
                "BINNING" in headarr[0]
            ):  # used in order to make the spec1d file generated by the pipeline usable with pypeit_sensfunc
                return headarr[0]["BINNING"]

        if meta_key == "rawshape":
            # Shape of the raw image written to disk, formatted as rows,columns.
            return parse.binning2string(headarr[0]["NAXIS2"], headarr[0]["NAXIS1"])

        if meta_key == "mjd":
            # Use custom scrubber + AstroPy to convert 'DATE-OBS' into a mjd.
            ttime = self.scrub_isot_dateobs(headarr[0]["DATE-BEG"])
            return ttime.mjd

        if meta_key == "exptime":
            # If older than mid-March 2026, fix bug by computing effective exptime
            if self.scrub_isot_dateobs(headarr[0]["DATE-BEG"]).mjd < 61114.0:
                # Total exposure time minus the time for the "pedestal" initial frame
                return np.round(headarr[0]["EXPTIME"] - headarr[0]["FRTIME"], 2)

            # Otherwise, return EXPTIMEE (correct effective exposure time)
            return np.round(headarr[0]["EXPTIMEE"], 2)

        if meta_key == "lampstat01":
            # The spectral comparison lamps turned on are listed in `LAMPCAL`, but
            #  if no lamps are on, then this string is blank.  Return either the
            #  populated `LAMPCAL` string, or 'off' to ensure a positive entry for
            #  `lampstat01`.
            lampcal = ""
            return "off" if lampcal == "" else lampcal

        if meta_key == "slitwid":
            # Convert the decker into a slitwidth in arcseconds
            decker = headarr[0]["FILTER4"].strip()
            if decker == "0.6''":
                return 0.6
            if decker == "1.0''":
                return 1.0
            if decker == "2.0''":
                return 2.0
            if decker == "1.2'' long":
                return 1.2
            return 0

        if meta_key == "dithpos":
            # Parse the dither type and location in pattern
            dithtype = headarr[0]["DITHTYP"].strip().upper()
            if dithtype == "ABBA":
                # For ABBA dithers, convert DITHNUM to A or B
                cyclepos = (headarr[0]["DITHNUM"] - 1) % 4
                return "A" if cyclepos in [0, 3] else "B"
            if dithtype == "ONOFF":
                # For OnOff dithers, convert DITHNUM to On or Off
                cyclepos = (headarr[0]["DITHNUM"] - 1) % 2
                return "On" if cyclepos == 0 else "Off"
            return "None"

        log.error("Not ready for compound meta %s for LDT/RIMAS", meta_key)

    def config_independent_frames(self) -> dict[str, str | list[str] | None]:
        """
        Define frame types that are independent of the fully defined
        instrument configuration.

        This method returns a dictionary where the keys of the dictionary are
        the list of configuration-independent frame types. The value of each
        dictionary element can be set to one or more metadata keys that can
        be used to assign each frame type to a given configuration group. See
        :func:`~pypeit.metadata.PypeItMetaData.set_configurations` and how it
        interprets the dictionary values, which can be None.

        Returns
        -------
        dict
            Dictionary where the keys are the frame types that are
            configuration-independent and the values are the metadata keywords
            that can be used to assign the frames to a configuration group.
        """
        # Dark frames are independent of the optical configuration, but they
        # must match the detector arm and raw readout geometry.
        return {"dark": ["arm", "rawshape"]}

    def configuration_keys(self) -> list[str]:
        """
        Return the metadata keys that define a unique instrument
        configuration.

        This list is used by :class:`~pypeit.metadata.PypeItMetaData` to
        identify the unique configurations among the list of frames read
        for a given reduction.

        Returns
        -------
        list
            List of keywords of data pulled from file headers and used to
            construct the :class:`~pypeit.metadata.PypeItMetaData` object.
        """
        return ["arm", "dispname", "decker"]

    def raw_header_cards(self) -> list[str]:
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

        Returns
        -------
        list
            List of keywords from the raw data files that should be propagated
            in output files.
        """
        return ["CAMNAME", "FILTER1", "FILTER2", "FILTER3", "FILTER4"]

    def pypeit_file_keys(self) -> list[str]:
        """
        Define the list of keys to be output into a standard PypeIt file.

        Returns
        -------
        list
            The list of keywords in the relevant
            :class:`~pypeit.metadata.PypeItMetaData` instance to print to the
            :ref:`pypeit_file`.
        """
        return super().pypeit_file_keys() + [
            "slitwid",
            "lampstat01",
            "filter1",
            "rawshape",
            "dither",
            "dithpat",
            "dithpos",
            "dithoff",
        ]

    @classmethod
    def default_pypeit_par(cls) -> pypeitpar.PypeItPar:
        """
        Return the default parameters to use for all of RIMAS.

        ..note ::

            Each of the child classes will have modifications on top of these,
            but some parameters are instrument-wide.

        Returns
        -------
        :class:`~pypeit.par.pypeitpar.PypeItPar`
            Parameters required by all PypeIt methods.
        """
        # Get the PypeIt default parameters
        par = super().default_pypeit_par()

        # No bias or overscan for IRFPAs.  Also disable flat-fielding globally;
        # only science/standard frames should apply the RIMAS flat products.
        turn_off = {
            "use_illumflat": False,
            "use_biasimage": False,
            "use_overscan": False,
            "use_darkimage": False,
            "use_pixelflat": False,
        }
        par.reset_all_processimages_par(**turn_off)

        # Use dark frames for ARC and TILT, which are commonly also science frames.
        par["calibrations"]["arcframe"]["process"]["use_darkimage"] = True
        par["calibrations"]["tiltframe"]["process"]["use_darkimage"] = True
        # Use dark frames and flats for SCIENCE and STANDARD frames.
        par["scienceframe"]["process"]["use_darkimage"] = True
        par["calibrations"]["standardframe"]["process"]["use_darkimage"] = True
        par["scienceframe"]["process"]["use_pixelflat"] = True
        par["scienceframe"]["process"]["use_illumflat"] = True
        par["calibrations"]["standardframe"]["process"]["use_pixelflat"] = True
        par["calibrations"]["standardframe"]["process"]["use_illumflat"] = True
        # Do not mask CRs in dark frames -- it actually removes the hot pixels!
        par["calibrations"]["darkframe"]["process"]["mask_cr"] = False

        # Everybody is going to use OH lines
        par["calibrations"]["wavelengths"]["lamps"] = ["OH_GNIRS"]

        # Shared RIMAS spectral reduction defaults.
        par["calibrations"]["slitedges"]["edge_thresh"] = 50.0
        par["calibrations"]["slitedges"]["fit_order"] = 2
        par["calibrations"]["slitedges"]["max_shift_adj"] = 0.5
        par["calibrations"]["slitedges"]["trace_thresh"] = 10.0
        par["calibrations"]["slitedges"]["fit_min_spec_length"] = 0.5
        par["calibrations"]["slitedges"]["left_right_pca"] = True
        par["calibrations"]["slitedges"]["length_range"] = 0.3

        par["calibrations"]["arcframe"]["process"]["clip"] = False
        par["calibrations"]["arcframe"]["process"]["combine"] = "mean"
        par["calibrations"]["tiltframe"]["process"]["clip"] = False
        par["calibrations"]["tiltframe"]["process"]["combine"] = "mean"

        par["calibrations"]["flatfield"]["slit_illum_finecorr"] = False
        par["calibrations"]["flatfield"]["spec_samp_coarse"] = 3.0
        par["calibrations"]["flatfield"]["spat_samp"] = 2.0
        par["calibrations"]["flatfield"]["tweak_slits"] = False

        par["scienceframe"]["process"]["sigclip"] = 5.0
        par["scienceframe"]["process"]["objlim"] = 2.0
        par["scienceframe"]["process"]["satpix"] = "nothing"

        assumed_seeing = 1.5  # arcsec
        par["reduce"]["findobj"]["trace_npoly"] = 3
        par["reduce"]["findobj"]["maxnumber_std"] = 1
        par["reduce"]["findobj"]["maxnumber_sci"] = 5
        par["reduce"]["findobj"]["find_fwhm"] = np.round(assumed_seeing / 0.34, 1)
        par["reduce"]["findobj"]["find_trim_edge"] = [0, 0]
        par["reduce"]["extraction"]["boxcar_radius"] = np.round(
            assumed_seeing * 1.28, 1
        )
        par["reduce"]["extraction"]["use_2dmodel_mask"] = False
        par["reduce"]["skysub"]["sky_sigrej"] = 4.0

        par["flexure"]["spec_method"] = "boxcar"
        par["flexure"]["spec_maxshift"] = 30

        par["sensfunc"]["algorithm"] = "IR"
        par["sensfunc"]["polyorder"] = 8
        par["sensfunc"]["extrap_red"] = 0.4
        par["sensfunc"]["extrap_blu"] = 0.4
        par["sensfunc"]["IR"]["telgridfile"] = "TellPCA_3000_26000_R10000.fits"
        par["sensfunc"]["IR"]["maxiter"] = 2
        par["sensfunc"]["IR"]["lower"] = 3.0
        par["sensfunc"]["IR"]["upper"] = 3.0

        # TODO tune up LA COSMICS parameters here for X-shooter as tellurics
        # are being excessively masked

        # # Tilt parameters
        # par["calibrations"]["arcframe"]["process"]["subtract_continuum"] = True
        # par["calibrations"]["tiltframe"]["process"]["subtract_continuum"] = True
        # par["calibrations"]["tilts"]["rm_continuum"] = True
        # par["calibrations"]["tilts"]["tracethresh"] = 25.0
        # par["calibrations"]["tilts"]["maxdev_tracefit"] = 0.04
        # par["calibrations"]["tilts"]["maxdev2d"] = 0.04
        # par["calibrations"]["tilts"]["spat_order"] = 3
        # par["calibrations"]["tilts"]["spec_order"] = 4

        return par

    def config_specific_par(
        self,
        inp: str | list | pathlib.Path | astropy.io.fits.Header | astropy.table.Table,
        inp_par: parset.ParSet = None,
    ) -> parset.ParSet:
        """
        Modify the PypeIt parameters to hard-wired values used for
        specific instrument configurations.

        Parameters
        ----------
        inp : str, list, pathlib.Path, `astropy.io.fits.Header`_, or `astropy.table.Table`_
            Input filename, an `astropy.io.fits.Header`_ object, a list of
            `astropy.io.fits.Header`_ objects, or a row from the metadata
            table.
        inp_par : :class:`~pypeit.par.parset.ParSet`, optional
            Parameter set used for the full run of PypeIt. If None, use
            :func:`default_pypeit_par`.

        Returns
        -------
        :class:`~pypeit.par.parset.ParSet`
            The PypeIt parameter set adjusted for configuration-specific
            parameter values.
        """
        # Start with instrument-wide parameters
        par = super().config_specific_par(inp, inp_par=inp_par)

        # Save the arm and decker as instance variables for later use
        self.arm = self.get_meta_value(inp, "arm")
        self.decker = self.get_meta_value(inp, "decker")

        # Set the detector number based on the arm
        par["rdx"]["detnum"] = 1 if self.arm == "YJ" else 2

        return par

    @staticmethod
    def adjust_config_specific_binning(
        par: parset.ParSet, binning: str
    ) -> parset.ParSet:
        """Adjust RIMAS parameters that are specified in binned pixels.

        Parameters
        ----------
        par : :class:`~pypeit.par.parset.ParSet`
            Parameter set to modify.
        binning : str
            Spectral and spatial binning string.

        Returns
        -------
        :class:`~pypeit.par.parset.ParSet`
            Modified parameter set.
        """
        binspec, binspat = parse.parse_binning(binning)
        par["reduce"]["findobj"][
            "find_fwhm"
        ] /= binspat  # Specified in pixels and not arcsec
        par["flexure"]["spec_maxshift"] //= binspec  # Must be an integer

        # SlitEdges Exclusion Regions (30 pixels at each edge) -- adjust based on binning
        excl_l, excl_r, last = np.array([30, 485, 515], dtype=int) // binspat
        par["calibrations"]["slitedges"][
            "exclude_regions"
        ] = f"1:0:{excl_l},1:{excl_r}:{last}"

        return par

    def check_frame_type(
        self,
        ftype: str,
        fitstbl: astropy.table.Table,
        exprng: list[float] | None = None,
    ) -> np.ndarray:
        """
        Check for frames of the provided type.

        Parameters
        ----------
        ftype : str
            Type of frame to check. Must be a valid frame type; see frame-type
            :ref:`frame_type_defs`.
        fitstbl : `astropy.table.Table`_
            The table with the metadata for one or more frames to check.
        exprng : list, optional
            Range in the allowed exposure time for a frame of type ``ftype``.
            See :func:`~pypeit.core.framematch.check_frame_exptime`.

        Returns
        -------
        `numpy.ndarray`_
            Boolean array with the flags selecting the exposures in
            ``fitstbl`` that are ``ftype`` type frames.
        """

        good_exp = framematch.check_frame_exptime(fitstbl["exptime"], exprng)
        if ftype in ["arc", "tilt"]:
            return (
                good_exp
                & (
                    (fitstbl["idname"] == "SCIENCE")
                    # | (fitstbl["idname"] == "SCIENCE_EXTENDED")  # Probably not good for arc lines
                    # | (fitstbl["idname"] == "TEST")  # Don't want test exposures
                    | (fitstbl["idname"] == "SCIENCE_ON")
                    | (fitstbl["idname"] == "SCIENCE_OFF")
                    | (fitstbl["idname"] == "COMPARISON_SKY")
                    | (fitstbl["idname"] == "COMPARISON_LAMP")
                )
                & (fitstbl["filter1"] != "blank")
            )

        if ftype in ["trace", "pixelflat", "illumflat"]:
            return (
                good_exp
                & self._match_flat_exptime(fitstbl)
                & (fitstbl["idname"] == "DOME_FLAT")
                & (fitstbl["filter1"] != "blank")
                # & (fitstbl["lampstat01"] == "off")
            )

        if ftype == "lampoffflats":
            return (
                good_exp
                & self._match_flat_exptime(fitstbl)
                & (fitstbl["idname"] == "DOME_BACKGROUND")
                & (fitstbl["filter1"] != "blank")
            )

        if ftype == "science":
            return (
                good_exp
                & (
                    (fitstbl["idname"] == "SCIENCE")
                    # | (fitstbl["idname"] == "TEST")  # Don't want test exposures
                    | (fitstbl["idname"] == "SCIENCE_EXTENDED")
                    | (fitstbl["idname"] == "SCIENCE_CROWDED")
                    | (fitstbl["idname"] == "SCIENCE_ON")
                    | (fitstbl["idname"] == "SCIENCE_OFF")
                )
                # & (fitstbl["lampstat01"] == "off")
            )

        if ftype == "standard":
            return (
                good_exp
                & ((fitstbl["idname"] == "STANDARD") | (fitstbl["idname"] == "A0V"))
                # & (fitstbl["lampstat01"] == "off")
            )

        if ftype == "dark":
            return (
                good_exp
                # It's a dark iff the AUX blank is in place
                & (fitstbl["filter1"] == "blank")
                # & (fitstbl["lampstat01"] == "off")
            )

        if ftype in [
            "bias",
            "pinhole",
            "align",
            "sky",
            "scattlight",
            "slitless_pixflat",
        ]:
            # DeVeny doesn't have any of these types of frames
            return np.zeros(len(fitstbl), dtype=bool)
        log.warning("Cannot determine if frames are of type %s", ftype)
        return np.zeros(len(fitstbl), dtype=bool)

    def get_comb_group(self, fitstbl: astropy.table.Table) -> astropy.table.Table:
        """
        Automatically assign combination groups and background images by parsing
        known dither patterns for LDT/RIMAS.

        This method is used in
        :func:`~pypeit.metadata.PypeItMetaData.set_combination_groups`, and
        directly modifies the ``calib``, ``comb_id``, and ``bkg_id`` columns in the
        provided table.

        Specifically here for RIMAS, science and standard frames are grouped
        together and assigned calibration groups according to their exposure times.
        Dark frames are assigned to the matching science/standard calibration
        group, while ``lampoffflats`` and ``pixelflat,illumflat,trace`` frames are
        assigned ``calib = "all"``. Any remaining frame types are assigned
        sequential calibration groups after the science/standard and dark frames.

        Moreover, this method parses from the header the dither pattern of the
        science/standard frames in a given calibration group and assigns to each of
        them a default ``comb_id`` and ``bkg_id``. The dither patterns used here
        are: ``ABBA``, ``OnOff``, and ``None``. For ``ABBA``, all A frames in
        one sequence share a ``comb_id`` and use the shared B ``comb_id`` as
        their ``bkg_id``; the B frames are assigned conversely. This produces
        the AA-BB and BB-AA nod-subtracted products. ``OnOff`` frames retain
        their own ``comb_id`` and use the nearest-in-time frame at the opposite
        dither position as their background. The ``comb_id`` and ``bkg_id``
        will *not* be assigned if the dither pattern recorded in the header is
        not recognized or is set to ``NONE`` or ``MANUAL``.

        If ``rimas_dither_corrections.ecsv`` is present in the setup working
        directory, it overrides ABBA or On/Off metadata for listed files. Raw-data
        directories are searched only when the working-directory sidecar is
        absent. Corrected and uncorrected sequences are never mixed.

        Parameters
        ----------
        fitstbl : `astropy.table.Table`_
            The table with the metadata for all the frames.

        Returns
        -------
        `astropy.table.Table`_
            Modified metadata table.
        """
        # Find index of fitstbl that contains science and standard frames
        sci_idx = np.array(
            ["science" in _tab for _tab in fitstbl["frametype"]]
        ) | np.array(["standard" in _tab for _tab in fitstbl["frametype"]])

        corrections = self._load_dither_corrections(fitstbl)
        corrected = np.zeros(len(fitstbl), dtype=bool)
        corrected_seq = np.full(len(fitstbl), -1, dtype=int)
        for i, filename in enumerate(fitstbl["filename"]):
            basename = pathlib.Path(str(filename)).name
            if basename not in corrections:
                continue

            pattern, dither_position, dithseq = corrections[basename]
            fitstbl["dithpat"][i] = pattern
            fitstbl["dithpos"][i] = dither_position
            corrected[i] = True
            corrected_seq[i] = dithseq

        if "calib" in fitstbl.keys():
            # Find index of fitstbl that contains lampoffflats, pixelflat,
            # illumflat, or trace frames
            flat_idx = (
                np.array(["lampoffflats" in _tab for _tab in fitstbl["frametype"]])
                | np.array(["pixelflat" in _tab for _tab in fitstbl["frametype"]])
                | np.array(["illumflat" in _tab for _tab in fitstbl["frametype"]])
                | np.array(["trace" in _tab for _tab in fitstbl["frametype"]])
            )
            # Set calib for those frames to "all"
            fitstbl["calib"][flat_idx] = "all"

            # Find index of fitstbl that contains dark frames
            dark_idx = np.array(["dark" in _tab for _tab in fitstbl["frametype"]])

            # Initialize calib_id / group by exposure times
            calib_id = 0
            exptime_to_calib = {}

            # Assign calib groups for science/standard frames by unique exptime
            sci_exptimes = np.unique(fitstbl[sci_idx]["exptime"])
            for exptime in sci_exptimes:
                this_sci = sci_idx & np.isclose(fitstbl["exptime"], exptime, atol=0.1)
                fitstbl["calib"][this_sci] = f"{calib_id}"
                exptime_to_calib[exptime] = calib_id
                calib_id += 1

            # Assign darks to the matching science/standard exptime calibration group
            dark_exptimes = np.unique(fitstbl[dark_idx]["exptime"])
            for exptime in dark_exptimes:
                this_dark = dark_idx & np.isclose(fitstbl["exptime"], exptime, atol=0.1)

                matched = None
                for sci_exptime, calib_id in exptime_to_calib.items():
                    if np.isclose(exptime, sci_exptime, atol=0.1):
                        matched = calib_id
                        break

                if matched is not None:
                    fitstbl["calib"][this_dark] = f"{matched}"
                else:
                    # Assign the next sequential ID to these frames
                    fitstbl["calib"][this_dark] = f"{calib_id}"
                    calib_id += 1

            # Assign calib groups to any remaining frame types
            assigned_idx = flat_idx | sci_idx | dark_idx
            remaining_frametypes = np.unique(fitstbl["frametype"][~assigned_idx])
            for ftype in remaining_frametypes:
                this_ftype = fitstbl["frametype"] == ftype
                fitstbl["calib"][this_ftype] = str(calib_id)
                calib_id += 1

        # Loop over the setups (generally there is only one setup, but we check anyway)
        for setup in np.unique(fitstbl[sci_idx]["setup"]):
            in_setup = sci_idx & (fitstbl["setup"] == setup)

            for target in np.unique(fitstbl[in_setup]["target"]):
                targ_idx = in_setup & (fitstbl["target"] == target)

                if np.sum(targ_idx) <= 1:
                    continue

                for dpat in np.unique(fitstbl[targ_idx]["dithpat"]):

                    if dpat.upper() in ["NONE", "MANUAL"]:
                        continue

                    dpat_idx = targ_idx & (fitstbl["dithpat"] == dpat)
                    if np.sum(dpat_idx) <= 1:
                        continue

                    combid = fitstbl["comb_id"][dpat_idx].data.copy()
                    bkgid = fitstbl["bkg_id"][dpat_idx].data.copy()
                    dpos = fitstbl["dithpos"][dpat_idx].data
                    mjd = fitstbl["mjd"][dpat_idx].data
                    is_corrected = corrected[dpat_idx]
                    dither_sequence = corrected_seq[dpat_idx]

                    if dpat == "ABBA":
                        # Sidecar metadata identifies distinct commanded
                        # sequences. Without it, use all matching header-based
                        # frames in this setup/target/pattern group.
                        sequence_masks = [~is_corrected]
                        sequence_masks += [
                            is_corrected & (dither_sequence == sequence)
                            for sequence in np.unique(dither_sequence[is_corrected])
                        ]
                        for sequence_mask in sequence_masks:
                            is_a = sequence_mask & (dpos == "A")
                            is_b = sequence_mask & (dpos == "B")
                            if not np.any(is_a) or not np.any(is_b):
                                continue

                            a_combid = np.min(combid[is_a])
                            b_combid = np.min(combid[is_b])
                            combid[is_a], bkgid[is_a] = a_combid, b_combid
                            combid[is_b], bkgid[is_b] = b_combid, a_combid

                    elif dpat == "OnOff":
                        for i in range(len(combid)):
                            opp = (
                                "Off"
                                if dpos[i] == "On"
                                else "On" if dpos[i] == "Off" else None
                            )
                            if opp is None:
                                continue

                            match = dpos == opp
                            if is_corrected[i]:
                                match &= is_corrected & (
                                    dither_sequence == dither_sequence[i]
                                )
                            else:
                                match &= ~is_corrected

                            if not np.any(match):
                                continue

                            j = np.argmin(np.abs(mjd[match] - mjd[i]))
                            bkgid[i] = combid[match][j]

                    fitstbl["bkg_id"][dpat_idx] = bkgid
                    fitstbl["comb_id"][dpat_idx] = combid

        return fitstbl

    def get_rawimage(
        self,
        raw_file: str | pathlib.Path,
        det: int,
        sec_includes_binning: bool = False,
    ) -> tuple[
        detector_container.DetectorContainer,
        np.ndarray,
        astropy.io.fits.HDUList,
        float,
        np.ndarray,
        np.ndarray,
    ]:
        """
        Read raw spectrograph image files and return data and relevant metadata
        needed for image processing.

        For LDT/RIMAS, we need to convert NaN pixels in the raw frames to
        finite saturated values.

        Parameters
        ----------
        raw_file : str or pathlib.Path
            File to read.
        det : int
            1-indexed detector number.
        sec_includes_binning : bool, optional
            Interpret data and overscan sections as already including the
            detector binning.

        Returns
        -------
        tuple
            Detector parameters, raw image, opened HDU list, exposure time, data
            section image, and overscan section image.
        """
        # Call the super()
        (
            detector_par,
            raw_img,
            hdu,
            exptime,
            rawdatasec_img,
            oscansec_img,
        ) = super().get_rawimage(
            raw_file, det, sec_includes_binning=sec_includes_binning
        )

        # Get the locations of NaN pixels & replace with saturated value
        nan_idx = np.logical_not(np.isfinite(raw_img))
        raw_img[nan_idx] = detector_par.saturation

        # Return the mess
        return detector_par, raw_img, hdu, exptime, rawdatasec_img, oscansec_img

    @staticmethod
    def scrub_isot_dateobs(dt_str: str) -> astropy.time.Time:
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
        dt_str : str
            Input datetime string from the ``DATE-OBS`` header keyword.

        Returns
        -------
        `astropy.time.Time`_
            The AstroPy Time object corresponding to the ``DATE-OBS`` input
            string.
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

    # Internal Helper Methods ============================#
    @staticmethod
    def _mode_compatible_darks(
        fitstbl: astropy.table.Table, mode_idx: np.ndarray
    ) -> np.ndarray:
        """
        Select dark frames compatible with data already assigned to a mode.

        A dark is compatible when its detector arm and raw readout geometry
        match at least one non-dark frame selected by ``mode_idx``.  This
        permits dark frames to be retained even though their optical metadata
        do not identify the science mode in which they will be used.

        Parameters
        ----------
        fitstbl : `astropy.table.Table`_
            Metadata table.  The table must contain the ``arm``, ``rawshape``,
            and ``filter1`` columns.
        mode_idx : `numpy.ndarray`_
            Boolean array selecting the non-dark frames assigned to the mode.
            Its length must match ``fitstbl``.

        Returns
        -------
        `numpy.ndarray`_
            Boolean array selecting dark frames whose ``(arm, rawshape)`` pair
            occurs among the mode frames.
        """
        # Optical metadata do not reliably identify a dark's mode, so compare
        # detector properties that remain meaningful for blank-filter frames.
        mode_geometry = set(
            zip(fitstbl["arm"][mode_idx], fitstbl["rawshape"][mode_idx])
        )
        is_dark = fitstbl["filter1"] == "blank"
        geometry_matches = np.array(
            [
                (arm, rawshape) in mode_geometry
                for arm, rawshape in zip(fitstbl["arm"], fitstbl["rawshape"])
            ],
            dtype=bool,
        )
        return is_dark & geometry_matches

    @staticmethod
    def _match_flat_exptime(fitstbl: astropy.table.Table) -> np.ndarray:
        """
        Select RIMAS dome flats with arm/mode-specific exposure times.

        RIMAS dome flats for the YJ channel are intentionally longer than
        those for the HK channel.  This selection has to happen during frame
        typing, before :meth:`config_specific_par` is called, so it uses the
        metadata table directly instead of the calibration frame ``exprng``
        parameters.

        Parameters
        ----------
        fitstbl : `astropy.table.Table`_
            Metadata table.  The table must contain the ``dispname``,
            ``decker``, ``arm``, and ``exptime`` columns.

        Returns
        -------
        `numpy.ndarray`_
            Boolean array selecting exposures within the configured flat-field
            exposure-time range for their disperser, slit, and detector arm.
        """
        good_exp = np.zeros(len(fitstbl), dtype=bool)

        # Evaluate each supported configuration independently because its
        # acceptable exposure-time interval can differ between detector arms.
        for (dispname, decker, arm), exprng in RIMAS_FLAT_EXPTIME_RANGES.items():
            in_config = (
                (fitstbl["dispname"] == dispname)
                & (fitstbl["decker"] == decker)
                & (fitstbl["arm"] == arm)
            )

            # Avoid passing empty table slices to the generic frame matcher.
            if np.any(in_config):
                good_exp[in_config] = framematch.check_frame_exptime(
                    fitstbl["exptime"][in_config], exprng
                )

        return good_exp

    @staticmethod
    def _load_dither_corrections(
        fitstbl: astropy.table.Table,
    ) -> dict[str, tuple[str, int, int]]:
        """Load optional observer corrections to RIMAS dither metadata.

        The setup working directory is searched first.  If it does not contain
        the conventionally named sidecar, each unique raw-data directory in
        ``fitstbl`` is searched instead.

        Parameters
        ----------
        fitstbl : `astropy.table.Table`_
            Metadata table containing at least ``filename`` and, for fallback
            discovery, ``directory``.

        Returns
        -------
        dict
            Mapping from raw-file basename to normalized ``(dithpat, dithpos,
            dithseq)`` tuples.  The returned pattern is either ``'ABBA'`` or
            ``'OnOff'``.

        Raises
        ------
        PypeItError
            Raised if a sidecar cannot be read, does not have the exact schema
            produced by :command:`rimas_dither_pattern`, contains duplicate
            filenames, or contains invalid or internally inconsistent dither
            metadata.

        Notes
        -----
        Only a sidecar in the current working directory is used when present;
        raw-data directories are searched only as a fallback.  Corrections for
        filenames absent from ``fitstbl`` are retained in the returned mapping
        and generate a warning.
        """
        # A setup-directory sidecar represents the observer's explicit choice
        # and therefore suppresses discovery beside the raw files.
        cwd_sidecar = pathlib.Path.cwd() / RIMAS_DITHER_CORRECTION_FILE
        if cwd_sidecar.is_file():
            sidecars = [cwd_sidecar]
        else:
            sidecars = []
            if "directory" in fitstbl.colnames:
                raw_directories = {
                    pathlib.Path(str(directory)).expanduser()
                    for directory in fitstbl["directory"]
                }
                sidecars = sorted(
                    path / RIMAS_DITHER_CORRECTION_FILE
                    for path in raw_directories
                    if (path / RIMAS_DITHER_CORRECTION_FILE).is_file()
                )

        # Keep the consumer synchronized with the fixed output contract of the
        # obstools ``rimas_dither_pattern`` command.
        expected_columns = [
            "filename",
            "dithpat",
            "dithnum",
            "dithpos",
            "dithseq",
            "source",
        ]
        patterns = {
            "ABBA": ("ABBA", ("A", "B", "B", "A")),
            "ONOFF": ("OnOff", ("On", "Off")),
        }

        corrections = {}
        for sidecar in sidecars:
            try:
                table = astropy.table.Table.read(sidecar, format="ascii.ecsv")
            except Exception as exc:
                raise PypeItError(
                    f"Unable to read RIMAS dither correction file {sidecar}: {exc}"
                ) from exc

            if table.colnames != expected_columns:
                raise PypeItError(
                    f"RIMAS dither correction file {sidecar} must contain exactly "
                    f"these columns: {', '.join(expected_columns)}"
                )

            for row in table:
                filename = pathlib.Path(str(row["filename"])).name
                # Basenames are the cross-directory identity used by fitstbl.
                if filename in corrections:
                    raise PypeItError(
                        f"Duplicate RIMAS dither correction for {filename}"
                    )

                # Derive the position from pattern membership and verify that
                # the redundant logged position agrees with it.
                try:
                    pattern, positions = patterns[str(row["dithpat"]).strip().upper()]
                    dithnum = int(row["dithnum"])
                    dithseq = int(row["dithseq"])
                    if dithnum < 1 or dithnum > len(positions) or dithseq < 1:
                        raise ValueError
                    expected_dithpos = positions[dithnum - 1]
                except (KeyError, TypeError, ValueError) as exc:
                    raise PypeItError(
                        f"Invalid RIMAS dither correction for {filename}"
                    ) from exc
                if str(row["dithpos"]).strip() != expected_dithpos:
                    raise PypeItError(
                        f"Invalid RIMAS dithpos {row['dithpos']!r} for {filename}; "
                        f"dithpat={pattern} and dithnum={dithnum} require {expected_dithpos}"
                    )

                corrections[filename] = (pattern, expected_dithpos, dithseq)

        table_filenames = {pathlib.Path(str(name)).name for name in fitstbl["filename"]}
        # Stale entries are harmless, but warn because they often indicate that
        # the wrong setup directory or observing log was selected.
        for filename in corrections.keys() - table_filenames:
            log.warning(
                "RIMAS dither correction for %s does not match a file in the metadata table",
                filename,
            )

        return corrections


# VPH Child Class for LDT/RIMAS ==============================================#
class LDTRIMASVphSpectrograph(LDTRIMASSpectrograph):
    """
    Child to handle common aspects of the LDT/RIMAS VPH (low-res) modes
    """

    name = "ldt_rimas_vph"
    comment = "LDT Rapid infrared IMAger Spectrometer, VPH Gratings"
    pypeline = "MultiSlit"

    def validate_fitstbl(self, fitstbl: astropy.table.Table) -> astropy.table.Table:
        """Validate the metadata table

        Because of the multiple arms and modes of RIMAS, this method removes
        from the metadata table frames not associated with this mode,
        namely removes frames that are not low-res (Vph) gratings with a slit.

        Parameters
        ----------
        fitstbl : `astropy.table.Table`_
            The metadata table to be validated.

        Returns
        -------
        `astropy.table.Table`_
            The validated metadata table.
        """
        is_dark = fitstbl["filter1"] == "blank"

        # Only keep non-dark frames with one of the VPH gratings -- no slitless mode!
        vph_idx = (
            ((fitstbl["dispname"] == "Vph300") | (fitstbl["dispname"] == "Vph30"))
            & (fitstbl["decker"] != "open")
            & np.logical_not(is_dark)
        )

        # Include darks taken in other optical configurations only when their
        # detector arm and raw readout geometry match data selected above.
        dark_idx = self._mode_compatible_darks(fitstbl, vph_idx)

        # Return the corrected table
        return fitstbl[vph_idx | dark_idx]

    @staticmethod
    def configuration_list() -> dict[str, dict[str, np.str_]]:
        """
        Return the default-ordered list of VPH configurations for this spectrograph

        Returns
        -------
        dict
            Dictionary of default-ordered spectrograph configurations.
        """
        return {
            # All of the VPH30 possible combinations first
            "30YJ12l": {
                "arm": np.str_("YJ"),
                "dispname": np.str_("Vph30"),
                "decker": np.str_("1.2'' long"),
            },
            "30HK12l": {
                "arm": np.str_("HK"),
                "dispname": np.str_("Vph30"),
                "decker": np.str_("1.2'' long"),
            },
            "30YJ06": {
                "arm": np.str_("YJ"),
                "dispname": np.str_("Vph30"),
                "decker": np.str_("0.6''"),
            },
            "30HK06": {
                "arm": np.str_("HK"),
                "dispname": np.str_("Vph30"),
                "decker": np.str_("0.6''"),
            },
            "30YJ10": {
                "arm": np.str_("YJ"),
                "dispname": np.str_("Vph30"),
                "decker": np.str_("1.0''"),
            },
            "30HK10": {
                "arm": np.str_("HK"),
                "dispname": np.str_("Vph30"),
                "decker": np.str_("1.0''"),
            },
            "30YJ20": {
                "arm": np.str_("YJ"),
                "dispname": np.str_("Vph30"),
                "decker": np.str_("2.0''"),
            },
            "30HK20": {
                "arm": np.str_("HK"),
                "dispname": np.str_("Vph30"),
                "decker": np.str_("2.0''"),
            },
            # All of the VPH300 possible combinations next
            "300YJ12l": {
                "arm": np.str_("YJ"),
                "dispname": np.str_("Vph300"),
                "decker": np.str_("1.2'' long"),
            },
            "300HK12l": {
                "arm": np.str_("HK"),
                "dispname": np.str_("Vph300"),
                "decker": np.str_("1.2'' long"),
            },
            "300YJ06": {
                "arm": np.str_("YJ"),
                "dispname": np.str_("Vph300"),
                "decker": np.str_("0.6''"),
            },
            "300HK06": {
                "arm": np.str_("HK"),
                "dispname": np.str_("Vph300"),
                "decker": np.str_("0.6''"),
            },
            "300YJ10": {
                "arm": np.str_("YJ"),
                "dispname": np.str_("Vph300"),
                "decker": np.str_("1.0''"),
            },
            "300HK10": {
                "arm": np.str_("HK"),
                "dispname": np.str_("Vph300"),
                "decker": np.str_("1.0''"),
            },
            "300YJ20": {
                "arm": np.str_("YJ"),
                "dispname": np.str_("Vph300"),
                "decker": np.str_("2.0''"),
            },
            "300HK20": {
                "arm": np.str_("HK"),
                "dispname": np.str_("Vph300"),
                "decker": np.str_("2.0''"),
            },
        }

    @classmethod
    def default_pypeit_par(cls) -> pypeitpar.PypeItPar:
        """
        Return the default parameters to use for this instrument.

        Returns
        -------
        :class:`~pypeit.par.pypeitpar.PypeItPar`
            Parameters required by all PypeIt methods.
        """
        # Get the PypeIt and RIMAS-wide default parameters
        par = super().default_pypeit_par()

        # For processing the arc frame, these settings allow for the combination of
        #   of frames from different lamps into a comprehensible Master
        par["calibrations"]["arcframe"]["process"]["subtract_continuum"] = False
        # par["calibrations"]["tiltframe"]["process"]["clip"] = False
        # par["calibrations"]["tiltframe"]["process"]["combine"] = "mean"
        par["calibrations"]["tiltframe"]["process"]["subtract_continuum"] = False

        # Wavelengths
        par["calibrations"]["wavelengths"]["rms_thresh_frac_fwhm"] = 0.4
        par["calibrations"]["wavelengths"]["sigdetect"] = 5.0
        par["calibrations"]["wavelengths"]["lamps"] = ["OH_GNIRS"]
        # par['calibrations']['wavelengths']['nonlinear_counts'] =
        #      self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par["calibrations"]["wavelengths"]["n_first"] = 1
        par["calibrations"]["wavelengths"]["n_final"] = 1

        # # Reidentification parameters
        # par["calibrations"]["wavelengths"]["method"] = "reidentify"
        # par["calibrations"]["wavelengths"]["cc_thresh"] = 0.6
        # par["calibrations"]["wavelengths"]["reid_arxiv"] = "gemini_gnirs.fits"

        # Tilts
        par["calibrations"]["tilts"]["tracethresh"] = 10
        par["calibrations"]["tilts"]["sig_neigh"] = 5.0
        par["calibrations"]["tilts"]["nfwhm_neigh"] = 2.0

        par["sensfunc"]["IR"]["pix_shift_bounds"] = (-5.0, 5.0)

        # # Wavelength Calibration Parameters
        # # Set this as default... but use `holy-grail` for DV4, DV8
        # par["calibrations"]["wavelengths"][
        #     "method"
        # ] = "full_template"  # Default: 'holy-grail'
        # # The DeVeny arc line FWHM varies based on slitwidth used
        # par["calibrations"]["wavelengths"]["fwhm"] = 3.0  # Default: 4.0
        # par["calibrations"]["wavelengths"]["nsnippet"] = 1  # Default: 2

        # # For the tilts, our lines are not as well-behaved as others',
        # #   possibly due to the Wynne version E camera.
        # par["calibrations"]["tilts"]["spat_order"] = 4  # Default: 3
        # par["calibrations"]["tilts"]["spec_order"] = 5  # Default: 4

        # # Cosmic ray rejection parameters for science frames
        # par["scienceframe"]["process"]["sigclip"] = 5.0  # Default: 4.5
        # par["scienceframe"]["process"]["objlim"] = 2.0  # Default: 3.0

        # # Object Finding, Extraction, and Sky Subtraction Parameters
        # assumed_seeing = 1.5  # arcsec
        # par["reduce"]["findobj"]["trace_npoly"] = 3  # Default: 5
        # par["reduce"]["findobj"]["snr_thresh"] = 50.0  # Default: 10.0
        # par["reduce"]["findobj"]["maxnumber_std"] = 1  # Default: 5
        # par["reduce"]["findobj"]["maxnumber_sci"] = 5  # Default: 10
        # par["reduce"]["findobj"]["find_fwhm"] = np.round(
        #     assumed_seeing / 0.34, 1
        # )  # Default: 5.0 pix
        # par["reduce"]["findobj"]["find_trim_edge"] = [0, 0]  # Default: [5, 5]
        # # Boxcar width = ±3σ of Gaussian profile = >99% enclosed flux; radius = 1.28 * seeing
        # par["reduce"]["extraction"]["boxcar_radius"] = np.round(
        #     assumed_seeing * 1.28, 1
        # )  # Default: 1.5"
        # par["reduce"]["extraction"]["use_2dmodel_mask"] = False  # Default: True
        # par["reduce"]["skysub"]["sky_sigrej"] = 4.0  # Default: 3.0

        # # Flexure Correction Parameters
        # par["flexure"]["spec_method"] = "boxcar"  # Default: 'skip'
        # par["flexure"]["spec_maxshift"] = 30  # Default: 20

        # # Slit-edge settings for long-slit data (DeVeny's slit is > 90" long)
        # par["calibrations"]["slitedges"]["bound_detector"] = True  # Defualt: False
        # par["calibrations"]["slitedges"]["sync_predict"] = "nearest"  # Default: 'pca'
        # par["calibrations"]["slitedges"]["minimum_slit_length"] = 170.0  # Default: None
        # par["calibrations"]["slitedges"]["max_nudge"] = 5  # Default: None

        # # Adjustments to slit and tilts for NIR
        # par["calibrations"]["slitedges"]["edge_thresh"] = 20.0  # Default: 20.0
        # par["calibrations"]["slitedges"]["fit_order"] = 2  # Default: 5
        # # par["calibrations"]["slitedges"]["max_nudge"] = 5  # Default: None
        # par["calibrations"]["slitedges"]["max_shift_adj"] = 0.5
        # par["calibrations"]["slitedges"]["trace_thresh"] = 20.0
        # # par["calibrations"]["slitedges"]["fit_min_spec_length"] = 0.5
        # # par["calibrations"]["slitedges"]["left_right_pca"] = True
        # # par["calibrations"]["slitedges"]["smash_range"] = [
        # #     0.25,
        # #     0.75,
        # # ]  # Default: [0.0, 1.0]
        # # par["calibrations"]["slitedges"]["sync_predict"] = "nearest"  # Default: 'pca'
        # par["calibrations"]["slitedges"]["sobel_enhance"] = 3  # Default: 0
        # par["calibrations"]["slitedges"]["trim_spec"] = [1024, 1024]

        # # Only use LONG arc frames
        # # par["calibrations"]["arcframe"]["exprng"] = [30, None]
        # # For processing the arc frame, these settings allow for the combination of
        # #   of frames from different lamps into a comprehensible Master
        # par["calibrations"]["arcframe"]["process"]["clip"] = False
        # par["calibrations"]["arcframe"]["process"]["combine"] = "mean"
        # # par['calibrations']['arcframe']['process']['subtract_continuum'] = True
        # par["calibrations"]["tiltframe"]["process"]["clip"] = False
        # par["calibrations"]["tiltframe"]["process"]["combine"] = "mean"
        # # par['calibrations']['tiltframe']['process']['subtract_continuum'] = True

        # # Wavelength Calibration Parameters
        # # Arc lamps list from header -- instead of defining the full list here
        # # Set this as default... but use `holy-grail` for DV4, DV8
        # par["calibrations"]["wavelengths"][
        #     "method"
        # ] = "holy-grail"  #'full_template'  # Default: 'holy-grail'
        # # Reidentification parameters
        # par["calibrations"]["wavelengths"]["reid_arxiv"] = "ldt_nihts.fits"
        # # The DeVeny arc line FWHM varies based on slitwidth used
        # par["calibrations"]["wavelengths"]["fwhm_fromlines"] = True  # Default: True
        # par["calibrations"]["wavelengths"]["nsnippet"] = 1  # Default: 2

        # # # For the tilts, our lines are not as well-behaved as others',
        # # #   possibly due to the Wynne version E camera.
        # # par["calibrations"]["tilts"]["spat_order"] = 4  # Default: 3
        # # par["calibrations"]["tilts"]["spec_order"] = 5  # Default: 4

        # # Flat-field parameter modification
        # par["calibrations"]["flatfield"]["pixelflat_min_wave"] = 3000.0  # Default: None
        # par["calibrations"]["flatfield"]["slit_illum_finecorr"] = False  # Default: True
        # par["calibrations"]["flatfield"]["spec_samp_fine"] = 30  # Default: 1.2
        # par["calibrations"]["flatfield"]["tweak_slits"] = False  # Default: True

        # # Cosmic ray rejection parameters for science frames
        # par["scienceframe"]["process"]["sigclip"] = 5.0  # Default: 4.5
        # par["scienceframe"]["process"]["objlim"] = 2.0  # Default: 3.0

        # # Object Finding, Extraction, and Sky Subtraction Parameters
        # assumed_seeing = 1.5  # arcsec
        # par["reduce"]["findobj"]["trace_npoly"] = 3  # Default: 5
        # par["reduce"]["findobj"]["snr_thresh"] = 50.0  # Default: 10.0
        # par["reduce"]["findobj"]["maxnumber_std"] = 1  # Default: 5
        # par["reduce"]["findobj"]["maxnumber_sci"] = 5  # Default: 10
        # par["reduce"]["findobj"]["find_fwhm"] = np.round(
        #     assumed_seeing / 0.34, 1
        # )  # Default: 5.0 pix
        # par["reduce"]["findobj"]["find_trim_edge"] = [0, 0]  # Default: [5, 5]
        # # Boxcar width = ±3σ of Gaussian profile = >99% enclosed flux; radius = 1.28 * seeing
        # par["reduce"]["extraction"]["boxcar_radius"] = np.round(
        #     assumed_seeing * 1.28, 1
        # )  # Default: 1.5"
        # par["reduce"]["extraction"]["use_2dmodel_mask"] = False  # Default: True
        # par["reduce"]["skysub"]["sky_sigrej"] = 4.0  # Default: 3.0

        # # Flexure Correction Parameters
        # par["flexure"]["spec_method"] = "boxcar"  # Default: 'skip'
        # par["flexure"]["spec_maxshift"] = 30  # Default: 20

        return par

    def config_specific_par(
        self,
        inp: str | list | pathlib.Path | astropy.io.fits.Header | astropy.table.Table,
        inp_par: parset.ParSet = None,
    ) -> parset.ParSet:
        """
        Modify the PypeIt parameters to hard-wired values used for
        specific instrument configurations.

        Parameters
        ----------
        inp : str, list, pathlib.Path, `astropy.io.fits.Header`_, or `astropy.table.Table`_
            Input filename, an `astropy.io.fits.Header`_ object, a list of
            `astropy.io.fits.Header`_ objects, or a row from the metadata
            table.
        inp_par : :class:`~pypeit.par.parset.ParSet`, optional
            Parameter set used for the full run of PypeIt. If None, use
            :func:`default_pypeit_par`.

        Returns
        -------
        :class:`~pypeit.par.parset.ParSet`
            The PypeIt parameter set adjusted for configuration-specific
            parameter values.
        """
        # Start with instrument-wide parameters
        par = super().config_specific_par(inp, inp_par=inp_par)

        # Adjust parameters based on instrument settings
        grating = self.get_meta_value(inp, "dispname")
        binning = self.get_meta_value(inp, "binning")

        # Check for the 1.2" long slit... edge tracing is the same for both arms and both gratings
        if "long" in self.decker:
            # Slit-edge settings for long-slit data (75" long)
            par["calibrations"]["slitedges"]["bound_detector"] = True
            par["calibrations"]["slitedges"]["sync_predict"] = "nearest"
            par["calibrations"]["slitedges"]["minimum_slit_length"] = 70.0  # arcsec
            par["calibrations"]["slitedges"]["max_nudge"] = 5  # Default: None
            par["calibrations"]["slitedges"]["sobel_enhance"] = 3
        else:
            log.warning("We have not set slit edge stuff for the short slits / VPH!")

        # Get the arm-specific parameters based on grating and decker
        if self.arm == "YJ":
            par = self.config_specific_par_vph_yj(par, grating, self.decker)
        else:
            par = self.config_specific_par_vph_hk(par, grating, self.decker)

        return self.adjust_config_specific_binning(par, binning)

    def config_specific_par_vph_yj(
        self, par: parset.ParSet, grating: str, decker: str
    ) -> parset.ParSet:
        """Set YJ arm configuration-specific parameters for VPH gratings

        Parameters
        ----------
        par : :class:`~pypeit.par.parset.ParSet`
            The instrument-wide parameter set to be modified.
        grating : str
            The grating used (from :meth:`get_meta_value`).
        decker : str
            The slit / decker used (from :meth:`get_meta_value`).

        Returns
        -------
        :class:`~pypeit.par.parset.ParSet`
            Modified parameter set for the YJ arm / Vph gratings.
        """
        par["calibrations"]["slitedges"]["edge_thresh"] = 20.0
        par["calibrations"]["slitedges"]["trace_thresh"] = 20.0
        par["calibrations"]["slitedges"]["sobel_enhance"] = 3
        par["calibrations"]["slitedges"]["trim_spec"] = [1024, 1024]
        par["calibrations"]["wavelengths"]["lamps"] = ["XeI"]
        par["calibrations"]["wavelengths"]["method"] = "holy-grail"
        par["calibrations"]["wavelengths"]["reid_arxiv"] = "ldt_nihts.fits"
        par["calibrations"]["wavelengths"]["fwhm_fromlines"] = True
        par["calibrations"]["wavelengths"]["nsnippet"] = 1
        par["reduce"]["findobj"]["snr_thresh"] = 50.0

        if grating == "Vph30":
            par["calibrations"]["wavelengths"]["xcorr_only"] = True
            par["calibrations"]["wavelengths"][
                "reid_arxiv"
            ] = "ldt_deveny_150_HgCdAr.fits"
            # Because of the wide wavelength range, split Vph30 arcs in half.
            par["calibrations"]["wavelengths"]["nsnippet"] = 2
            par["calibrations"]["wavelengths"]["n_first"] = 3
            par["calibrations"]["wavelengths"]["n_final"] = 5
            par["sensfunc"]["IR"]["resln_guess"] = 30

        elif grating == "Vph300":
            par["calibrations"]["wavelengths"][
                "reid_arxiv"
            ] = "ldt_deveny_300_HgCdAr.fits"
            par["calibrations"]["wavelengths"]["n_first"] = 3
            par["calibrations"]["wavelengths"]["n_final"] = 5
            par["sensfunc"]["IR"]["resln_guess"] = 300

        else:
            raise ValueError(f"Grating {grating} not recognized for RIMAS VPH modes")

        return par

    def config_specific_par_vph_hk(
        self, par: parset.ParSet, grating: str, decker: str
    ) -> parset.ParSet:
        """Set HK arm configuration-specific parameters for VPH gratings

        Parameters
        ----------
        par : :class:`~pypeit.par.parset.ParSet`
            The instrument-wide parameter set to be modified.
        grating : str
            The grating used (from :meth:`get_meta_value`).
        decker : str
            The slit / decker used (from :meth:`get_meta_value`).

        Returns
        -------
        :class:`~pypeit.par.parset.ParSet`
            Modified parameter set for the HK arm / Vph gratings.
        """
        par["calibrations"]["slitedges"]["edge_thresh"] = 20.0
        par["calibrations"]["slitedges"]["fit_min_spec_length"] = 0.4
        par["calibrations"]["slitedges"]["det_min_spec_length"] = 0.4
        par["calibrations"]["slitedges"]["sync_predict"] = "nearest"
        par["calibrations"]["bpm_usebias"] = True
        par["calibrations"]["wavelengths"]["method"] = "full_template"
        par["calibrations"]["wavelengths"]["fwhm"] = 5.0
        par["calibrations"]["wavelengths"]["nsnippet"] = 1
        par["calibrations"]["wavelengths"]["sigdetect"] = 1.0
        par["calibrations"]["wavelengths"]["echelle"] = False
        par["calibrations"]["tilts"]["spat_order"] = 4
        par["calibrations"]["tilts"]["spec_order"] = 5
        par["reduce"]["findobj"]["snr_thresh"] = 10.0

        if grating == "Vph30":
            par["calibrations"]["wavelengths"]["xcorr_only"] = True
            par["calibrations"]["slitedges"]["edge_thresh"] = 50.0
            par["calibrations"]["slitedges"]["fit_min_spec_length"] = 0.05
            par["calibrations"]["slitedges"]["det_min_spec_length"] = 0.05
            par["calibrations"]["wavelengths"]["reid_arxiv"] = "ldt_rimas_HK_30_OH.fits"
            # Because of the wide wavelength range, split Vph30 arcs in half.
            par["calibrations"]["wavelengths"]["nsnippet"] = 2
            par["calibrations"]["wavelengths"]["n_first"] = 3
            par["calibrations"]["wavelengths"]["n_final"] = 5
            par["calibrations"]["wavelengths"]["lamps"] = ["OH_RIMAS_HK"]

            par["reduce"]["findobj"]["find_fwhm"] = 5
            par["reduce"]["findobj"]["snr_thresh"] = 7
            par["reduce"]["findobj"]["find_min_max"] = [65, 103]
            par["sensfunc"]["IR"]["resln_guess"] = 30

            # par["reduce"]["skysub"]["no_local_sky"] = True
            # par["reduce"]["skysub"]["no_poly"] = True

        elif grating == "Vph300":
            par["calibrations"]["wavelengths"][
                "reid_arxiv"
            ] = "ldt_rimas_HK_300_OH.fits"
            par["calibrations"]["wavelengths"]["n_first"] = 3
            par["calibrations"]["wavelengths"]["n_final"] = 5
            par["calibrations"]["wavelengths"]["sigdetect"] = 5
            par["calibrations"]["wavelengths"]["lamps"] = ["OH_RIMAS_HK"]
            par["reduce"]["findobj"]["find_fwhm"] = 7
            par["reduce"]["findobj"]["snr_thresh"] = 10
            par["sensfunc"]["IR"]["resln_guess"] = 300

        else:
            raise ValueError(f"Grating {grating} not recognized for RIMAS VPH modes")

        return par


# jump
# GRISM Child Class for LDT/RIMAS ============================================#
class LDTRIMASGrismSpectrograph(LDTRIMASSpectrograph):
    """
    Child to handle common aspects of the LDT/RIMAS Grism (echelle) modes
    """

    name = "ldt_rimas_grism"
    comment = "LDT Rapid infrared IMAger Spectrometer, Med-Res Echelle Grism"
    pypeline = "Echelle"
    ech_fixed_format = True

    # TODO: Measure these values for RIMAS GRISM spectra!!!
    _grism_geometry = {
        "YJ": EchelleProps(
            norders=15,
            orders=np.arange(44, 29, -1, dtype=int),  # orders 30-44
            order_spat_pos=np.array(
                [
                    0.151,
                    0.21075,
                    0.2665,
                    0.319,
                    0.368,
                    0.41525,
                    0.4595,
                    0.50175,
                    0.54225,
                    0.58075,
                    0.61725,
                    0.6525,
                    0.68625,
                    0.71875,
                    0.74975,
                ]
            ),
            spec_min=np.array(
                [
                    707,
                    664,
                    642,
                    632,
                    592,
                    607,
                    603,
                    606,
                    621,
                    614,
                    628,
                    631,
                    649,
                    649,
                    674,
                ]
            ),
            spec_max=np.array(
                [
                    1348,
                    1344,
                    1352,
                    1370,
                    1362,
                    1359,
                    1366,
                    1362,
                    1334,
                    1330,
                    1309,
                    1295,
                    1281,
                    1280,
                    1241,
                ]
            ),
            platescale=0.19,
            dloglam=2.8188583080448797e-5,
            loglam_minmax=(np.log10(8597.0), np.log10(14040.0)),
        ),
        "HK": EchelleProps(
            norders=17,
            orders=np.arange(39, 22, -1, dtype=int),  # orders 23-39
            order_spat_pos=np.array(
                [
                    0.2475,
                    0.3135,
                    0.373,
                    0.4285,
                    0.48,
                    0.5275,
                    0.572,
                    0.613,
                    0.652,
                    0.689,
                    0.724,
                    0.756,
                    0.78725,
                    0.8165,
                    0.8445,
                    0.871,
                    0.896,
                ]
            ),
            spec_min=np.array(
                [
                    338,
                    284,
                    316,
                    302,
                    313,
                    317,
                    370,
                    388,
                    377,
                    402,
                    428,
                    420,
                    445,
                    470,
                    500,
                    506,
                    510,
                ]
            ),
            spec_max=np.array(
                [
                    1607,
                    1660,
                    1631,
                    1617,
                    1613,
                    1595,
                    1581,
                    1545,
                    1542,
                    1495,
                    1484,
                    1495,
                    1459,
                    1445,
                    1420,
                    1373,
                    1327,
                ]
            ),
            platescale=0.19,
            dloglam=2.3904436313733664e-5,
            loglam_minmax=(np.log10(13423.0), np.log10(25075.0)),
        ),
    }

    def _geometry(self) -> EchelleProps:
        """Return the echelle geometry for the active RIMAS arm.

        Returns
        -------
        EchelleProps
            Fixed-format echelle geometry for ``self.arm``.
        """
        if getattr(self, "arm", None) is None:
            raise PypeItError(
                "RIMAS grism arm is undefined; call config_specific_par with an example frame."
            )
        return self._grism_geometry[self.arm]

    def validate_fitstbl(self, fitstbl: astropy.table.Table) -> astropy.table.Table:
        """Validate the metadata table

        Because of the multiple arms and modes of RIMAS, this method removes
        from the metadata table frames not associated with this mode,
        namely removes frames that are not the med-res echelle grism.

        Parameters
        ----------
        fitstbl : `astropy.table.Table`_
            The metadata table to be validated.

        Returns
        -------
        `astropy.table.Table`_
            The validated metadata table.
        """
        is_dark = fitstbl["filter1"] == "blank"

        # Only keep non-dark frames with the echelle grism -- no IFU mode!
        grism_idx = (
            (fitstbl["dispname"] == "Grism")
            & (fitstbl["decker"] != "open")
            & np.logical_not(is_dark)
        )

        # Include darks taken in other optical configurations only when their
        # detector arm and raw readout geometry match data selected above.
        dark_idx = self._mode_compatible_darks(fitstbl, grism_idx)

        # Return the corrected table
        return fitstbl[grism_idx | dark_idx]

    @staticmethod
    def configuration_list() -> dict[str, dict[str, np.str_]]:
        """
        Return the default-ordered list of GRISM configurations for this spectrograph

        Returns
        -------
        dict
            Dictionary of default-ordered spectrograph configurations.
        """
        return {
            # All of the GRISM possible combinations
            "YJ12l": {
                "arm": np.str_("YJ"),
                "dispname": np.str_("Grism"),
                "decker": np.str_("1.2'' long"),
            },
            "HK12l": {
                "arm": np.str_("HK"),
                "dispname": np.str_("Grism"),
                "decker": np.str_("1.2'' long"),
            },
            "YJ06": {
                "arm": np.str_("YJ"),
                "dispname": np.str_("Grism"),
                "decker": np.str_("0.6''"),
            },
            "HK06": {
                "arm": np.str_("HK"),
                "dispname": np.str_("Grism"),
                "decker": np.str_("0.6''"),
            },
            "YJ10": {
                "arm": np.str_("YJ"),
                "dispname": np.str_("Grism"),
                "decker": np.str_("1.0''"),
            },
            "HK10": {
                "arm": np.str_("HK"),
                "dispname": np.str_("Grism"),
                "decker": np.str_("1.0''"),
            },
            "YJ20": {
                "arm": np.str_("YJ"),
                "dispname": np.str_("Grism"),
                "decker": np.str_("2.0''"),
            },
            "HK20": {
                "arm": np.str_("HK"),
                "dispname": np.str_("Grism"),
                "decker": np.str_("2.0''"),
            },
        }

    @classmethod
    def default_pypeit_par(cls) -> pypeitpar.PypeItPar:
        """
        Return the default parameters to use for this instrument.

        Returns
        -------
        :class:`~pypeit.par.pypeitpar.PypeItPar`
            Parameters required by all PypeIt methods.
        """
        # Get the PypeIt and RIMAS-wide default parameters
        par = super().default_pypeit_par()

        # Grism-specific slit tracing parameters.
        par["calibrations"]["slitedges"]["det_min_spec_length"] = 0.10
        par["calibrations"]["slitedges"]["sync_predict"] = "nearest"
        par["calibrations"]["slitedges"]["add_missed_orders"] = False

        # Make a bad pixel mask
        par["calibrations"]["bpm_usebias"] = True

        # 1D wavelength solution
        par["calibrations"]["wavelengths"]["rms_thresh_frac_fwhm"] = 0.15
        par["calibrations"]["wavelengths"]["sigdetect"] = 2
        par["calibrations"]["wavelengths"]["fwhm"] = 4.0
        par["calibrations"]["wavelengths"]["n_final"] = 4
        # Wavelength-identification parameters
        par["calibrations"]["wavelengths"]["method"] = "holy-grail"
        par["calibrations"]["wavelengths"]["reid_arxiv"] = None
        par["calibrations"]["wavelengths"]["cc_thresh"] = 0.50
        par["calibrations"]["wavelengths"]["cc_local_thresh"] = 0.50
        # # Echelle parameters
        par["calibrations"]["wavelengths"]["echelle"] = True
        par["calibrations"]["wavelengths"]["ech_nspec_coeff"] = 5
        par["calibrations"]["wavelengths"]["ech_norder_coeff"] = 5
        par["calibrations"]["wavelengths"]["ech_sigrej"] = 3.0
        par["calibrations"]["wavelengths"]["qa_log"] = False
        # Measured FWHM is correct, but resulting wavelength solution is poor.
        # This should be explored further, but for now, turning off fwhm_fromlines helps.
        par["calibrations"]["wavelengths"]["fwhm_fromlines"] = False

        # Tilts
        par["calibrations"]["tilts"]["rm_continuum"] = True
        par["calibrations"]["tilts"]["tracethresh"] = 25.0
        par["calibrations"]["tilts"]["maxdev_tracefit"] = 0.04
        par["calibrations"]["tilts"]["maxdev2d"] = 0.04
        par["calibrations"]["tilts"]["spat_order"] = 3
        par["calibrations"]["tilts"]["spec_order"] = 4

        # Flats
        par["calibrations"]["flatfield"]["tweak_slits_thresh"] = 0.90
        par["calibrations"]["flatfield"]["tweak_slits_maxfrac"] = 0.10

        # Standards
        par["calibrations"]["standardframe"]["process"]["mask_cr"] = False

        # Extraction
        par["reduce"]["skysub"]["bspline_spacing"] = 0.8
        par["reduce"]["skysub"][
            "global_sky_std"
        ] = False  # Do not perform global sky subtraction for standard stars
        par["reduce"]["extraction"][
            "model_full_slit"
        ] = True  # local sky subtraction operates on entire slit
        par["sensfunc"]["IR"]["pix_shift_bounds"] = (-5.0, 5.0)
        par["sensfunc"]["IR"]["resln_guess"] = 4000

        # Telluric parameters
        par["telluric"]["pix_shift_bounds"] = (-10.0, 10.0)
        par["telluric"]["resln_frac_bounds"] = (0.4, 2.0)

        # Coadding
        par["coadd1d"]["wave_method"] = "log10"

        # Split-arm grism defaults from the original RIMAS classes. These
        # intentionally follow the broad echelle defaults above so that the
        # original child-class overrides are preserved in the consolidated class.
        par["calibrations"]["wavelengths"]["fwhm"] = 3.0
        par["calibrations"]["wavelengths"]["nsnippet"] = 1
        par["calibrations"]["tilts"]["spat_order"] = 4
        par["calibrations"]["tilts"]["spec_order"] = 5
        par["reduce"]["findobj"]["snr_thresh"] = 50.0

        return par

    def config_specific_par(
        self,
        inp: str | list | pathlib.Path | astropy.io.fits.Header | astropy.table.Table,
        inp_par: parset.ParSet = None,
    ) -> parset.ParSet:
        """
        Modify the PypeIt parameters to hard-wired values used for
        specific instrument configurations.

        Parameters
        ----------
        inp : str, list, pathlib.Path, `astropy.io.fits.Header`_, or `astropy.table.Table`_
            Input filename, an `astropy.io.fits.Header`_ object, a list of
            `astropy.io.fits.Header`_ objects, or a row from the metadata
            table.
        inp_par : :class:`~pypeit.par.parset.ParSet`, optional
            Parameter set used for the full run of PypeIt. If None, use
            :func:`default_pypeit_par`.

        Returns
        -------
        :class:`~pypeit.par.parset.ParSet`
            The PypeIt parameter set adjusted for configuration-specific
            parameter values.
        """
        # Start with instrument-wide parameters
        par = super().config_specific_par(inp, inp_par=inp_par)

        # Adjust parameters based on instrument settings
        binning = self.get_meta_value(inp, "binning")

        # Get the arm-specific parameters based on grating and decker
        if self.arm == "YJ":
            par = self.config_specific_par_grism_yj(par, self.decker)
        else:
            par = self.config_specific_par_grism_hk(par, self.decker)

        return self.adjust_config_specific_binning(par, binning)

    def config_specific_par_grism_yj(
        self, par: parset.ParSet, decker: str
    ) -> parset.ParSet:
        """Set YJ arm configuration-specific parameters for grism

        Parameters
        ----------
        par : :class:`~pypeit.par.parset.ParSet`
            The instrument-wide parameter set to be modified.
        decker : str
            The slit / decker used (from :meth:`get_meta_value`).

        Returns
        -------
        :class:`~pypeit.par.parset.ParSet`
            Modified parameter set for the YJ arm / grism.
        """
        # The YJ detector is noisy and makes tracing challenging.  Use the
        #   following "process" settings for tracing only.
        flat_process = {
            "mask_cr": True,
            "rmcompact": True,
            "sigclip": 4.5,
            "objlim": 3.0,
        }
        for frame in ["traceframe"]:  # , "pixelflatframe", "illumflatframe"]:
            for key, value in flat_process.items():
                par["calibrations"][frame]["process"][key] = value

        par["calibrations"]["slitedges"]["fit_min_spec_length"] = 0.15
        par["calibrations"]["slitedges"]["smash_range"] = [0.35, 0.70]

        match decker:
            case "0.6''":
                pass
            case "1.0''":
                pass
            case "2.0''":
                pass
            case _:
                raise PypeItError("Long Slits are not compatible with the grism!")

        return par

    def config_specific_par_grism_hk(
        self, par: parset.ParSet, decker: str
    ) -> parset.ParSet:
        """Set HK arm configuration-specific parameters for grism

        Parameters
        ----------
        par : :class:`~pypeit.par.parset.ParSet`
            The instrument-wide parameter set to be modified.
        decker : str
            The slit / decker used (from :meth:`get_meta_value`).

        Returns
        -------
        :class:`~pypeit.par.parset.ParSet`
            Modified parameter set for the HK arm / grism.
        """
        par["calibrations"]["slitedges"]["fit_min_spec_length"] = 0.20
        par["calibrations"]["slitedges"]["smash_range"] = [0.25, 0.98]

        match decker:
            case "0.6''":
                pass
            case "1.0''":
                pass
            case "2.0''":
                pass
            case _:
                raise PypeItError("Long Slits are not compatible with the grism!")

        return par

    # Fixed-format Echelle properties consumed by the larger PypeIt code =====#
    # These properties select the correct instrument arm to return based on the
    # `self.arm` instance variable.
    @property
    def norders(self) -> int:
        """
        Number of orders for this spectrograph. Should only be defined for
        echelle spectrographs, and it is undefined for the base class.
        """
        return self._geometry().norders

    @property
    def order_spat_pos(self) -> np.ndarray:
        """
        Return the expected spatial position of each echelle order.
        """
        return self._geometry().order_spat_pos

    @property
    def orders(self) -> np.ndarray:
        """
        Return the order number for each echelle order.
        """
        return self._geometry().orders

    @property
    def spec_min_max(self) -> np.ndarray:
        """
        Return the minimum and maximum spectral pixel expected for the
        spectral range of each order.
        """
        return np.vstack((self._geometry().spec_min, self._geometry().spec_max))

    def order_platescale(
        self, order_vec: np.ndarray, binning: str | None = None
    ) -> np.ndarray:
        """
        Return the platescale for each echelle order.

        This routine is only defined for echelle spectrographs, and it is
        undefined in the base class.

        Parameters
        ----------
        order_vec : `numpy.ndarray`_
            The vector providing the order numbers.
        binning : str, optional
            The string defining the spectral and spatial binning.

        Returns
        -------
        `numpy.ndarray`_
            An array with the platescale for each order provided by
            ``order_vec``.
        """
        return np.full(order_vec.size, self._geometry().platescale)

    @property
    def dloglam(self) -> float:
        """
        Return the logarithmic step in wavelength for output spectra (in Angstroms).
        """
        return self._geometry().dloglam

    @property
    def loglam_minmax(self) -> tuple[float, float]:
        """
        Return the base-10 logarithm of the first and last wavelength for
        output spectra (in Angstroms).
        """
        return self._geometry().loglam_minmax
