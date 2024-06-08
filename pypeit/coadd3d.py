"""
Module containing routines used by 3D datacubes.

.. include:: ../include/links.rst
"""

import os
import copy
import inspect

from astropy import wcs, units
from astropy.io import fits
import erfa
from scipy.interpolate import interp1d
import numpy as np

from pypeit import msgs
from pypeit import alignframe, datamodel, flatfield, io, spec2dobj, utils
from pypeit.core.flexure import calculate_image_phase
from pypeit.core import datacube, extract, flux_calib, parse
from pypeit.spectrographs.util import load_spectrograph

from IPython import embed


class DataCube(datamodel.DataContainer):
    """
    DataContainer to hold the products of a datacube

    The datamodel attributes are:

    .. include:: ../include/class_datamodel_datacube.rst

    Args:
        flux (`numpy.ndarray`_):
            The science datacube (nwave, nspaxel_y, nspaxel_x)
        sig (`numpy.ndarray`_):
            The error datacube (nwave, nspaxel_y, nspaxel_x)
        bpm (`numpy.ndarray`_):
            The bad pixel mask of the datacube (nwave, nspaxel_y, nspaxel_x).
            True values indicate a bad pixel
        wave (`numpy.ndarray`_):
            A 1D numpy array containing the wavelength array for convenience (nwave)
        blaze_wave (`numpy.ndarray`_):
            Wavelength array of the spectral blaze function
        blaze_spec (`numpy.ndarray`_):
            The spectral blaze function
        sensfunc (`numpy.ndarray`_, None):
            Sensitivity function (nwave,). Only saved if the data are fluxed.
        PYP_SPEC (str):
            Name of the PypeIt Spectrograph
        fluxed (bool):
            If the cube has been flux calibrated, this will be set to "True"

    Attributes:
        head0 (`astropy.io.fits.Header`_):
            Primary header
        filename (str):
            Filename to use when loading from file
        spect_meta (:obj:`dict`):
            Parsed meta from the header
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            Build from PYP_SPEC
        _ivar (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            Build from PYP_SPEC
        _wcs (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            Build from PYP_SPEC

    """
    version = '1.2.0'

    datamodel = {'flux': dict(otype=np.ndarray, atype=np.floating,
                              descr='Flux datacube in units of counts/s/Ang/arcsec^2 or '
                                    '10^-17 erg/s/cm^2/Ang/arcsec^2'),
                 'sig': dict(otype=np.ndarray, atype=np.floating,
                             descr='Error datacube (matches units of flux)'),
                 'bpm': dict(otype=np.ndarray, atype=np.uint8,
                             descr='Bad pixel mask of the datacube (0=good, 1=bad)'),
                 'wave': dict(otype=np.ndarray, atype=np.floating,
                              descr='Wavelength of each slice in the spectral direction. '
                                    'The units are Angstroms.'),
                 'blaze_wave': dict(otype=np.ndarray, atype=np.floating,
                                    descr='Wavelength array of the spectral blaze function'),
                 'blaze_spec': dict(otype=np.ndarray, atype=np.floating,
                                    descr='The spectral blaze function'),
                 'sensfunc': dict(otype=np.ndarray, atype=np.floating,
                                  descr='Sensitivity function 10^-17 erg/(counts/cm^2)'),
                 'PYP_SPEC': dict(otype=str, descr='PypeIt: Spectrograph name'),
                 'fluxed': dict(otype=bool, descr='Boolean indicating if the datacube is fluxed.')}

    internals = ['head0',
                 'filename',
                 'spectrograph',
                 'spect_meta',
                 '_ivar',  # This is set internally, and should be accessed with self.ivar
                 '_wcs'  # This is set internally, and should be accessed with self.wcs
                ]

    def __init__(self, flux, sig, bpm, wave, PYP_SPEC, blaze_wave, blaze_spec, sensfunc=None,
                 fluxed=None):

        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        _d = dict([(k, values[k]) for k in args[1:]])
        # Setup the DataContainer
        datamodel.DataContainer.__init__(self, d=_d)
        # Initialise the internals
        self._ivar = None
        self._wcs = None

    def _bundle(self):
        """
        Over-write default _bundle() method to separate the DetectorContainer
        into its own HDU

        Returns:
            :obj:`list`: A list of dictionaries, each list element is
            written to its own fits extension. See the description
            above.
        """
        d = []
        # Rest of the datamodel
        for key in self.keys():
            # Skip Nones
            if self[key] is None:
                continue
            # Array?
            if self.datamodel[key]['otype'] == np.ndarray:
                tmp = {}
                if self.datamodel[key]['atype'] == np.floating:
                    tmp[key] = self[key].astype(np.float32)
                else:
                    tmp[key] = self[key]
                d.append(tmp)
            else:
                # Add to header of the primary image
                d[0][key] = self[key]
        # Return
        return d

    def to_file(self, ofile, primary_hdr=None, hdr=None, **kwargs):
        """
        Over-load :func:`~pypeit.datamodel.DataContainer.to_file`
        to deal with the header

        Args:
            ofile (:obj:`str`):
                Filename
            primary_hdr (`astropy.io.fits.Header`_, optional):
                Base primary header.  Updated with new subheader keywords.  If
                None, initialized using :func:`~pypeit.io.initialize_header`.
            wcs (`astropy.io.fits.Header`_, optional):
                The World Coordinate System, represented by a fits header
            kwargs (dict):
                Keywords passed directly to parent ``to_file`` function.

        """
        if primary_hdr is None:
            primary_hdr = io.initialize_header()
        # Build the header
        if self.head0 is not None and self.PYP_SPEC is not None:
            spectrograph = load_spectrograph(self.PYP_SPEC)
            subheader = spectrograph.subheader_for_spec(self.head0, self.head0)
        else:
            subheader = {}
        # Add em in
        for key in subheader:
            primary_hdr[key] = subheader[key]
        # Do it
        super(DataCube, self).to_file(ofile, primary_hdr=primary_hdr, hdr=hdr, **kwargs)

    @classmethod
    def from_file(cls, ifile, verbose=True, chk_version=True, **kwargs):
        """
        Instantiate the object from an extension in the specified fits file.

        Over-load :func:`~pypeit.datamodel.DataContainer.from_file`
        to deal with the header
        
        Args:
            ifile (:obj:`str`, `Path`_):
                Fits file with the data to read
            verbose (:obj:`bool`, optional):
                Print informational messages (not currently used)
            chk_version (:obj:`bool`, optional):
                Passed to :func:`from_hdu`.
            kwargs (:obj:`dict`, optional):
                Arguments passed directly to :func:`from_hdu`.
        """
        with io.fits_open(ifile) as hdu:
            # Read using the base class
            self = cls.from_hdu(hdu, chk_version=chk_version, **kwargs)
            # Internals
            self.filename = ifile
            self.head0 = hdu[1].header  # Actually use the first extension here, since it contains the WCS
            # Meta
            self.spectrograph = load_spectrograph(self.PYP_SPEC)
            self.spect_meta = self.spectrograph.parse_spec_header(hdu[0].header)
            self._ivar = None
            self._wcs = None
        return self

    @property
    def ivar(self):
        """
        Utility function to compute the inverse variance cube

        Returns
        -------
        self._ivar : `numpy.ndarray`_
            The inverse variance of the datacube. Note that self._ivar should
            not be accessed directly, and you should only call self.ivar
        """
        if self._ivar is None:
            self._ivar = utils.inverse(self.sig**2)
        return self._ivar

    @property
    def wcs(self):
        """
        Utility function to provide the world coordinate system of the datacube

        Returns
        -------
        self._wcs : `astropy.wcs.WCS`_
            The WCS based on the stored header information. Note that self._wcs should
            not be accessed directly, and you should only call self.wcs
        """
        if self._wcs is None:
            self._wcs = wcs.WCS(self.head0)
        return self._wcs


class DARcorrection:
    """
    This class holds all of the functions needed to quickly compute the differential atmospheric refraction correction.
    """
    def __init__(self, airmass, parangle, pressure, temperature, humidity, cosdec, wave_ref=4500.0):
        """
        Args:
            airmass (:obj:`float`):
                The airmass of the observations (unitless)
            parangle (:obj:`float`):
                The parallactic angle of the observations (units=radians, relative to North, towards East is postive)
            pressure (:obj:`float`):
                The atmospheric pressure during the observations in Pascal. Valid range is from 10kPa - 140 kPa.
            temperature (:obj:`float`):
                Temperature in degree Celsius. Valid temperate range is -40 to
                100 degree Celsius.
            humidity (:obj:`float`):
                The humidity during the observations (Expressed as a percentage, not a fraction!).
                Valid range is 0 to 100.
            cosdec (:obj:`float`):
                Cosine of the target declination.
            wave_ref (:obj:`float`, optional):
                Reference wavelength (The DAR correction will be performed relative to this wavelength)
        """
        msgs.info("Preparing the parameters for the DAR correction")

        # Get DAR parameters
        self.airmass = airmass  # unitless
        self.parangle = parangle
        self.pressure = pressure * units.mbar
        self.temperature = temperature * units.deg_C
        self.humidity = humidity/100.0
        self.wave_ref = wave_ref*units.Angstrom
        self.cosdec = cosdec

        # Calculate the coefficients of the correction
        self.refa, self.refb = erfa.refco(self.pressure.to_value(units.hPa), self.temperature.to_value(units.deg_C),
                                          self.humidity, self.wave_ref.to_value(units.micron))

        # Print out the DAR parameters
        msgs.info("DAR correction parameters:" + msgs.newline() +
                  "   Airmass = {0:.2f}".format(self.airmass) + msgs.newline() +
                  "   Pressure = {0:.2f} mbar".format(self.pressure.to_value(units.mbar)) + msgs.newline() +
                  "   Humidity = {0:.2f} %".format(self.humidity*100.0) + msgs.newline() +
                  "   Temperature = {0:.2f} deg C".format(self.temperature.to_value(units.deg_C)) + msgs.newline() +
                  "   Reference wavelength = {0:.2f} Angstrom".format(self.wave_ref.to_value(units.Angstrom)))

    def calculate_dispersion(self, waves):
        """ Calculate the total atmospheric dispersion relative to the reference wavelength

        Parameters
        ----------
        waves : `numpy.ndarray`_
            1D array of wavelengths (units must be Angstroms)

        Returns
        -------
        full_dispersion : :obj:`float`
            The atmospheric dispersion (in degrees) for each wavelength input
        """

        # Calculate the zenith angle
        z = np.arccos(1.0/self.airmass)

        # Calculate the coefficients of the correction
        # self.refa, self.refb = erfa.refco(self.pressure.to_value(units.hPa), self.temperature.to_value(units.deg_C),
        #                                   self.humidity, self.wave_ref.to_value(units.micron))
        cnsa, cnsb = erfa.refco(self.pressure.to_value(units.hPa), self.temperature.to_value(units.deg_C),
                                self.humidity, (waves*units.Angstrom).to_value(units.micron))
        dar_full = np.rad2deg((self.refa-cnsa) * np.tan(z) + (self.refb-cnsb) * np.tan(z)**3)
        return dar_full

    def correction(self, waves):
        """
        Main routine that computes the DAR correction for both right ascension and declination.

        Parameters
        ----------
        waves : `numpy.ndarray`_
            1D array of wavelengths (units must be Angstroms)

        Returns
        -------
        ra_corr : `numpy.ndarray`_
            The RA component of the atmospheric dispersion correction (in degrees) for each wavelength input.
        dec_corr : `numpy.ndarray`_
            The Dec component of the atmospheric dispersion correction (in degrees) for each wavelength input.
        """
        # Determine the correction angle
        corr_ang = self.parangle - np.pi/2
        # Calculate the full amount of refraction
        dar_full = self.calculate_dispersion(waves)

        # Calculate the correction in dec and RA for each detector pixel
        # These numbers should be ADDED to the original RA and Dec values
        ra_corr = (dar_full/self.cosdec)*np.cos(corr_ang)
        dec_corr = -dar_full*np.sin(corr_ang)
        return ra_corr, dec_corr


class CoAdd3D:
    """
    Main routine to convert processed PypeIt spec2d frames into
    DataCube (spec3d) files. This routine is only used for IFU
    data reduction.

    Algorithm steps are detailed in the coadd routine.
    """
    # Superclass factory method generates the subclass instance
    @classmethod
    def get_instance(cls, spec2dfiles, par, skysub_frame=None, scale_corr=None, ra_offsets=None,
                     dec_offsets=None, spectrograph=None, det=1, overwrite=False, show=False,
                     debug=False):
        """
        Instantiate the subclass appropriate for the provided spectrograph.

        The class to instantiate must match the ``pypeline``
        attribute of the provided ``spectrograph``, and must be a
        subclass of :class:`CoAdd3D`; see the parent class
        instantiation for parameter descriptions.

        Returns:
            :class:`CoAdd3D`: One of the subclasses with
            :class:`CoAdd3D` as its base.
        """
        return next(c for c in cls.__subclasses__()
                    if c.__name__ == (spectrograph.pypeline + 'CoAdd3D'))(
                        spec2dfiles, par, skysub_frame=skysub_frame, scale_corr=scale_corr,
                        ra_offsets=ra_offsets, dec_offsets=dec_offsets, spectrograph=spectrograph,
                        det=det, overwrite=overwrite, show=show, debug=debug)

    def __init__(self, spec2dfiles, par, skysub_frame=None, scale_corr=None, ra_offsets=None,
                 dec_offsets=None, spectrograph=None, det=None, overwrite=False, show=False,
                 debug=False):
        """

        Args:
            spec2dfiles (:obj:`list`):
                List of all spec2D files
            par (:class:`~pypeit.par.pypeitpar.PypeItPar`):
                An instance of the parameter set.  If None, assumes that detector 1
                is the one reduced and uses the default reduction parameters for the
                spectrograph (see
                :func:`~pypeit.spectrographs.spectrograph.Spectrograph.default_pypeit_par`
                for the relevant spectrograph class).
            skysub_frame (:obj:`list`, optional):
                If not None, this should be a list of frames to use for the sky subtraction of each individual
                entry of spec2dfiles. It should be the same length as spec2dfiles.
            scale_corr (:obj:`list`, optional):
                If not None, this should be a list of relative scale correction options. It should be the
                same length as spec2dfiles.
            ra_offsets (:obj:`list`, optional):
                If not None, this should be a list of relative RA offsets of each frame. It should be the
                same length as spec2dfiles. The units should be degrees.
            dec_offsets (:obj:`list`, optional):
                If not None, this should be a list of relative Dec offsets of each frame. It should be the
                same length as spec2dfiles. The units should be degrees.
            spectrograph (:obj:`str`, :class:`~pypeit.spectrographs.spectrograph.Spectrograph`, optional):
                The name or instance of the spectrograph used to obtain the data.
                If None, this is pulled from the file header.
            det (:obj:`int`_, optional):
                Detector index
            overwrite (:obj:`bool`, optional):
                Overwrite the output file, if it exists?
            show (:obj:`bool`, optional):
                Show results in ginga
            debug (:obj:`bool`, optional):
                Show QA for debugging.
        """
        # TODO :: Consider loading all calibrations into a single variable within the main CoAdd3D parent class.
        # Set the variables
        self.spec2d = spec2dfiles
        self.numfiles = len(spec2dfiles)
        self.par = par
        self.overwrite = overwrite
        self.chk_version = self.par['rdx']['chk_version']
        # Extract some parsets for simplicity
        self.cubepar = self.par['reduce']['cube']
        self.flatpar = self.par['calibrations']['flatfield']
        self.senspar = self.par['sensfunc']
        # Extract some commonly used variables
        self.method = self.cubepar['method']
        self.combine = self.cubepar['combine']
        self.align = self.cubepar['align']
        self.correct_dar = self.cubepar['correct_dar']
        # Do some quick checks on the input options
        if skysub_frame is not None and len(skysub_frame) != self.numfiles:
            msgs.error("The skysub_frame list should be identical length to the spec2dfiles list")
        if scale_corr is not None and len(scale_corr) != self.numfiles:
            msgs.error("The scale_corr list should be identical length to the spec2dfiles list")
        if ra_offsets is not None and len(ra_offsets) != self.numfiles:
            msgs.error("The ra_offsets list should be identical length to the spec2dfiles list")
        if dec_offsets is not None and len(dec_offsets) != self.numfiles:
            msgs.error("The dec_offsets list should be identical length to the spec2dfiles list")
        # Make sure both ra_offsets and dec_offsets are either both None or both lists
        if ra_offsets is None and dec_offsets is not None:
            msgs.error("If you provide dec_offsets, you must also provide ra_offsets")
        if ra_offsets is not None and dec_offsets is None:
            msgs.error("If you provide ra_offsets, you must also provide dec_offsets")
        # Set the frame specific options
        self.skysub_frame = skysub_frame
        self.scale_corr = scale_corr
        self.ra_offsets = list(ra_offsets) if isinstance(ra_offsets, np.ndarray) else ra_offsets
        self.dec_offsets = list(dec_offsets) if isinstance(dec_offsets, np.ndarray) else dec_offsets
        # If there is only one frame being "combined" AND there's no reference image, then don't compute the translation.
        if self.numfiles == 1 and self.cubepar["reference_image"] is None:
            if self.align:
                msgs.warn("Parameter 'align' should be False when there is only one frame and no reference image")
                msgs.info("Setting 'align' to False")
            self.align = False
        if self.ra_offsets is not None:
            if not self.align:
                msgs.warn("When 'ra_offset' and 'dec_offset' are set, 'align' must be True.")
                msgs.info("Setting 'align' to True")
            self.align = True
        # If no ra_offsets or dec_offsets have been provided, initialise the lists
        self.user_alignment = True
        if self.ra_offsets is None and self.dec_offsets is None:
            msgs.info("No RA or Dec offsets have been provided.")
            if self.align:
                msgs.info("An automatic alignment will be performed using WCS information from the headers.")
            # User offsets are not provided, so turn off the user_alignment
            self.user_alignment = False
            # Initialise the lists of ra_offsets and dec_offsets
            self.ra_offsets = [0.0]*self.numfiles
            self.dec_offsets = [0.0]*self.numfiles

        # Check on Spectrograph input
        if spectrograph is None:
            with fits.open(spec2dfiles[0]) as hdu:
                spectrograph = hdu[0].header['PYP_SPEC']

        self.spec = load_spectrograph(spectrograph)
        self.specname = self.spec.name

        # Initialise arrays for storage
        self.ifu_ra, self.ifu_dec = np.array([]), np.array([])  # The RA and Dec at the centre of the IFU, as stored in the header

        self.all_sci, self.all_ivar, self.all_wave, self.all_slitid, self.all_wghts = [], [], [], [], []
        self.all_tilts, self.all_slits, self.all_align = [], [], []
        self.all_wcs, self.all_ra, self.all_dec, self.all_dar = [], [], [], []
        self.weights = np.ones(self.numfiles)  # Weights to use when combining cubes

        self._dspat = None if self.cubepar['spatial_delta'] is None else self.cubepar['spatial_delta'] / 3600.0  # binning size on the sky (/3600 to convert to degrees)
        self._dwv = self.cubepar['wave_delta']  # linear binning size in wavelength direction (in Angstroms)

        # TODO :: The default behaviour (combine=False, align=False) produces a datacube that uses the instrument WCS
        #  It should be possible (and perhaps desirable) to do a spatial alignment (i.e. align=True), apply this to the
        #  RA,Dec values of each pixel, and then use the instrument WCS to save the output (or, just adjust the crval).
        #  At the moment, if the user wishes to spatially align the frames, a different WCS is generated.

        # Determine what method is requested
        self.spec_subpixel, self.spat_subpixel, self.slice_subpixel = 1, 1, 1
        self.skip_subpix_weights = True
        if self.method == "subpixel":
            self.spec_subpixel, self.spat_subpixel, self.slice_subpixel = self.cubepar['spec_subpixel'], self.cubepar['spat_subpixel'], self.cubepar['slice_subpixel']
            self.skip_subpix_weights = False
            msgs.info("Adopting the subpixel algorithm to generate the datacube, with subpixellation scales:" + msgs.newline() +
                      f"  Spectral: {self.spec_subpixel}" + msgs.newline() +
                      f"  Spatial: {self.spat_subpixel}" + msgs.newline() +
                      f"  Slices: {self.slice_subpixel}")
        elif self.method == "ngp":
            msgs.info("Adopting the nearest grid point (NGP) algorithm to generate the datacube.")
            self.skip_subpix_weights = True
        else:
            msgs.error(f"The following datacube method is not allowed: {self.method}")

        # Get the detector number and string representation
        if det is None:
            det = 1 if self.par['rdx']['detnum'] is None else self.par['rdx']['detnum']
        self.detname = self.spec.get_det_name(det)

        # Check if the output file exists
        self.check_outputs()

        # Check the reference cube and image exist, if requested
        self.fluxcal = False
        self.blaze_wave, self.blaze_spec = None, None
        self.blaze_spline, self.flux_spline = None, None
        self.flat_splines = dict()  # A dictionary containing the splines of the flatfield
        if self.cubepar['standard_cube'] is not None:
            self.make_sensfunc()

        # If a reference image has been set, check that it exists
        if self.cubepar['reference_image'] is not None:
            if not os.path.exists(self.cubepar['reference_image']):
                msgs.error("Reference image does not exist:" + msgs.newline() + self.cubepar['reference_image'])

        # Load the default scaleimg frame for the scale correction
        self.scalecorr_default = "none"
        self.relScaleImgDef = np.array([1])
        self.set_default_scalecorr()

        # Load the default sky frame to be used for sky subtraction
        self.skysub_default = "image"
        self.skyImgDef, self.skySclDef = None, None  # This is the default behaviour (i.e. to use the "image" for the sky subtraction)
        self.set_default_skysub()

    def check_outputs(self):
        """
        Check if any of the intended output files already exist. This check should be done near the
        beginning of the coaddition, to avoid any computation that won't be saved in the event that
        files won't be overwritten.
        """
        if self.combine:
            outfile = datacube.get_output_filename("", self.cubepar['output_filename'], self.combine)
            out_whitelight = datacube.get_output_whitelight_filename(outfile)
            if os.path.exists(outfile) and not self.overwrite:
                msgs.error("Output filename already exists:"+msgs.newline()+outfile)
            if os.path.exists(out_whitelight) and self.cubepar['save_whitelight'] and not self.overwrite:
                msgs.error("Output filename already exists:"+msgs.newline()+out_whitelight)
        else:
            # Finally, if there's just one file, check if the output filename is given
            if self.numfiles == 1 and self.cubepar['output_filename'] != "":
                outfile = datacube.get_output_filename("", self.cubepar['output_filename'], True, -1)
                out_whitelight = datacube.get_output_whitelight_filename(outfile)
                if os.path.exists(outfile) and not self.overwrite:
                    msgs.error("Output filename already exists:" + msgs.newline() + outfile)
                if os.path.exists(out_whitelight) and self.cubepar['save_whitelight'] and not self.overwrite:
                    msgs.error("Output filename already exists:" + msgs.newline() + out_whitelight)
            else:
                for ff in range(self.numfiles):
                    outfile = datacube.get_output_filename(self.spec2d[ff], self.cubepar['output_filename'], self.combine, ff+1)
                    out_whitelight = datacube.get_output_whitelight_filename(outfile)
                    if os.path.exists(outfile) and not self.overwrite:
                        msgs.error("Output filename already exists:" + msgs.newline() + outfile)
                    if os.path.exists(out_whitelight) and self.cubepar['save_whitelight'] and not self.overwrite:
                        msgs.error("Output filename already exists:" + msgs.newline() + out_whitelight)

    def set_blaze_spline(self, wave_spl, spec_spl):
        """
        Generate a spline that represents the blaze function. This only needs to be done once,
        because it is used as the reference blaze. It is only important if you are combining
        frames that require a grating correction (i.e. have slightly different grating angles).

        Args:
            wave_spl (`numpy.ndarray`_):
                1D wavelength array where the blaze has been evaluated
            spec_spl (`numpy.ndarray`_):
                1D array (same size as wave_spl), that represents the blaze function for each wavelength.
        """
        # Check if a reference blaze spline exists (either from a standard star if fluxing or from a previous
        # exposure in this for loop)
        if self.blaze_spline is None:
            self.blaze_wave, self.blaze_spec = wave_spl, spec_spl
            self.blaze_spline = interp1d(wave_spl, spec_spl, kind='linear',
                                         bounds_error=False, fill_value="extrapolate")

    def make_sensfunc(self):
        """
        Generate the sensitivity function to be used for the flux calibration.
        """
        self.fluxcal = True
        # The first standard star cube is used as the reference blaze spline
        if self.cubepar['grating_corr']:
            # Load the blaze information
            stdcube = fits.open(self.cubepar['standard_cube'])
            # If a reference blaze spline has not been set, do that now.
            self.set_blaze_spline(stdcube['BLAZE_WAVE'].data, stdcube['BLAZE_SPEC'].data)
        # Generate a spline representation of the sensitivity function
        self.flux_spline = datacube.make_sensfunc(self.cubepar['standard_cube'], self.senspar,
                                                  blaze_wave=self.blaze_wave, blaze_spline=self.blaze_spline,
                                                  grating_corr=self.cubepar['grating_corr'])

    def set_default_scalecorr(self):
        """
        Set the default mode to use for relative spectral scale correction.
        """
        if self.cubepar['scale_corr'] is not None:
            if self.cubepar['scale_corr'] == "image":
                msgs.info("The default relative spectral illumination correction will use the science image")
                self.scalecorr_default = "image"
            else:
                msgs.info("Loading default scale image for relative spectral illumination correction:" +
                          msgs.newline() + self.cubepar['scale_corr'])
                try:
                    spec2DObj = spec2dobj.Spec2DObj.from_file(self.cubepar['scale_corr'],
                                                              self.detname,
                                                              chk_version=self.chk_version)
                except Exception as e:
                    msgs.warn(f'Loading spec2d file raised {type(e).__name__}:\n{str(e)}')
                    msgs.warn("Could not load scaleimg from spec2d file:" + msgs.newline() +
                              self.cubepar['scale_corr'] + msgs.newline() +
                              "scale correction will not be performed unless you have specified the correct" + msgs.newline() +
                              "scale_corr file in the spec2d block")
                    self.cubepar['scale_corr'] = None
                    self.scalecorr_default = "none"
                else:
                    self.relScaleImgDef = spec2DObj.scaleimg
                    self.scalecorr_default = self.cubepar['scale_corr']

    def get_current_scalecorr(self, spec2DObj, scalecorr=None):
        """
        Determine the scale correction that should be used to correct
        for the relative spectral scaling of the science frame

        Args:
            spec2DObj (:class:`~pypeit.spec2dobj.Spec2DObj`):
                2D PypeIt spectra object.

            scalecorr (:obj:`str`, optional):
                A string that describes what mode should be used for the sky
                subtraction. The allowed values are:

                    * default: Use the default value, as defined in
                      :func:`set_default_scalecorr`.

                    * image: Use the relative scale that was derived from the
                      science frame

                    * none: Do not perform relative scale correction

        Returns:
            :obj:`tuple`: Contains (this_scalecorr, relScaleImg) where
            this_scalecorr is a :obj:`str` that describes the scale correction
            mode to be used (see scalecorr description) and relScaleImg is a
            `numpy.ndarray`_ (2D, same shape as science frame) containing the
            relative spectral scaling to apply to the science frame.
        """
        this_scalecorr = self.scalecorr_default
        relScaleImg = self.relScaleImgDef.copy()
        if scalecorr is not None:
            if scalecorr.lower() == 'default':
                if self.scalecorr_default == "image":
                    relScaleImg = spec2DObj.scaleimg
                    this_scalecorr = "image"  # Use the current spec2d for the relative spectral illumination scaling
                else:
                    this_scalecorr = self.scalecorr_default  # Use the default value for the scale correction
            elif scalecorr.lower() == 'image':
                relScaleImg = spec2DObj.scaleimg
                this_scalecorr = "image"  # Use the current spec2d for the relative spectral illumination scaling
            elif scalecorr.lower() == 'none':
                relScaleImg = np.array([1])
                this_scalecorr = "none"  # Don't do relative spectral illumination scaling
            else:
                # Load a user specified frame for sky subtraction
                msgs.info("Loading the following frame for the relative spectral illumination correction:" +
                          msgs.newline() + scalecorr)
                try:
                    spec2DObj_scl = spec2dobj.Spec2DObj.from_file(scalecorr, self.detname,
                                                                  chk_version=self.chk_version)
                except Exception as e:
                    msgs.warn(f'Loading spec2d file raised {type(e).__name__}:\n{str(e)}')
                    msgs.error("Could not load skysub image from spec2d file:" + msgs.newline() + scalecorr)
                else:
                    relScaleImg = spec2DObj_scl.scaleimg
                    this_scalecorr = scalecorr
        if this_scalecorr == "none":
            msgs.info("Relative spectral illumination correction will not be performed.")
        else:
            msgs.info("Using the following frame for the relative spectral illumination correction:" +
                      msgs.newline() + this_scalecorr)
        # Return the scaling correction for this frame
        return this_scalecorr, relScaleImg

    def set_default_skysub(self):
        """
        Set the default mode to use for sky subtraction.
        """
        if self.cubepar['skysub_frame'] in [None, 'none', '', 'None']:
            self.skysub_default = "none"
            self.skyImgDef = np.array([0.0])  # Do not perform sky subtraction
            self.skySclDef = np.array([0.0])  # Do not perform sky subtraction
        elif self.cubepar['skysub_frame'] == "image":
            msgs.info("The sky model in the spec2d science frames will be used for sky subtraction" + msgs.newline() +
                      "(unless specific skysub frames have been specified)")
            self.skysub_default = "image"
        else:
            msgs.info("Loading default image for sky subtraction:" +
                      msgs.newline() + self.cubepar['skysub_frame'])
            try:
                spec2DObj = spec2dobj.Spec2DObj.from_file(self.cubepar['skysub_frame'],
                                                          self.detname,
                                                          chk_version=self.chk_version)
                skysub_exptime = self.spec.get_meta_value([spec2DObj.head0], 'exptime')
            except:
                msgs.error("Could not load skysub image from spec2d file:" + msgs.newline() + self.cubepar['skysub_frame'])
            else:
                self.skysub_default = self.cubepar['skysub_frame']
                self.skyImgDef = spec2DObj.sciimg / skysub_exptime  # Sky counts/second
                # self.skyImgDef = spec2DObj.skymodel/skysub_exptime  # Sky counts/second
                self.skySclDef = spec2DObj.scaleimg

    def get_current_skysub(self, spec2DObj, exptime, opts_skysub=None):
        """
        Determine the sky frame that should be used to subtract from the science frame

        Args:
            spec2DObj (:class:`~pypeit.spec2dobj.Spec2DObj`):
                2D PypeIt spectra object.
            exptime (:obj:`float`):
                The exposure time of the science frame (in seconds)
            opts_skysub (:obj:`str`, optional):
                A string that describes what mode should be used for the sky
                subtraction. The allowed values are:

                    * default: Use the default value, as defined in
                      :func:`set_default_skysub`

                    * image: Use the sky model derived from the science frame

                    * none: Do not perform sky subtraction

        Returns:
            :obj:`tuple`: Contains (this_skysub, skyImg, skyScl) where
            this_skysub is a :obj:`str` that describes the sky subtration mode
            to be used (see opts_skysub description), skyImg is a
            `numpy.ndarray`_ (2D, same shape as science frame) containing the
            sky frame to be subtracted from the science frame, and skyScl is a
            `numpy.ndarray`_ (2D, same shape as science frame) containing the
            relative spectral scaling that has been applied to the returned sky
            frame.
        """
        this_skysub = self.skysub_default
        if self.skysub_default == "image":
            skyImg = spec2DObj.skymodel
            skyScl = spec2DObj.scaleimg
        else:
            skyImg = self.skyImgDef.copy() * exptime
            skyScl = self.skySclDef.copy()
        # See if there's any changes from the default behaviour
        if opts_skysub is not None:
            if opts_skysub.lower() == 'default':
                if self.skysub_default == "image":
                    skyImg = spec2DObj.skymodel
                    skyScl = spec2DObj.scaleimg
                    this_skysub = "image"  # Use the current spec2d for sky subtraction
                else:
                    skyImg = self.skyImgDef.copy() * exptime
                    skyScl = self.skySclDef.copy() * exptime
                    this_skysub = self.skysub_default  # Use the global value for sky subtraction
            elif opts_skysub.lower() == 'image':
                skyImg = spec2DObj.skymodel
                skyScl = spec2DObj.scaleimg
                this_skysub = "image"  # Use the current spec2d for sky subtraction
            elif opts_skysub.lower() == 'none':
                skyImg = np.array([0.0])
                skyScl = np.array([1.0])
                this_skysub = "none"  # Don't do sky subtraction
            else:
                # Load a user specified frame for sky subtraction
                msgs.info("Loading skysub frame:" + msgs.newline() + opts_skysub)
                try:
                    spec2DObj_sky = spec2dobj.Spec2DObj.from_file(opts_skysub, self.detname,
                                                                  chk_version=self.chk_version)
                    skysub_exptime = self.spec.get_meta_value([spec2DObj_sky.head0], 'exptime')
                except:
                    msgs.error("Could not load skysub image from spec2d file:" + msgs.newline() + opts_skysub)
                skyImg = spec2DObj_sky.sciimg * exptime / skysub_exptime  # Sky counts
                skyScl = spec2DObj_sky.scaleimg
                this_skysub = opts_skysub  # User specified spec2d for sky subtraction
        if this_skysub == "none":
            msgs.info("Sky subtraction will not be performed.")
        else:
            msgs.info("Using the following frame for sky subtraction:" + msgs.newline() + this_skysub)
        # Return the skysub params for this frame
        return this_skysub, skyImg, skyScl

    def add_grating_corr(self, flatfile, waveimg, slits, spat_flexure=None):
        """
        Calculate the relative spectral sensitivity correction due to grating
        shifts with the input frames.

        Parameters
        ----------
        flatfile : :obj:`str`
            Unique path of a flatfield frame used to calculate the relative
            spectral sensitivity of the corresponding science frame.
        waveimg : `numpy.ndarray`_
            2D image (same shape as the science frame) indicating the wavelength
            of each detector pixel.
        slits : :class:`~pypeit.slittrace.SlitTraceSet`
            Class containing information about the slits
        spat_flexure : :obj:`float`, optional:
            Spatial flexure in pixels
        """
        # Check if the Flat file exists
        if not os.path.exists(flatfile):
            msgs.warn("Grating correction requested, but the following file does not exist:" + msgs.newline() + flatfile)
            return
        if flatfile not in self.flat_splines.keys():
            msgs.info("Calculating relative sensitivity for grating correction")
            # Load the Flat file
            flatimages = flatfield.FlatImages.from_file(flatfile, chk_version=self.chk_version)
            total_illum = flatimages.fit2illumflat(slits, finecorr=False, frametype='illum', initial=True, spat_flexure=spat_flexure) * \
                          flatimages.fit2illumflat(slits, finecorr=True, frametype='illum', initial=True, spat_flexure=spat_flexure)
            flatframe = flatimages.pixelflat_raw / total_illum
            if flatimages.pixelflat_spec_illum is None:
                # Calculate the relative scale
                scale_model = flatfield.illum_profile_spectral(flatframe, waveimg, slits,
                                                               slit_illum_ref_idx=self.flatpar['slit_illum_ref_idx'],
                                                               model=None, trim=self.flatpar['slit_trim'],
                                                               flexure=spat_flexure,
                                                               smooth_npix=self.flatpar['slit_illum_smooth_npix'])
            else:
                msgs.info("Using relative spectral illumination from FlatImages")
                scale_model = flatimages.pixelflat_spec_illum
            # Extract a quick spectrum of the flatfield
            wave_spl, spec_spl = extract.extract_hist_spectrum(waveimg, flatframe*utils.inverse(scale_model),
                                                               gpm=waveimg != 0, bins=slits.nspec)
            # Store the result
            self.flat_splines[flatfile] = interp1d(wave_spl, spec_spl, kind='linear', bounds_error=False, fill_value="extrapolate")
            self.flat_splines[flatfile + "_wave"] = wave_spl.copy()
            # Finally, if a reference blaze spline has not been set, do that now.
            self.set_blaze_spline(wave_spl, spec_spl)

    def run(self):
        """
        Main entry routine to set the order of operations to coadd the data. For specific
        details of this procedure, see the child routines.
        """
        msgs.bug("This routine should be overridden by child classes.")
        msgs.error("Cannot proceed without coding the run() routine.")


class SlicerIFUCoAdd3D(CoAdd3D):
    """
    Child of CoAdd3D for SlicerIFU data reduction. For documentation, see CoAdd3d parent class above.

    This child class of the IFU datacube creation performs the series of steps that are specific to
    slicer-based IFUs, including the following steps

    Data preparation:

    * Loads individual spec2d files
    * If requested, subtract the sky (either from a dedicated sky frame, or use the sky model stored in the science spec2d file)
    * The sky regions near the spectral edges of the slits are masked
    * Apply a relative spectral illumination correction (scalecorr) that registers all input frames to the scale illumination.
    * Generate a WCS of each individual frame, and calculate the RA and DEC of each individual detector pixel
    * Calculate the astrometric correction that is needed to align spatial positions along the slices
    * Compute the differential atmospheric refraction correction
    * Apply the extinction correction
    * Apply a grating correction (gratcorr) - This corrects for the relative spectral efficiency of combining data taken with multiple different grating angles
    * Flux calibrate

    Data cube generation:

    * If frames are not being combined, individual data cubes are generated and saved as a DataCube object. A white light image is also produced, if requested
    * If frames are being aligned and/or combined, the following steps are followed:
        - The output voxel sampling is computed (this must be consistent for all frames)
        - Frames are aligned (either by user-specified offsets, or by a fancy cross-correlation)
        - The relative weights to each for each detector pixel is computed
        - If frames are not being combined, individual DataCube's will be generated for each frame
        - If frames are being combined, a single DataCube will be generated.
        - White light images are also produced, if requested.

    """
    def __init__(self, spec2dfiles, par, skysub_frame=None, scale_corr=None, ra_offsets=None,
                 dec_offsets=None, spectrograph=None, det=1, overwrite=False, show=False,
                 debug=False):
        super().__init__(spec2dfiles, par, skysub_frame=skysub_frame, scale_corr=scale_corr,
                         ra_offsets=ra_offsets, dec_offsets=dec_offsets, spectrograph=spectrograph,
                         det=det, overwrite=overwrite, show=show, debug=debug)
        self.mnmx_wv = None  # Will be used to store the minimum and maximum wavelengths of every slit and frame.
        self._spatscale = np.zeros((self.numfiles, 2))  # index 0, 1 = pixel scale, slicer scale
        self._specscale = np.zeros(self.numfiles)
        # Loop through all of the frames, load the data, and save datacubes if no combining is required
        self.load()

    def get_alignments(self, spec2DObj, slits, spat_flexure=None):
        """
        Generate and return the spline interpolation fitting functions to be used for
        the alignment frames, as part of the astrometric correction.

        Parameters
        ----------
        spec2DObj : :class:`~pypeit.spec2dobj.Spec2DObj`
            2D PypeIt spectra object.
        slits : :class:`~pypeit.slittrace.SlitTraceSet`
            Class containing information about the slits
        spat_flexure: :obj:`float`, optional
            Spatial flexure in pixels

        Returns
        -------
        alignSplines : :class:`~pypeit.alignframe.AlignmentSplines`
            Alignment splines used for the astrometric correction
        """
        # Loading the alignments frame for these data
        alignments = None
        if self.cubepar['astrometric']:
            key = alignframe.Alignments.calib_type.upper()
            if key in spec2DObj.calibs:
                alignfile = os.path.join(spec2DObj.calibs['DIR'], spec2DObj.calibs[key])
                if os.path.exists(alignfile) and self.cubepar['astrometric']:
                    msgs.info("Loading alignments")
                    alignments = alignframe.Alignments.from_file(alignfile,
                                                                 chk_version=self.chk_version)
            else:
                msgs.warn(f'Processed alignment frame not recorded or not found!')
                msgs.info("Using slit edges for astrometric transform")
        else:
            msgs.info("Using slit edges for astrometric transform")
        # If nothing better was provided, use the slit edges
        if alignments is None:
            left, right, _ = slits.select_edges(initial=True, flexure=spat_flexure)
            locations = [0.0, 1.0]
            traces = np.append(left[:, None, :], right[:, None, :], axis=1)
        else:
            locations = self.par['calibrations']['alignment']['locations']
            traces = alignments.traces
        msgs.info("Generating alignment splines")
        return alignframe.AlignmentSplines(traces, locations, spec2DObj.tilts)

    def load(self):
        """
        This is the main function that loads in the data, and performs several frame-specific corrections.
        If the user does not wish to align or combine the individual datacubes, then this routine will also
        produce a spec3d file, which is a DataCube representation of a PypeIt spec2d frame for SlicerIFU data.

        This function should be called in the __init__ method, and initialises multiple variables. The variables
        initialised by this function include:

        * self.ifu_ra  -  The RA of the IFU pointing
        * self.ifu_dec  -  The Dec of the IFU pointing
        * self.mnmx_wv  -  The minimum and maximum wavelengths of every slit and frame.
        * self._spatscale  -  The native spatial scales of all spec2d frames.
        * self._specscale  -  The native spectral scales of all spec2d frames.
        * self.weights  -  Weights to use when combining cubes
        * self.flat_splines  -  Spline representations of the blaze function (based on the illumflat).
        * self.blaze_spline  -  Spline representation of the reference blaze function
        * self.blaze_wave  -  Wavelength array used to construct the reference blaze function
        * self.blaze_spec  -  Spectrum used to construct the reference blaze function

        As well as the primary arrays that store the pixel information for multiple spec2d frames, including:

        * self.all_sci
        * self.all_ivar
        * self.all_wave
        * self.all_slitid
        * self.all_wghts
        * self.all_tilts
        * self.all_slits
        * self.all_align
        * self.all_wcs
        * self.all_ra
        * self.all_dec
        * self.all_dar
        """
        # Load all spec2d files and prepare the data for making a datacube
        for ff, fil in enumerate(self.spec2d):
            # Load it up
            msgs.info("Loading PypeIt spec2d frame:" + msgs.newline() + fil)
            spec2DObj = spec2dobj.Spec2DObj.from_file(fil, self.detname,
                                                      chk_version=self.chk_version)
            detector = spec2DObj.detector
            spat_flexure = None  # spec2DObj.sci_spat_flexure

            # Load the header
            hdr0 = spec2DObj.head0
            self.ifu_ra = np.append(self.ifu_ra, self.spec.compound_meta([hdr0], 'ra'))
            self.ifu_dec = np.append(self.ifu_dec, self.spec.compound_meta([hdr0], 'dec'))

            # Get the exposure time
            exptime = self.spec.compound_meta([hdr0], 'exptime')

            # Initialise the slit edges
            msgs.info("Constructing slit image")
            slits = spec2DObj.slits
            slitid_img_init = slits.slit_img(pad=0, initial=True, flexure=spat_flexure)
            slits_left, slits_right, _ = slits.select_edges(initial=True, flexure=spat_flexure)

            # The order of operations below proceeds as follows:
            #  (1) Get science image
            #  (2) Subtract sky (note, if a joint fit has been performed, the relative scale correction is applied in the reduction!)
            #  (3) Apply relative scale correction to both science and ivar

            # Set the default behaviour if a global skysub frame has been specified
            this_skysub, skyImg, skyScl = self.get_current_skysub(spec2DObj, exptime,
                                                                  opts_skysub=self.skysub_frame[ff])

            # Load the relative scale image, if something other than the default has been provided
            this_scalecorr, relScaleImg = self.get_current_scalecorr(spec2DObj,
                                                                     scalecorr=self.scale_corr[ff])
            # Prepare the relative scaling factors
            relSclSky = skyScl / spec2DObj.scaleimg  # This factor ensures the sky has the same relative scaling as the science frame
            relScale = spec2DObj.scaleimg / relScaleImg  # This factor is applied to the sky subtracted science frame

            # Extract the relevant information from the spec2d file
            sciImg = (spec2DObj.sciimg - skyImg * relSclSky) * relScale  # Subtract sky and apply relative illumination
            ivar = spec2DObj.ivarraw / relScale ** 2
            waveimg = spec2DObj.waveimg
            bpmmask = spec2DObj.bpmmask

            # Mask the edges of the spectrum where the sky model is bad
            sky_is_good = datacube.make_good_skymask(slitid_img_init, spec2DObj.tilts)

            # TODO :: Really need to write some detailed information in the docs about all of the various corrections that can optionally be applied

            # TODO :: Include a flexure correction from the sky frame? Note, you cannot use the waveimg from a sky frame,
            #  since the heliocentric correction may have been applied to the sky frame. Need to recalculate waveimg using
            #  the slitshifts from a skyimage, and then apply the vel_corr from the science image.

            wnonzero = (waveimg != 0.0)
            if not np.any(wnonzero):
                msgs.error("The wavelength image contains only zeros - You need to check the data reduction.")
            wave0 = waveimg[wnonzero].min()
            # Calculate the delta wave in every pixel on the slit
            waveimp = np.roll(waveimg, 1, axis=0)
            waveimn = np.roll(waveimg, -1, axis=0)
            dwaveimg = np.zeros_like(waveimg)
            # All good pixels
            wnz = np.where((waveimg != 0) & (waveimp != 0))
            dwaveimg[wnz] = np.abs(waveimg[wnz] - waveimp[wnz])
            # All bad pixels
            wnz = np.where((waveimg != 0) & (waveimp == 0))
            dwaveimg[wnz] = np.abs(waveimg[wnz] - waveimn[wnz])
            # All endpoint pixels
            dwaveimg[0, :] = np.abs(waveimg[0, :] - waveimn[0, :])
            dwaveimg[-1, :] = np.abs(waveimg[-1, :] - waveimp[-1, :])
            dwv = np.median(dwaveimg[dwaveimg != 0.0]) if self.cubepar['wave_delta'] is None else self.cubepar['wave_delta']

            msgs.info("Using wavelength solution: wave0={0:.3f}, dispersion={1:.3f} Angstrom/pixel".format(wave0, dwv))

            # Obtain the minimum and maximum wavelength of all slits
            if self.mnmx_wv is None:
                self.mnmx_wv = np.zeros((len(self.spec2d), slits.nslits, 2))
            for slit_idx, slit_spat in enumerate(slits.spat_id):
                onslit_init = (slitid_img_init == slit_spat)
                self.mnmx_wv[ff, slit_idx, 0] = np.min(waveimg[onslit_init])
                self.mnmx_wv[ff, slit_idx, 1] = np.max(waveimg[onslit_init])

            # Find the largest spatial scale of all images being combined
            # TODO :: probably need to put this in the DetectorContainer
            pxscl = detector.platescale * parse.parse_binning(detector.binning)[1] / 3600.0  # This is degrees/pixel
            slscl = self.spec.get_meta_value([spec2DObj.head0], 'slitwid')
            self._spatscale[ff, 0] = pxscl
            self._spatscale[ff, 1] = slscl
            self._specscale[ff] = dwv

            # If the spatial scale has been set by the user, check that it doesn't exceed the pixel or slicer scales
            if self._dspat is not None:
                if pxscl > self._dspat:
                    msgs.warn("Spatial scale requested ({0:f} arcsec) is less than the pixel scale ({1:f} arcsec)".format(
                        3600.0 * self._dspat, 3600.0 * pxscl))
                if slscl > self._dspat:
                    msgs.warn("Spatial scale requested ({0:f} arcsec) is less than the slicer scale ({1:f} arcsec)".format(
                        3600.0 * self._dspat, 3600.0 * slscl))

            # Construct a good pixel mask
            # TODO: This should use the mask function to figure out which elements are masked.
            onslit_gpm = (slitid_img_init > 0) & (bpmmask.mask == 0) & sky_is_good

            # Generate the alignment splines, and then retrieve images of the RA and Dec of every pixel,
            # and the number of spatial pixels in each slit
            alignSplines = self.get_alignments(spec2DObj, slits, spat_flexure=spat_flexure)

            # Grab the WCS of this frame, and generate the RA and Dec images
            # NOTE :: These RA and Dec images are only used to setup the WCS of the datacube. The actual RA and Dec
            #         of each pixel in the datacube is calculated in the datacube.subpixellate() method.
            crval_wv = self.cubepar['wave_min'] if self.cubepar['wave_min'] is not None else wave0
            cd_wv = self.cubepar['wave_delta'] if self.cubepar['wave_delta'] is not None else dwv
            self.all_wcs.append(self.spec.get_wcs(spec2DObj.head0, slits, detector.platescale, crval_wv, cd_wv))
            ra_img, dec_img, minmax = slits.get_radec_image(self.all_wcs[ff], alignSplines, spec2DObj.tilts, initial=True, flexure=spat_flexure)

            # Extract wavelength and delta wavelength arrays from the images
            wave_ext = waveimg[onslit_gpm]
            dwav_ext = dwaveimg[onslit_gpm]

            # For now, work in sorted wavelengths
            wvsrt = np.argsort(wave_ext)
            wave_sort = wave_ext[wvsrt]
            dwav_sort = dwav_ext[wvsrt]
            # Here's an array to get back to the original ordering
            resrt = np.argsort(wvsrt)

            # Compute the DAR correction
            cosdec = np.cos(self.ifu_dec[ff] * np.pi / 180.0)
            airmass = self.spec.get_meta_value([spec2DObj.head0], 'airmass')  # unitless
            parangle = self.spec.get_meta_value([spec2DObj.head0], 'parangle')
            pressure = self.spec.get_meta_value([spec2DObj.head0], 'pressure')  # units are pascals
            temperature = self.spec.get_meta_value([spec2DObj.head0], 'temperature')  # units are degrees C
            humidity = self.spec.get_meta_value([spec2DObj.head0], 'humidity')  # Expressed as a percentage (not a fraction!)
            darcorr = DARcorrection(airmass, parangle, pressure, temperature, humidity, cosdec)

            # Compute the extinction correction
            msgs.info("Applying extinction correction")
            extinct = flux_calib.load_extinction_data(self.spec.telescope['longitude'],
                                                      self.spec.telescope['latitude'],
                                                      self.senspar['UVIS']['extinct_file'])
            # extinction_correction requires the wavelength is sorted
            extcorr_sort = flux_calib.extinction_correction(wave_sort * units.AA, airmass, extinct)

            # Correct for sensitivity as a function of grating angle
            # (this assumes the spectrum of the flatfield lamp has the same shape for all setups)
            gratcorr_sort = 1.0
            if self.cubepar['grating_corr']:
                # Load the flatfield file
                key = flatfield.FlatImages.calib_type.upper()
                if key not in spec2DObj.calibs:
                    msgs.error('Processed flat calibration file not recorded by spec2d file!')
                flatfile = os.path.join(spec2DObj.calibs['DIR'], spec2DObj.calibs[key])
                # Setup the grating correction
                self.add_grating_corr(flatfile, waveimg, slits, spat_flexure=spat_flexure)
                # Calculate the grating correction
                gratcorr_sort = datacube.correct_grating_shift(wave_sort, self.flat_splines[flatfile + "_wave"],
                                                               self.flat_splines[flatfile],
                                                               self.blaze_wave, self.blaze_spline)
            # Sensitivity function
            sensfunc_sort = 1.0
            if self.fluxcal:
                msgs.info("Calculating the sensitivity function")
                sensfunc_sort = self.flux_spline(wave_sort)
            # Convert the flux_sav to counts/s,  correct for the relative sensitivity of different setups
            extcorr_sort *= sensfunc_sort / (exptime * gratcorr_sort)
            # Correct for extinction
            sciImg[onslit_gpm] *= extcorr_sort[resrt]
            ivar[onslit_gpm] /= extcorr_sort[resrt] ** 2

            # Convert units to Counts/s/Ang/arcsec2
            # Slicer sampling * spatial pixel sampling
            sl_deg = np.sqrt(self.all_wcs[ff].wcs.cd[0, 0] ** 2 + self.all_wcs[ff].wcs.cd[1, 0] ** 2)
            px_deg = np.sqrt(self.all_wcs[ff].wcs.cd[1, 1] ** 2 + self.all_wcs[ff].wcs.cd[0, 1] ** 2)
            scl_units = dwav_sort * (3600.0 * sl_deg) * (3600.0 * px_deg)
            sciImg[onslit_gpm] /= scl_units[resrt]
            ivar[onslit_gpm] *= scl_units[resrt] ** 2

            # Calculate the weights relative to the zeroth cube
            self.weights[ff] = 1.0  # exptime  #np.median(flux_sav[resrt]*np.sqrt(ivar_sav[resrt]))**2
            wghts = self.weights[ff] * np.ones(sciImg.shape)

            # Get the slit image and then unset pixels in the slit image that are bad
            slitid_img_gpm = slitid_img_init * onslit_gpm.astype(int)

            # If individual frames are to be output without aligning them,
            # there's no need to store information, just make the cubes now
            if not self.combine and not self.align:
                # Get the output filename
                if self.numfiles == 1 and self.cubepar['output_filename'] != "":
                    outfile = datacube.get_output_filename("", self.cubepar['output_filename'], True, -1)
                else:
                    outfile = datacube.get_output_filename(fil, self.cubepar['output_filename'], self.combine, ff + 1)
                # Get the coordinate bounds
                slitlength = int(np.round(np.median(slits.get_slitlengths(initial=True, median=True))))
                numwav = int((np.max(waveimg) - wave0) / dwv)
                bins = self.spec.get_datacube_bins(slitlength, minmax, numwav)
                # Set the wavelength range of the white light image.
                wl_wvrng = None
                if self.cubepar['save_whitelight']:
                    wl_wvrng = datacube.get_whitelight_range(np.max(self.mnmx_wv[ff, :, 0]),
                                                             np.min(self.mnmx_wv[ff, :, 1]),
                                                             self.cubepar['whitelight_range'])
                # Make the datacube
                if self.method in ['subpixel', 'ngp']:
                    # Generate the datacube
                    flxcube, sigcube, bpmcube, wave = \
                        datacube.generate_cube_subpixel(self.all_wcs[ff], bins, sciImg, ivar, waveimg, slitid_img_gpm, wghts,
                                                        self.all_wcs[ff], spec2DObj.tilts, slits, alignSplines, darcorr,
                                                        self.ra_offsets[ff], self.dec_offsets[ff],
                                                        overwrite=self.overwrite, whitelight_range=wl_wvrng, outfile=outfile,
                                                        spec_subpixel=self.spec_subpixel,
                                                        spat_subpixel=self.spat_subpixel,
                                                        slice_subpixel=self.slice_subpixel,
                                                        skip_subpix_weights=self.skip_subpix_weights,
                                                        correct_dar=self.correct_dar)
                    # Prepare the header
                    hdr = self.all_wcs[ff].to_header()
                    if self.fluxcal:
                        hdr['FLUXUNIT'] = (flux_calib.PYPEIT_FLUX_SCALE, "Flux units -- erg/s/cm^2/Angstrom/arcsec^2")
                    else:
                        hdr['FLUXUNIT'] = (1, "Flux units -- counts/s/Angstrom/arcsec^2")
                    # Write out the datacube
                    msgs.info("Saving datacube as: {0:s}".format(outfile))
                    final_cube = DataCube(flxcube, sigcube, bpmcube, wave, self.specname, self.blaze_wave, self.blaze_spec,
                                          sensfunc=None, fluxed=self.fluxcal)
                    final_cube.to_file(outfile, hdr=hdr, overwrite=self.overwrite)
                # No need to proceed and store arrays - we are writing individual datacubes
                continue

            # Store the information if we are combining multiple frames
            self.all_sci.append(sciImg.copy())
            self.all_ivar.append(ivar.copy())
            self.all_wave.append(waveimg.copy())
            self.all_ra.append(ra_img.copy())
            self.all_dec.append(dec_img.copy())
            self.all_slitid.append(slitid_img_gpm.copy())
            self.all_wghts.append(wghts.copy())
            self.all_tilts.append(spec2DObj.tilts)
            self.all_slits.append(slits)
            self.all_align.append(alignSplines)
            self.all_dar.append(darcorr)

    def run_align(self):
        """
        This routine aligns multiple cubes by using manual input offsets or by cross-correlating white light images.

        Returns:
            `numpy.ndarray`_: A new set of RA values that have been aligned
            `numpy.ndarray`_: A new set of Dec values that has been aligned
        """
        # Grab cos(dec) for convenience
        cosdec = np.cos(np.mean(self.ifu_dec[0]) * np.pi / 180.0)
        # Initialize the RA and Dec offset arrays
        ra_offsets, dec_offsets = [0.0]*self.numfiles, [0.0]*self.numfiles
        # Register spatial offsets between all frames
        if self.user_alignment:
            # The user has specified offsets - update these values accounting for the difference in header RA/DEC
            ra_offsets, dec_offsets = datacube.align_user_offsets(self.ifu_ra, self.ifu_dec,
                                                                  self.ra_offsets, self.dec_offsets)
        else:
            # Find the wavelength range where all frames overlap
            min_wl, max_wl = datacube.get_whitelight_range(np.max(self.mnmx_wv[:, :, 0]),  # The max blue wavelength
                                                           np.min(self.mnmx_wv[:, :, 1]),  # The min red wavelength
                                                           self.cubepar['whitelight_range'])  # The user-specified values (if any)
            # Get the good white light pixels
            slitid_img_gpm, wavediff = datacube.get_whitelight_pixels(self.all_wave, self.all_slitid, min_wl, max_wl)
            # Iterate over white light image generation and spatial shifting
            numiter = 2
            for dd in range(numiter):
                msgs.info(f"Iterating on spatial translation - ITERATION #{dd+1}/{numiter}")
                # Generate the WCS
                image_wcs, voxedge, reference_image = \
                    datacube.create_wcs(self.all_ra, self.all_dec, self.all_wave, slitid_img_gpm, self._dspat, wavediff,
                                        ra_offsets=ra_offsets, dec_offsets=dec_offsets,
                                        ra_min=self.cubepar['ra_min'], ra_max=self.cubepar['ra_max'],
                                        dec_min=self.cubepar['dec_min'], dec_max=self.cubepar['dec_max'],
                                        wave_min=self.cubepar['wave_min'], wave_max=self.cubepar['wave_max'],
                                        reference=self.cubepar['reference_image'], collapse=True, equinox=2000.0,
                                        specname=self.specname)
                if voxedge[2].size != 2:
                    msgs.error("Spectral range for WCS is incorrect for white light image")

                wl_imgs = datacube.generate_image_subpixel(image_wcs, voxedge, self.all_sci, self.all_ivar, self.all_wave,
                                                           slitid_img_gpm, self.all_wghts, self.all_wcs,
                                                           self.all_tilts, self.all_slits, self.all_align, self.all_dar,
                                                           ra_offsets, dec_offsets,
                                                           spec_subpixel=self.spec_subpixel,
                                                           spat_subpixel=self.spat_subpixel,
                                                           slice_subpixel=self.slice_subpixel)
                if reference_image is None:
                    # ref_idx will be the index of the cube with the highest S/N
                    ref_idx = np.argmax(self.weights)
                    reference_image = wl_imgs[:, :, ref_idx].copy()
                    msgs.info("Calculating spatial translation of each cube relative to cube #{0:d})".format(ref_idx+1))
                else:
                    msgs.info("Calculating the spatial translation of each cube relative to user-defined 'reference_image'")

                # Calculate the image offsets relative to the reference image
                for ff in range(self.numfiles):
                    # Calculate the shift
                    ra_shift, dec_shift = calculate_image_phase(reference_image.copy(), wl_imgs[:, :, ff], maskval=0.0)
                    # Convert pixel shift to degrees shift
                    ra_shift *= self._dspat/cosdec
                    dec_shift *= self._dspat
                    msgs.info("Spatial shift of cube #{0:d}:".format(ff + 1) + msgs.newline() +
                              "RA, DEC (arcsec) = {0:+0.3f} E, {1:+0.3f} N".format(ra_shift*3600.0, dec_shift*3600.0))
                    # Store the shift in the RA and DEC offsets in degrees
                    ra_offsets[ff] += ra_shift
                    dec_offsets[ff] += dec_shift
        return ra_offsets, dec_offsets

    def compute_weights(self):
        """
        Compute the relative weights to apply to pixels that are collected into the voxels of the output DataCubes

        Returns:
            `numpy.ndarray`_: The individual pixel weights for each detector pixel, and every frame.
        """
        # If there is only one file, then all pixels have the same weight
        if self.numfiles == 1:
            return np.ones_like(self.all_sci)

        # Calculate the relative spectral weights of all pixels
        return datacube.compute_weights_frompix(self.all_ra, self.all_dec, self.all_wave, self.all_sci, self.all_ivar,
                                                self.all_slitid, self._dspat, self._dwv, self.mnmx_wv, self.all_wghts,
                                                self.all_wcs, self.all_tilts, self.all_slits, self.all_align, self.all_dar,
                                                self.ra_offsets, self.dec_offsets,
                                                ra_min=self.cubepar['ra_min'], ra_max=self.cubepar['ra_max'],
                                                dec_min=self.cubepar['dec_min'], dec_max=self.cubepar['dec_max'],
                                                wave_min=self.cubepar['wave_min'], wave_max=self.cubepar['wave_max'],
                                                weight_method=self.cubepar['weight_method'],
                                                whitelight_range=self.cubepar['whitelight_range'],
                                                reference_image=self.cubepar['reference_image'],
                                                correct_dar=self.correct_dar,
                                                specname=self.specname)

    def run(self):
        """
        This is the main routine called to convert PypeIt spec2d files into PypeIt DataCube objects. It is specific
        to the SlicerIFU data.

        First the data are loaded and several corrections are made. These include:

        * A sky frame or model is subtracted from the science data, and the relative spectral illumination
          of different slices is corrected.
        * A mask of good pixels is identified
        * A common spaxel scale is determined, and the astrometric correction is derived
        * An RA and Dec image is created for each pixel.
        * Based on atmospheric conditions, a differential atmospheric refraction correction is applied.
        * Extinction correction
        * Flux calibration (optional - this calibration is only applied if a standard star cube is supplied)

        If the input frames will not be combined (combine=False) if they won't be aligned (align=False), then
        each individual spec2d file is converted into a spec3d file (i.e. a PypeIt DataCube object). These fits
        files can be loaded/viewed in other software packages to display or combine multiple datacubes into a
        single datacube. However, note that different software packages use combination algorithms that may not
        conserve flux, or may produce covariance between adjacent voxels.

        If the user wishes to either spatially align multiple exposures (align=True) or combine multiple
        exposures (combine=True), then the next set of operations include:

        * Generate white light images of each individual cube (according to a user-specified wavelength range)
        * Align multiple frames if align=True (either manually by user input, or automatically by cross-correlation)
        * Create the output WCS, and apply the flux calibration to the data
        * Generate individual datacubes (combine=False) or one master datacube containing all exposures (combine=True).
          Note, there are several algorithms used to combine multiple frames. Refer to the subpixellate() routine for
          more details about the combination options.
        """
        # No need to continue if we are not combining nor aligning frames
        if not self.combine and not self.align:
            return

        # If the user is aligning or combining, the spatial scale of the output cubes needs to be consistent.
        # Set the spatial and spectral scales of the output datacube
        self._dspat, self._dwv = datacube.set_voxel_sampling(self._spatscale, self._specscale,
                                                             dspat=self._dspat, dwv=self._dwv)

        # Align the frames
        if self.align:
            self.ra_offsets, self.dec_offsets = self.run_align()

        # Compute the relative weights on the spectra
        self.all_wghts = self.compute_weights()

        # Generate the WCS, and the voxel edges
        cube_wcs, vox_edges, _ = \
            datacube.create_wcs(self.all_ra, self.all_dec, self.all_wave, self.all_slitid, self._dspat, self._dwv,
                                ra_offsets=self.ra_offsets, dec_offsets=self.dec_offsets,
                                ra_min=self.cubepar['ra_min'], ra_max=self.cubepar['ra_max'],
                                dec_min=self.cubepar['dec_min'], dec_max=self.cubepar['dec_max'],
                                wave_min=self.cubepar['wave_min'], wave_max=self.cubepar['wave_max'],
                                reference=self.cubepar['reference_image'], collapse=False, equinox=2000.0,
                                specname=self.specname)

        sensfunc = None
        if self.flux_spline is not None:
            # Get wavelength of each pixel
            numwav = vox_edges[2].size - 1
            wcs_scale = (1.0 * cube_wcs.spectral.wcs.cunit[0]).to(units.Angstrom).value  # Ensures the WCS is in Angstroms
            senswave = wcs_scale * cube_wcs.spectral.wcs_pix2world(np.arange(numwav), 0)[0]
            sensfunc = self.flux_spline(senswave)

        # Generate a datacube
        if self.method in ['subpixel', 'ngp']:
            # Generate the datacube
            wl_wvrng = None
            if self.cubepar['save_whitelight']:
                wl_wvrng = datacube.get_whitelight_range(np.max(self.mnmx_wv[:, :, 0]),
                                                np.min(self.mnmx_wv[:, :, 1]),
                                                self.cubepar['whitelight_range'])
            if self.combine:
                outfile = datacube.get_output_filename("", self.cubepar['output_filename'], True, -1)
                # Generate the datacube
                flxcube, sigcube, bpmcube, wave = \
                    datacube.generate_cube_subpixel(cube_wcs, vox_edges, self.all_sci, self.all_ivar, self.all_wave,
                                                    self.all_slitid, self.all_wghts, self.all_wcs,
                                                    self.all_tilts, self.all_slits, self.all_align, self.all_dar,
                                                    self.ra_offsets, self.dec_offsets,
                                                    outfile=outfile, overwrite=self.overwrite, whitelight_range=wl_wvrng,
                                                    spec_subpixel=self.spec_subpixel,
                                                    spat_subpixel=self.spat_subpixel,
                                                    slice_subpixel=self.slice_subpixel,
                                                    skip_subpix_weights=self.skip_subpix_weights,
                                                    correct_dar=self.correct_dar)
                # Prepare the header
                hdr = cube_wcs.to_header()
                if self.fluxcal:
                    hdr['FLUXUNIT'] = (flux_calib.PYPEIT_FLUX_SCALE, "Flux units -- erg/s/cm^2/Angstrom/arcsec^2")
                else:
                    hdr['FLUXUNIT'] = (1, "Flux units -- counts/s/Angstrom/arcsec^2")
                # Write out the datacube
                msgs.info("Saving datacube as: {0:s}".format(outfile))
                final_cube = DataCube(flxcube, sigcube, bpmcube, wave, self.specname, self.blaze_wave, self.blaze_spec,
                                      sensfunc=sensfunc, fluxed=self.fluxcal)
                final_cube.to_file(outfile, hdr=hdr, overwrite=self.overwrite)
            else:
                for ff in range(self.numfiles):
                    outfile = datacube.get_output_filename("", self.cubepar['output_filename'], False, ff)
                    # Generate the datacube
                    flxcube, sigcube, bpmcube, wave = \
                        datacube.generate_cube_subpixel(cube_wcs, vox_edges,
                                                        self.all_sci[ff], self.all_ivar[ff], self.all_wave[ff],
                                                        self.all_slitid[ff], self.all_wghts[ff], self.all_wcs[ff],
                                                        self.all_tilts[ff], self.all_slits[ff], self.all_align[ff], self.all_dar[ff],
                                                        self.ra_offsets[ff], self.dec_offsets[ff],
                                                        overwrite=self.overwrite, whitelight_range=wl_wvrng,
                                                        outfile=outfile, spec_subpixel=self.spec_subpixel,
                                                        spat_subpixel=self.spat_subpixel,
                                                        slice_subpixel=self.slice_subpixel,
                                                        skip_subpix_weights=self.skip_subpix_weights,
                                                        correct_dar=self.correct_dar)
                    # Prepare the header
                    hdr = cube_wcs.to_header()
                    if self.fluxcal:
                        hdr['FLUXUNIT'] = (flux_calib.PYPEIT_FLUX_SCALE, "Flux units -- erg/s/cm^2/Angstrom/arcsec^2")
                    else:
                        hdr['FLUXUNIT'] = (1, "Flux units -- counts/s/Angstrom/arcsec^2")
                    # Write out the datacube
                    msgs.info("Saving datacube as: {0:s}".format(outfile))
                    final_cube = DataCube(flxcube, sigcube, bpmcube, wave, self.specname, self.blaze_wave, self.blaze_spec,
                                          sensfunc=sensfunc, fluxed=self.fluxcal)
                    final_cube.to_file(outfile, hdr=hdr, overwrite=self.overwrite)
