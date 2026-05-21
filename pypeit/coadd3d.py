"""
Module containing routines used by 3D datacubes.

.. include:: ../include/links.rst
"""

from pathlib import Path
import inspect

from astropy import wcs, units
from astropy.io import fits
import erfa
from IPython import embed
import numpy as np
from scipy.interpolate import interp1d

from pypeit import log
from pypeit import PypeItError
from pypeit import alignframe, datamodel, flatfield, io, sensfunc, spec2dobj, utils
from pypeit.core.flexure import calculate_image_phase
from pypeit.core import datacube, extract, flux_calib, parse, combine 
from pypeit.spectrographs.util import load_spectrograph
from pypeit.manual_extract import ManualCubeExtractionObj


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
                 '_wcs'
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
        self.head0 = None  # This contains the primary header of the spec2d used to make the datacube

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
        # Add them in
        for key in subheader:
            primary_hdr[key] = subheader[key]
        # Set the exposure time to 1, since datacubes are counts/second.
        # This is needed for the sensitivity function calculation
        primary_hdr['EXPTIME'] = 1.0
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
            self.head0 = hdu[0].header
            # Meta
            self.spectrograph = load_spectrograph(self.PYP_SPEC)
            self.spect_meta = self.spectrograph.parse_spec_header(hdu[0].header)
            self._ivar = None
            self._wcs = wcs.WCS(hdu[1].header)
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

    def extract_spec(self, parset, output_dir=None, overwrite=False, debug=False):
        """
        Extract a spectrum from the datacube

        Parameters
        ----------
        parset : dict
            A dictionary containing the :class:`~pypeit.par.pypeitpar.ReducePar` parameters.
        boxcar_radius : float, optional
            Radius of the circular boxcar (in arcseconds) to use for the extraction
        overwrite : bool, optional
            Overwrite any existing files
        output_dir : str, optional
            The directory for the output files. If None, the output files are written to the
            directory in which the class is run.
        """
        _output_dir = Path().absolute() if output_dir is None else Path(output_dir).absolute()
        
        # Check if the files exist, if so crash out
        file_suffix = (
            Path(self.filename).name if parset['cube']['extraction']['output_filename'] is None
            else parset['cube']['extraction']['output_filename']
        )
        spec1d_filename = _output_dir / f'spec1d_{file_suffix}'
        if spec1d_filename.suffix != '.fits':
            spec1d_filename = spec1d_filename.with_suffix('.fits')
        if spec1d_filename.is_file() and not overwrite:
            raise PypeItError(f"{spec1d_filename} exists!  Set overwrite=True to overwrite.")

        spec2d_filename = _output_dir / f'spec2d_{file_suffix}'
        if spec2d_filename.suffix != '.fits':
            spec2d_filename = spec2d_filename.with_suffix('.fits')
        if spec2d_filename.is_file() and not overwrite:
            raise PypeItError(f"{spec2d_filename} exists!  Set overwrite=True to overwrite.")

        out_whitelight = Path(
            datacube.get_output_whitelight_filename(str(_output_dir), file_suffix)
        ).absolute()
        if out_whitelight.is_file() and not overwrite:
            raise PypeItError(f"{out_whitelight} exists!  Set overwrite=True to overwrite.")

        if parset['cube']['extraction']['manual'] is not None and len(parset['cube']['extraction']['manual']) > 0:
            manual_dict= ManualCubeExtractionObj.parse(parset['cube']['extraction']['manual']).to_dict()
            manual_position = (manual_dict['spatx'][0], manual_dict['spaty'][0])
        else:
            manual_position = None

        # Datacube's are counts/second, so set the exposure time to 1
        exptime = 1.0
        sobjs, spec2d, wl_img, wl_ivar, wl_gpm = datacube.extract_point_source(
            self.wave, self.flux, self.ivar, self.bpm, self._wcs, exptime,
            fluxed=self.fluxed, min_frac_use=parset['extraction']['min_frac_prof'],
            whitelight_range=parset['cube']['extraction']['whitelight_range'],
            fwhm=parset['cube']['extraction']['fwhm'],
            skysub_resid=parset['cube']['extraction']['skysub_resid'],
            snr_thresh=parset['cube']['extraction']['snr_thresh'], manual_position=manual_position,
            boxcar_radius=parset['cube']['extraction']['boxcar_radius'], 
            opt_prof_method=parset['cube']['extraction']['opt_prof_method'],
            spectrograph=self.spectrograph, show_qa=debug
        )

        # Save the extracted spectrum
        sobjs.write_to_fits(self.head0, str(spec1d_filename), overwrite=overwrite)
        # Save the psuedo spec2d images
        all_spec2d = spec2dobj.AllSpec2DObj()
        all_spec2d[spec2d.detector.name] = spec2d
        # Build header for spec2d
        all_spec2d.write_to_fits(str(spec2d_filename), pri_hdr=fits.Header(), overwrite=overwrite)
        # Write out the white light image
        # TODO This is replicated code from datacube.make_whitelight, clean this up.
        # TODO : Note that this overwrites the whitelight generated from the main code.
        log.info(f"Saving white light image as: {out_whitelight}")
        primary_hdu = fits.PrimaryHDU(wl_img, header=self._wcs.celestial.to_header())
        ivar_hdu = fits.ImageHDU(wl_ivar, name='IVAR')
        gpm_hdu = fits.ImageHDU(wl_gpm.astype(np.uint8), name='GPM')
        hdul = fits.HDUList([primary_hdu, ivar_hdu, gpm_hdu])
        hdul.writeto(out_whitelight, overwrite=overwrite)


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
                The atmospheric pressure during the observations in mbar. Valid range is from 100 mbar - 1400 mbar.
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
        log.info("Preparing the parameters for the DAR correction")

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
        log.info(
            "DAR correction parameters:\n"
            f"   Airmass = {self.airmass:.2f}\n"
            f"   Pressure = {self.pressure.to_value(units.mbar):.2f} mbar\n"
            f"   Humidity = {self.humidity*100.0:.2f} %\n"
            f"   Temperature = {self.temperature.to_value(units.deg_C):.2f} deg C\n"
            f"   Reference wavelength = {self.wave_ref.to_value(units.Angstrom):.2f} Angstrom"
        )

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

    # TODO Make parangle and cosdec arguments to this function rather than class attributes
    # required upon init, since they are not required to calculate the DAR dispersion, only 
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

    Parameters
    ----------
    spec2dfiles : :obj:`list`
        List of all spec2D files
    par : :class:`~pypeit.par.pypeitpar.PypeItPar`
        An instance of the parameter set.  If None, assumes that detector 1 is
        the one reduced and uses the default reduction parameters for the
        spectrograph (see
        :func:`~pypeit.spectrographs.spectrograph.Spectrograph.default_pypeit_par`
        for the relevant spectrograph class).
    redux_path : :obj: str, optional
        The top-level directory for all output files.  If None, this is set to
        the current working directory.  Output files may also be put into
        subdirectories within this top-level directory, like Science_cube and
        QA_cube.
    skysub_frame : :obj:`list`, optional
        A list of frames to use for the sky subtraction of each spec2D file; the
        list length should match ``spec2dfiles``.  Ignored if None.
    sensfile : :obj:`list`, optional
        A list of frames to use for the sensitivity function of each spec2D
        file; the list length should match ``spec2dfiles``.  Ignored if None.
    scale_corr : :obj:`list`, optional
        A list of relative scale corrections to use for each spec2D file; the
        list length should match ``spec2dfiles``.  Ignored if None.
    grating_corr : :obj:`list`, optional
        A list of strings with the relative path to the flat calibration file
        used to reduce each spec2D file; the list length should match
        ``spec2dfiles``.  Ignored if None.
    ra_offsets : :obj:`list`, optional
        A list of relative RA offsets of each spec2D file; the list length
        should match ``spec2dfiles``.  Ignored if None.
    dec_offsets : :obj:`list`, optional
        A list of relative Dec offsets of each spec2D file; the list length
        should match ``spec2dfiles``.  Ignored if None.
    spectrograph : :obj:`str`, :class:`~pypeit.spectrographs.spectrograph.Spectrograph`, optional
        The name or instance of the spectrograph used to obtain the data.  If
        None, this is pulled from the header of the first spec2D file provided.
    det : :obj:`int`_, optional
        Detector index
    overwrite : :obj:`bool`, optional
        Overwrite existing output files
    show : :obj:`bool`, optional
        Show results in ginga
    debug : :obj:`bool`, optional
        Run in debugging mode
    """

    @classmethod
    def get_instance(
        cls, spec2dfiles, par, redux_path=None, skysub_frame=None, sensfile=None, scale_corr=None,
        grating_corr=None, ra_offsets=None, dec_offsets=None, spectrograph=None, det=1,
        overwrite=False, show=False, debug=False
    ):
        """
        Instantiate the subclass appropriate for the provided spectrograph.

        The class to instantiate must match the ``pypeline``
        attribute of the provided ``spectrograph``, and must be a
        subclass of :class:`CoAdd3D`; see the parent class
        instantiation for parameter descriptions.

        Returns
        -------
        :class:`CoAdd3D`
            The :class:`CoAdd3D` subclass instance relevant to the spectrograph
            used to obtain the data.
        """
        return next(
            c for c in cls.__subclasses__() if c.__name__ == (spectrograph.pypeline + 'CoAdd3D')
        )(
            spec2dfiles, par, redux_path=redux_path, skysub_frame=skysub_frame,
            sensfile=sensfile, scale_corr=scale_corr, grating_corr=grating_corr,
            ra_offsets=ra_offsets, dec_offsets=dec_offsets, spectrograph=spectrograph,
            det=det, overwrite=overwrite, show=show, debug=debug
        )

    # TODO: det is needed to load the spec2d files. Although, what if there are
    # multiple detectors?  Do we need a det for each spec2d file?
    def __init__(
        self, spec2dfiles, par, redux_path=None, skysub_frame=None, sensfile=None, scale_corr=None,
        grating_corr=None, ra_offsets=None, dec_offsets=None, spectrograph=None, det=None,
        overwrite=False, show=False, debug=False
    ):
        # TODO: Consider loading all calibrations into a single variable within
        # the main CoAdd3D parent class.

        # Set the variables
        self.spec2d = spec2dfiles
        self.numfiles = len(spec2dfiles)
        self.par = par
        self.redux_path = Path().absolute() if redux_path is None else Path(redux_path).absolute()
        self.overwrite = overwrite
        self.chk_version = self.par['rdx']['chk_version']

        # Get the paths
        self.scidir, self.qadir = map(
            lambda x : Path(x).absolute(), CoAdd3D.output_paths(
                spec2dfiles, par['rdx']['scidir'], par['rdx']['qadir'], coadd_dir=self.redux_path
            )
        )

        # Write the par to disk
        # TODO: Add an argument that sets whether or not this should be saved?
        par_outfile = (
            self.scidir.parent / f"{par['reduce']['cube']['output_filename']}_datacube.par"
        )
        log.info(f'Writing full parameter set to {par_outfile}.')
        par.to_config(par_outfile, exclude_defaults=True, include_descr=False)

        # Extract some parsets for simplicity
        self.cubepar = self.par['reduce']['cube']
        self.flatpar = self.par['calibrations']['flatfield']
        self.senspar = self.par['sensfunc']

        # Extract some commonly used variables
        self.method = self.cubepar['method']
        self.combine = self.cubepar['combine']
        self.alignment_method = None if self.cubepar['alignment_method'].lower() in ["none"] else self.cubepar['alignment_method']
        self.align = False if self.alignment_method is None else True
        self.native = self.cubepar['save_native']
        self.correct_dar = self.cubepar['correct_dar']

        # TODO: Only need one of show or debug probably
        self.show = show
        self.debug = debug 

        # Check the input
        for name, var in zip(
            [
                'sky-subtraction frames', 'sensitivity function files',
                'relative scale correction options', 'flat calibration files', 'RA offsets',
                'Dec offsets'
            ], [skysub_frame, sensfile, scale_corr, grating_corr, ra_offsets, dec_offsets]
        ):
            if var is not None and len(var) != self.numfiles:
                raise PypeItError(
                    f'The list of {name} should have the same length as the spec2dfiles list'
                )

        # Make sure both ra_offsets and dec_offsets are either both None or both lists
        if (
            (ra_offsets is None and dec_offsets is not None)
            or (ra_offsets is not None and dec_offsets is None)
        ):
            raise PypeItError(
                'You must provide both ra_offsets and dec_offsets if you provide either.'
            )

        # Set the frame specific options
        self.sensfile = None
        if sensfile is None:
            # User didn't provide a sensfile for each frame. Check if they provided a single one.
            if self.cubepar['sensfile'] is not None:
                # User provided a single sensfile. Use this for all frames.
                self.sensfile = self.numfiles*[self.cubepar['sensfile']]
        else:
            # User provided a sensfile for each frame. Use these.
            self.sensfile = sensfile

        self.skysub_frame = skysub_frame
        self.scale_corr = scale_corr
        self.grating_corr = grating_corr
        self.ra_offsets = list(ra_offsets) if isinstance(ra_offsets, np.ndarray) else ra_offsets
        self.dec_offsets = list(dec_offsets) if isinstance(dec_offsets, np.ndarray) else dec_offsets

        # If there is only one frame being "combined" AND there's no reference
        # image, then don't compute the translation.
        if self.numfiles == 1 and self.cubepar["reference_image"] is None:
            if self.align:
                log.warning(
                    "Parameter 'alignment_method' should be 'none' when there is only one frame and no reference image"
                )
                log.info("Setting 'alignment_method' to 'none'")
            self.align = False
            self.alignment_method = None
        if self.ra_offsets is not None:
            if self.alignment_method is None:
                log.warning(
                    "You have provided ra_offsets and dec_offsets, but set 'alignment_method' to 'none'. "
                    "The datacubes will not be aligned. If you wish to align the datacubes using the "
                    "user-specified offsets, please set 'alignment_method' to 'user'.")
                self.align = False
                self.alignment_method = None
            elif self.alignment_method != "user":
                log.warning("When 'ra_offset' and 'dec_offset' are set, 'alignment_method' must be 'user'.")
                log.info("Setting 'alignment_method' to 'user'")
                self.align = True
                self.alignment_method = "user"
        # If no ra_offsets or dec_offsets have been provided, initialize the lists
        if self.ra_offsets is None and self.dec_offsets is None:
            log.info("No RA or Dec offsets have been provided.")
            if self.align:
                log.info(f"An automatic alignment will be performed the {self.alignment_method} method.")
            # Initialise the lists of ra_offsets and dec_offsets
            self.ra_offsets = [0.0]*self.numfiles
            self.dec_offsets = [0.0]*self.numfiles
        if self.grating_corr is None:
            self.grating_corr = [None] * self.numfiles

        # Initialize the spectrograph
        if spectrograph is None:
            with fits.open(spec2dfiles[0]) as hdu:
                spectrograph = hdu[0].header['PYP_SPEC']
        self.spectrograph = load_spectrograph(spectrograph)
        self.specname = self.spectrograph.name

        # Initialise arrays for storage.
        # - The RA and Dec at the centre of the IFU, as stored in the header
        self.ifu_ra, self.ifu_dec = np.array([]), np.array([])
        # - Others
        self.all_sci, self.all_ivar, self.all_wave, self.all_slitid, self.all_wghts = [], [], [], [], []
        self.all_tilts, self.all_slits, self.all_align, self.all_header = [], [], [], []
        self.all_wcs, self.all_ra, self.all_dec, self.all_deltapix, self.all_dar = [], [], [], [], []
        # - Weights to use when combining cubes
        self.weights = np.ones(self.numfiles)

        # Set the spatial sampling.  Use "/3600.0" to convert to degrees.
        self._dspat = (
            None if self.cubepar['spatial_delta'] is None
            else self.cubepar['spatial_delta'] / 3600.0
        )
        # Linear wavelength sampling in Angstroms
        self._dwv = self.cubepar['wave_delta']

        # TODO: The default behaviour (combine=False, align=False) produces a
        # datacube that uses the instrument WCS.  It should be possible (and
        # perhaps desirable) to do a spatial alignment (i.e. align=True), apply
        # this to the RA,Dec values of each pixel, and then use the instrument
        # WCS to save the output (or, just adjust the crval).  At the moment, if
        # the user wishes to spatially align the frames, a different WCS is
        # generated.

        # Determine what method is requested
        self.spec_subpixel, self.spat_subpixel, self.slice_subpixel = 1, 1, 1
        self.skip_subpix_weights = True
        if self.method == "subpixel":
            self.spec_subpixel = self.cubepar['spec_subpixel']
            self.spat_subpixel = self.cubepar['spat_subpixel']
            self.slice_subpixel = self.cubepar['slice_subpixel']
            self.skip_subpix_weights = False
            log.info(
                'Adopting the subpixel algorithm to generate the datacube, with the following '
                f'subpixellation scales:  Spectral: {self.spec_subpixel}; '
                f'Spatial: {self.spat_subpixel}; Slices: {self.slice_subpixel}.'
            )
        elif self.method == "ngp":
            log.info("Adopting the nearest grid point (NGP) algorithm to generate the datacube.")
            self.skip_subpix_weights = True
        else:
            raise PypeItError(f"The following datacube method is not allowed: {self.method}")

        # Get the detector number and string representation
        if det is None:
            det = 1 if self.par['rdx']['detnum'] is None else self.par['rdx']['detnum']
        self.detname = self.spectrograph.get_det_name(det)

        # Check if the output file exists
        self.check_outputs()

        # Check the reference cube and image exist, if requested
        self.fluxcal = self.sensfile is not None
        self.blaze_wave, self.blaze_spec = None, None
        self.blaze_spline, self.flux_spline = None, None
        self.flat_splines = dict()  # A dictionary containing the splines of the flatfield

        # If a reference image has been set, check that it exists
        if self.cubepar['reference_image'] is not None:
            if not Path(self.cubepar['reference_image']).is_file():
                raise PypeItError(
                    f"Reference image does not exist: {self.cubepar['reference_image']}"
                )

        # Load the default scaleimg frame for the scale correction
        self.scalecorr_default = "none"
        self.relScaleImgDef = np.array([1])
        self.set_default_scalecorr()

        # Load the default sky frame to be used for sky subtraction
        self.skysub_default = "image"
        # This is the default behaviour (i.e. to use the "image" for the sky subtraction)
        self.skyImgDef, self.skySclDef = None, None
        self.set_default_skysub()

    @staticmethod
    def output_paths(spec2d_files, science_dir, qa_dir, coadd_dir=None):
        """
        Construct the names and ensure the existence of the science and QA output directories.

        Args:
            spec2d_files (:obj:`list`):
                The list of PypeIt spec2d files to be coadded.  The top-level
                directory for the coadd3d output directories is assumed to be
                same as used by the basic reductions.  For example, if one of
                the spec2d files is
                ``/path/to/reductions/Science/spec2d_file.fits``, the parent
                directory for the coadd2d directories is
                ``/path/to/reductions/``.
            science_dir (:obj:`str`):
                The name of the science directory to use for the coadd3d output.
                 For example, if scidir is "Science", the science output directory will be
                ``/path/to/reductions/Science_cube/``.
            qa_dir (:obj:`str`):
                The name of the QA directory to use for the coadd3d output.  For
                example, if qadir is "QA", the QA output directory will be
                ``/path/to/reductions/QA_cube/``.
            coadd_dir (:obj:`str`, optional):
                Path to the directory to use for the coadd3d output.
                If None, the parent of the science directory is used.

        Returns:
            :obj:`tuple`: Two strings with the names of (1) the science output
            directory and (2) the QA output directory.  The function also
            creates both directories if they do not exist.
        """
        # Science output directory
        if coadd_dir is not None:
            pypeit_scidir = Path(coadd_dir).absolute() / 'Science'
        else:
            pypeit_scidir = Path(spec2d_files[0]).parent
        coadd_scidir = pypeit_scidir.parent / f'{science_dir}_cube'
        if not coadd_scidir.exists():
            coadd_scidir.mkdir(parents=True)
        # QA directory
        qa_path = pypeit_scidir.parent / f'{qa_dir}_cube' / 'PNGs'
        if not qa_path.exists():
            qa_path.mkdir(parents=True)
        return str(coadd_scidir), str(qa_path)

    def check_outputs(self):
        """
        Check if any of the intended output files already exist.

        This check should be done near the beginning of the coaddition, to avoid
        any computation that won't be saved in the event that files won't be
        overwritten.
        """
        if self.combine:
            outfile = datacube.get_output_filename(
                str(self.scidir), "", self.cubepar['output_filename'], self.combine
            )
            if Path(outfile).is_file() and not self.overwrite:
                raise PypeItError(
                    f"{outfile} exists!  Use overwrite flag or parameter to overwrite."
                )

            out_whitelight = datacube.get_output_whitelight_filename(str(self.scidir), outfile)
            if (
                Path(out_whitelight).is_file() and self.cubepar['save_whitelight']
                and not self.overwrite
            ):
                raise PypeItError(
                    f"{out_whitelight} exists!  Use overwrite flag or parameter to overwrite, or "
                    "use the save_whitelight parameter to skip saving the image."
                )
            return

        # Finally, if there's just one file, check if the output filename is given
        if self.numfiles == 1 and self.cubepar['output_filename'] != "":
            outfile = datacube.get_output_filename(
                str(self.scidir), "", self.cubepar['output_filename'], True
            )
            if Path(outfile).is_file() and not self.overwrite:
                raise PypeItError(
                    f"{outfile} exists!  Use overwrite flag or parameter to overwrite."
                )

            out_whitelight = datacube.get_output_whitelight_filename(
                str(self.scidir), outfile
            )
            if (
                Path(out_whitelight).is_file() and self.cubepar['save_whitelight']
                and not self.overwrite
            ):
                raise PypeItError(
                    f"{out_whitelight} exists!  Use overwrite flag or parameter to overwrite, "
                    "or use the save_whitelight parameter to skip saving the image."
                )
            return

        for ff in range(self.numfiles):
            # Check native first
            if self.native:
                outfile = datacube.get_output_filename(
                    str(self.scidir), self.spec2d[ff], self.cubepar['output_filename'],
                    self.combine, native=True, idx=ff+1
                )
                if Path(outfile).is_file() and not self.overwrite:
                    raise PypeItError(
                        f"{outfile} exists!  Use overwrite flag or parameter to overwrite."
                    )
                # Now check the whitelight of the native sampling files
                out_whitelight = datacube.get_output_whitelight_filename(str(self.scidir), outfile)
                if (
                    Path(out_whitelight).is_file() and self.cubepar['save_whitelight']
                    and not self.overwrite
                ):
                    raise PypeItError(
                        f"{out_whitelight} exists!  Use overwrite flag or parameter to overwrite, "
                        "or use the save_whitelight parameter to skip saving the image."
                    )
            # Check not native
            outfile = datacube.get_output_filename(
                str(self.scidir), self.spec2d[ff], self.cubepar['output_filename'],
                self.combine, native=False, idx=ff+1
            )
            if Path(outfile).is_file() and not self.overwrite:
                raise PypeItError(
                    f"{outfile} exists!  Use overwrite flag or parameter to overwrite."
                )

            out_whitelight = datacube.get_output_whitelight_filename(str(self.scidir), outfile)
            if (
                Path(out_whitelight).is_file() and self.cubepar['save_whitelight']
                and not self.overwrite
            ):
                raise PypeItError(
                    f"{out_whitelight} exists!  Use overwrite flag or parameter to overwrite, "
                    "or use the save_whitelight parameter to skip saving the image."
                )

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

    def set_default_scalecorr(self):
        """
        Set the default mode to use for relative spectral scale correction.
        """
        if self.cubepar['scale_corr'] is not None:
            if self.cubepar['scale_corr'] == "image":
                log.info("The default relative spectral illumination correction will use the science image")
                self.scalecorr_default = "image"
            else:
                log.info(
                    "Loading default scale image for relative spectral illumination correction:\n"
                    +self.cubepar['scale_corr']
                )
                try:
                    spec2DObj = spec2dobj.Spec2DObj.from_file(self.cubepar['scale_corr'],
                                                              self.detname,
                                                              chk_version=self.chk_version)
                except Exception as e:
                    log.warning(f'Loading spec2d file raised {type(e).__name__}:\n{str(e)}')
                    log.warning(
                        "Could not load scaleimg from spec2d file:\n" +
                        self.cubepar['scale_corr'] +
                        "\nscale correction will not be performed unless you have specified the "
                        "correct\nscale_corr file in the spec2d block")
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
                log.info(
                    "Loading the following frame for the relative spectral illumination "
                    "correction:\n" + scalecorr
                )
                try:
                    spec2DObj_scl = spec2dobj.Spec2DObj.from_file(scalecorr, self.detname,
                                                                  chk_version=self.chk_version)
                except Exception as e:
                    log.warning(f'Loading spec2d file raised {type(e).__name__}:\n{str(e)}')
                    raise PypeItError("Could not load skysub image from spec2d file:\n" + scalecorr)
                else:
                    relScaleImg = spec2DObj_scl.scaleimg
                    this_scalecorr = scalecorr
        if this_scalecorr == "none":
            log.info("Relative spectral illumination correction will not be performed.")
        else:
            log.info(
                "Using the following frame for the relative spectral illumination correction:\n"
                + this_scalecorr
            )
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
            log.info("The sky model in the spec2d science frames will be used for sky "
                      "subtraction\n(unless specific skysub frames have been specified)")
            self.skysub_default = "image"
        else:
            log.info("Loading default image for sky subtraction:\n"
                      + self.cubepar['skysub_frame'])
            try:
                spec2DObj = spec2dobj.Spec2DObj.from_file(self.cubepar['skysub_frame'],
                                                          self.detname,
                                                          chk_version=self.chk_version)
                skysub_exptime = self.spectrograph.get_meta_value([spec2DObj.head0], 'exptime')
            except:
                raise PypeItError("Could not load skysub image from spec2d file:\n"
                                  + self.cubepar['skysub_frame'])
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
                log.info("Loading skysub frame:\n" + opts_skysub)
                try:
                    spec2DObj_sky = spec2dobj.Spec2DObj.from_file(opts_skysub, self.detname,
                                                                  chk_version=self.chk_version)
                    skysub_exptime = self.spectrograph.get_meta_value([spec2DObj_sky.head0], 'exptime')
                except:
                    raise PypeItError(
                        "Could not load skysub image from spec2d file:\n" + opts_skysub
                    )
                skyImg = spec2DObj_sky.sciimg * exptime / skysub_exptime  # Sky counts
                skyScl = spec2DObj_sky.scaleimg
                this_skysub = opts_skysub  # User specified spec2d for sky subtraction
        if this_skysub == "none":
            log.info("Sky subtraction will not be performed.")
        else:
            log.info("Using the following frame for sky subtraction:\n" + this_skysub)
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
        if not Path(flatfile).is_file():
            log.warning(f"Grating correction requested, but {flatfile} does not exist!")
            return

        if flatfile not in self.flat_splines.keys():
            log.info("Calculating relative sensitivity for grating correction")
            # Load the Flat file
            flatimages = flatfield.FlatImages.from_file(flatfile, chk_version=self.chk_version)
            total_illum = flatimages.fit2illumflat(slits, finecorr=False, frametype='illum', spat_flexure=spat_flexure) * \
                          flatimages.fit2illumflat(slits, finecorr=True, frametype='illum', spat_flexure=spat_flexure)
            flatframe = flatimages.pixelflat_raw / total_illum
            if flatimages.pixelflat_spec_illum is None:
                # Calculate the relative scale
                scale_model = flatfield.illum_profile_spectral(flatframe, waveimg, slits,
                                                               slit_illum_ref_idx=self.flatpar['slit_illum_ref_idx'],
                                                               model=None, trim=self.flatpar['slit_trim'],
                                                               flexure=spat_flexure,
                                                               smooth_npix=self.flatpar['slit_illum_smooth_npix'])
            else:
                log.info("Using relative spectral illumination from FlatImages")
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
        raise NotImplementedError('Parent class run() function must be overwritten in child class.')


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
    def __init__(self, spec2dfiles, par,
                 redux_path=None, skysub_frame=None, sensfile=None, scale_corr=None, grating_corr=None,
                 ra_offsets=None, dec_offsets=None, spectrograph=None, det=1,
                 overwrite=False, show=False, debug=False):
        super().__init__(spec2dfiles, par,
                         redux_path=redux_path, skysub_frame=skysub_frame, sensfile=sensfile,
                         scale_corr=scale_corr, grating_corr=grating_corr,
                         ra_offsets=ra_offsets, dec_offsets=dec_offsets, spectrograph=spectrograph, det=det,
                         overwrite=overwrite, show=show, debug=debug)
        self.mnmx_wv = None  # Will be used to store the minimum and maximum wavelengths of every slit and frame.
        self._spatscale = np.zeros((self.numfiles, 2))  # index 0, 1 = pixel scale, slicer scale
        self._specscale = np.zeros(self.numfiles)

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
                alignfile = Path(spec2DObj.calibs['DIR']).absolute() / spec2DObj.calibs[key]
                if alignfile.is_file() and self.cubepar['astrometric']:
                    log.info("Loading alignments")
                    alignments = alignframe.Alignments.from_file(
                        str(alignfile), chk_version=self.chk_version
                    )
            else:
                log.warning(f'Processed alignment frame not recorded or not found!')
                log.info("Using slit edges for astrometric transform")
        else:
            log.info("Using slit edges for astrometric transform")
        # If nothing better was provided, use the slit edges
        if alignments is None:
            left, right, _ = slits.select_edges(flexure=spat_flexure)
            locations = [0.0, 1.0]
            traces = np.append(left[:, None, :], right[:, None, :], axis=1)
        else:
            locations = self.par['calibrations']['alignment']['locations']
            traces = alignments.traces
        log.info("Generating alignment splines")
        return alignframe.AlignmentSplines(traces, locations, spec2DObj.tilts)

    def load_data(self):
        """
        This is the main function that loads in the data, and performs several frame-specific corrections.

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
            log.info(f"Loading PypeIt spec2d frame ({ff+1}/{len(self.spec2d)}):\n" + fil)
            spec2DObj = spec2dobj.Spec2DObj.from_file(fil, self.detname,
                                                      chk_version=self.chk_version)
            detector = spec2DObj.detector
            spat_flexure = None  # spec2DObj.sci_spat_flexure

            # Load the header
            hdr0 = spec2DObj.head0
            self.all_header.append(hdr0)
            self.ifu_ra = np.append(self.ifu_ra, self.spectrograph.compound_meta([hdr0], 'ra'))
            self.ifu_dec = np.append(self.ifu_dec, self.spectrograph.compound_meta([hdr0], 'dec'))

            # Get the exposure time
            exptime = self.spectrograph.compound_meta([hdr0], 'exptime')

            # Initialise the slit edges
            log.info("Constructing slit image")
            slits = spec2DObj.slits
            slitid_img = slits.slit_img(pad=0, flexure=spat_flexure)

            # The order of operations below proceeds as follows:
            #  (1) Get science image
            #  (2) Subtract sky (note, if a joint fit has been performed, the
            #      relative scale correction is applied in the reduction!)
            #  (3) Apply relative scale correction to both science and ivar

            # Set the default behaviour if a global skysub frame has been specified
            this_skysub, skyImg, skyScl = self.get_current_skysub(
                spec2DObj, exptime, opts_skysub=self.skysub_frame[ff]
            )

            # Load the relative scale image, if something other than the default has been provided
            this_scalecorr, relScaleImg = self.get_current_scalecorr(
                spec2DObj, scalecorr=self.scale_corr[ff]
            )

            # Prepare the relative scaling factors
            # - This factor ensures the sky has the same relative scaling as the
            # science frame
            relSclSky = skyScl / spec2DObj.scaleimg
            # - This factor is applied to the sky subtracted science frame
            relScale = spec2DObj.scaleimg / relScaleImg

            # Extract the relevant information from the spec2d file
            # - Subtract sky and apply relative illumination
            sciImg = (spec2DObj.sciimg - skyImg * relSclSky) * relScale
            ivar = spec2DObj.ivarraw / relScale ** 2
            waveimg = spec2DObj.waveimg
            bpmmask = spec2DObj.bpmmask

            # Mask the edges of the spectrum where the sky model is bad
            sky_is_good = datacube.make_good_skymask(slitid_img, spec2DObj.tilts)

            # TODO: Really need to write some detailed information in the docs
            # about all of the various corrections that can optionally be
            # applied

            # TODO: Include a flexure correction from the sky frame? Note, you
            # cannot use the waveimg from a sky frame, since the heliocentric
            # correction may have been applied to the sky frame. Need to
            # recalculate waveimg using the slitshifts from a skyimage, and then
            # apply the vel_corr from the science image.

            wnonzero = (waveimg != 0.0)
            if not np.any(wnonzero):
                raise PypeItError(
                    "The wavelength image contains only zeros - You need to check the data "
                    "reduction."
                )
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
            dwv = (
                np.median(dwaveimg[dwaveimg != 0.0]) if self.cubepar['wave_delta'] is None
                else self.cubepar['wave_delta']
            )

            log.info(
                f"Using wavelength solution: wave0={wave0:.3f}, dispersion={dwv:.3f} "
                "Angstrom/pixel"
            )

            # Obtain the minimum and maximum wavelength of all slits
            if self.mnmx_wv is None:
                self.mnmx_wv = np.zeros((len(self.spec2d), slits.nslits, 2))
            for slit_idx, slit_spat in enumerate(slits.spat_id):
                onslit_init = (slitid_img == slit_spat)
                self.mnmx_wv[ff, slit_idx, 0] = np.min(waveimg[onslit_init])
                self.mnmx_wv[ff, slit_idx, 1] = np.max(waveimg[onslit_init])

            # Find the largest spatial scale of all images being combined
            # TODO: probably need to put this in the DetectorContainer
            # pxscl is in degrees/pixel
            pxscl = detector.platescale * parse.parse_binning(detector.binning)[1] / 3600.0
            slscl = self.spectrograph.get_meta_value([spec2DObj.head0], 'slitwid')
            self._spatscale[ff, 0] = pxscl
            self._spatscale[ff, 1] = slscl
            self._specscale[ff] = dwv

            # If the spatial scale has been set by the user, check that it
            # doesn't exceed the pixel or slicer scales
            if self._dspat is not None:
                if pxscl > self._dspat:
                    log.warning(
                        f"Spatial scale requested ({3600.0 * self._dspat:f} arcsec) is less than "
                        f"the pixel scale ({3600.0 * pxscl:f} arcsec)"
                    )
                if slscl > self._dspat:
                    log.warning(
                        f"Spatial scale requested ({3600.0 * self._dspat:f} arcsec) is less than "
                        f"the slicer scale ({3600.0 * slscl:f} arcsec)"
                    )

            # Construct a good pixel mask
            # TODO: This should use the mask function to figure out which elements are masked.
            onslit_gpm = (slitid_img > 0) & (bpmmask.mask == 0) & sky_is_good

            # Generate the alignment splines, and then retrieve images of the RA
            # and Dec of every pixel, and the number of spatial pixels in each
            # slit
            alignSplines = self.get_alignments(spec2DObj, slits, spat_flexure=spat_flexure)

            # Grab the WCS of this frame, and generate the RA and Dec images
            # NOTE: These RA and Dec images are only used to setup the WCS of
            # the datacube. The actual RA and Dec of each pixel in the datacube
            # is calculated in the datacube.subpixellate() method.
            crval_wv = self.cubepar['wave_min'] if self.cubepar['wave_min'] is not None else wave0
            cd_wv = self.cubepar['wave_delta'] if self.cubepar['wave_delta'] is not None else dwv
            self.all_wcs.append(self.spectrograph.get_wcs(
                spec2DObj.head0, slits, detector.platescale, crval_wv, cd_wv
            ))
            ra_img, dec_img, delta_pix = slits.get_radec_image(
                self.all_wcs[ff], alignSplines, spec2DObj.tilts, flexure=spat_flexure
            )

            # Extract wavelength and delta wavelength arrays from the images
            wave_ext = waveimg[onslit_gpm]
            dwav_ext = dwaveimg[onslit_gpm]

            # For now, work in sorted wavelengths
            wvsrt = np.argsort(wave_ext, kind='stable')
            wave_sort = wave_ext[wvsrt]
            dwav_sort = dwav_ext[wvsrt]
            # Here's an array to get back to the original ordering
            resrt = np.argsort(wvsrt, kind='stable')

            # Compute the DAR correction
            cosdec = np.cos(self.ifu_dec[ff] * np.pi / 180.0)
            airmass = self.spectrograph.get_meta_value([spec2DObj.head0], 'airmass')  # unitless
            parangle = self.spectrograph.get_meta_value([spec2DObj.head0], 'parangle')
            pressure = self.spectrograph.get_meta_value([spec2DObj.head0], 'pressure')  # units are pascals
            temperature = self.spectrograph.get_meta_value([spec2DObj.head0], 'temperature')  # units are degrees C
            humidity = self.spectrograph.get_meta_value([spec2DObj.head0], 'humidity')  # Expressed as a percentage (not a fraction!)
            darcorr = DARcorrection(airmass, parangle, pressure, temperature, humidity, cosdec)

            # TODO: Need to make a note somewhere that the extinction correction
            # cannot currently be done in the datacube because the sensitivity
            # function algorithms correct the standard star for extinction when
            # generating the sensitivity function. Including the extinction
            # correction in the datacube would result in a double correction of
            # the standard star for extinction.  This could be wrong when
            # combining multiple standard star exposures if the airmass of the
            # standard star exposures is significantly different. For now, to
            # stay consistent with the current pipeline, the extinction
            # correction is done in the sensitivity function algorithms, with
            # the caveat that the standard star exposures are assumed to have
            # similar airmasses.  For now, comment out the extinction
            # correction, and reapply this later when the sensitivity function
            # algorithms are unified.
            extcorr_sort = 1.0
            if False:
                # Compute the extinction correction
                log.info("Applying extinction correction")
                atmext = self.spectrograph.get_atmospheric_extinction(self.senspar['UVIS']['extinct_file'])
                extcorr_sort = atmext.correction_factor(wave_sort, airmass=airmass)

            # Correct for sensitivity as a function of grating angle
            # (this assumes the spectrum of the flatfield lamp has the same shape for all setups)
            gratcorr_sort = 1.0
            if self.grating_corr[ff] is not None:
                # Setup the grating correction
                flatfile = self.grating_corr[ff]
                self.add_grating_corr(flatfile, waveimg, slits, spat_flexure=spat_flexure)
                # Calculate the grating correction
                gratcorr_sort = datacube.correct_grating_shift(
                    wave_sort, self.flat_splines[flatfile + "_wave"], self.flat_splines[flatfile],
                    self.blaze_wave, self.blaze_spline
                )

            # Sensitivity function - note that the sensitivity function factors
            # in the exposure time and the wavelength sampling, so if the flux
            # calibration will not be applied, the sens_factor needs to be
            # scaled by the exposure time and the wavelength sampling
            sens_sort = 1.0/(exptime * dwav_sort)  # If no sensitivity function is provided
            if self.fluxcal:
                log.info("Calculating the sensitivity function")
                # Load the sensitivity function
                sens = sensfunc.SensFunc.from_file(
                    self.sensfile[ff], chk_version=self.par['rdx']['chk_version']
                )
                # Interpolate the sensitivity function onto the wavelength grid of the data
                # TODO: Change the ['UVIS']['extinct_file'] here when the
                # sensitivity function calculation is unified.
                atmext = self.spectrograph.get_atmospheric_extinction(self.senspar['UVIS']['extinct_file'])
                sens_sort = flux_calib.get_sensfunc_factor(
                    wave_sort, sens.wave[:, 0], sens.zeropoint[:, 0], exptime,
                    delta_wave=dwav_sort, atmext=atmext,
                    airmass=airmass, extrap_sens=self.par['fluxcalib']['extrap_sens']
                )

            # Convert the flux units to counts/s, and correct for the relative
            # sensitivity of different setups
            sens_sort *= extcorr_sort/gratcorr_sort
            # Correct for extinction
            sciImg[onslit_gpm] *= sens_sort[resrt]
            ivar[onslit_gpm] /= sens_sort[resrt] ** 2

            # Convert units to Counts/s/Ang/arcsec2
            # Slicer sampling * spatial pixel sampling
            unitscale = (
                self.all_wcs[ff].wcs.cunit[0].to(units.arcsec)
                * self.all_wcs[ff].wcs.cunit[1].to(units.arcsec)
            )
            sl_deg = np.sqrt(
                self.all_wcs[ff].wcs.cd[0, 0] ** 2 + self.all_wcs[ff].wcs.cd[1, 0] ** 2
            )
            px_deg = np.sqrt(
                self.all_wcs[ff].wcs.cd[1, 1] ** 2 + self.all_wcs[ff].wcs.cd[0, 1] ** 2
            )
            scl_units = unitscale * sl_deg * px_deg
            sciImg[onslit_gpm] /= scl_units
            ivar[onslit_gpm] *= scl_units ** 2

            # Calculate the weights relative to the zeroth cube
            # TODO: Weights are always uniform now, right?  Doesn't this mess up
            # the selection of the reference image when using a
            # cross-correlation to get the registration offsets?
            self.weights[ff] = 1.0  # exptime  #np.median(flux_sav[resrt]*np.sqrt(ivar_sav[resrt]))**2
            wghts = self.weights[ff] * np.ones(sciImg.shape)

            # Get the slit image and then unset pixels in the slit image that are bad
            slitid_img_gpm = slitid_img * onslit_gpm.astype(int)

            # Store the information if we are combining multiple frames
            self.all_sci.append(sciImg.copy())
            self.all_ivar.append(ivar.copy())
            self.all_wave.append(waveimg.copy())
            self.all_ra.append(ra_img.copy())
            self.all_dec.append(dec_img.copy())
            self.all_deltapix.append(delta_pix.copy())
            self.all_slitid.append(slitid_img_gpm.copy())
            self.all_wghts.append(wghts.copy())
            self.all_tilts.append(spec2DObj.tilts)
            self.all_slits.append(slits)
            self.all_align.append(alignSplines)
            self.all_dar.append(darcorr)

    def run_align(self, fwhm=1.5, show_qa=False):
        """
        This routine aligns multiple cubes by using manual input offsets or 
        by cross-correlating white light images.
        
        Parameters
        ----------
        fwhm : float, optional
            The full-width half-maximum of the PSF in arcseconds. This is used
            only if the offsets are computed from point source positions. 
        show_qa : bool, optional
            If True, show QA plots for point source alignment.

        Returns
        -------
        ra_offsets : :class:`numpy.ndarray`
            A new set of RA values that have been aligned
        dec_offsets : :class:`numpy.ndarray`
            A new set of Dec values that have been aligned            
        """
        # Grab cos(dec) for convenience
        cosdec = np.cos(np.mean(self.ifu_dec[0]) * np.pi / 180.0)
        dspat = (self._dspat*units.deg).to(units.arcsec).value
        # Initialize the RA and Dec offset arrays
        ra_offsets, dec_offsets = [0.0]*self.numfiles, [0.0]*self.numfiles
        # Register spatial offsets between all frames
        if self.alignment_method == "user":
            # The user has specified offsets - update these values accounting for the difference in header RA/DEC
            return datacube.align_user_offsets(
                self.ifu_ra, self.ifu_dec, self.ra_offsets, self.dec_offsets
            )

        # Find the wavelength range where all frames overlap
        min_wl, max_wl = datacube.get_whitelight_range(np.max(self.mnmx_wv[:, :, 0]),  # The max blue wavelength
                                                       np.min(self.mnmx_wv[:, :, 1]),  # The min red wavelength
                                                       self.cubepar['whitelight_range'])  # The user-specified values (if any)
        # Get the good white light pixels
        slitid_img_gpm, wavediff = datacube.get_whitelight_pixels(self.all_wave, self.all_slitid, min_wl, max_wl)
        # Iterate over white light image generation and spatial shifting
        # TODO : Should we add this as a parameter the user can update?
        numiter = 5
        for dd in range(numiter):
            log.info(f"Iterating on spatial translation - ITERATION #{dd+1}/{numiter}")
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
                raise PypeItError("Spectral range for WCS is incorrect for white light image")

            wl_imgs, sig_imgs, bpm_imgs = datacube.generate_image_subpixel(
                image_wcs, voxedge, self.all_sci, self.all_ivar, self.all_wave,
                slitid_img_gpm, self.all_wghts, self.all_wcs,
                self.all_tilts, self.all_slits, self.all_align, self.all_dar,
                ra_offsets, dec_offsets, spec_subpixel=self.spec_subpixel,
                spat_subpixel=self.spat_subpixel, slice_subpixel=self.slice_subpixel)
            if reference_image is None:
                # ref_idx will be the index of the cube with the highest S/N
                ref_idx = np.argmax(self.weights)
                reference_image = wl_imgs[:, :, ref_idx].copy()
                log.info("Calculating spatial translation of each cube relative to cube #{0:d})".format(ref_idx+1))
            else:
                log.info("Calculating the spatial translation of each cube relative to user-defined 'reference_image'")

            # Calculate the image offsets relative to the reference image
            if self.alignment_method == 'phase':
                for ff in range(self.numfiles):
                    # Calculate the shift
                        ra_shift, dec_shift = calculate_image_phase(
                            reference_image.copy(), wl_imgs[:, :, ff], maskval=0.0)
                        # Convert pixel shift to degrees shift
                        ra_shift *= self._dspat/cosdec
                        dec_shift *= self._dspat
                        log.info(
                            f"Spatial shift of cube #{ff+1}:\n"
                            f"RA, DEC (arcsec) = {ra_shift*3600.0:+0.3f} E, "
                            f"{dec_shift*3600.0:+0.3f} N"
                        )
                        # Store the shift in the RA and DEC offsets in degrees
                        ra_offsets[ff] += ra_shift
                        dec_offsets[ff] += dec_shift
            elif self.alignment_method == 'fit':
                ra_pix_star = np.zeros(self.numfiles)
                dec_pix_star = np.zeros(self.numfiles)
                for ff in range(self.numfiles):
                    popt, pcov, model, init_obj_position, flux_opt, sigma_opt = \
                        datacube.fitGaussian2D(
                            wl_imgs[:, :, ff], ivar=utils.inverse(np.square(sig_imgs[:,:, ff])),
                            gpm=np.logical_not(bpm_imgs[:, :, ff]), fwhm=fwhm/dspat, norm=False
                        )
                    gaussian_position = popt[1], popt[2]
                    if show_qa and dd == numiter-1:
                        datacube.whitelight_objfind_qa(
                            wl_imgs[:, :, ff], utils.inverse(np.square(sig_imgs[:, :, ff])),
                            np.logical_not(bpm_imgs[:, :, ff]), model, gaussian_position,
                            init_obj_position, channel_prefix = f'Img_{ff}'
                        )
                    ra_pix_star[ff], dec_pix_star[ff] = gaussian_position

                ra_shifts = (ra_pix_star - ra_pix_star[ref_idx]) * self._dspat / cosdec
                dec_shifts = (dec_pix_star - dec_pix_star[ref_idx]) * self._dspat
                ra_offsets =[ra_offsets[ff] + ra_shifts[ff] for ff in range(self.numfiles)]
                dec_offsets =[dec_offsets[ff] + dec_shifts[ff] for ff in range(self.numfiles)]
            else:
                raise PypeItError(f"self.alignment_method method '{self.alignment_method}' is not supported.")

            for ff in range(self.numfiles):
                log.info(
                    f"Spatial shift of cube #{ff + 1}:\n"
                    f"RA, DEC (arcsec) = {ra_offsets[ff]*3600.0:+0.3f} E, "
                    f"{dec_offsets[ff]*3600.0:+0.3f} N"
                )

        return ra_offsets, dec_offsets

    def save_native(self):
        """
        If the user requests to store datacubes at the native sampling of the spectrograph,
        loop through all files and write them to disk.
        """
        # Loop through all input spec2d files and save datacubes at a spectrographs native sampling.
        for ff, fil in enumerate(self.spec2d):
            # Get the output filename
            if self.numfiles == 1 and self.cubepar['output_filename'] != "":
                outfile = datacube.get_output_filename(self.scidir, "", self.cubepar['output_filename'],
                                                       False, native=True, idx=-1)
            else:
                outfile = datacube.get_output_filename(self.scidir, fil, self.cubepar['output_filename'],
                                                       False, native=True, idx=ff+1)
            # Get the coordinate bounds
            wave0 = self.all_wave[ff][self.all_wave[ff] != 0.0].min()
            slitlength = int(np.round(np.median(self.all_slits[ff].get_slitlengths(median=True))))
            numwav = int((np.max(self.all_wave[ff]) - wave0) / self._specscale[ff])
            bins = self.spectrograph.get_datacube_bins(slitlength, self.all_deltapix[ff], numwav)
            # Set the wavelength range of the white light image.
            wl_wvrng = None
            if self.cubepar['save_whitelight']:
                wl_wvrng = datacube.get_whitelight_range(np.max(self.mnmx_wv[ff, :, 0]),
                                                         np.min(self.mnmx_wv[ff, :, 1]),
                                                         self.cubepar['whitelight_range'])
            # Make the datacube
            if self.method in ['subpixel', 'ngp']:
                # Generate the datacube
                wghts = np.ones(self.all_sci[ff].shape)
                flxcube, sigcube, bpmcube, normcube, wave = \
                    datacube.generate_cube_subpixel(self.all_wcs[ff], bins, self.all_sci[ff], self.all_ivar[ff],
                                                    self.all_wave[ff], self.all_slitid[ff], wghts,
                                                    self.all_wcs[ff], self.all_tilts[ff], self.all_slits[ff],
                                                    self.all_align[ff], self.all_dar[ff],
                                                    self.ra_offsets[ff], self.dec_offsets[ff],
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
                log.info(f"Saving datacube at the native sampling of {self.specname.replace("_", " ")}: {outfile}")
                final_cube = DataCube(
                    flxcube, sigcube, bpmcube.astype(np.uint8),
                    wave, self.specname, self.blaze_wave, self.blaze_spec,
                    sensfunc=None, fluxed=self.fluxcal
                )
                final_cube.to_file(
                    str(self.scidir / outfile), primary_hdr=self.all_header[ff],
                    hdr=hdr, overwrite=self.overwrite
                )

                ivarcube = utils.inverse(np.square(sigcube))
                if self.cubepar['save_whitelight']:
                    datacube.make_whitelight(
                        self.all_wcs[ff], flxcube, ivarcube, np.logical_not(bpmcube),
                        wave, self.scidir, outfile, whitelight_range=wl_wvrng, overwrite=self.overwrite)

    def compute_weights(self, show_qa=False):
        """
        Compute the relative weights to apply to pixels that are collected into the voxels of the output DataCubes

        Parameters
        ----------
        show_qa : bool, optional
            If True, show QA plots for point source alignment.

        Returns
        -------
        :class:`numpy.ndarray`
            The individual pixel weights for each detector pixel, and every frame.
        """
        # If there is only one file, then all pixels have the same weight
        if self.numfiles == 1:
            return np.ones_like(self.all_sci)
        elif self.cubepar['weight_method'] == 'uniform':
            # No need to compute weights if they are uniform
            return [np.ones_like(sci) for sci in self.all_sci]
        else:
            # Calculate the relative spectral weights of all pixels
            
            if self.cubepar['weights_init_obj_pos'] is not None and len(self.cubepar['weights_init_obj_pos']) > 0:
                manual_dict= ManualCubeExtractionObj.parse(self.cubepar['weights_init_obj_pos']).to_dict()
                init_obj_position = (manual_dict['spatx'][0], manual_dict['spaty'][0])
            else: 
                init_obj_position = None
            
            return datacube.compute_weights_frompix(
                self.all_ra, self.all_dec, self.all_wave, self.all_sci, self.all_ivar,
                self.all_slitid, self._dspat, self._dwv, self.mnmx_wv, self.all_wghts,
                self.all_wcs, self.all_tilts, self.all_slits, self.all_align, self.all_dar,
                self.ra_offsets, self.dec_offsets,
                ra_min=self.cubepar['ra_min'], ra_max=self.cubepar['ra_max'],
                dec_min=self.cubepar['dec_min'], dec_max=self.cubepar['dec_max'],
                wave_min=self.cubepar['wave_min'], wave_max=self.cubepar['wave_max'],
                weight_method=self.cubepar['weight_method'],
                sn_smooth_npix = self.cubepar['sn_smooth_npix'],
                whitelight_range=self.cubepar['whitelight_range'],
                reference_image=self.cubepar['reference_image'],
                correct_dar=self.correct_dar,
                specname=self.specname, init_obj_position=init_obj_position,
                show_qa=show_qa)

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
        # Loop through all frames and load the data
        self.load_data()

        # If we are combining frames, check that alignment has been requested. 
        # If not, then print out a warning. 
        if self.combine and not self.align:
            log.warning(
                "Combining frames without aligning them.  Make sure that you know what you are "
                "doing!  Even if your frames are taken at the same position, alignment is still "
                "recommended because of differential atmospheric refraction."
            )

        # If the user is aligning or combining, the spatial scale of the output cubes needs to be consistent.
        # Set the spatial and spectral scales of the output datacube
        self._dspat, self._dwv = datacube.set_voxel_sampling(self._spatscale, self._specscale,
                                                             dspat=self._dspat, dwv=self._dwv)

        # Align the frames
        if self.align:
            self.ra_offsets, self.dec_offsets = self.run_align(show_qa=self.debug)
            log.info('Alignment offsets are:')
            for i, (rao, dco) in enumerate(zip(self.ra_offsets, self.dec_offsets)):
                log.info(
                    f"Cube {i + 1}: RA, DEC (arcsec) = {rao*3600.0:+0.3f} E, {dco*3600.0:+0.3f} N"
                )

        # If individual frames are to be output with the native resolution of the instrument, write those out now.
        if self.native:
            self.save_native()

        # TODO There should be an if self.combine here, as we only need these weights now if we are going to
        # combine the cubes.  Furthermore, since the images are aligned, we should be using the full cube to 
        # compute the whitelight image since we do that anyway below. So basically the weight computation 
        # should be moved just before the final combined datacube generation below. Specifically, we should be: 
        # 1. Performing an intiial sigma clipping of the cubes. 
        # 2. Computing an initial preliminary stacked cube. 
        # 3. Generate a whitelight image from the preliminary stacked cube.
        # 4. Perform object finding on this stacked cube. 
        # 5. Compute the weights at the location of the object by extracting spectra from the individual cubes, 
        #    probably there should be an option to use the optimal extraction method  (extract_point_source) 
        #    or one can use the single pixel computation (extended sources) in compute_weights. 
        # 6. Re-combined the cubes using these weights, again performing the final round of sigma clipping
        # 7. Write out the individual cubes with their sigma clipped pixels masked (?)
        # 8. Write out the final combined cube.

        # Compute the relative weights on the spectra
        self.all_wghts = self.compute_weights(show_qa=self.debug)

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
            numwav = vox_edges[0].size - 1
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

            for ff in range(self.numfiles):
                outfile = datacube.get_output_filename(
                    self.scidir, "", self.cubepar['output_filename'], False, idx=ff+1
                )
                # Generate the datacube       
                
                # TODO Put in a self.native flag to allow for the datacube to be generated in the native resolution
                # of the data as it is currently being done in the load method?    
                flxcube, sigcube, bpmcube, normcube, wave = \
                    datacube.generate_cube_subpixel(cube_wcs, vox_edges,
                                                    self.all_sci[ff], self.all_ivar[ff], self.all_wave[ff],
                                                    self.all_slitid[ff], self.all_wghts[ff], self.all_wcs[ff],
                                                    self.all_tilts[ff], self.all_slits[ff], self.all_align[ff], 
                                                    self.all_dar[ff],
                                                    self.ra_offsets[ff], self.dec_offsets[ff],
                                                    spec_subpixel=self.spec_subpixel,
                                                    spat_subpixel=self.spat_subpixel,
                                                    slice_subpixel=self.slice_subpixel,
                                                    skip_subpix_weights=self.skip_subpix_weights,
                                                    correct_dar=self.correct_dar)
                if self.combine: #& self.align:                   
                    # If we are combining cubes, then we need to save these for the final combination
                    # with sigma clipping below, otherwise no need to store these and use more memory
                    if ff == 0: 
                        stack_shape = (self.numfiles,) + flxcube.shape
                        flxcube_stack = np.zeros(stack_shape)
                        varcube_stack = np.zeros(stack_shape)
                        bpmcube_stack = np.zeros(stack_shape)
                        normcube_stack = np.zeros(stack_shape)
                        # TODO Add proper weights
                        weightcube_stack = np.ones(stack_shape)

                    flxcube_stack[ff, :] = flxcube
                    varcube_stack[ff, :] = np.square(sigcube)
                    bpmcube_stack[ff, :] = bpmcube
                    normcube_stack[ff, :] = normcube
                    
                # Prepare the header
                hdr = cube_wcs.to_header()
                if self.fluxcal:
                    hdr['FLUXUNIT'] = (flux_calib.PYPEIT_FLUX_SCALE, "Flux units -- erg/s/cm^2/Angstrom/arcsec^2")
                else:
                    hdr['FLUXUNIT'] = (1, "Flux units -- counts/s/Angstrom/arcsec^2")
                # Write out the datacube
                log.info("Saving datacube as: {0:s}".format(str(outfile)))
                final_cube = DataCube(
                    flxcube, sigcube, bpmcube.astype(np.uint8), wave, self.specname,
                    self.blaze_wave, self.blaze_spec, sensfunc=sensfunc, fluxed=self.fluxcal
                )
                final_cube.to_file(
                    str(self.scidir / outfile), primary_hdr=self.all_header[ff],
                    hdr=hdr, overwrite=self.overwrite
                )
                ivarcube = final_cube.ivar
                if self.cubepar['save_whitelight']:
                    datacube.make_whitelight(
                        cube_wcs, flxcube, ivarcube, np.logical_not(bpmcube), wave,
                        self.scidir, outfile, whitelight_range=wl_wvrng,
                        overwrite=self.overwrite
                    )

            if self.combine:
                sigrej = 3.0
                maxiters = 10                
                sci_list_out, var_list_out, combined_gpm, nused = combine.weighted_combine(
                    weightcube_stack, [flxcube_stack], [varcube_stack],
                    np.logical_not(bpmcube_stack), sigma_clip=True, sigma_clip_stack=flxcube_stack,
                    sigrej=sigrej, maxiters=maxiters
                )
                combined_cube = sci_list_out[0]
                combined_sigma = np.sqrt(var_list_out[0])
                combined_ivar = utils.inverse(var_list_out[0])
                combined_bpm = np.logical_not(combined_gpm)
                combined_outfile = datacube.get_output_filename(
                    self.scidir, "", self.cubepar['output_filename'], True, idx=-1
                )
                log.info(f"Saving combined datacube as: {str(combined_outfile)}")
                final_combined_cube = DataCube(
                    combined_cube, combined_sigma, combined_bpm.astype(np.uint8), wave,
                    self.specname, self.blaze_wave, self.blaze_spec, sensfunc=sensfunc,
                    fluxed=self.fluxcal
                )
                final_combined_cube.to_file(
                    str(self.scidir / combined_outfile), primary_hdr=self.all_header[ff],
                    hdr=hdr, overwrite=self.overwrite
                )
                # Make combined white light image if whitelight is requested
                if self.cubepar['save_whitelight']:                
                    datacube.make_whitelight(
                        cube_wcs, combined_cube, combined_ivar, combined_gpm, wave,
                        self.scidir, combined_outfile, whitelight_range=wl_wvrng,
                        overwrite=self.overwrite
                    )
