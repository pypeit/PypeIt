"""
Module containing routines used by 3D datacubes.

.. include:: ../include/links.rst
"""

import os
import copy
import inspect

from astropy import wcs, units
from astropy.coordinates import AltAz, SkyCoord
from astropy.io import fits
import scipy.optimize as opt
from scipy.interpolate import interp1d
import numpy as np

from pypeit import msgs
from pypeit import alignframe, datamodel, flatfield, io, specobj, spec2dobj, utils
from pypeit.core.flexure import calculate_image_phase
from pypeit.core import coadd, datacube, extract, findobj_skymask, flux_calib, parse, skysub
from pypeit.core.procimg import grow_mask
from pypeit.spectrographs.util import load_spectrograph

# Use a fast histogram for speed!
try:
    from fast_histogram import histogramdd
except ImportError:
    histogramdd = None

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

    """
    version = '1.1.0'

    datamodel = {'flux': dict(otype=np.ndarray, atype=np.floating,
                              descr='Flux datacube in units of counts/s/Ang/arcsec^2 or '
                                    '10^-17 erg/s/cm^2/Ang/arcsec^2'),
                 'sig': dict(otype=np.ndarray, atype=np.floating,
                             descr='Error datacube (matches units of flux)'),
                 'bpm': dict(otype=np.ndarray, atype=np.uint8,
                             descr='Bad pixel mask of the datacube (0=good, 1=bad)'),
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
                 'spect_meta'
                ]

    def __init__(self, flux, sig, bpm, PYP_SPEC, blaze_wave, blaze_spec, sensfunc=None,
                 fluxed=None):

        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        _d = dict([(k, values[k]) for k in args[1:]])
        # Setup the DataContainer
        datamodel.DataContainer.__init__(self, d=_d)

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
    def from_file(cls, ifile):
        """
        Over-load :func:`~pypeit.datamodel.DataContainer.from_file`
        to deal with the header

        Args:
            ifile (str):  Filename holding the object
        """
        with io.fits_open(ifile) as hdu:
            # Read using the base class
            self = super().from_hdu(hdu)
            # Internals
            self.filename = ifile
            self.head0 = hdu[1].header  # Actually use the first extension here, since it contains the WCS
            # Meta
            self.spectrograph = load_spectrograph(self.PYP_SPEC)
            self.spect_meta = self.spectrograph.parse_spec_header(hdu[0].header)
        return self

    @property
    def ivar(self):
        return utils.inverse(self.sig**2)

    @property
    def wcs(self):
        return wcs.WCS(self.head0)


class CoAdd3D:
    """
    Main routine to convert processed PypeIt spec2d frames into
    DataCube (spec3d) files. This routine is only used for IFU
    data reduction.

    Algorithm steps are as follows:
        - TODO :: Fill this in.

    """
    # Superclass factory method generates the subclass instance
    @classmethod
    def get_instance(cls, spec2dfiles, opts, spectrograph=None, par=None, det=1, overwrite=False,
                     show=False, debug=False):
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
                        spec2dfiles, spectrograph=spectrograph, par=par, det=det, overwrite=overwrite,
                        show=show, debug=debug)

    def __init__(self, files, opts, spectrograph=None, par=None, det=1, overwrite=False,
                 show=False, debug=False):
        """

        Args:
            files (:obj:`list`):
                List of all spec2D files
            opts (:obj:`dict`):
                Options associated with each spec2d file
            spectrograph (:obj:`str`, :class:`~pypeit.spectrographs.spectrograph.Spectrograph`, optional):
                The name or instance of the spectrograph used to obtain the data.
                If None, this is pulled from the file header.
            par (:class:`~pypeit.par.pypeitpar.PypeItPar`, optional):
                An instance of the parameter set.  If None, assumes that detector 1
                is the one reduced and uses the default reduction parameters for the
                spectrograph (see
                :func:`~pypeit.spectrographs.spectrograph.Spectrograph.default_pypeit_par`
                for the relevant spectrograph class).
            det (int):
                Detector index
            overwrite (:obj:`bool`, optional):
                Overwrite the output file, if it exists?
            show (:obj:`bool`, optional):
                Show results in ginga
            debug (:obj:`bool`, optional):
                Show QA for debugging.

        """
        self.spec2d = files
        self.numfiles = len(files)
        self.opts = opts
        self.overwrite = overwrite

        # Check on Spectrograph input
        if spectrograph is None:
            with fits.open(files[0]) as hdu:
                spectrograph = hdu[0].header['PYP_SPEC']

        if isinstance(spectrograph, str):
            self.spec = load_spectrograph(spectrograph)
            self.specname = spectrograph
        else:
            # Assume it's a Spectrograph instance
            self.spec = spectrograph
            self.specname = spectrograph.name

        # Grab the parset, if not provided
        if par is None:
            # TODO :: Use config_specific_par instead?
            par = self.spec.default_pypeit_par()
        self.par = par
        # Extract some parsets for simplicity
        self.cubepar = self.par['reduce']['cube']
        self.flatpar = self.par['calibrations']['flatfield']
        self.senspar = self.par['sensfunc']

        # Initialise arrays for storage
        self.ifu_ra, self.ifu_dec = np.array([]), np.array([])  # The RA and Dec at the centre of the IFU, as stored in the header
        self.all_ra, self.all_dec, self.all_wave = np.array([]), np.array([]), np.array([])
        self.all_sci, self.all_ivar, self.all_idx, self.all_wghts = np.array([]), np.array([]), np.array([]), np.array([])
        self.all_spatpos, self.all_specpos, self.all_spatid = np.array([], dtype=int), np.array([], dtype=int), np.array([], dtype=int)
        self.all_tilts, self.all_slits, self.all_align = [], [], []
        self.all_wcs = []
        self.weights = np.ones(self.numfiles)  # Weights to use when combining cubes




        # TODO :: need to sort out what to do with these - make them self. as well?
        assert False
        dspat = None if self.cubepar['spatial_delta'] is None else self.cubepar['spatial_delta'] / 3600.0  # binning size on the sky (/3600 to convert to degrees)
        dwv = self.cubepar['wave_delta']  # binning size in wavelength direction (in Angstroms)
        flat_splines = dict()  # A dictionary containing the splines of the flatfield





        # Extract some commonly used variables
        self.method = self.cubepar['method'].lower()
        self.combine = self.cubepar['combine']
        self.align = self.cubepar['align']
        # If there is only one frame being "combined" AND there's no reference image, then don't compute the translation.
        if self.numfiles == 1 and self.cubepar["reference_image"] is None:
            if not self.align:
                msgs.warn("Parameter 'align' should be False when there is only one frame and no reference image")
                msgs.info("Setting 'align' to False")
            self.align = False
        if self.opts['ra_offset'] is not None:
            if not self.align:
                msgs.warn("When 'ra_offset' and 'dec_offset' are set, 'align' must be True.")
                msgs.info("Setting 'align' to True")
            self.align = True
        # TODO :: The default behaviour (combine=False, align=False) produces a datacube that uses the instrument WCS
        #  It should be possible (and perhaps desirable) to do a spatial alignment (i.e. align=True), apply this to the
        #  RA,Dec values of each pixel, and then use the instrument WCS to save the output (or, just adjust the crval).
        #  At the moment, if the user wishes to spatially align the frames, a different WCS is generated.
        # Check if fast-histogram exists
        if histogramdd is None:
            msgs.warn("Generating a datacube is faster if you install fast-histogram:"+msgs.newline()+
                      "https://pypi.org/project/fast-histogram/")
            if self.method != 'ngp':
                msgs.warn("Forcing NGP algorithm, because fast-histogram is not installed")
                self.method = 'ngp'

        # Determine what method is requested
        self.spec_subpixel, self.spat_subpixel = 1, 1
        if self.method == "subpixel":
            msgs.info("Adopting the subpixel algorithm to generate the datacube.")
            spec_subpixel, spat_subpixel = self.cubepar['spec_subpixel'], self.cubepar['spat_subpixel']
        elif self.method == "ngp":
            msgs.info("Adopting the nearest grid point (NGP) algorithm to generate the datacube.")
        else:
            msgs.error(f"The following datacube method is not allowed: {self.method}")

        # Get the detector number and string representation
        det = 1 if self.par['rdx']['detnum'] is None else self.par['rdx']['detnum']
        self.detname = self.spec.get_det_name(det)

        # Check if the output file exists
        self.check_outputs()

        # Check the reference cube and image exist, if requested
        self.fluxcal = False
        self.blaze_wave, self.blaze_spec = None, None
        self.blaze_spline, self.flux_spline = None, None
        if self.cubepar['standard_cube'] is not None:
            self.make_sensfunc()

        # If a reference image has been set, check that it exists
        if self.cubepar['reference_image'] is not None:
            if not os.path.exists(self.cubepar['reference_image']):
                msgs.error("Reference image does not exist:" + msgs.newline() + self.cubepar['reference_image'])

    def check_outputs(self):
        """
        TODO :: docstring
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

    def make_sensfunc(self):
        """
        TODO :: docstring
        """
        self.fluxcal = True
        ss_file = self.cubepar['standard_cube']
        if not os.path.exists(ss_file):
            msgs.error("Standard cube does not exist:" + msgs.newline() + ss_file)
        msgs.info(f"Loading standard star cube: {ss_file:s}")
        # Load the standard star cube and retrieve its RA + DEC
        stdcube = fits.open(ss_file)
        star_ra, star_dec = stdcube[1].header['CRVAL1'], stdcube[1].header['CRVAL2']

        # Extract a spectrum of the standard star
        wave, Nlam_star, Nlam_ivar_star, gpm_star = datacube.extract_standard_spec(stdcube)

        # Extract the information about the blaze
        if self.cubepar['grating_corr']:
            blaze_wave_curr, blaze_spec_curr = stdcube['BLAZE_WAVE'].data, stdcube['BLAZE_SPEC'].data
            blaze_spline_curr = interp1d(blaze_wave_curr, blaze_spec_curr,
                                         kind='linear', bounds_error=False, fill_value="extrapolate")
            # The first standard star cube is used as the reference blaze spline
            if self.blaze_spline is None:
                self.blaze_wave, self.blaze_spec = stdcube['BLAZE_WAVE'].data, stdcube['BLAZE_SPEC'].data
                self.blaze_spline = interp1d(self.blaze_wave, self.blaze_spec,
                                             kind='linear', bounds_error=False, fill_value="extrapolate")
            # Perform a grating correction
            grat_corr = datacube.correct_grating_shift(wave.value, blaze_wave_curr, blaze_spline_curr, self.blaze_wave,
                                              self.blaze_spline)
            # Apply the grating correction to the standard star spectrum
            Nlam_star /= grat_corr
            Nlam_ivar_star *= grat_corr ** 2

        # Read in some information above the standard star
        std_dict = flux_calib.get_standard_spectrum(star_type=self.senspar['star_type'],
                                                    star_mag=self.senspar['star_mag'],
                                                    ra=star_ra, dec=star_dec)
        # Calculate the sensitivity curve
        # TODO :: This needs to be addressed... unify flux calibration into the main PypeIt routines.
        msgs.warn("Datacubes are currently flux-calibrated using the UVIS algorithm... this will be deprecated soon")
        zeropoint_data, zeropoint_data_gpm, zeropoint_fit, zeropoint_fit_gpm = \
            flux_calib.fit_zeropoint(wave.value, Nlam_star, Nlam_ivar_star, gpm_star, std_dict,
                                     mask_hydrogen_lines=self.senspar['mask_hydrogen_lines'],
                                     mask_helium_lines=self.senspar['mask_helium_lines'],
                                     hydrogen_mask_wid=self.senspar['hydrogen_mask_wid'],
                                     nresln=self.senspar['UVIS']['nresln'],
                                     resolution=self.senspar['UVIS']['resolution'],
                                     trans_thresh=self.senspar['UVIS']['trans_thresh'],
                                     polyorder=self.senspar['polyorder'],
                                     polycorrect=self.senspar['UVIS']['polycorrect'],
                                     polyfunc=self.senspar['UVIS']['polyfunc'])
        wgd = np.where(zeropoint_fit_gpm)
        sens = np.power(10.0, -0.4 * (zeropoint_fit[wgd] - flux_calib.ZP_UNIT_CONST)) / np.square(wave[wgd])
        self.flux_spline = interp1d(wave[wgd], sens, kind='linear', bounds_error=False, fill_value="extrapolate")

        # Load the default scaleimg frame for the scale correction
        self.scalecorr_default = "none"
        self.relScaleImgDef = np.array([1])
        self.set_default_scalecorr()

        # Load the default sky frame to be used for sky subtraction
        self.skysub_default = "image"
        self.skyImgDef, self.skySclDef = None, None  # This is the default behaviour (i.e. to use the "image" for the sky subtraction)
        self.set_default_skysub()

    def set_default_scalecorr(self):
        """
        TODO :: docstring
        """
        if self.cubepar['scale_corr'] is not None:
            if self.cubepar['scale_corr'] == "image":
                msgs.info("The default relative spectral illumination correction will use the science image")
                self.scalecorr_default = "image"
            else:
                msgs.info("Loading default scale image for relative spectral illumination correction:" +
                          msgs.newline() + self.cubepar['scale_corr'])
                try:
                    spec2DObj = spec2dobj.Spec2DObj.from_file(self.cubepar['scale_corr'], self.detname)
                    self.relScaleImgDef = spec2DObj.scaleimg
                    self.scalecorr_default = self.cubepar['scale_corr']
                except:
                    msgs.warn("Could not load scaleimg from spec2d file:" + msgs.newline() +
                              self.cubepar['scale_corr'] + msgs.newline() +
                              "scale correction will not be performed unless you have specified the correct" + msgs.newline() +
                              "scale_corr file in the spec2d block")
                    self.cubepar['scale_corr'] = None
                    self.scalecorr_default = "none"

    def get_current_scalecorr(self, spec2DObj, opts_scalecorr=None):
        """
        TODO :: docstring
        """
        this_scalecorr = self.scalecorr_default
        relScaleImg = self.relScaleImgDef.copy()
        if opts_scalecorr is not None:
            if opts_scalecorr.lower() == 'default':
                if self.scalecorr_default == "image":
                    relScaleImg = spec2DObj.scaleimg
                    this_scalecorr = "image"  # Use the current spec2d for the relative spectral illumination scaling
                else:
                    this_scalecorr = self.scalecorr_default  # Use the default value for the scale correction
            elif opts_scalecorr.lower() == 'image':
                relScaleImg = spec2DObj.scaleimg
                this_scalecorr = "image"  # Use the current spec2d for the relative spectral illumination scaling
            elif opts_scalecorr.lower() == 'none':
                relScaleImg = np.array([1])
                this_scalecorr = "none"  # Don't do relative spectral illumination scaling
            else:
                # Load a user specified frame for sky subtraction
                msgs.info("Loading the following frame for the relative spectral illumination correction:" +
                          msgs.newline() + opts_scalecorr)
                try:
                    spec2DObj_scl = spec2dobj.Spec2DObj.from_file(opts_scalecorr, self.detname)
                except:
                    msgs.error(
                        "Could not load skysub image from spec2d file:" + msgs.newline() + opts_scalecorr)
                relScaleImg = spec2DObj_scl.scaleimg
                this_scalecorr = opts_scalecorr
        if this_scalecorr == "none":
            msgs.info("Relative spectral illumination correction will not be performed.")
        else:
            msgs.info("Using the following frame for the relative spectral illumination correction:" +
                      msgs.newline() + this_scalecorr)
        # Return the scaling correction for this frame
        return this_scalecorr, relScaleImg

    def set_default_skysub(self):
        """
        TODO :: Add docstring
        """
        if self.cubepar['skysub_frame'] in [None, 'none', '', 'None']:
            self.skysub_default = "none"
            self.skyImgDef = np.array([0.0])  # Do not perform sky subtraction
            self.skySclDef = np.array([0.0])  # Do not perform sky subtraction
        elif self.cubepar['skysub_frame'].lower() == "image":
            msgs.info("The sky model in the spec2d science frames will be used for sky subtraction" + msgs.newline() +
                      "(unless specific skysub frames have been specified)")
            self.skysub_default = "image"
        else:
            msgs.info("Loading default image for sky subtraction:" +
                      msgs.newline() + self.cubepar['skysub_frame'])
            try:
                spec2DObj = spec2dobj.Spec2DObj.from_file(self.cubepar['skysub_frame'], self.detname)
                skysub_exptime = fits.open(self.cubepar['skysub_frame'])[0].header['EXPTIME']
                self.skysub_default = self.cubepar['skysub_frame']
                self.skyImgDef = spec2DObj.sciimg / skysub_exptime  # Sky counts/second
                # self.skyImgDef = spec2DObj.skymodel/skysub_exptime  # Sky counts/second
                self.skySclDef = spec2DObj.scaleimg
            except:
                msgs.error("Could not load skysub image from spec2d file:" + msgs.newline() + self.cubepar['skysub_frame'])

    def get_current_skysub(self, spec2DObj, exptime, opts_skysub=None):
        """
        TODO :: docstring
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
                    spec2DObj_sky = spec2dobj.Spec2DObj.from_file(opts_skysub, self.detname)
                    skysub_exptime = fits.open(opts_skysub)[0].header['EXPTIME']
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

    def compute_DAR(self, hdr0, raimg, decimg, waveimg, onslit_gpm, wave_ref=None):
        """
        TODO :: docstring
        """
        if wave_ref is None:
            wave_ref = 0.5 * (np.min(waveimg[onslit_gpm]) + np.max(waveimg[onslit_gpm]))
        # Get DAR parameters
        raval = self.spec.get_meta_value([hdr0], 'ra')
        decval = self.spec.get_meta_value([hdr0], 'dec')
        obstime = self.spec.get_meta_value([hdr0], 'obstime')
        pressure = self.spec.get_meta_value([hdr0], 'pressure')
        temperature = self.spec.get_meta_value([hdr0], 'temperature')
        rel_humidity = self.spec.get_meta_value([hdr0], 'humidity')
        coord = SkyCoord(raval, decval, unit=(units.deg, units.deg))
        location = self.spec.location  # TODO :: spec.location should probably end up in the TelescopePar (spec.telescope.location)
        if pressure == 0.0:
            msgs.warn("Pressure is set to zero - DAR correction will not be performed")
        else:
            msgs.info("DAR correction parameters:" + msgs.newline() +
                      "   Pressure = {0:f} bar".format(pressure) + msgs.newline() +
                      "   Temperature = {0:f} deg C".format(temperature) + msgs.newline() +
                      "   Humidity = {0:f}".format(rel_humidity))
            ra_corr, dec_corr = datacube.correct_dar(waveimg[onslit_gpm], coord, obstime, location,
                                                     pressure * units.bar, temperature * units.deg_C, rel_humidity,
                                                     wave_ref=wave_ref)
            raimg[onslit_gpm] += ra_corr * np.cos(np.mean(decimg[onslit_gpm]) * np.pi / 180.0)
            decimg[onslit_gpm] += dec_corr

    def coadd(self):
        """
        TODO :: Docstring
        """
        msgs.bug("This routine should be overridden by child classes.")
        msgs.error("Cannot proceed without coding the coadd routine.")
        return


class SlicerIFUCoAdd3D(CoAdd3D):
    """
    Child of CoAdd3D for SlicerIFU data reduction. For documentation, see CoAdd3d parent class above.
    spec2dfiles, opts, spectrograph=None, par=None, det=1, overwrite=False,
                     show=False, debug=False

    """
    def __init__(self, spec2dfiles, opts, spectrograph=None, par=None, det=1, overwrite=False,
                 show=False, debug=False):
        super().__init__(spec2dfiles, opts, spectrograph=spectrograph, par=par, det=det, overwrite=overwrite,
                         show=show, debug=debug)

    def get_alignments(self, spec2DObj, slits, frame_wcs, spat_flexure=None):
        """
        TODO :: docstring
        """
        # Loading the alignments frame for these data
        alignments = None
        if self.cubepar['astrometric']:
            key = alignframe.Alignments.calib_type.upper()
            if key in spec2DObj.calibs:
                alignfile = os.path.join(spec2DObj.calibs['DIR'], spec2DObj.calibs[key])
                if os.path.exists(alignfile) and self.cubepar['astrometric']:
                    msgs.info("Loading alignments")
                    alignments = alignframe.Alignments.from_file(alignfile)
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
        # Generate an RA/DEC image
        msgs.info("Generating RA/DEC image")
        alignSplines = alignframe.AlignmentSplines(traces, locations, spec2DObj.tilts)
        # Return the alignment splines
        return alignSplines

    def load(self):
        """
        TODO :: docstring
        """
        # Initialise variables
        wave_ref = None
        mnmx_wv = None  # Will be used to store the minimum and maximum wavelengths of every slit and frame.
        # Load all spec2d files and prepare the data for making a datacube
        for ff, fil in enumerate(self.spec2d):
            # Load it up
            msgs.info("Loading PypeIt spec2d frame:" + msgs.newline() + fil)
            spec2DObj = spec2dobj.Spec2DObj.from_file(fil, self.detname)
            detector = spec2DObj.detector
            spat_flexure = None  # spec2DObj.sci_spat_flexure

            # Load the header
            hdr0 = spec2DObj.head0
            self.ifu_ra = np.append(self.ifu_ra, self.spec.compound_meta([hdr0], 'ra'))
            self.ifu_dec = np.append(self.ifu_dec, self.spec.compound_meta([hdr0], 'dec'))

            # Get the exposure time
            # TODO :: Surely this should be retrieved from metadata...
            exptime = hdr0['EXPTIME']

            # Setup for PypeIt imports
            msgs.reset(verbosity=2)

            # TODO :: Consider loading all calibrations into a single variable within the main CoAdd3D parent class.

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
                                                                  opts_skysub=self.opts['skysub_frame'][ff])

            # Load the relative scale image, if something other than the default has been provided
            this_scalecorr, relScaleImg = self.get_current_scalecorr(spec2DObj,
                                                                     opts_scalecorr=self.opts['scale_corr'][ff])

            # Prepare the relative scaling factors
            relSclSky = skyScl / spec2DObj.scaleimg  # This factor ensures the sky has the same relative scaling as the science frame
            relScale = spec2DObj.scaleimg / relScaleImg  # This factor is applied to the sky subtracted science frame

            # Extract the relevant information from the spec2d file
            sciImg = (spec2DObj.sciimg - skyImg * relSclSky) * relScale  # Subtract sky and apply relative illumination
            ivar = spec2DObj.ivarraw / relScale ** 2
            waveimg = spec2DObj.waveimg
            bpmmask = spec2DObj.bpmmask

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
            if mnmx_wv is None:
                mnmx_wv = np.zeros((len(self.spec2d), slits.nslits, 2))
            for slit_idx, slit_spat in enumerate(slits.spat_id):
                onslit_init = (slitid_img_init == slit_spat)
                mnmx_wv[ff, slit_idx, 0] = np.min(waveimg[onslit_init])
                mnmx_wv[ff, slit_idx, 1] = np.max(waveimg[onslit_init])

            # Remove edges of the spectrum where the sky model is bad
            sky_is_good = datacube.make_good_skymask(slitid_img_init, spec2DObj.tilts)

            # Construct a good pixel mask
            # TODO: This should use the mask function to figure out which elements are masked.
            onslit_gpm = (slitid_img_init > 0) & (bpmmask.mask == 0) & sky_is_good

            # Grab the WCS of this frame
            frame_wcs = self.spec.get_wcs(spec2DObj.head0, slits, detector.platescale, wave0, dwv)
            self.all_wcs.append(copy.deepcopy(frame_wcs))

            # Find the largest spatial scale of all images being combined
            # TODO :: probably need to put this in the DetectorContainer
            pxscl = detector.platescale * parse.parse_binning(detector.binning)[
                1] / 3600.0  # This should be degrees/pixel
            slscl = self.spec.get_meta_value([spec2DObj.head0], 'slitwid')
            if dspat is None:
                dspat = max(pxscl, slscl)
            if pxscl > dspat:
                msgs.warn("Spatial scale requested ({0:f} arcsec) is less than the pixel scale ({1:f} arcsec)".format(
                    3600.0 * dspat, 3600.0 * pxscl))
            if slscl > dspat:
                msgs.warn("Spatial scale requested ({0:f} arcsec) is less than the slicer scale ({1:f} arcsec)".format(
                    3600.0 * dspat, 3600.0 * slscl))

            # Generate the alignment splines, and then
            # retrieve images of the RA and Dec of every pixel,
            # and the number of spatial pixels in each slit
            alignSplines = self.get_alignments(spec2DObj, slits, frame_wcs, spat_flexure=spat_flexure)
            raimg, decimg, minmax = slits.get_radec_image(frame_wcs, alignSplines, spec2DObj.tilts,
                                                          initial=True, flexure=spat_flexure)

            # Perform the DAR correction
            self.compute_DAR(spec2DObj.head0, raimg, decimg, waveimg, onslit_gpm, wave_ref=wave_ref)

            # Get copies of arrays to be saved
            wave_ext = waveimg[onslit_gpm].copy()
            flux_ext = sciImg[onslit_gpm].copy()
            ivar_ext = ivar[onslit_gpm].copy()
            dwav_ext = dwaveimg[onslit_gpm].copy()

            # Correct for sensitivity as a function of grating angle
            # (this assumes the spectrum of the flatfield lamp has the same shape for all setups)
            key = flatfield.FlatImages.calib_type.upper()
            if key not in spec2DObj.calibs:
                msgs.error('Processed flat calibration file not recorded by spec2d file!')
            flatfile = os.path.join(spec2DObj.calibs['DIR'], spec2DObj.calibs[key])
            if cubepar['grating_corr'] and flatfile not in flat_splines.keys():
                msgs.info("Calculating relative sensitivity for grating correction")
                # Check if the Flat file exists
                if not os.path.exists(flatfile):
                    msgs.error("Grating correction requested, but the following file does not exist:" +
                               msgs.newline() + flatfile)
                # Load the Flat file
                flatimages = flatfield.FlatImages.from_file(flatfile)
                total_illum = flatimages.fit2illumflat(slits, finecorr=False, frametype='illum', initial=True,
                                                       spat_flexure=spat_flexure) * \
                              flatimages.fit2illumflat(slits, finecorr=True, frametype='illum', initial=True,
                                                       spat_flexure=spat_flexure)
                flatframe = flatimages.pixelflat_raw / total_illum
                if flatimages.pixelflat_spec_illum is None:
                    # Calculate the relative scale
                    scale_model = flatfield.illum_profile_spectral(flatframe, waveimg, slits,
                                                                   slit_illum_ref_idx=flatpar['slit_illum_ref_idx'],
                                                                   model=None,
                                                                   skymask=None, trim=flatpar['slit_trim'],
                                                                   flexure=spat_flexure,
                                                                   smooth_npix=flatpar['slit_illum_smooth_npix'])
                else:
                    msgs.info("Using relative spectral illumination from FlatImages")
                    scale_model = flatimages.pixelflat_spec_illum
                # Apply the relative scale and generate a 1D "spectrum"
                onslit = waveimg != 0
                wavebins = np.linspace(np.min(waveimg[onslit]), np.max(waveimg[onslit]), slits.nspec)
                hist, edge = np.histogram(waveimg[onslit], bins=wavebins,
                                          weights=flatframe[onslit] / scale_model[onslit])
                cntr, edge = np.histogram(waveimg[onslit], bins=wavebins)
                cntr = cntr.astype(float)
                norm = (cntr != 0) / (cntr + (cntr == 0))
                spec_spl = hist * norm
                wave_spl = 0.5 * (wavebins[1:] + wavebins[:-1])
                flat_splines[flatfile] = interp1d(wave_spl, spec_spl, kind='linear',
                                                  bounds_error=False, fill_value="extrapolate")
                flat_splines[flatfile + "_wave"] = wave_spl.copy()
                # Check if a reference blaze spline exists (either from a standard star if fluxing or from a previous
                # exposure in this for loop)
                if blaze_spline is None:
                    blaze_wave, blaze_spec = wave_spl, spec_spl
                    blaze_spline = interp1d(wave_spl, spec_spl, kind='linear',
                                            bounds_error=False, fill_value="extrapolate")

            # Perform extinction correction
            msgs.info("Applying extinction correction")
            longitude = self.spec.telescope['longitude']
            latitude = self.spec.telescope['latitude']
            airmass = spec2DObj.head0[self.spec.meta['airmass']['card']]
            extinct = flux_calib.load_extinction_data(longitude, latitude, self.senspar['UVIS']['extinct_file'])
            # extinction_correction requires the wavelength is sorted
            wvsrt = np.argsort(wave_ext)
            ext_corr = flux_calib.extinction_correction(wave_ext[wvsrt] * units.AA, airmass, extinct)
            # Grating correction
            grat_corr = 1.0
            if self.cubepar['grating_corr']:
                grat_corr = correct_grating_shift(wave_ext[wvsrt], flat_splines[flatfile + "_wave"],
                                                  flat_splines[flatfile],
                                                  blaze_wave, blaze_spline)
            # Sensitivity function
            sens_func = 1.0
            if self.fluxcal:
                msgs.info("Calculating the sensitivity function")
                sens_func = flux_spline(wave_ext[wvsrt])
            # Convert the flux_sav to counts/s,  correct for the relative sensitivity of different setups
            ext_corr *= sens_func / (exptime * grat_corr)
            # Correct for extinction
            flux_sav = flux_ext[wvsrt] * ext_corr
            ivar_sav = ivar_ext[wvsrt] / ext_corr ** 2

            # Convert units to Counts/s/Ang/arcsec2
            # Slicer sampling * spatial pixel sampling
            sl_deg = np.sqrt(frame_wcs.wcs.cd[0, 0] ** 2 + frame_wcs.wcs.cd[1, 0] ** 2)
            px_deg = np.sqrt(frame_wcs.wcs.cd[1, 1] ** 2 + frame_wcs.wcs.cd[0, 1] ** 2)
            scl_units = dwav_ext[wvsrt] * (3600.0 * sl_deg) * (3600.0 * px_deg)
            flux_sav /= scl_units
            ivar_sav *= scl_units ** 2

            # sort back to the original ordering
            resrt = np.argsort(wvsrt)
            numpix = raimg[onslit_gpm].size

            # Calculate the weights relative to the zeroth cube
            weights[ff] = 1.0  # exptime  #np.median(flux_sav[resrt]*np.sqrt(ivar_sav[resrt]))**2

            # Get the slit image and then unset pixels in the slit image that are bad
            this_specpos, this_spatpos = np.where(onslit_gpm)
            this_spatid = slitid_img_init[onslit_gpm]

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
                # Generate the output WCS for the datacube
                crval_wv = self.cubepar['wave_min'] if self.cubepar['wave_min'] is not None else 1.0E10 * frame_wcs.wcs.crval[2]
                cd_wv = self.cubepar['wave_delta'] if self.cubepar['wave_delta'] is not None else 1.0E10 * frame_wcs.wcs.cd[2, 2]
                output_wcs = self.spec.get_wcs(spec2DObj.head0, slits, detector.platescale, crval_wv, cd_wv)
                # Set the wavelength range of the white light image.
                wl_wvrng = None
                if self.cubepar['save_whitelight']:
                    wl_wvrng = datacube.get_whitelight_range(np.max(mnmx_wv[ff, :, 0]),
                                                    np.min(mnmx_wv[ff, :, 1]),
                                                    self.cubepar['whitelight_range'])
                # Make the datacube
                if self.method in ['subpixel', 'ngp']:
                    # Generate the datacube
                    generate_cube_subpixel(outfile, output_wcs, raimg[onslit_gpm], decimg[onslit_gpm], wave_ext,
                                           flux_sav[resrt], ivar_sav[resrt], np.ones(numpix),
                                           this_spatpos, this_specpos, this_spatid,
                                           spec2DObj.tilts, slits, alignSplines, bins,
                                           all_idx=None, overwrite=self.overwrite, blaze_wave=blaze_wave,
                                           blaze_spec=blaze_spec,
                                           fluxcal=self.fluxcal, specname=self.specname, whitelight_range=wl_wvrng,
                                           spec_subpixel=spec_subpixel, spat_subpixel=spat_subpixel)
                continue

            # Store the information if we are combining multiple frames
            self.all_ra = np.append(self.all_ra, raimg[onslit_gpm].copy())
            self.all_dec = np.append(self.all_dec, decimg[onslit_gpm].copy())
            self.all_wave = np.append(self.all_wave, wave_ext.copy())
            self.all_sci = np.append(self.all_sci, flux_sav[resrt].copy())
            self.all_ivar = np.append(self.all_ivar, ivar_sav[resrt].copy())
            self.all_idx = np.append(self.all_idx, ff * np.ones(numpix))
            self.all_wghts = np.append(self.all_wghts, weights[ff] * np.ones(numpix) / weights[0])
            self.all_spatpos = np.append(self.all_spatpos, this_spatpos)
            self.all_specpos = np.append(self.all_specpos, this_specpos)
            self.all_spatid = np.append(self.all_spatid, this_spatid)
            self.all_tilts.append(spec2DObj.tilts)
            self.all_slits.append(slits)
            self.all_align.append(alignSplines)

    def run_align(self):
        """
        TODO :: Add docstring
        """
        # Grab cos(dec) for convenience
        cosdec = np.cos(np.mean(self.all_dec) * np.pi / 180.0)

        # Register spatial offsets between all frames
        if self.opts['ra_offset'] is not None:
            # First, translate all coordinates to the coordinates of the first frame
            # Note :: Don't need cosdec here, this just overrides the IFU coordinate centre of each frame
            ref_shift_ra = self.ifu_ra[0] - self.ifu_ra
            ref_shift_dec = self.ifu_dec[0] - self.ifu_dec
            for ff in range(self.numfiles):
                # Apply the shift
                self.all_ra[self.all_idx == ff] += ref_shift_ra[ff] + self.opts['ra_offset'][ff]/3600.0
                self.all_dec[self.all_idx == ff] += ref_shift_dec[ff] + self.opts['dec_offset'][ff]/3600.0
                msgs.info("Spatial shift of cube #{0:d}: RA, DEC (arcsec) = {1:+0.3f} E, {2:+0.3f} N".format(ff + 1, opts['ra_offset'][ff], opts['dec_offset'][ff]))
        else:
            # Find the wavelength range where all frames overlap
            min_wl, max_wl = datacube.get_whitelight_range(np.max(mnmx_wv[:, :, 0]),  # The max blue wavelength
                                                           np.min(mnmx_wv[:, :, 1]),  # The min red wavelength
                                                           self.cubepar['whitelight_range'])  # The user-specified values (if any)
            # Get the good whitelight pixels
            ww, wavediff = get_whitelight_pixels(self.all_wave, min_wl, max_wl)
            # Iterate over white light image generation and spatial shifting
            numiter = 2
            for dd in range(numiter):
                msgs.info(f"Iterating on spatial translation - ITERATION #{dd+1}/{numiter}")
                # Setup the WCS to use for all white light images
                ref_idx = None  # Don't use an index - This is the default behaviour when a reference image is supplied
                image_wcs, voxedge, reference_image = create_wcs(cubepar, all_ra[ww], all_dec[ww], all_wave[ww],
                                                                 dspat, wavediff, collapse=True)
                if voxedge[2].size != 2:
                    msgs.error("Spectral range for WCS is incorrect for white light image")

                wl_imgs = generate_image_subpixel(image_wcs, all_ra[ww], all_dec[ww], all_wave[ww],
                                                  all_sci[ww], all_ivar[ww], all_wghts[ww],
                                                  all_spatpos[ww], all_specpos[ww], all_spatid[ww],
                                                  all_tilts, all_slits, all_align, voxedge, all_idx=all_idx[ww],
                                                  spec_subpixel=spec_subpixel, spat_subpixel=spat_subpixel)
                if reference_image is None:
                    # ref_idx will be the index of the cube with the highest S/N
                    ref_idx = np.argmax(weights)
                    reference_image = wl_imgs[:, :, ref_idx].copy()
                    msgs.info("Calculating spatial translation of each cube relative to cube #{0:d})".format(ref_idx+1))
                else:
                    msgs.info("Calculating the spatial translation of each cube relative to user-defined 'reference_image'")

                # Calculate the image offsets relative to the reference image
                for ff in range(self.numfiles):
                    # Calculate the shift
                    ra_shift, dec_shift = calculate_image_phase(reference_image.copy(), wl_imgs[:, :, ff], maskval=0.0)
                    # Convert pixel shift to degrees shift
                    ra_shift *= dspat/cosdec
                    dec_shift *= dspat
                    msgs.info("Spatial shift of cube #{0:d}: RA, DEC (arcsec) = {1:+0.3f} E, {2:+0.3f} N".format(ff+1, ra_shift*3600.0, dec_shift*3600.0))
                    # Apply the shift
                    all_ra[all_idx == ff] += ra_shift
                    all_dec[all_idx == ff] += dec_shift

    def compute_weights(self):
        # Calculate the relative spectral weights of all pixels
        if self.numfiles == 1:
            # No need to calculate weights if there's just one frame
            self.all_wghts = np.ones_like(self.all_sci)
        else:
            # Find the wavelength range where all frames overlap
            min_wl, max_wl = datacube.get_whitelight_range(np.max(mnmx_wv[:, :, 0]),  # The max blue wavelength
                                                  np.min(mnmx_wv[:, :, 1]),  # The min red wavelength
                                                  self.cubepar['whitelight_range'])  # The user-specified values (if any)
            # Get the good white light pixels
            ww, wavediff = datacube.get_whitelight_pixels(all_wave, min_wl, max_wl)
            # Get a suitable WCS
            image_wcs, voxedge, reference_image = create_wcs(cubepar, all_ra, all_dec, all_wave, dspat, wavediff,
                                                             collapse=True)
            # Generate the white light image (note: hard-coding subpixel=1 in both directions, and combining into a single image)
            wl_full = generate_image_subpixel(image_wcs, all_ra, all_dec, all_wave,
                                              all_sci, all_ivar, all_wghts,
                                              all_spatpos, all_specpos, all_spatid,
                                              all_tilts, all_slits, all_align, voxedge, all_idx=all_idx,
                                              spec_subpixel=1, spat_subpixel=1, combine=True)
            # Compute the weights
            all_wghts = datacube.compute_weights(all_ra, all_dec, all_wave, all_sci, all_ivar, all_idx, wl_full[:, :, 0],
                                                 dspat, dwv, relative_weights=self.cubepar['relative_weights'])

    def coadd(self):
        """
        TODO :: Add docstring
        """
        # First loop through all of the frames, load the data, and save datacubes if no combining is required
        self.load()

        # No need to continue if we are not combining nor aligning frames
        if not self.combine and not self.align:
            return

        # Align the frames
        if self.align:
            self.run_align()

        # Compute the relative weights on the spectra
        self.compute_weights()

        # Generate the WCS, and the voxel edges
        cube_wcs, vox_edges, _ = datacube.create_wcs(self.cubepar, self.all_ra, self.all_dec, self.all_wave, dspat, dwv)

        sensfunc = None
        if self.flux_spline is not None:
            # Get wavelength of each pixel, and note that the WCS gives this in m, so convert to Angstroms (x 1E10)
            numwav = vox_edges[2].size - 1
            senswave = cube_wcs.spectral.wcs_pix2world(np.arange(numwav), 0)[0] * 1.0E10
            sensfunc = self.flux_spline(senswave)

        # Generate a datacube
        outfile = datacube.get_output_filename("", self.cubepar['output_filename'], True, -1)
        if self.method in ['subpixel', 'ngp']:
            # Generate the datacube
            wl_wvrng = None
            if self.cubepar['save_whitelight']:
                wl_wvrng = datacube.get_whitelight_range(np.max(mnmx_wv[:, :, 0]),
                                                np.min(mnmx_wv[:, :, 1]),
                                                self.cubepar['whitelight_range'])
            if self.combine:
                generate_cube_subpixel(outfile, cube_wcs, all_ra, all_dec, all_wave, all_sci, all_ivar,
                                       np.ones(all_wghts.size),  # all_wghts,
                                       all_spatpos, all_specpos, all_spatid, all_tilts, all_slits, all_align, vox_edges,
                                       all_idx=all_idx, overwrite=overwrite, blaze_wave=blaze_wave,
                                       blaze_spec=blaze_spec,
                                       fluxcal=fluxcal, sensfunc=sensfunc, specname=specname, whitelight_range=wl_wvrng,
                                       spec_subpixel=spec_subpixel, spat_subpixel=spat_subpixel)
            else:
                for ff in range(self.numfiles):
                    outfile = datacube.get_output_filename("", self.cubepar['output_filename'], False, ff)
                    ww = np.where(self.all_idx == ff)
                    generate_cube_subpixel(outfile, cube_wcs, all_ra[ww], all_dec[ww], all_wave[ww], all_sci[ww],
                                           all_ivar[ww], np.ones(all_wghts[ww].size),
                                           all_spatpos[ww], all_specpos[ww], all_spatid[ww], all_tilts[ff],
                                           all_slits[ff], all_align[ff], vox_edges,
                                           all_idx=all_idx[ww], overwrite=overwrite, blaze_wave=blaze_wave,
                                           blaze_spec=blaze_spec,
                                           fluxcal=fluxcal, sensfunc=sensfunc, specname=specname,
                                           whitelight_range=wl_wvrng,
                                           spec_subpixel=spec_subpixel, spat_subpixel=spat_subpixel)
