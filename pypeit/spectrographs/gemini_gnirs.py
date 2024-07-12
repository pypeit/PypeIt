"""
Module for Gemini/GNIRS specific methods.

.. include:: ../include/links.rst
"""
import numpy as np
from astropy import wcs, units
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.time import Time

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch, parse
from pypeit.images import detector_container
from pypeit.spectrographs import spectrograph


class GeminiGNIRSSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Gemini/GNIRS specific code
    """
    ndet = 1
    camera = 'GNIRS'
    url = 'https://www.gemini.edu/instrumentation/gnirs'
    header_name = 'GNIRS'
    telescope = telescopes.GeminiNTelescopePar()

    def __init__(self):
        super().__init__()

        # TODO :: Might consider changing TelescopePar to use the astropy EarthLocation.
        self.location = EarthLocation.of_site('Gemini North')

    def get_detector_par(self, det, hdu=None):
        """
        Return metadata for the selected detector.

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
        # Detector 1
        detector_dict = dict(
            binning         = '1,1',
            det             = 1,
            dataext         = 1,
            specaxis        = 0,
            specflip=True,
            spatflip=True,
            platescale      = 0.15,
            darkcurr        = 540.0,  # e-/hour/pixel  (=0.15 e-/pixel/s)
            saturation      = 150000.,
            nonlinear       = 0.71,
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d(13.5),
            ronoise         = np.atleast_1d(7.0),
            datasec         = np.atleast_1d('[:,:]'),
            oscansec        = None,
        )
        return detector_container.DetectorContainer(**detector_dict)

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=0, card='RA')
        self.meta['dec'] = dict(ext=0, card='DEC')
        self.meta['target'] = dict(ext=0, card='OBJECT')
        self.meta['decker'] = dict(ext=0, card='SLIT')

        self.meta['binning'] = dict(ext=0, card=None, default='1,1')
        self.meta['mjd'] = dict(ext=0, card='MJD_OBS')
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')

        # Dithering
        self.meta['dithpos'] = dict(ext=0, card='QOFFSET')
        self.meta['dithoff'] = dict(card=None, compound=True)

        # Extras for config and frametyping
        self.meta['filter1'] = dict(ext=0, card='FILTER2')
        self.meta['slitwid'] = dict(ext=0, compound=True, card=None)
        self.meta['dispname'] = dict(ext=0, card='GRATING')
        self.meta['hatch'] = dict(ext=0, card='COVER')
        self.meta['dispangle'] = dict(ext=0, card='GRATTILT', rtol=1e-4)
        self.meta['idname'] = dict(ext=0, card='OBSTYPE')
        self.meta['instrument'] = dict(ext=0, card='INSTRUME')

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
        if meta_key == 'dithoff':
            if headarr[0].get('OBSTYPE') == 'OBJECT':
                return headarr[0].get('QOFFSET')
            else:
                return 0.0
        elif meta_key == 'slitwid':
            deckname = headarr[0].get('DECKER')
            if 'LR-IFU' in deckname:
                return 0.15/3600.0  # divide by 3600 for degrees
            elif 'HR-IFU' in deckname:
                return 0.05/3600.0  # divide by 3600 for degrees
            else:
                # TODO :: Need to provide a more complete set of options here
                return None
        elif meta_key == 'obstime':
            try:
                return Time(headarr[0]['DATE-OBS'] + "T" + headarr[0]['TIME-OBS'])
            except KeyError:
                msgs.warn("Time of observation is not in header")
                return 0.0
        elif meta_key == 'pressure':
            try:
                return headarr[0]['PRESSUR2']/100.0  # Must be in astropy.units.mbar
            except KeyError:
                msgs.warn("Pressure is not in header - The default pressure (611 mbar) will be assumed")
                return 611.0
        elif meta_key == 'temperature':
            try:
                return headarr[0]['TAMBIENT']  # Must be in astropy.units.deg_C
            except KeyError:
                msgs.warn("Temperature is not in header - The default temperature (1.5 deg C) will be assumed")
                return 1.5  # van Kooten & Izett, arXiv:2208.11794
        elif meta_key == 'humidity':
            try:
                # Humidity expressed as a percentage, not a fraction
                return headarr[0]['HUMIDITY']
            except KeyError:
                msgs.warn("Humidity is not in header - The default relative humidity (20 %) will be assumed")
                return 20.0  # van Kooten & Izett, arXiv:2208.11794
        elif meta_key == 'parangle':
            try:
                # Humidity expressed as a percentage, not a fraction
                msgs.warn("Parallactic angle is not available for GNIRS - DAR correction may be incorrect")
                return headarr[0]['PARANGLE']  # Must be expressed in radians
            except KeyError:
                msgs.warn("Parallactic angle is not in header - The default parallactic angle (0 degrees) will be assumed")
                return 0.0
        else:
            msgs.error("Not ready for this compound meta")

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
        return ['decker', 'dispname', 'dispangle']

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
        return ['SLIT', 'GRATING', 'GRATTILT']

    def pypeit_file_keys(self):
        """
        Define the list of keys to be output into a standard ``PypeIt`` file.

        Returns:
            :obj:`list`: The list of keywords in the relevant
            :class:`~pypeit.metadata.PypeItMetaData` instance to print to the
            :ref:`pypeit_file`.
        """
        return super().pypeit_file_keys() + ['dithoff']

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
        if ftype == 'science':
            return good_exp & (fitstbl['idname'] == 'OBJECT')
        if ftype == 'standard':
            return good_exp & (fitstbl['idname'] == 'OBJECT')
        if ftype in ['pixelflat', 'trace']:
            # Flats and trace frames are typed together
            return good_exp & (fitstbl['idname'] == 'FLAT')
        if ftype in ['pinhole', 'dark', 'bias']:
            # Don't type pinhole, dark, or bias frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc', 'tilt']:
            ## FW ToDo: self.dispname does not work yet. need to replace the following later.
            if '32/mm' in fitstbl['dispname'][0]:
                return good_exp & (fitstbl['idname'] == 'OBJECT')
            elif '10/mmLBSX' in fitstbl['dispname'][0]:
                return good_exp & (fitstbl['idname'] == 'ARC')

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.

        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        # Image processing steps
        turn_off = dict(use_illumflat=False, use_biasimage=False, use_overscan=False,
                        use_darkimage=False)
        par.reset_all_processimages_par(**turn_off)

        # Flats
        par['calibrations']['flatfield']['tweak_slits_thresh'] = 0.90
        par['calibrations']['flatfield']['tweak_slits_maxfrac'] = 0.10

        # Relatively short slit, so keep the spatial tilt order low
        par['calibrations']['tilts']['spat_order'] = 1

        # Reduce parameters
        # par['reduce']['findobj']['snr_thresh'] = 5.0          # Object finding threshold
        par['reduce']['findobj']['find_trim_edge'] = [2, 2]  # Slit is too short to trim 5,5 especially
        par['reduce']['skysub']['bspline_spacing'] = 0.8
        par['reduce']['skysub']['global_sky_std'] = False  # Do not perform global sky subtraction for standard stars
        par['reduce']['skysub']['no_poly'] = True  # Do not use polynomial degree of freedom for global skysub
        par['reduce']['extraction']['model_full_slit'] = True  # local sky subtraction operates on entire slit
        par['reduce']['findobj']['maxnumber_sci'] = 2  # Slit is narrow so allow two objects per order
        par['reduce']['findobj']['maxnumber_std'] = 1  # Slit is narrow so allow one object per order
        # Standards
        par['calibrations']['standardframe']['process']['mask_cr'] = False  # Do not mask_cr standards

        # Do not correct for flexure
        par['flexure']['spec_method'] = 'skip'

        # Set the default exposure time ranges for the frame typing
        par['calibrations']['pixelflatframe']['exprng'] = [None, 30]
        par['calibrations']['traceframe']['exprng'] = [None, 30]
        par['calibrations']['standardframe']['exprng'] = [None, 30]
        par['scienceframe']['exprng'] = [30, None]

        # Sensitivity function parameters
        par['sensfunc']['algorithm'] = 'IR'
        par['sensfunc']['polyorder'] = 6
        par['sensfunc']['IR']['telgridfile'] = 'TellPCA_3000_26000_R10000.fits'
        return par

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
        par = super().config_specific_par(scifile, inp_par=inp_par)
        # TODO This is a hack for now until we figure out how to set dispname
        # and other meta information in the spectrograph class itself
        self.dispname = self.get_meta_value(scifile, 'dispname')
        # 32/mmSB_G5533 setup, covering XYJHK with short blue camera
        if '32/mm' in self.dispname:
            # Edges
            par['calibrations']['slitedges']['edge_thresh'] = 20.
            par['calibrations']['slitedges']['trace_thresh'] = 10.
            par['calibrations']['slitedges']['fit_order'] = 5
            par['calibrations']['slitedges']['max_shift_adj'] = 0.5
            par['calibrations']['slitedges']['fit_min_spec_length'] = 0.5

            # Wavelengths
            par['calibrations']['wavelengths']['rms_thresh_frac_fwhm'] = 0.4
            par['calibrations']['wavelengths']['sigdetect'] = 10.0
            par['calibrations']['wavelengths']['lamps'] = ['OH_GNIRS']
            # par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
            par['calibrations']['wavelengths']['n_first'] = 2
            par['calibrations']['wavelengths']['n_final'] = 3

            # Reidentification parameters
            par['calibrations']['wavelengths']['method'] = 'reidentify'
            par['calibrations']['wavelengths']['cc_thresh'] = 0.6
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gemini_gnirs.fits'

            # Tilts
            par['calibrations']['tilts']['tracethresh'] = 10
            par['calibrations']['tilts']['sig_neigh'] = 5.0
            par['calibrations']['tilts']['nfwhm_neigh'] = 2.0

            # Coadding. Not for longslit data this might be problematic but that is not yet supported.
            par['coadd1d']['wave_method'] = 'log10'

        # 10/mmLBSX_G5532 setup, covering YJHK with the long blue camera and SXD prism
        elif '10/mmLBSX' in self.dispname:
            # Edges
            par['calibrations']['slitedges']['edge_thresh'] = 20.
            par['calibrations']['slitedges']['trace_thresh'] = 10.
            par['calibrations']['slitedges']['fit_order'] = 2
            par['calibrations']['slitedges']['max_shift_adj'] = 0.5
            par['calibrations']['slitedges']['det_min_spec_length'] = 0.20
            par['calibrations']['slitedges']['fit_min_spec_length'] = 0.20
            par['calibrations']['slitedges']['sync_predict'] = 'nearest'

            # Wavelengths
            par['calibrations']['wavelengths']['rms_thresh_frac_fwhm'] = 0.05
            par['calibrations']['wavelengths']['sigdetect'] = 5.0
            par['calibrations']['wavelengths']['lamps'] = ['Ar_IR_GNIRS']
            # par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
            par['calibrations']['wavelengths']['n_first'] = 2
            par['calibrations']['wavelengths']['n_final'] = 3
            # Reidentification parameters
            par['calibrations']['wavelengths']['method'] = 'reidentify'
            par['calibrations']['wavelengths']['cc_thresh'] = 0.6
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gemini_gnirs_10mm_LBSX.fits'

            # Tilts
            par['calibrations']['tilts']['tracethresh'] = 10
            par['calibrations']['tilts']['sig_neigh'] = 5.0
            par['calibrations']['tilts']['nfwhm_neigh'] = 2.0
        elif '10/mmLBHR_G5532' in self.dispname:
            # TODO :: Need to fill this in
            pass
        else:
            msgs.error(f'Unrecognized GNIRS dispname: {self.dispname}')

        return par

    def bpm(self, filename, det, shape=None, msbias=None):
        """
        Generate a default bad-pixel mask.

        Even though they are both optional, either the precise shape for
        the image (``shape``) or an example file that can be read to get
        the shape (``filename`` using :func:`get_image_shape`) *must* be
        provided.

        Args:
            filename (:obj:`str` or None):
                An example file to use to get the image shape.
            det (:obj:`int`):
                1-indexed detector number to use when getting the image
                shape from the example file.
            shape (tuple, optional):
                Processed image shape
                Required if filename is None
                Ignored if filename is not None
            msbias (`numpy.ndarray`_, optional):
                Processed bias frame used to identify bad pixels

        Returns:
            `numpy.ndarray`_: An integer array with a masked value set
            to 1 and an unmasked value set to 0.  All values are set to
            0.
        """
        msgs.info("Custom bad pixel mask for GNIRS")
        # Call the base-class method to generate the empty bpm
        bpm_img = super().bpm(filename, det, shape=shape, msbias=msbias)

        # JFH Changed. Dealing with detector scratch
        if det == 1:
            bpm_img[687:765,12:16] = 1.
            bpm_img[671:687,8:13] = 1.
        #    bpm_img[:, 1000:] = 1.

        return bpm_img


class GeminiGNIRSEchelleSpectrograph(GeminiGNIRSSpectrograph):
    """
    Child to handle Gemini/GNIRS echelle specific code
    """
    name = 'gemini_gnirs_echelle'
    pypeline = 'Echelle'
    ech_fixed_format = True

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
        par = super().config_specific_par(scifile, inp_par=inp_par)
        # TODO This is a hack for now until we figure out how to set dispname
        # and other meta information in the spectrograph class itself
        self.dispname = self.get_meta_value(scifile, 'dispname')
        # 32/mmSB_G5533 setup, covering XYJHK with short blue camera
        if '32/mm' in self.dispname:
            # Edges
            par['calibrations']['slitedges']['left_right_pca'] = True
            par['calibrations']['slitedges']['pca_order'] = 3

            # Wavelengths
            par['calibrations']['wavelengths']['sigdetect'] = [4.0, 5.0, 5.0, 5.0, 5.0, 5.0] #5.0
            par['calibrations']['wavelengths']['n_final'] = [1, 3, 3, 3, 3, 3]

            # Echelle parameters
            # JFH This is provisional these IDs should be checked.
            par['calibrations']['wavelengths']['echelle'] = True
            par['calibrations']['wavelengths']['ech_nspec_coeff'] = 3
            par['calibrations']['wavelengths']['ech_norder_coeff'] = 5
            par['calibrations']['wavelengths']['ech_sigrej'] = 3.0

            # Tilts
            par['calibrations']['tilts']['tracethresh'] = [5.0, 10, 10, 10, 10, 10]
        # 10/mmLBSX_G5532 setup, covering YJHK with the long blue camera and SXD prism
        elif '10/mmLBSX' in self.dispname:
            # Edges
            par['calibrations']['slitedges']['left_right_pca'] = True # Actually we need a parameter to disable PCA entirely

            # Wavelengths
            par['calibrations']['wavelengths']['n_final'] = [3, 3, 3, 3]
            # Echelle parameters
            par['calibrations']['wavelengths']['echelle'] = True
            par['calibrations']['wavelengths']['ech_nspec_coeff'] = 3
            par['calibrations']['wavelengths']['ech_norder_coeff'] = 3
            par['calibrations']['wavelengths']['ech_sigrej'] = 3.0

            # Tilts
            par['calibrations']['tilts']['tracethresh'] = [10, 10, 10, 10]
        else:
            msgs.error('Unrecognized GNIRS dispname')

        return par

    def order_platescale(self, order_vec, binning=None):
        """
        Return the platescale for each echelle order.

        Args:
            order_vec (`numpy.ndarray`_):
                The vector providing the order numbers.
            binning (:obj:`str`, optional):
                The string defining the spectral and spatial binning.

        Returns:
            `numpy.ndarray`_: An array with the platescale for each order
            provided by ``order``.
        """
        # TODO: Binning is ignored.  Should it be?
        self.check_disperser()
        if '10/mmLBSX' in self.dispname:
            return np.full(order_vec.size, 0.05)
        elif '32/mm' in self.dispname:
            return np.full(order_vec.size, 0.15)
        else:
            msgs.error('Unrecognized disperser')

    @property
    def norders(self):
        """
        Number of orders for this spectograph.
        """
        self.check_disperser()
        if '10/mmLBSX' in self.dispname:
            return 4
        elif '32/mm' in self.dispname:
            return 6
        else:
            msgs.error('Unrecognized disperser')

    @property
    def order_spat_pos(self):
        """
        Return the expected spatial position of each echelle order.
        """
        self.check_disperser()
        if '10/mmLBSX' in self.dispname:
            return np.array([0.050, 0.215, 0.442, 0.759])
        elif '32/mm' in self.dispname:
            #ToDo: create self.date similar to self.dispname and use that to decide which numbers to use
            ## Old data, i.e. before 2011
            #return np.array([0.241211 , 0.3173828, 0.387695, 0.456054, 0.530273, 0.640625])
            ##New data
            return np.array([0.2955097 , 0.37635756, 0.44952223, 0.51935601, 0.59489503, 0.70210309])
        else:
            msgs.error('Unrecognized disperser')

    @property
    def orders(self):
        """
        Return the order number for each echelle order.
        """
        self.check_disperser()
        if '10/mmLBSX' in self.dispname:
            return np.arange(6,2,-1, dtype=int)
        elif '32/mm' in self.dispname:
            return np.arange(8,2,-1,dtype=int)
        else:
            msgs.error('Unrecognized disperser')

    @property
    def spec_min_max(self):
        """
        Return the minimum and maximum spectral pixel expected for the
        spectral range of each order.
        """
        # TODO: Why aren't these numbers in fraction of the detector
        # size instead of the pixel number?
        self.check_disperser()
        if '10/mmLBSX' in self.dispname:
            spec_max = np.asarray([1022, 1022, 1022, 1022])
            spec_min = np.asarray([450, 0, 0, 0])
            return np.vstack((spec_min, spec_max))
        elif '32/mm' in self.dispname:
            spec_max = np.asarray([1022, 1022, 1022, 1022, 1022, 1022])
            spec_min = np.asarray([512, 280, 0, 0, 0, 0])
            return np.vstack((spec_min, spec_max))
        else:
            msgs.error('Unrecognized disperser')


class GNIRSIFUSpectrograph(GeminiGNIRSSpectrograph):
    # TODO :: A list of steps that could improve the reduction
    #   * Have a high threshold for detecting slit edges (par['calibrations']['slitedges']['edge_thresh'] = 100.), and have an option when inserting new traces to be the median of all other slit lengths (or a fit to the slit lengths).
    #   * Need to store a wavelength solution for different grating options (Note, the Holy Grail algorithm works pretty well, most of the time)
    name = 'gemini_gnirs_ifu'
    pypeline = 'SlicerIFU'

    def init_meta(self):
        super().init_meta()
        self.meta['obstime'] = dict(card=None, compound=True, required=False)
        self.meta['pressure'] = dict(card=None, compound=True, required=False)
        self.meta['temperature'] = dict(card=None, compound=True, required=False)
        self.meta['humidity'] = dict(card=None, compound=True, required=False)
        self.meta['parangle'] = dict(card=None, compound=True, required=False)

    @classmethod
    def default_pypeit_par(cls):
        par = super().default_pypeit_par()

        # LACosmics parameters
        par['scienceframe']['process']['sigclip'] = 4.0
        par['scienceframe']['process']['objlim'] = 1.5
        par['scienceframe']['process']['use_illumflat'] = False  # illumflat is applied when building the relative scale image in reduce.py, so should be applied to scienceframe too.
        par['scienceframe']['process']['use_specillum'] = False  # apply relative spectral illumination
        par['scienceframe']['process']['spat_flexure_correct'] = False  # don't correct for spatial flexure - varying spatial illumination profile could throw this correction off. Also, there's no way to do astrometric correction if we can't correct for spatial flexure of the contbars frames
        par['scienceframe']['process']['use_biasimage'] = False
        par['scienceframe']['process']['use_darkimage'] = False
        par['calibrations']['flatfield']['slit_illum_finecorr'] = False
        # Don't do 1D extraction for 3D data - it's meaningless because the DAR correction must be performed on the 3D data.
        par['reduce']['extraction']['skip_extraction'] = True  # Because extraction occurs before the DAR correction, don't extract

        #par['calibrations']['flatfield']['tweak_slits'] = False  # Do not tweak the slit edges (we want to use the full slit)
        par['calibrations']['flatfield']['tweak_slits_thresh'] = 0.0  # Make sure the full slit is used (i.e. when the illumination fraction is > 0.5)
        par['calibrations']['flatfield']['tweak_slits_maxfrac'] = 0.0  # Make sure the full slit is used (i.e. no padding)
        par['calibrations']['flatfield']['slit_trim'] = 2  # Trim the slit edges
        par['calibrations']['slitedges']['pad'] = 2  # Need to pad out the tilts for the astrometric transform when creating a datacube.

        # Decrease the wave tilts order, given the shorter slits of the IFU
        par['calibrations']['tilts']['spat_order'] = 1
        par['calibrations']['tilts']['spec_order'] = 1

        # Make sure that this is reduced as a slit (as opposed to fiber) spectrograph
        par['reduce']['cube']['slit_spec'] = True
        par['reduce']['cube']['grating_corr'] = False
        par['reduce']['cube']['combine'] = False  # Make separate spec3d files from the input spec2d files

        # Sky subtraction parameters
        par['reduce']['skysub']['no_poly'] = True
        par['reduce']['skysub']['bspline_spacing'] = 0.6
        par['reduce']['skysub']['joint_fit'] = False

        # Don't correct flexure by default since the OH lines are used for wavelength calibration
        # If someone does want to do a spectral flexure correction, you should use slitcen,
        # because this is a slit-based IFU where no objects are extracted.
        par['flexure']['spec_method'] = 'skip'
        par['flexure']['spec_maxshift'] = 0  # The sky lines are used for calibration - don't allow flexure

        # Flux calibration parameters
        par['sensfunc']['UVIS']['extinct_correct'] = False  # This must be False - the extinction correction is performed when making the datacube

        return par

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
        par = super().config_specific_par(scifile, inp_par=inp_par)
        # Obtain a header keyword to determine which range is being used
        filter = self.get_meta_value(scifile, 'filter1')
        par['calibrations']['slitedges']['edge_thresh'] = 30.
        # TODO :: The following wavelength solutions are not general enough - need to implement a solution for each setup+grating
        # TODO BEFORE PR MERGE :: The full_template solutions below were generated (quickly!) from holy-grail... might want to redo this...
        if filter == 'X_G0518':  # H band
            par['calibrations']['wavelengths']['method'] = 'holy-grail'
        elif filter == 'J_G0517':  # K band
            par['calibrations']['wavelengths']['method'] = 'holy-grail'
        elif filter == 'H_G0516':  # H band
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gemini_gnirs_lrifu_H.fits'
        elif filter == 'K_G0515':  # K band
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gemini_gnirs_lrifu_K.fits'
        else:
            par['calibrations']['wavelengths']['method'] = 'holy-grail'

        return par

    def get_wcs(self, hdr, slits, platescale, wave0, dwv, spatial_scale=None):
        """
        Construct/Read a World-Coordinate System for a frame.

        Args:
            hdr (`astropy.io.fits.Header`_):
                The header of the raw frame. The information in this
                header will be extracted and returned as a WCS.
            slits (:class:`~pypeit.slittrace.SlitTraceSet`):
                Slit traces.
            platescale (:obj:`float`):
                The platescale of an unbinned pixel in arcsec/pixel (e.g.
                detector.platescale).
            wave0 (:obj:`float`):
                The wavelength zeropoint.
            dwv (:obj:`float`):
                Change in wavelength per spectral pixel.

        Returns:
            `astropy.wcs.WCS`_: The world-coordinate system.
        """
        msgs.info("Calculating the WCS")
        # Get the x and y binning factors, and the typical slit length
        binspec, binspat = parse.parse_binning(self.get_meta_value([hdr], 'binning'))

        # Get the pixel and slice scales
        pxscl = platescale * binspat / 3600.0  # Need to convert arcsec to degrees
        msgs.work("NEED TO WORK OUT SLICER SCALE AND PIXEL SCALE")
        slscl = self.get_meta_value([hdr], 'slitwid')
        if spatial_scale is not None:
            if pxscl > spatial_scale / 3600.0:
                msgs.warn("Spatial scale requested ({0:f}'') is less than the pixel scale ({1:f}'')".format(spatial_scale, pxscl*3600.0))
            # Update the pixel scale
            pxscl = spatial_scale / 3600.0  # 3600 is to convert arcsec to degrees

        # Get the typical slit length (this changes by ~0.3% over all slits, so a constant is fine for now)
        slitlength = int(np.round(np.median(slits.get_slitlengths(initial=True, median=True))))

        # Get RA/DEC
        raval = self.get_meta_value([hdr], 'ra')
        decval = self.get_meta_value([hdr], 'dec')

        # Create a coordinate
        coord = SkyCoord(raval, decval, unit=(units.deg, units.deg))

        # Get rotator position
        msgs.warn("CURRENTLY A HACK --- NEED TO FIGURE OUT RPOS and RREF FOR HRIFU FROM HEADER INFO")
        if 'ROTPOSN' in hdr:
            rpos = hdr['ROTPOSN']
        else:
            rpos = 0.
        if 'ROTREFAN' in hdr:
            rref = hdr['ROTREFAN']
        else:
            rref = 0.
        # Get the offset and PA
        rotoff = 0.0  # IFU-SKYPA offset (degrees)
        skypa = rpos + rref  # IFU position angle (degrees)
        crota = np.radians(-(skypa + rotoff))

        # Calculate the fits coordinates
        cdelt1 = -slscl
        cdelt2 = pxscl
        if coord is None:
            ra = 0.
            dec = 0.
            crota = 1
        else:
            ra = coord.ra.degree
            dec = coord.dec.degree
        # Calculate the CD Matrix
        cd11 = cdelt1 * np.cos(crota)                          # RA degrees per column
        cd12 = abs(cdelt2) * np.sign(cdelt1) * np.sin(crota)   # RA degrees per row
        cd21 = -abs(cdelt1) * np.sign(cdelt2) * np.sin(crota)  # DEC degress per column
        cd22 = cdelt2 * np.cos(crota)                          # DEC degrees per row
        # Get reference pixels (set these to the middle of the FOV)
        crpix1 = 11   # i.e. see get_datacube_bins (11 is used as the reference point - somewhere in the middle of the FOV)
        crpix2 = slitlength / 2.
        crpix3 = 1.
        # Get the offset
        msgs.warn("HACK FOR HRIFU --- Need to obtain offset from header?")
        off1 = 0.
        off2 = 0.
        off1 /= binspec
        off2 /= binspat
        crpix1 += off1
        crpix2 += off2

        # Create a new WCS object.
        msgs.info("Generating GNIRS IFU WCS")
        w = wcs.WCS(naxis=3)
        w.wcs.equinox = hdr['EQUINOX']
        w.wcs.name = 'GNIRS IFU'
        w.wcs.radesys = 'FK5'
        # Insert the coordinate frame
        w.wcs.cname = ['RA', 'DEC', 'Wavelength']
        w.wcs.cunit = [units.degree, units.degree, units.Angstrom]
        w.wcs.ctype = ["RA---TAN", "DEC--TAN", "WAVE"]  # Note, WAVE is in vacuum
        w.wcs.crval = [ra, dec, wave0]  # RA, DEC, and wavelength zeropoints
        w.wcs.crpix = [crpix1, crpix2, crpix3]  # RA, DEC, and wavelength reference pixels
        w.wcs.cd = np.array([[cd11, cd12, 0.0], [cd21, cd22, 0.0], [0.0, 0.0, dwv]])
        w.wcs.lonpole = 180.0  # Native longitude of the Celestial pole
        w.wcs.latpole = 0.0  # Native latitude of the Celestial pole

        return w

    def get_datacube_bins(self, slitlength, minmax, num_wave):
        r"""
        Calculate the bin edges to be used when making a datacube.

        Args:
            slitlength (:obj:`int`):
                Length of the slit in pixels
            minmax (`numpy.ndarray`_):
                An array with the minimum and maximum pixel locations on each
                slit relative to the reference location (usually the centre
                of the slit). Shape must be :math:`(N_{\rm slits},2)`, and is
                typically the array returned by
                :func:`~pypeit.slittrace.SlitTraceSet.get_radec_image`.
            num_wave (:obj:`int`):
                Number of wavelength steps.  Given by::
                    int(round((wavemax-wavemin)/delta_wave))

        Args:
            :obj:`tuple`: Three 1D `numpy.ndarray`_ providing the bins to use
            when constructing a histogram of the spec2d files. The elements
            are :math:`(x,y,\lambda)`.
        """
        # TODO :: The HRIFU might have 25 slits with 13 being the reference
        xbins = np.arange(1 + 21) - 11.0 - 0.5  # 21 is for 21 slices, and 11 is the reference slit
        ybins = np.linspace(np.min(minmax[:, 0]), np.max(minmax[:, 1]), 1+slitlength) - 0.5
        spec_bins = np.arange(1+num_wave) - 0.5
        return xbins, ybins, spec_bins

    def pypeit_file_keys(self):
        """
        Define the list of keys to be output into a standard ``PypeIt`` file.

        Returns:
            :obj:`list`: The list of keywords in the relevant
            :class:`~pypeit.metadata.PypeItMetaData` instance to print to the
            :ref:`pypeit_file`.
        """
        return super().pypeit_file_keys() + ['filter']
