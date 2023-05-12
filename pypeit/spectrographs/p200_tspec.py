"""
Module for P200/Triplespec specific methods.

.. include:: ../include/links.rst
"""
import numpy as np

from astropy.time import Time

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container


class P200TSPECSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle P200/TripleSpec specific code
    """
    ndet = 1
    name = 'p200_tspec'
    telescope = telescopes.P200TelescopePar()
    camera = 'TSPEC'
    url = 'https://sites.astro.caltech.edu/palomar/observer/200inchResources/tspeccookbook.html'
    header_name = 'TSPEC_SPEC'
    pypeline = 'Echelle'
    ech_fixed_format = True
    supported = True
    comment = 'TripleSpec spectrograph'

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
        self.meta['decker'] = dict(ext=0, card=None, default='default')
        self.meta['binning'] = dict(ext=0, card=None, default='1,1')

        self.meta['mjd'] = dict(ext=0, card=None, compound=True)
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Extras for config and frametyping
        self.meta['dispname'] = dict(ext=0, card='FPA')
        self.meta['idname'] = dict(ext=0, card='OBSTYPE')
        self.meta['instrument'] = dict(ext=0, card='FPA')

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
        if meta_key == 'mjd':
            time = headarr[0]['UTSHUT']
            ttime = Time(time, format='isot')
            return ttime.mjd
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
        return ['dispname']

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
        return ['FPA']

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
            binning='1,1',
            det=1,
            dataext         = 0,
            specaxis        = 1,
            specflip        = True,
            spatflip=False,
            platescale      = 0.37,
            darkcurr        = 0.085,
            saturation      = 28000,
            nonlinear       = 0.9,
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d(3.8),
            ronoise         = np.atleast_1d(3.5),
            datasec         = np.atleast_1d('[:,:]'),
            oscansec        = None #np.atleast_1d('[:,:]')
            )
        return detector_container.DetectorContainer(**detector_dict)

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.
        
        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        # Wavelengths
        # 1D wavelength solution
        par['calibrations']['wavelengths']['rms_threshold'] = 0.3
        par['calibrations']['wavelengths']['sigdetect']=5.0
        par['calibrations']['wavelengths']['fwhm']= 5.0
        par['calibrations']['wavelengths']['n_final']= [3,4,4,4,4]
        par['calibrations']['wavelengths']['lamps'] = ['OH_NIRES']
        par['calibrations']['wavelengths']['method'] = 'reidentify'

        # Reidentification parameters
        par['calibrations']['wavelengths']['reid_arxiv'] = 'p200_triplespec.fits'
#        par['calibrations']['wavelengths']['ech_fix_format'] = True
        # Echelle parameters
        par['calibrations']['wavelengths']['echelle'] = True
        par['calibrations']['wavelengths']['ech_nspec_coeff'] = 4
        par['calibrations']['wavelengths']['ech_norder_coeff'] = 6
        par['calibrations']['wavelengths']['ech_sigrej'] = 3.0

        #par['calibrations']['slitedges']['edge_thresh'] = 15.
        par['calibrations']['slitedges']['trace_thresh'] = 5.
        par['calibrations']['slitedges']['fit_min_spec_length'] = 0.3
        par['calibrations']['slitedges']['left_right_pca'] = True
        par['calibrations']['slitedges']['fwhm_gaussian'] = 4.0

        # Tilt parameters
        par['calibrations']['tilts']['tracethresh'] =  10.0

        # Processing steps
        turn_off = dict(use_illumflat=False, use_biasimage=False, use_overscan=False,
                        use_darkimage=False)
        par.reset_all_processimages_par(**turn_off)

        # Extraction
        par['reduce']['skysub']['bspline_spacing'] = 0.8
        par['reduce']['extraction']['sn_gauss'] = 4.0

        # Model entire slit
        par['reduce']['extraction']['model_full_slit'] = True  # local sky subtraction operates on entire slit
        par['reduce']['findobj']['maxnumber_sci'] = 2  # Slit is narrow so allow one object per order
        par['reduce']['findobj']['maxnumber_std'] = 1  # Slit is narrow so allow one object per order

        # Flexure
        par['flexure']['spec_method'] = 'skip'

        par['scienceframe']['process']['sigclip'] = 20.0
        par['scienceframe']['process']['satpix'] ='nothing'
        par['reduce']['extraction']['boxcar_radius'] = 0.75  # arcsec


        # Set the default exposure time ranges for the frame typing
        par['calibrations']['standardframe']['exprng'] = [None, 60]
        par['calibrations']['arcframe']['exprng'] = [100, None]
        par['calibrations']['tiltframe']['exprng'] = [100, None]
        par['calibrations']['darkframe']['exprng'] = [0, None]
        par['scienceframe']['exprng'] = [60, None]

        # Sensitivity function parameters
        par['sensfunc']['algorithm'] = 'IR'
        par['sensfunc']['polyorder'] = 8
        par['sensfunc']['IR']['telgridfile'] = 'TelFit_MaunaKea_3100_26100_R20000.fits'

        return par

    def pypeit_file_keys(self):
        """
        Define the list of keys to be output into a standard PypeIt file.

        Returns:
            :obj:`list`: The list of keywords in the relevant
            :class:`~pypeit.metadata.PypeItMetaData` instance to print to the
            :ref:`pypeit_file`.
        """
        pypeit_keys = super().pypeit_file_keys()
        # TODO: Why are these added here? See
        # pypeit.metadata.PypeItMetaData.set_pypeit_cols
        pypeit_keys += ['calib', 'comb_id', 'bkg_id']
        return pypeit_keys

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

        if ftype in ['pinhole','bias']:
            # No pinhole frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype == 'dark':
            return good_exp & (fitstbl['target'] == 'lamp_off')
        if ftype == 'standard':
            return good_exp & ((fitstbl['idname'] == 'object') | (fitstbl['idname'] == 'Object'))
        if ftype in ['pixelflat', 'trace']:
            return good_exp & (fitstbl['target'] == 'lamp_on')
        if ftype in 'science':
            return good_exp & ((fitstbl['idname'] == 'object') | (fitstbl['idname'] == 'Object'))
        if ftype in ['arc', 'tilt']:
            return good_exp & ((fitstbl['idname'] == 'object') | (fitstbl['idname'] == 'Object'))
        return np.zeros(len(fitstbl), dtype=bool)


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
                Processed bias frame used to identify bad pixels. **This is
                always ignored.**

        Returns:
            `numpy.ndarray`_: An integer array with a masked value set
            to 1 and an unmasked value set to 0.  All values are set to
            0.
        """
        # Call the base-class method to generate the empty bpm
        msgs.info("Custom bad pixel mask for TSPEC")
        return super().bpm(filename, det, shape=shape, msbias=None)

    @property
    def norders(self):
        """
        Number of orders for this spectograph. Should only defined for
        echelle spectrographs, and it is undefined for the base class.
        """
        return 5

    @property
    def order_spat_pos(self):
        """
        Return the expected spatial position of each echelle order.
        """
        return np.array([0.3096, 0.4863, 0.6406, 0.7813, 0.9424])

    @property
    def orders(self):
        """
        Return the order number for each echelle order.
        """
        return np.arange(7, 2, -1, dtype=int)

    @property
    def spec_min_max(self):
        """
        Return the minimum and maximum spectral pixel expected for the
        spectral range of each order.
        """
        spec_max = np.asarray([np.inf]*self.norders)
        spec_min = np.asarray([1024, -np.inf, -np.inf, -np.inf, -np.inf])
        return np.vstack((spec_min, spec_max))

    def order_platescale(self, order_vec, binning=None):
        """
        Return the platescale for each echelle order.

        This routine is only defined for echelle spectrographs, and it is
        undefined in the base class.

        Args:
            order_vec (`numpy.ndarray`_):
                The vector providing the order numbers.
            binning (:obj:`str`, optional):
                The string defining the spectral and spatial binning.

        Returns:
            `numpy.ndarray`_: An array with the platescale for each order
            provided by ``order``.
        """
        return np.full(order_vec.size, 0.37)


