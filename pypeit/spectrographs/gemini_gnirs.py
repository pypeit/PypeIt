"""
Module for Gemini/GNIRS specific methods.

.. include:: ../include/links.rst
"""
import numpy as np

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.images import detector_container
from pypeit.spectrographs import spectrograph


class GeminiGNIRSSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Gemini/GNIRS specific code
    """
    ndet = 1
    name = 'gemini_gnirs'
    camera = 'GNIRS'
    url = 'https://www.gemini.edu/instrumentation/gnirs'
    header_name = 'GNIRS'
    telescope = telescopes.GeminiNTelescopePar()
    pypeline = 'Echelle'
    ech_fixed_format = True
    supported = True

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
            darkcurr        = 0.15,
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

        # Reduce parameters
        #par['reduce']['findobj']['snr_thresh'] = 5.0          # Object finding threshold
        par['reduce']['findobj']['find_trim_edge'] = [2,2]    # Slit is too short to trim 5,5 especially
        par['reduce']['skysub']['bspline_spacing'] = 0.8
        par['reduce']['skysub']['global_sky_std']  = False    # Do not perform global sky subtraction for standard stars
        # TODO: JFH: Is this the correct behavior?  (Is why we have sky-subtraction problems for GNIRS?)
        par['reduce']['skysub']['no_poly'] = True             # Do not use polynomial degree of freedom for global skysub
        par['reduce']['extraction']['model_full_slit'] = True  # local sky subtraction operates on entire slit
        par['reduce']['findobj']['maxnumber_sci'] = 2  # Slit is narrow so allow one object per order
        par['reduce']['findobj']['maxnumber_std'] = 1  # Slit is narrow so allow one object per order
        # Standards
        par['calibrations']['standardframe']['process']['mask_cr'] = False # Do not mask_cr standards

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
        par['sensfunc']['IR']['telgridfile'] = 'TelFit_MaunaKea_3100_26100_R20000.fits'
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
            par['calibrations']['slitedges']['left_right_pca'] = True
            par['calibrations']['slitedges']['pca_order'] = 3

            # Wavelengths
            par['calibrations']['wavelengths']['rms_threshold'] = 1.0  # Might be grating dependent..
            par['calibrations']['wavelengths']['sigdetect'] =  [4.0, 5.0, 5.0, 5.0, 5.0, 5.0] #5.0
            par['calibrations']['wavelengths']['lamps'] = ['OH_GNIRS']
            #par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
            par['calibrations']['wavelengths']['n_first'] = 2
            par['calibrations']['wavelengths']['n_final'] = [1, 3, 3, 3, 3, 3]

            # Reidentification parameters
            par['calibrations']['wavelengths']['method'] = 'reidentify'
            par['calibrations']['wavelengths']['cc_thresh'] = 0.6
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gemini_gnirs.fits'
#            par['calibrations']['wavelengths']['ech_fix_format'] = True
            # Echelle parameters
            # JFH This is provisional these IDs should be checked.
            par['calibrations']['wavelengths']['echelle'] = True
            par['calibrations']['wavelengths']['ech_nspec_coeff'] = 3
            par['calibrations']['wavelengths']['ech_norder_coeff'] = 5
            par['calibrations']['wavelengths']['ech_sigrej'] = 3.0

            # Tilts
            par['calibrations']['tilts']['tracethresh'] = [5.0, 10, 10, 10, 10, 10]
            par['calibrations']['tilts']['sig_neigh'] = 5.0
            par['calibrations']['tilts']['nfwhm_neigh'] = 2.0
        # 10/mmLBSX_G5532 setup, covering YJHK with the long blue camera and SXD prism
        elif '10/mmLBSX' in self.dispname:
            # Edges
            par['calibrations']['slitedges']['edge_thresh'] = 20.
            par['calibrations']['slitedges']['trace_thresh'] = 10.
            par['calibrations']['slitedges']['fit_order'] = 2
            par['calibrations']['slitedges']['max_shift_adj'] = 0.5
            par['calibrations']['slitedges']['det_min_spec_length'] = 0.20
            par['calibrations']['slitedges']['fit_min_spec_length'] = 0.20
            par['calibrations']['slitedges']['left_right_pca'] = True # Actually we need a parameter to disable PCA entirely
            par['calibrations']['slitedges']['sync_predict'] = 'nearest'

            # Wavelengths
            par['calibrations']['wavelengths']['rms_threshold'] = 1.0  # Might be grating dependent..
            par['calibrations']['wavelengths']['sigdetect'] = 5.0
            par['calibrations']['wavelengths']['lamps'] = ['Ar_IR_GNIRS']
            #par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
            par['calibrations']['wavelengths']['n_first'] = 2
            par['calibrations']['wavelengths']['n_final'] = [3, 3, 3, 3]
            # Reidentification parameters
            par['calibrations']['wavelengths']['method'] = 'reidentify'
            par['calibrations']['wavelengths']['cc_thresh'] = 0.6
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gemini_gnirs_10mm_LBSX.fits'
#            par['calibrations']['wavelengths']['ech_fix_format'] = True
            # Echelle parameters
            par['calibrations']['wavelengths']['echelle'] = True
            par['calibrations']['wavelengths']['ech_nspec_coeff'] = 3
            par['calibrations']['wavelengths']['ech_norder_coeff'] = 3
            par['calibrations']['wavelengths']['ech_sigrej'] = 3.0

            # Tilts
            par['calibrations']['tilts']['tracethresh'] = [10, 10, 10, 10]
            par['calibrations']['tilts']['sig_neigh'] = 5.0
            par['calibrations']['tilts']['nfwhm_neigh'] = 2.0
        else:
            msgs.error('Unrecognized GNIRS dispname')

        return par

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
        Define the list of keys to be output into a standard PypeIt file.

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



