"""
Module for VLT X-Shooter

.. include:: ../include/links.rst
"""
import numpy as np

from astropy.coordinates import SkyCoord
from astropy import units

from pypeit import msgs
from pypeit import telescopes
from pypeit import io
from pypeit.core import parse
from pypeit.core import framematch
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container
from pypeit import data

from IPython import embed

class VLTXShooterSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle VLT/XSHOOTER specific code
    """
    ndet = 1
    telescope = telescopes.VLTTelescopePar()
    pypeline = 'Echelle'
    url = 'https://www.eso.org/sci/facilities/paranal/instruments/xshooter.html'
    ech_fixed_format = True
    header_name = 'XSHOOTER'

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=0, card='RA', 
            required_ftypes=['science', 'standard'])  # Need to convert to : separated
        self.meta['dec'] = dict(ext=0, card='DEC', required_ftypes=['science', 'standard'])
        self.meta['target'] = dict(ext=0, card='OBJECT')
        self.meta['binning'] = dict(card=None, compound=True)

        self.meta['mjd'] = dict(ext=0, card='MJD-OBS')
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='HIERARCH ESO TEL AIRM START', required_ftypes=['science', 'standard'])
        # Extras for config and frametyping
        self.meta['dispname'] = dict(ext=0, card=None, default='default')
        self.meta['idname'] = dict(ext=0, card='HIERARCH ESO DPR CATG')
        self.meta['arm'] = dict(ext=0, card='HIERARCH ESO SEQ ARM')
        self.meta['instrument'] = dict(ext=0, card='INSTRUME')
        # Dithering -- Not required for redux
        self.meta['dither'] = dict(ext=0, card='HIERARCH ESO SEQ CUMOFF Y',
            required=False)  # This header card is *not* always present in science/standard frames

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
        if meta_key == 'binning':
            if 'HIERARCH ESO DET WIN1 BINX' in headarr[0]:
                binspatial = headarr[0]['HIERARCH ESO DET WIN1 BINX']
            else:
                binspatial = 1
            if 'HIERARCH ESO DET WIN1 BINY' in headarr[0]:
                binspec = headarr[0]['HIERARCH ESO DET WIN1 BINY']
            else:
                binspec = 1
            return parse.binning2string(binspec, binspatial)
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
        return ['arm']

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
        return ['HIERARCH ESO SEQ ARM']

    def pypeit_file_keys(self):
        """
        Define the list of keys to be output into a standard PypeIt file.

        Returns:
            :obj:`list`: The list of keywords in the relevant
            :class:`~pypeit.metadata.PypeItMetaData` instance to print to the
            :ref:`pypeit_file`.
        """
        return super().pypeit_file_keys() + ['dither']

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
        # TODO: Allow for 'sky' frame type, for now include sky in
        # 'science' category
        if ftype == 'science':
            return good_exp & ((fitstbl['idname'] == 'SCIENCE')
                                | (fitstbl['target'] == 'STD,TELLURIC')
                                | (fitstbl['target'] == 'STD,SKY'))
        if ftype == 'standard':
            return good_exp & (fitstbl['target'] == 'STD,FLUX')
        if ftype == 'bias':
            return good_exp & (fitstbl['target'] == 'BIAS')
        if ftype == 'dark':
            return good_exp & (fitstbl['target'] == 'DARK')
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            # Flats and trace frames are typed together
            return good_exp & ((fitstbl['target'] == 'LAMP,DFLAT')
                               | (fitstbl['target'] == 'LAMP,QFLAT')
                               | (fitstbl['target'] == 'LAMP,FLAT'))
        if ftype == 'pinhole':
            # Don't type pinhole
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['target'] == 'LAMP,WAVE')

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)


class VLTXShooterNIRSpectrograph(VLTXShooterSpectrograph):
    """
    Child to handle VLT/XSHOOTER specific code
    """

    name = 'vlt_xshooter_nir'
    camera = 'XShooter_NIR'
    supported = True
    comment = 'See :doc:`xshooter`'

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
            binning         = '1,1',  # No binning in near-IR
            det             = 1,
            dataext         = 0,
            specaxis        = 1,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.197, # average between order 11 & 30, see manual
            darkcurr        = 0.0,
            saturation      = 2.0e5, # I think saturation may never be a problem here since there are many DITs
            nonlinear       = 0.86,
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d(2.12), #
            ronoise         = np.atleast_1d(8.0), # ?? more precise value? #TODO the read noise is exposure time  dependent and should be grabbed from header
            datasec         = np.atleast_1d('[4:2044,4:]'), # These are all unbinned pixels
            # EMA: No real overscan for XSHOOTER-NIR:
            # See Table 6 in http://www.eso.org/sci/facilities/paranal/instruments/xshooter/doc/VLT-MAN-ESO-14650-4942_P103v1.pdf
            # The overscan region below contains only zeros
            # ToDo should we just set it as empty?
            #  JXP says yes
            #oscansec        = np.atleast_1d('[4:2044,1:3]'), # These are all unbinned pixels.
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

        # Turn off illumflat
        turn_off = dict(use_illumflat=False, use_biasimage=False, use_overscan=False,
                        use_darkimage=False)
        par.reset_all_processimages_par(**turn_off)
        # Require dark images to be subtracted from the flat images used for
        # tracing, pixelflats, and illumflats
        # par['calibrations']['traceframe']['process']['use_darkimage'] = True
        # par['calibrations']['pixelflatframe']['process']['use_darkimage'] = True
        # par['calibrations']['illumflatframe']['process']['use_darkimage'] = True
        # TODO: `mask_cr` now defaults to True for darks.  Should this be turned off?

        # Is this needed below?
        par['scienceframe']['process']['sigclip'] = 20.0
        par['scienceframe']['process']['satpix'] = 'nothing'
        # TODO tune up LA COSMICS parameters here for X-shooter as tellurics are being excessively masked


        # Adjustments to slit and tilts for NIR
        par['calibrations']['slitedges']['edge_thresh'] = 50.
        par['calibrations']['slitedges']['fit_order'] = 8
        par['calibrations']['slitedges']['max_shift_adj'] = 0.5
        par['calibrations']['slitedges']['trace_thresh'] = 10.
        par['calibrations']['slitedges']['fit_min_spec_length'] = 0.5
        par['calibrations']['slitedges']['left_right_pca'] = True
        par['calibrations']['slitedges']['length_range'] = 0.3

        # Tilt parameters
        par['calibrations']['tilts']['rm_continuum'] = True
        par['calibrations']['tilts']['tracethresh'] =  25.0
        par['calibrations']['tilts']['maxdev_tracefit'] =  0.04
        par['calibrations']['tilts']['maxdev2d'] =  0.04
        par['calibrations']['tilts']['spat_order'] =  3
        par['calibrations']['tilts']['spec_order'] =  4

        # 1D wavelength solution
        par['calibrations']['wavelengths']['lamps'] = ['OH_XSHOOTER']
        par['calibrations']['wavelengths']['rms_threshold'] = 0.25
        par['calibrations']['wavelengths']['sigdetect'] = 10.0
        par['calibrations']['wavelengths']['fwhm'] = 5.0
        par['calibrations']['wavelengths']['n_final'] = 4
        # Reidentification parameters
        par['calibrations']['wavelengths']['method'] = 'reidentify'
        par['calibrations']['wavelengths']['reid_arxiv'] = 'vlt_xshooter_nir.fits'
        par['calibrations']['wavelengths']['cc_thresh'] = 0.50
        par['calibrations']['wavelengths']['cc_local_thresh'] = 0.50
#        par['calibrations']['wavelengths']['ech_fix_format'] = True
        # Echelle parameters
        par['calibrations']['wavelengths']['echelle'] = True
        par['calibrations']['wavelengths']['ech_nspec_coeff'] = 5
        par['calibrations']['wavelengths']['ech_norder_coeff'] = 5
        par['calibrations']['wavelengths']['ech_sigrej'] = 3.0

        # Flats
        #par['calibrations']['standardframe']['process']['illumflatten'] = False
        par['calibrations']['flatfield']['tweak_slits_thresh'] = 0.90
        par['calibrations']['flatfield']['tweak_slits_maxfrac'] = 0.10

        # Standards
        par['calibrations']['standardframe']['process']['mask_cr'] = False

        # Extraction
        par['reduce']['skysub']['bspline_spacing'] = 0.8
        par['reduce']['skysub']['global_sky_std']  = False # Do not perform global sky subtraction for standard stars
        par['reduce']['extraction']['model_full_slit'] = True  # local sky subtraction operates on entire slit
        par['reduce']['findobj']['trace_npoly'] = 8
        par['reduce']['findobj']['maxnumber_sci'] = 2  # Assume that there is only one object on the slit.
        par['reduce']['findobj']['maxnumber_std'] = 1  # Assume that there is only one object on the slit.


        # The settings below enable X-shooter dark subtraction from the traceframe and pixelflatframe, but enforce
        # that this bias won't be subtracted from other images. It is a hack for now, because eventually we want to
        # perform this operation with the dark frame class, and we want to attach individual sets of darks to specific
        # images.
        #par['calibrations']['biasframe']['useframe'] = 'bias'
        #par['calibrations']['traceframe']['process']['bias'] = 'force'
        #par['calibrations']['pixelflatframe']['process']['bias'] = 'force'
        #par['calibrations']['arcframe']['process']['bias'] = 'skip'
        #par['calibrations']['tiltframe']['process']['bias'] = 'skip'
        #par['calibrations']['standardframe']['process']['bias'] = 'skip'
        #par['scienceframe']['process']['bias'] = 'skip'

        # Sensitivity function parameters
        par['sensfunc']['algorithm'] = 'IR'
        par['sensfunc']['polyorder'] = 8
        par['sensfunc']['IR']['telgridfile'] = 'TelFit_Paranal_NIR_9800_25000_R25000.fits'

        return par


    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        super().init_meta()
        # No binning in the NIR
        self.meta['binning'] = dict(card=None, default='1,1')

        # Required
        self.meta['decker'] = dict(ext=0, card='HIERARCH ESO INS OPTI5 NAME')
    
        # Dark-flat identification via exposure number
        self.meta['seq_expno'] = dict(ext=0, card='HIERARCH ESO TPL EXPNO')

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
        pypeit_keys += ['comb_id', 'bkg_id']
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

        # Default NIR calibration behavior is to take flat/darks in sequence
        #  These are marked by the seq_expno column
        good_flat_seq = np.array([seq is not None and int(seq) % 2 == 1 for seq in fitstbl['seq_expno']])
        good_dark_seq = np.array([seq is not None and int(seq) % 2 == 0 for seq in fitstbl['seq_expno']])

        # TODO: Allow for 'sky' frame type, for now include sky in
        # 'science' category
        if ftype == 'science':
            return good_exp & ((fitstbl['idname'] == 'SCIENCE')
                                | (fitstbl['target'] == 'STD,TELLURIC')
                                | (fitstbl['target'] == 'STD,SKY'))
        if ftype == 'standard':
            return good_exp & (fitstbl['target'] == 'STD,FLUX')
        if ftype == 'bias':
            return good_exp & (fitstbl['target'] == 'BIAS')
        if ftype == 'sky':
            return good_exp & (fitstbl['target'] == 'DARK')
        
        if ftype in ['pixelflat', 'trace']:
            # Flats and trace frames are typed together
            # Lamp on flats are taken first (odd exposure number)
            return good_exp & (((fitstbl['target'] == 'LAMP,DFLAT')
                               | (fitstbl['target'] == 'LAMP,QFLAT')
                               | (fitstbl['target'] == 'LAMP,FLAT'))
                               & good_flat_seq)
        
        if ftype in ['dark']:
            # Lamp off flats are taken second (even exposure number)
            return good_exp & (((fitstbl['target'] == 'LAMP,DFLAT')
                                | (fitstbl['target'] == 'LAMP,QFLAT')
                                | (fitstbl['target'] == 'LAMP,FLAT'))
                               & good_dark_seq)
        
        if ftype == 'pinhole':
            # Don't type pinhole
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc', 'tilt']:
            return good_exp & ((fitstbl['target'] == 'LAMP,WAVE') | (fitstbl['target'] == 'SCIENCE'))

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
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
                Processed bias frame used to identify bad pixels

        Returns:
            `numpy.ndarray`_: An integer array with a masked value set
            to 1 and an unmasked value set to 0.  All values are set to
            0.
        """
        # Call the base-class method to generate the empty bpm
        bpm_img = super().bpm(filename, det, shape=shape, msbias=msbias)

        if det == 1:
            bpm_dir = data.Paths.static_calibs / 'vlt_xshoooter'
            try :
                bpm_loc = np.loadtxt(bpm_dir / 'BP_MAP_RP_NIR.dat', usecols=(0,1))
            except IOError :
                msgs.warn('BP_MAP_RP_NIR.dat not present in the static database')
                bpm_fits = io.fits_open(bpm_dir / 'BP_MAP_RP_NIR.fits.gz')
                # ToDo: this depends on datasec, biassec, specflip, and specaxis
                #       and should become able to adapt to these parameters.
                # Flipping and shifting BPM to match the PypeIt format
                y_shift = -2
                x_shift = 18
                bpm_data = np.flipud(bpm_fits[0].data)
                y_len = len(bpm_data[:,0])
                x_len = len(bpm_data[0,:])
                bpm_data_pypeit = np.full( ((y_len+abs(y_shift)),(x_len+abs(x_shift))) , 0)
                bpm_data_pypeit[:-abs(y_shift),:-abs(x_shift)] = bpm_data_pypeit[:-abs(y_shift),:-abs(x_shift)] + bpm_data
                bpm_data_pypeit = np.roll(bpm_data_pypeit,-y_shift,axis=0)
                bpm_data_pypeit = np.roll(bpm_data_pypeit,x_shift,axis=1)
                filt_bpm = bpm_data_pypeit[1:y_len,1:x_len]>100.
                y_bpm, x_bpm = np.where(filt_bpm)
                bpm_loc = np.array([y_bpm,x_bpm]).T
                np.savetxt(bpm_dir / 'BP_MAP_RP_NIR.dat', bpm_loc, fmt=['%d','%d'])
            finally :
                bpm_img[bpm_loc[:,0].astype(int),bpm_loc[:,1].astype(int)] = 1.

        return bpm_img

    @property
    def norders(self):
        """
        Number of orders for this spectograph. Should only defined for
        echelle spectrographs, and it is undefined for the base class.
        """
        return 16

    @property
    def order_spat_pos(self):
        """
        Return the expected spatial position of each echelle order.
        """
        return np.array([0.08284662, 0.1483813 , 0.21158701, 0.27261607,
                         0.33141317, 0.38813936, 0.44310197, 0.49637422,
                         0.54839496, 0.59948157, 0.65005956, 0.70074477,
                         0.75240745, 0.80622583, 0.86391259, 0.9280528 ])

    @property
    def orders(self):
        """
        Return the order number for each echelle order.
        """
        return np.arange(26, 10, -1, dtype=int)

    @property
    def spec_min_max(self):
        """
        Return the minimum and maximum spectral pixel expected for the
        spectral range of each order.
        """
        spec_max = np.asarray([1467,1502,1540, 1580,1620,1665,1720, 1770,1825,1895, 1966, 2000,2000,2000,2000,2000])
        spec_min = np.asarray([420 ,390 , 370,  345, 315, 285, 248,  210, 165, 115,   63,   10,   0,   0,   0,   0])
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
        # TODO: Either assume a linear trend or measure this
        # X-shooter manual says, but gives no exact numbers per order.
        # NIR: 52.4 pixels (0.210"/pix) at order 11 to 59.9 pixels (0.184"/pix) at order 26.

        # Right now I just assume a simple linear trend
        plate_scale = 0.184 + (order_vec - 26)*(0.184-0.210)/(26 - 11)
        return plate_scale

    @property
    def dloglam(self):
        """
        Return the logarithmic step in wavelength for output spectra.
        """
        # This number was computed by taking the mean of the dloglam for all
        # the X-shooter orders. The specific loglam across the orders deviates
        # from this value by +-6% from this first to final order
        return 1.93724e-5

    @property
    def loglam_minmax(self):
        """
        Return the base-10 logarithm of the first and last wavelength for
        ouput spectra.
        """
        return np.log10(9500.0), np.log10(26000)


class VLTXShooterVISSpectrograph(VLTXShooterSpectrograph):
    """
    Child to handle VLT/XSHOOTER specific code
    """

    name = 'vlt_xshooter_vis'
    camera = 'XShooter_VIS'
    supported = True
    comment = 'See :doc:`xshooter`'

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
        # Binning
        # TODO: Could this be detector dependent??
        binning = '1,1' if hdu is None else self.get_meta_value(self.get_headarr(hdu), 'binning')

        # Detector 1
        detector_dict = dict(
            binning         = binning,
            det              =1,
            dataext         = 0,
            specaxis        = 0,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.16, # average from order 17 and order 30, see manual
            darkcurr        = 0.0,
            saturation      = 65535.,
            nonlinear       = 0.86,
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d(0.595), # FITS format is flipped: PrimaryHDU  (2106, 4000) w/respect to Python
            ronoise         = np.atleast_1d(3.1), # raw unbinned images are (4000,2106) (spec, spat)
            datasec=np.atleast_1d('[:,11:2058]'),  # pre and oscan are in the spatial direction
            oscansec=np.atleast_1d('[:,2059:2106]'),
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

        # Adjustments to parameters for VIS
        turn_on = dict(use_biasimage=False, use_overscan=True, overscan_method='median',
                       use_darkimage=False, use_illumflat=False, use_pixelflat=False,
                       use_specillum=False)
        par.reset_all_processimages_par(**turn_on)

        # X-SHOOTER arcs/tilts are also have different binning with bias
        # frames, so don't use bias frames. Don't use the biases for any
        # calibrations since it appears to be a different amplifier readout
        par['calibrations']['traceframe']['process']['overscan_method'] = 'median'

        par['scienceframe']['process']['use_biasimage']=True
        par['scienceframe']['process']['use_illumflat']=True
        par['scienceframe']['process']['use_pixelflat']=True
        # Right now we are using the overscan and not biases becuase the
        # standards are read with a different read mode and we don't yet have
        # the option to use different sets of biases for different standards,
        # or use the overscan for standards but not for science frames
        par['calibrations']['standardframe']['process']['use_illumflat']=True
        par['calibrations']['standardframe']['process']['use_pixelflat']=True

        par['calibrations']['slitedges']['edge_thresh'] = 8.0
        par['calibrations']['slitedges']['fit_order'] = 8
        par['calibrations']['slitedges']['max_shift_adj'] = 0.5
        par['calibrations']['slitedges']['trace_thresh'] = 10.
        par['calibrations']['slitedges']['left_right_pca'] = True
        par['calibrations']['slitedges']['length_range'] = 0.3

        # These are the defaults
        par['calibrations']['tilts']['tracethresh'] = 15
        par['calibrations']['tilts']['spat_order'] =  3
        par['calibrations']['tilts']['spec_order'] =  5 

        # 1D wavelength solution
        par['calibrations']['wavelengths']['lamps'] = ['ThAr_XSHOOTER_VIS']
        # The following is for 1x1 binning. TODO GET BINNING SORTED OUT!!
        par['calibrations']['wavelengths']['rms_threshold'] = 0.50
        par['calibrations']['wavelengths']['sigdetect'] = 5.0
        par['calibrations']['wavelengths']['n_final'] = [3] + 13*[4] + [3]
        # This is for 1x1 binning. Needs to be divided by binning for binned data!!
        par['calibrations']['wavelengths']['fwhm'] = 11.0
        # Reidentification parameters
        par['calibrations']['wavelengths']['method'] = 'reidentify'
        # TODO: the arxived solution is for 1x1 binning. It needs to be
        # generalized for different binning!
        par['calibrations']['wavelengths']['reid_arxiv'] = 'vlt_xshooter_vis1x1.fits'
        par['calibrations']['wavelengths']['cc_thresh'] = 0.50
        par['calibrations']['wavelengths']['cc_local_thresh'] = 0.50
#        par['calibrations']['wavelengths']['ech_fix_format'] = True
        # Echelle parameters
        par['calibrations']['wavelengths']['echelle'] = True
        par['calibrations']['wavelengths']['ech_nspec_coeff'] = 4
        par['calibrations']['wavelengths']['ech_norder_coeff'] = 4
        par['calibrations']['wavelengths']['ech_sigrej'] = 3.0

        # Flats
        par['calibrations']['flatfield']['tweak_slits_thresh'] = 0.90
        par['calibrations']['flatfield']['tweak_slits_maxfrac'] = 0.10

        # Extraction
        par['reduce']['skysub']['bspline_spacing'] = 0.5
        par['reduce']['skysub']['global_sky_std'] = False
        # local sky subtraction operates on entire slit
        par['reduce']['extraction']['model_full_slit'] = True
        # Mask 3 edges pixels since the slit is short, insted of default (5,5)
        par['reduce']['findobj']['find_trim_edge'] = [3,3]
        par['reduce']['findobj']['maxnumber_sci'] = 2  # Assume that there is only one object on the slit.
        par['reduce']['findobj']['maxnumber_std'] = 1  # Assume that there is only one object on the slit.
        # Continnum order for determining thresholds

        # Sensitivity function parameters
        par['sensfunc']['algorithm'] = 'IR'
        par['sensfunc']['polyorder'] = [9, 11, 11, 9, 9, 8, 8, 7, 7, 7, 7, 7, 7, 7, 7]
        par['sensfunc']['IR']['telgridfile'] = 'TelFit_Paranal_VIS_4900_11100_R25000.fits'
        return par

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        super().init_meta()
        # Add the name of the dispersing element
        # dispangle and filter1 are not defined for Shane Kast Blue

        # Required
        self.meta['decker'] = dict(ext=0, card='HIERARCH ESO INS OPTI4 NAME')
    

    @property
    def norders(self):
        """
        Number of orders observed for this spectograph.
        """
        return 15

    @property
    def order_spat_pos(self):
        """
        Return the expected spatial position of each echelle order.
        """
        return np.array([0.13540436, 0.21055672, 0.2817009, 0.34907542, 0.41289127, 0.4733839 ,
                         0.53072208, 0.58509916, 0.63671413, 0.685754, 0.73236772, 0.77676367,
                         0.8191196 , 0.85968302, 0.89877932])

    @property
    def orders(self):
        """
        Return the order number for each echelle order.
        """
        return np.arange(30, 15, -1, dtype=int)

    @property
    def spec_min_max(self):
        """
        Return the minimum and maximum spectral pixel expected for the
        spectral range of each order.
        """
        spec_max = np.asarray([4000]*14 + [3000])
        spec_min = np.asarray([2000,1000] + [0]*13)
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
        # VIS has no binning, but for an instrument with binning we would do this
        binspectral, binspatial = parse.parse_binning(binning)

        # ToDO Either assume a linear trend or measure this
        # X-shooter manual says, but gives no exact numbers per order.
        # VIS: 65.9 pixels (0.167"/pix) at order 17 to 72.0 pixels (0.153"/pix) at order 30.

        # Right now I just assume a simple linear trend
        plate_scale = 0.153 + (order_vec - 30)*(0.153-0.167)/(30 - 17)
        return plate_scale*binspatial

    @property
    def dloglam(self):
        """
        Return the logarithmic step in wavelength for output spectra.
        """
        # This number was computed by taking the mean of the dloglam for all
        # the X-shooter orders. The specific loglam across the orders deviates
        # from this value by +-7% from this first to final order. This is the
        # unbinned value. It was actually measured to be 1.69207e-5 from a 2x1
        # data and then divided by two.
        return 8.46035e-06

    @property
    def loglam_minmax(self):
        """
        Return the base-10 logarithm of the first and last wavelength for
        ouput spectra.
        """
        return np.log10(5000.0), np.log10(11000)

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
        # Call the base-class method to generate the empty bpm
        bpm_img = super().bpm(filename, det, shape=shape, msbias=msbias)

        shape = bpm_img.shape
        #
        # TODO: Ema: This is just a workaround to deal with different binning.
        # I guess binspatial and binspectral should be passed in.
        if shape[0]<3000.:
            binspectral_bpm=2
        else:
            binspectral_bpm=1
        if shape[1]<1500.:
            binspatial_bpm=2
        else:
            binspatial_bpm=1

        if det == 1:
            bpm_img[2912//binspectral_bpm:,842//binspatial_bpm:844//binspatial_bpm] = 1.
            bpm_img[3548//binspectral_bpm:,1249//binspatial_bpm:1252//binspatial_bpm] = 1.
        return bpm_img


class VLTXShooterUVBSpectrograph(VLTXShooterSpectrograph):
    """
    Child to handle VLT/XSHOOTER specific code for the UVB arm
    """

    name = 'vlt_xshooter_uvb'
    camera = 'XShooter_UVB'
    supported = True
    comment = 'See :doc:`xshooter`'
    
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
        # Binning
        binning = '1,1' if hdu is None else self.get_meta_value(self.get_headarr(hdu), 'binning')

        # Detector 1
        detector_dict = dict(
            binning         = binning,
            det             = 1,
            dataext         = 0,
            specaxis        = 0,
            specflip        = True,
            spatflip        = True,
            platescale      = 0.161, # average from order 14 and order 24, see manual
            darkcurr        = 0.0,
            saturation      = 65000.,
            nonlinear       = 0.86,  
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d(1.61),
            ronoise         = np.atleast_1d(2.60),
            datasec         = np.atleast_1d('[:,49:2096]'), # '[49:2000,1:2999]',
            oscansec        = np.atleast_1d('[:,1:48]'), # '[1:48, 1:2999]',
            )
        # Return
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

        # Adjustments to parameters for UVB (following VIS)
        turn_on = dict(use_biasimage=False, use_overscan=True, overscan_method='median',
                       use_darkimage=False, use_illumflat=False, use_pixelflat=False,
                       use_specillum=False)
        par.reset_all_processimages_par(**turn_on)

        # X-SHOOTER arcs/tilts are also have different binning with bias
        # frames, so don't use bias frames. Don't use the biases for any
        # calibrations since it appears to be a different amplifier readout

        # Adjustments to slit and tilts for UVB
        par['calibrations']['slitedges']['edge_thresh'] = 8.
        par['calibrations']['slitedges']['max_shift_adj'] = 0.5
        par['calibrations']['slitedges']['trace_thresh'] = 10.
        par['calibrations']['slitedges']['left_right_pca'] = True
        par['calibrations']['slitedges']['length_range'] = 0.3

        # 1D wavelength solution
        par['calibrations']['wavelengths']['lamps'] = ['ThAr_XSHOOTER_UVB']
        par['calibrations']['wavelengths']['n_final'] = [3] + 10*[4] 
        par['calibrations']['wavelengths']['rms_threshold'] = 0.60 
        par['calibrations']['wavelengths']['sigdetect'] = 3.0 # Pretty faint lines in places
        # Reidentification parameters
        par['calibrations']['wavelengths']['method'] = 'reidentify'
        par['calibrations']['wavelengths']['reid_arxiv'] = 'vlt_xshooter_uvb1x1.fits'
#        par['calibrations']['wavelengths']['ech_fix_format'] = True
        # Echelle parameters
        par['calibrations']['wavelengths']['echelle'] = True
        par['calibrations']['wavelengths']['ech_nspec_coeff'] = 4
        par['calibrations']['wavelengths']['ech_norder_coeff'] = 4
        par['calibrations']['wavelengths']['ech_sigrej'] = 3.0
        
        par['calibrations']['wavelengths']['cc_thresh'] = 0.50
        par['calibrations']['wavelengths']['cc_local_thresh'] = 0.50

        # Right now we are using the overscan and not biases becuase the
        # standards are read with a different read mode and we don't yet have
        # the option to use different sets of biases for different standards,
        # or use the overscan for standards but not for science frames
        par['scienceframe']['process']['use_biasimage']=True
        par['scienceframe']['process']['use_illumflat']=True
        par['scienceframe']['process']['use_pixelflat']=True
        par['calibrations']['standardframe']['process']['use_illumflat']=True
        par['calibrations']['standardframe']['process']['use_pixelflat']=True


        # Extraction
        par['reduce']['skysub']['bspline_spacing'] = 0.5
        par['reduce']['skysub']['global_sky_std'] = False
        par['reduce']['extraction']['model_full_slit'] = True
        # Mask 3 edges pixels since the slit is short, insted of default (5,5)
        par['reduce']['findobj']['find_trim_edge'] = [3,3]
        # Continnum order for determining thresholds
        #par['reduce']['findobj']['find_npoly_cont'] = 0
        # Don't attempt to fit a continuum to the trace rectified image
        #par['reduce']['findobj']['find_cont_fit'] = False
        par['reduce']['findobj']['maxnumber_sci'] = 2  # Assume that there is a max of 2 objects on the slit
        par['reduce']['findobj']['maxnumber_std'] = 1  # Assume that there is only one object on the slit.

        return par

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        super().init_meta()
        # Add the name of the dispersing element

        # Required
        self.meta['decker'] = dict(ext=0, card='HIERARCH ESO INS OPTI3 NAME')

    @property
    def norders(self):
        """
        Number of orders observed for this spectograph.
        """
        return 11

    @property
    def order_spat_pos(self):
        """
        Return the expected spatial position of each echelle order.

        The following lines generated the values below:

        .. code-block:: python

            from pypeit import edgetrace
            edges = edgetrace.EdgeTraceSet.from_file('Edges_A_1_DET01.fits.gz')

            nrm_edges = edges.edge_fit[edges.nspec//2,:] / edges.nspat
            slit_cen = ((nrm_edges + np.roll(nrm_edges,1))/2)[np.arange(nrm_edges.size//2)*2+1]

        """
        # This starts by ignoring the first, partial order (25?)
        #  Order 24 is very faint and not included here
        #  Order 12 is very faint and not included here (the flat crashes out with issues..)
        return np.array([0.32671887, 0.39553878, 0.45989826, 0.52009878, 0.5764598,
            0.62917188, 0.67859507, 0.72482729, 0.76815531, 0.80879042,
            0.84700373])#, 0.88317493])

    @property
    def orders(self):
        """
        Return the order number for each echelle order.
        """
        return np.arange(23, 12, -1, dtype=int)  # 11 orders; the reddest is too faint to use
        #return np.arange(23, 11, -1, dtype=int)   # 12 orders

    @property
    def spec_min_max(self):
        """
        Return the minimum and maximum spectral pixel expected for the
        spectral range of each order.
        """
        spec_max = np.asarray([4000]*13)# + [3000])
        spec_min = np.asarray([0]*13)
        return np.vstack((spec_min, spec_max))

    def order_platescale(self, order_vec, binning = None):
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
        binspectral, binspatial = parse.parse_binning(binning)

        # ToDO Either assume a linear trend or measure this
        # X-shooter manual says, but gives no exact numbers per order.
        # UVB: 65.9 pixels (0.167“/pix) at order 14 to 70.8 pixels (0.155”/pix) at order 24

        # Assume a simple linear trend
        plate_scale = 0.155 + (order_vec - 24)*(0.155-0.167)/(24 - 14)

        # Right now I just took the average
        return np.full(self.norders, 0.161)*binspatial

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
        # Call the base-class method to generate the empty bpm
        bpm_img = super().bpm(filename, det, shape=shape, msbias=msbias)

        # TODO -- Mask bad column if it is problematic (it isn't so far)

        return bpm_img



