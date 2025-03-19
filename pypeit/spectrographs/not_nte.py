"""
Module for NOT NTE

Started by copying VLT/XSHOOTER, so watch out for inherited stuff
Also copying quite a lot from the not_alfosc.py file.
In particular can remove ESO type things

.. include:: ../include/links.rst
"""
import numpy as np

from astropy.coordinates import SkyCoord
from astropy import units
from astropy.time import Time

from pypeit import msgs
from pypeit import telescopes
from pypeit import io
from pypeit.core import parse
from pypeit.core import framematch
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container
from pypeit.core import parse
from pypeit import data

from IPython import embed

class NOTNTESpectrograph(spectrograph.Spectrograph):
    """
    Child to handle NOT/NTE specific code
    Note that, like for Xshooter, we have 3 seperate classes for the arms
    """
    ndet = 1 # number of detectors
    name = 'not_nte'
    telescope = telescopes.NOTTelescopePar() # already have NOTTelescopePar from ALFOSC dev
    pypeline = 'Echelle'
    url = 'https://nte.nbi.ku.dk/'
    ech_fixed_format = True
    header_name = 'NTE'

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the ``PypeIt``-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)
        # Have used the options from not_alfosc, should be consistent
        # Dithering stuff from keck_nires
        self.meta['ra'] = dict(ext=0, card="RA")
        self.meta['dec'] = dict(ext=0, card='DEC')
        self.meta['target'] = dict(ext=0, card='OBJECT')
        self.meta['binning'] = dict(card=None, compound=True) # uses NOT style keys, see compound_meta
        self.meta['mjd'] = dict(ext=0, card=None, compound=True) # uses NOT style keys, see compound_meta
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        self.meta['decker'] = dict(ext=0, card='SLIT') # This is maybe not the best choice, but the decker header is needed

        # Extras for config and frametyping
        self.meta['dispname'] = dict(ext=0, card='MODES')
        self.meta['idname'] = dict(ext=0, card='IMAGETYP')
        self.meta['arm'] = dict(ext=0, card="ARM")
        self.meta['instrument'] = dict(ext=0, card='INSTRUME')

        # Dithering
        # Need to edit this for NIR
##        self.meta['dithpat'] = dict(ext=0, card='DPATNAME')
##        self.meta['dithpos'] = dict(card=None, compound=True)
##        self.meta['dithoff'] = dict(ext=0, card='YOFFSET')


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
        # This comes from not_alfosc
        if meta_key == 'binning':
            # PypeIt frame
            binspatial = headarr[0]['DETXBIN']
            binspec = headarr[0]['DETYBIN']
            return parse.binning2string(binspec, binspatial)
        elif meta_key == 'mjd':
            time = headarr[0]['DATE-AVG']
            ttime = Time(time, format='isot')
            return ttime.mjd
        #elif meta_key == 'ra':
        #    objra = headarr[0]['OBJRA'] # Given in hours, not deg
        #    return objra*15.
        else :
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
        # Makes sense to use arm here
        # later should have more?
        return ['arm',"decker"]

    def pypeit_file_keys(self):
        """
        Define the list of keys to be output into a standard ``PypeIt`` file.

        Returns:
            :obj:`list`: The list of keywords in the relevant
            :class:`~pypeit.metadata.PypeItMetaData` instance to print to the
            :ref:`pypeit_file`.
        """
        # Not sure about this one
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

        # Copying not_alfosc, might need to add more
        # Also using the header decisions made in Lise's python script
        # etalon? NIR? 
        # Might need to include telluric and/or sky in the science classified ftype

        if ftype == 'science':
            return good_exp & (fitstbl['idname'] == 'OBJECT') #| (fitstbl['target'] == 'STD,TELLURIC')  | (fitstbl['target'] == 'STD,SKY'))
        if ftype == 'standard':
            return good_exp & ((fitstbl['idname'] == 'STD') | (fitstbl['target'] == 'STD') | (fitstbl['target'] == 'STD,SLIT'))
        if ftype == 'bias':
            return good_exp & (fitstbl['idname'] == 'BIAS')
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            return good_exp & ((fitstbl['idname'] == 'LAMP,FLAT') | (fitstbl['target'] == 'LAMP,TRACE'))
        if ftype == 'dark':
            return good_exp & (fitstbl['idname'] == 'DARK')
        if ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['idname'] == 'LAMP,WAVE')
        if ftype == 'pinhole':
            return good_exp & (fitstbl['idname'] == 'LAMP,TRACE')
        if ftype == "sky":
            return good_exp & (fitstbl['idname'] == 'SKY')

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

class NOTNTEVISSpectrograph(NOTNTESpectrograph):
    """
    Child to handle NOT/NTE specific code
    """

    name = 'not_nte_vis'
    camera = 'nte_vis'
    supported = False
    comment = 'See :doc:`nte`'

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

        # Detector 1 (this all needs checking for a final version)
        detector_dict = dict(
            binning         = binning,
            det              =1,
            dataext         = 0,
            specaxis        = 1, #check this, opposite of xshooter?
            specflip        = False, # check this
            spatflip        = False, # check this
            platescale      = 0.23, # taken from NTE_NOT 2022 presentation slides
            darkcurr        = 0.0,
            saturation      = 65535., # check
            nonlinear       = 0.86, # check
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d(0.595), # Get this value
            ronoise         = np.atleast_1d(3.0), # Get this value
            datasec=np.atleast_1d('[{}:{},:]'.format(1,1024)),  # Just trying something here
            oscansec= None ,  # Overscan not actually required
        )

        return detector_container.DetectorContainer(**detector_dict)

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.

        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of ``PypeIt`` methods.
        """
        par = super().default_pypeit_par()

        # Adjustments to parameters for VIS
        #turn_on = dict(use_biasimage=True, use_overscan=False, overscan_method='median',
        #               use_darkimage=False, use_illumflat=False, use_pixelflat=False,
        #               use_specillum=False)
        turn_off = dict(use_overscan=False)
        par.reset_all_processimages_par(**turn_off)

        # The below is sufficient for OK edge tracing, probably not a necessary set
        par['calibrations']['slitedges']['edge_thresh'] = 5.0
        par['calibrations']['slitedges']['fit_order'] = 5
        par['calibrations']['slitedges']['max_shift_adj'] = 0.5
        par['calibrations']['slitedges']['trace_thresh'] = 10
        par['calibrations']['slitedges']['fit_min_spec_length'] = 0.1
        par['calibrations']['slitedges']['length_range'] = 0.3

        par['calibrations']['slitedges']['det_buffer'] = 1
        par['calibrations']['slitedges']['max_nudge'] = 1
        #par['calibrations']['slitedges']['left_right_pca'] = False
        #par['calibrations']['slitedges']['add_slits'] = "1:2280:35:124"
        #par['calibrations']['slitedges']['sync_predict'] = "nearest"
        par['calibrations']['slitedges']['smash_range'] = [0.3,0.7]
        #par['calibrations']['slitedges']['sobel_mode'] = "constant"


        # Start on wl calib
        par['calibrations']['wavelengths']['lamps'] = ["HgAr_NTE_VIS"]
        par['calibrations']['wavelengths']['rms_threshold'] = 0.4
        par['calibrations']['wavelengths']['sigdetect'] = 2.0
        par['calibrations']['wavelengths']['fwhm'] = 4.0
        par['calibrations']['wavelengths']['n_final'] = 4# [2, 4, 4, 4, 4, 4, 4, 4]
        par['calibrations']['wavelengths']['nreid_min'] = 1 # important
        
        par['calibrations']['wavelengths']['reference'] = 'arc'
        par['calibrations']['wavelengths']['reid_arxiv'] = 'not_nte_vis.fits'
        par['calibrations']['wavelengths']['method'] = 'full_template'
        par['calibrations']['wavelengths']['nsnippet'] = 1 # important

        # Echelle parameters
        par['calibrations']['wavelengths']['echelle'] = True
        par['calibrations']['wavelengths']['ech_nspec_coeff'] = 5
        par['calibrations']['wavelengths']['ech_norder_coeff'] = 5
        par['calibrations']['wavelengths']['ech_sigrej'] = 3.0

        # tilts
        #par['calibrations']['tilts']['spat_order'] =  3
        
        # Flat
        par['calibrations']['flatfield']['slit_illum_finecorr'] = False # turn off for now

        # skysub
        par['reduce']['skysub']['bspline_spacing'] = 1

        # extraction
        par['reduce']['findobj']['maxnumber_sci'] = 1
        par['reduce']['findobj']['maxnumber_std'] = 1


        # Sensitivity function parameters
        par['sensfunc']['algorithm'] = 'IR'
        #par['sensfunc']['polyorder'] = [9, 11, 11, 9, 9, 8, 8, 7, 7, 7, 7, 7, 7, 7, 7]
        #par['sensfunc']['IR']['telgridfile'] = 'TelFit_Paranal_VIS_4900_11100_R25000.fits'

        return par

    @property
    def norders(self):
        """
        Number of orders observed for this spectograph.
        """
        return 8

    @property
    def order_spat_pos(self):
        """
        Return the expected spatial position of each echelle order.
        """

        return np.array([0.17285156, 0.31542969, 0.43554688, 0.54980469,
                         0.64941406,0.75195312, 0.84960938, 0.9375])


        #np.array([73, 207, 332, 449, 558, 660 ,757,854]) were the positions used
        #np.array([177, 323, 446, 563, 665, 770 ,870,960]) are the new positions used

        # normalised by the detector height

    @property
    def orders(self):
        """
        Return the order number for each echelle order.
        """
        return np.arange(15, 7, -1, dtype=int) # orders 15-8, from NTE_NOT_2022 slides

    @property
    def spec_min_max(self):
        """
        Return the minimum and maximum spectral pixel expected for the
        spectral range of each order.
        """
        spec_min = np.asarray([0]*8)
        spec_max = np.asarray([2000,3000,4000,4000,4000,4000,4000,4000])
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

        # ToDO Work this out
        
        # Right now I just assume constant
        plate_scale = np.ones(8) * 0.23 
        return plate_scale*binspatial

        # Not sure about this, commenting out
##    @property
##    def dloglam(self):
##        """
##        Return the logarithmic step in wavelength for output spectra.
##        """
##        # This number was computed by taking the mean of the dloglam for all
##        # the X-shooter orders. The specific loglam across the orders deviates
##        # from this value by +-7% from this first to final order. This is the
##        # unbinned value. It was actually measured to be 1.69207e-5 from a 2x1
##        # data and then divided by two.
##        return 8.46035e-06

    @property
    def loglam_minmax(self):
        """
        Return the base-10 logarithm of the first and last wavelength for
        ouput spectra.
        """
        return np.log10(4200), np.log10(9150)

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
                Master bias frame used to identify bad pixels

        Returns:
            `numpy.ndarray`_: An integer array with a masked value set
            to 1 and an unmasked value set to 0.  All values are set to
            0.
        """
        # Call the base-class method to generate the empty bpm

        # THE BELOW MASKS OUT THE +/- 1 ORDERS, NOT A REAL BPM AND MAY NOT BE IN THE FINAL CODE
        
        bpm_img = super().bpm(filename, det, shape=shape, msbias=msbias)
        bpm_dir = data.Paths.static_calibs / 'not_nte'
        bpm_loc = np.loadtxt(bpm_dir / 'mask_VIS.dat', usecols=(0,1))
        
        for i in range(0,4096):
            bpm_img[i] = np.append(np.ones(int(bpm_loc[i][0])),np.zeros(int(1024-bpm_loc[i][0])))

        
        return bpm_img


class NOTNTENIRSpectrograph(NOTNTESpectrograph):
    """
    Child to handle NTE/NIR specific code
    """

    name = 'not_nte_nir'
    camera = 'nte_nir'
    supported = False
    comment = 'See :doc:`nte`'

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
            platescale      = 0.46, # taken from NTE_NOT 2022 presentation slides, requires checking
            darkcurr        = 0.0, # CHECK
            saturation      = 2.0e5, # CHECK, although saturation should never be a problem in IR if you are observing properly
            nonlinear       = 0.86, # CHECK,
            mincounts       = -1e10, 
            numamplifiers   = 1,
            gain            = np.atleast_1d(2.16), # e/ADU, from the MPIA test report
            ronoise         = np.atleast_1d(10.6), # e-, from the MPIA test report
            datasec=np.atleast_1d('[{}:{},:]'.format(1,2048)),  # Just trying something here
            #datasec         = np.atleast_1d('[4:2044,4:]'), # These are all unbinned pixels
            oscansec= None , # Should check if we want an overscan
            )
        return detector_container.DetectorContainer(**detector_dict)

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.

        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of ``PypeIt`` methods.
        """

        par = super().default_pypeit_par()
        # Turn off illumflat, bias, oversscan and dark (?)
        turn_off = dict(use_illumflat=False, use_biasimage=False, use_overscan=False,
                        use_darkimage=False)
        par.reset_all_processimages_par(**turn_off)


    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the ``PypeIt``-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        super().init_meta()
        # No binning in the NIR (true also for nte?)
        self.meta['binning'] = dict(card=None, default='1,1')

        # Required
        #self.meta['decker'] = dict(ext=0, card='HIERARCH ESO INS OPTI5 NAME')

        # Dark-flat identification via exposure number
        #self.meta['seq_expno'] = dict(ext=0, card='HIERARCH ESO TPL EXPNO')

    def pypeit_file_keys(self):
        """
        Define the list of keys to be output into a standard ``PypeIt`` file.

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


    # Use this if we need different frame type rules here than for VIS
##    def check_frame_type(self, ftype, fitstbl, exprng=None):
##        """
##        Check for frames of the provided type.
##
##        Args:
##            ftype (:obj:`str`):
##                Type of frame to check. Must be a valid frame type; see
##                frame-type :ref:`frame_type_defs`.
##            fitstbl (`astropy.table.Table`_):
##                The table with the metadata for one or more frames to check.
##            exprng (:obj:`list`, optional):
##                Range in the allowed exposure time for a frame of type
##                ``ftype``. See
##                :func:`pypeit.core.framematch.check_frame_exptime`.
##
##        Returns:
##            `numpy.ndarray`_: Boolean array with the flags selecting the
##            exposures in ``fitstbl`` that are ``ftype`` type frames.
##        """
##        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
##
##        # Default NIR calibration behavior is to take flat/darks in sequence
##        #  These are marked by the seq_expno column
##        good_flat_seq = np.array([seq is not None and int(seq) % 2 == 1 for seq in fitstbl['seq_expno']])
##        good_dark_seq = np.array([seq is not None and int(seq) % 2 == 0 for seq in fitstbl['seq_expno']])
##
##        # TODO: Allow for 'sky' frame type, for now include sky in
##        # 'science' category
##        if ftype == 'science':
##            return good_exp & ((fitstbl['idname'] == 'SCIENCE')
##                                | (fitstbl['target'] == 'STD,TELLURIC')
##                                | (fitstbl['target'] == 'STD,SKY'))
##        if ftype == 'standard':
##            return good_exp & (fitstbl['target'] == 'STD,FLUX')
##        if ftype == 'bias':
##            return good_exp & (fitstbl['target'] == 'BIAS')
##        if ftype == 'sky':
##            return good_exp & (fitstbl['target'] == 'DARK')
##
##        if ftype in ['pixelflat', 'trace']:
##            # Flats and trace frames are typed together
##            # Lamp on flats are taken first (odd exposure number)
##            return good_exp & (((fitstbl['target'] == 'LAMP,DFLAT')
##                               | (fitstbl['target'] == 'LAMP,QFLAT')
##                               | (fitstbl['target'] == 'LAMP,FLAT'))
##                               & good_flat_seq)
##
##        if ftype in ['dark']:
##            # Lamp off flats are taken second (even exposure number)
##            return good_exp & (((fitstbl['target'] == 'LAMP,DFLAT')
##                                | (fitstbl['target'] == 'LAMP,QFLAT')
##                                | (fitstbl['target'] == 'LAMP,FLAT'))
##                               & good_dark_seq)
##
##        if ftype == 'pinhole':
##            # Don't type pinhole
##            return np.zeros(len(fitstbl), dtype=bool)
##        if ftype in ['arc', 'tilt']:
##            return good_exp & ((fitstbl['target'] == 'LAMP,WAVE') | (fitstbl['target'] == 'SCIENCE'))
##
##        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
##        return np.zeros(len(fitstbl), dtype=bool)

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
                Master bias frame used to identify bad pixels

        Returns:
            `numpy.ndarray`_: An integer array with a masked value set
            to 1 and an unmasked value set to 0.  All values are set to
            0.
        """
        # Call the base-class method to generate the empty bpm
        bpm_img = super().bpm(filename, det, shape=shape, msbias=msbias)
        # Should return an empty bpm
        return bpm_img



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
        return np.array([0.13720703, 0.12792969, 0.22949219, 0.30761719, 0.45898438])
    
        # using 281,262,470,630,940 / 2048

    @property
    def orders(self):
        """
        Return the order number for each echelle order.
        """
        return np.arange(7, 2, -1, dtype=int) # orders 7-3, from NTE_NOT_2022 slides

    @property
    def spec_min_max(self):
        """
        Return the minimum and maximum spectral pixel expected for the
        spectral range of each order.
        """
        spec_max = np.asarray([2048]*5)
        spec_min = np.asarray([0]*5)
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
        # Should work this out properly

        # Right now I just assume constant, this is not correct
        plate_scale = np.ones(5) * 0.46
        return plate_scale


    @property
    def loglam_minmax(self):
        """
        Return the base-10 logarithm of the first and last wavelength for
        ouput spectra.
        """
        return np.log10(8500.0), np.log10(25000)


class NOTNTEUVSpectrograph(NOTNTESpectrograph):
    """
    Child to handle NOT/NTE specific code
    """

    name = 'not_nte_uv'
    camera = 'nte_uv'
    supported = False
    comment = 'See :doc:`nte`'

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

        # Detector 1 (this all needs checking for a final version)
        detector_dict = dict(
            binning         = binning,
            det              =1,
            dataext         = 0,
            specaxis        = 1, #check this, opposite of xshooter?
            specflip        = False, # check this
            spatflip        = False, # check this
            platescale      = 0.23, # taken from NTE_NOT 2022 presentation slides
            darkcurr        = 0.0,
            saturation      = 65535., # check
            nonlinear       = 0.86, # check
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d(0.595), # Get this value
            ronoise         = np.atleast_1d(3.0), # Get this value
            datasec=np.atleast_1d('[{}:{},:]'.format(1,1024)),  # Just trying something here
            oscansec= None ,  # Overscan not actually required
        )

        return detector_container.DetectorContainer(**detector_dict)

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.

        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of ``PypeIt`` methods.
        """
        par = super().default_pypeit_par()

        # Adjustments to parameters for VIS
        #turn_on = dict(use_biasimage=True, use_overscan=False, overscan_method='median',
        #               use_darkimage=False, use_illumflat=False, use_pixelflat=False,
        #               use_specillum=False)
        turn_off = dict(use_overscan=False , use_illumflat=False , use_darkimage=False)
        par.reset_all_processimages_par(**turn_off)


        par['calibrations']['standardframe']['process']['mask_cr'] = False
        par['scienceframe']['process']['mask_cr'] = False

        # The below is sufficient for OK edge tracing, probably not a necessary set
        par['calibrations']['slitedges']['edge_thresh'] = 2.0
        par['calibrations']['slitedges']['fit_order'] = 5
        par['calibrations']['slitedges']['max_shift_adj'] = 0.8
        par['calibrations']['slitedges']['trace_thresh'] = 3.0
        par['calibrations']['slitedges']['fit_min_spec_length'] = 0.1
        par['calibrations']['slitedges']['length_range'] = 0.3

        par['calibrations']['slitedges']['det_buffer'] = 1
        par['calibrations']['slitedges']['max_nudge'] = 1
        #par['calibrations']['slitedges']['left_right_pca'] = False
        #par['calibrations']['slitedges']['add_slits'] = "1:2280:35:124"
        #par['calibrations']['slitedges']['sync_predict'] = "nearest"
        par['calibrations']['slitedges']['smash_range'] = [0.3,0.7]
        #par['calibrations']['slitedges']['sobel_mode'] = "constant"


        # Start on wl calib
        par['calibrations']['wavelengths']['lamps'] = ["HgAr_NTE_UV"]
        #par['calibrations']['wavelengths']['rms_threshold'] = 0.8 
        par['calibrations']['wavelengths']['sigdetect'] = 2.0
        par['calibrations']['wavelengths']['fwhm'] = 4.0
        par['calibrations']['wavelengths']['n_final'] =  6
        par['calibrations']['wavelengths']['n_first'] =  2
        par['calibrations']['wavelengths']['nreid_min'] = 1 # important
        par['calibrations']['wavelengths']['reference'] = 'arc'
        par['calibrations']['wavelengths']['reid_arxiv'] = 'not_nte_uv_5.fits'
        par['calibrations']['wavelengths']['method'] = 'full_template'
        par['calibrations']['wavelengths']['nsnippet'] = 1 # important

        # Echelle parameters
        par['calibrations']['wavelengths']['echelle'] = True
        par['calibrations']['wavelengths']['ech_nspec_coeff'] = 5
        par['calibrations']['wavelengths']['ech_norder_coeff'] = 5
        par['calibrations']['wavelengths']['ech_sigrej'] = 3.0

        # tilts
        #par['calibrations']['tilts']['spat_order'] =  2
        
        # Flat
        par['calibrations']['flatfield']['slit_illum_finecorr'] = False # turn off for now

        # skysub
        par['reduce']['skysub']['user_regions'] = ':30,70:'
        #par['reduce']['skysub']['global_sky_std'] = False
        #par['reduce']['skysub']['no_local_sky'] = True

        # extraction
        par['reduce']['findobj']['maxnumber_sci'] = 1
        par['reduce']['findobj']['maxnumber_std'] = 1
        par['reduce']['extraction']['use_2dmodel_mask'] = False


        # Sensitivity function parameters
        par['sensfunc']['algorithm'] = 'UVIS'

        return par

    @property
    def norders(self):
        """
        Number of orders observed for this spectograph.
        """
        return 4

    @property
    def order_spat_pos(self):
        """
        Return the expected spatial position of each echelle order.
        """

        return np.array(#[0.14648438, 0.29296875, 
            [0.44238312, 0.5878906,
                         0.72558594, 0.852539])


        #np.array([150, 300, 453, 601, 743, 873]) were positions used for UV
        #np.array([177, 323, 446, 563, 665, 770 ,870,960])are positions used for VIS

        # normalised by the detector height

    @property
    def orders(self):
        """
        Return the order number for each echelle order.
        """
        return np.arange(19, 15, -1, dtype=int) # orders 15-8, from NTE_NOT_2022 slides

    @property
    def spec_min_max(self):
        """
        Return the minimum and maximum spectral pixel expected for the
        spectral range of each order.
        """
        spec_min = np.asarray([0]*4)
        spec_max = np.asarray([4400,4400,4400,4400])
        return np.vstack((spec_min, spec_max))
    
#Really don't know what this section will do

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

        # ToDO Work this out
        
        # Right now I just assume constant
        plate_scale = np.ones(4) * 0.23 
        return plate_scale*binspatial

        # Not sure about this, commenting out
##    @property
##    def dloglam(self):
##        """
##        Return the logarithmic step in wavelength for output spectra.
##        """
##        # This number was computed by taking the mean of the dloglam for all
##        # the X-shooter orders. The specific loglam across the orders deviates
##        # from this value by +-7% from this first to final order. This is the
##        # unbinned value. It was actually measured to be 1.69207e-5 from a 2x1
##        # data and then divided by two.
##        return 8.46035e-06

    @property
    def loglam_minmax(self):
        """
        Return the base-10 logarithm of the first and last wavelength for
        ouput spectra.
        """
        return np.log10(3300), np.log10(4300)

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
                Master bias frame used to identify bad pixels

        Returns:
            `numpy.ndarray`_: An integer array with a masked value set
            to 1 and an unmasked value set to 0.  All values are set to
            0.
        """
        # Call the base-class method to generate the empty bpm

        # THE BELOW MASKS OUT THE +/- 1 ORDERS, NOT A REAL BPM AND MAY NOT BE IN THE FINAL CODE
        
        bpm_img = super().bpm(filename, det, shape=shape, msbias=msbias)
        # Should return an empty bpm
        return bpm_img