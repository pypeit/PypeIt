"""
Module for ARC/ARCES

.. include:: ../include/links.rst
"""
import numpy as np
from astropy.time import Time

from pypeit import log
from pypeit import PypeItError
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.spectrographs import spectrograph
from pypeit.core import parse
from pypeit.images import detector_container


class ARCARCESSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle ARC KOSMOS instrument+detector
    """
    ndet = 1
    name = 'arc_arces'
    telescope = telescopes.ARCTelescopePar()
    camera = 'ARCARCES'
    url = 'https://www.apo.nmsu.edu/arc35m/Instruments/ARCES'
    header_name = 'ARCES'
    supported = False
    comment = 'ARC ARCES spectrometer'
    pypeline = 'Echelle'
    ech_fixed_format = True


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
            binning         = '1,1' if hdu is None 
                                    else self.get_meta_value(self.get_headarr(hdu), 'binning'),
            det=1,
            dataext         = 0,
            specaxis        = 1,
            specflip        = True,
            spatflip        = True,
            xgap            = 0.,
            ygap            = 0.,
            ysize           = 1.,
            platescale      = 0.52,
            mincounts       = -1e10,
            darkcurr        = 0.0,  # e-/pixel/hour
            saturation      = 65535.,
            nonlinear       = 0.86,
            numamplifiers   = 1,
            gain            = np.atleast_1d([3.8]),
            ronoise         = np.atleast_1d([7.0]),
            datasec         = np.atleast_1d(['[1:2048,23:2048]']),
            oscansec        = np.atleast_1d(['[1:2048,2075:2125]']),
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


        # Ignore PCA
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'

        # Set pixel flat combination method
        par['calibrations']['pixelflatframe']['process']['combine'] = 'mean'
        # Wavelength calibration methods
        par['calibrations']['wavelengths']['method'] = 'reidentify'
        par['calibrations']['wavelengths']['lamps'] = ['ThAr']
        par['calibrations']['wavelengths']['reid_arxiv'] = 'arc_arces.fits'
        par['calibrations']['wavelengths']['sigdetect'] = 10.0
        par['calibrations']['wavelengths']['echelle'] = True

        # allow for multiple wavecals of different lamps and/or exptimes
        par['calibrations']['arcframe']['process']['clip'] = False
        par['calibrations']['tiltframe']['process']['clip'] = False

        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, None]
        par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['arcframe']['exprng'] = [None, None]  # Long arc exposures on this telescope
        par['calibrations']['standardframe']['exprng'] = [None, 120]
        par['scienceframe']['exprng'] = [1, None]

        #Debora from here

        # edge tracing
        par['calibrations']['traceframe']['process']['scale_to_mean'] = True
        par['calibrations']['slitedges']['edge_thresh'] = 15.
        par['calibrations']['slitedges']['smash_range'] = [0.3,0.7]
        par['calibrations']['slitedges']['order_match'] = 0.005
        par['calibrations']['slitedges']['fwhm_gaussian'] = 1.5
        par['calibrations']['slitedges']['fwhm_uniform'] = 1.5
        par['calibrations']['slitedges']['pad'] = 5

        # Wavelength
        # 1D wavelength solution
        par['calibrations']['wavelengths']['n_final'] = 4
        par['calibrations']['wavelengths']['cc_thresh'] = 0.4
        par['calibrations']['wavelengths']['sigdetect'] = 3.
        par['calibrations']['wavelengths']['fwhm'] = 3.5
        par['calibrations']['wavelengths']['fwhm_fromlines'] = False
        par['calibrations']['wavelengths']['match_toler'] = 2.
        par['calibrations']['wavelengths']['rms_thresh_frac_fwhm'] = 0.5
        par['calibrations']['wavelengths']['bad_orders_maxfrac'] = 0.5
        
        # Echelle parameters
        par['calibrations']['wavelengths']['ech_nspec_coeff'] = 4
        par['calibrations']['wavelengths']['ech_norder_coeff'] = 6
        par['calibrations']['wavelengths']['ech_sigrej'] = 3.0

        # wave tilts calibration
        par['calibrations']['tilts']['tracethresh'] = 10.

        # flat fileding
        # Set pixel flat combination method
        par['calibrations']['pixelflatframe']['process']['scale_to_mean'] = True
        par['calibrations']['illumflatframe']['process']['scale_to_mean'] = True
        par['calibrations']['pixelflatframe']['process']['combine'] = 'mean'
        par['calibrations']['flatfield']['slit_illum_finecorr'] = False
        par['calibrations']['flatfield']['tweak_slits'] = False
        par['calibrations']['flatfield']['spat_samp'] = 1  # default is 5
        par['calibrations']['flatfield']['slit_trim'] = 1

        # no sky subtraction on standard stars
        par['reduce']['skysub']['global_sky_std'] = False
        # skip sky subtraction when searching for objects
        par['reduce']['findobj']['skip_skysub'] = True
        # no local sky subtraction
        par['reduce']['skysub']['no_local_sky'] = True
        par['reduce']['skysub']['mask_by_boxcar'] = True

        # find objects
        par['reduce']['findobj']['find_trim_edge'] = [0, 0]
        # extraction
        par['reduce']['extraction']['boxcar_radius'] = 1.56
        par['reduce']['extraction']['model_full_slit'] = True
        par['reduce']['extraction']['sn_gauss'] = 4000 # basically always use the Gaussian model for optimal extraction

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

        par['calibrations']['wavelengths']['reid_arxiv'] = "arc_arces.fits"
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
        self.meta['target'] = dict(ext=0, card='OBJNAME')
        self.meta['decker'] = dict(ext=0, card='INSTRUME')
        self.meta['binning'] = dict(card=None, compound=True)

        self.meta['mjd'] = dict(card=None,compound=True)
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Extras for config and frametyping
        self.meta['dispname'] = dict(ext=0, card='INSTRUME')
        self.meta['idname'] = dict(ext=0, card='IMAGETYP')
        # Lamps
        self.meta['mirror'] = dict(ext=0, card='MIRROR')
        self.meta['lampstat01'] = dict(ext=0, card='LAMPW')
        self.meta['lampstat02'] = dict(ext=0, card='LAMPT')
        self.meta['instrument'] = dict(ext=0, card='INSTRUME')
        # Mirror
        #self.meta['mirror'] = dict(card=None)

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
            binspatial = headarr[0]['CCDBIN1']
            binspec = headarr[0]['CCDBIN2']
            return parse.binning2string(binspec, binspatial)
        elif meta_key == 'mjd':
            time = headarr[0]['DATE-OBS']
            ttime = Time(time, format='isot')
            return ttime.mjd
        else:
            raise PypeItError("Not ready for this compound meta")

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
        return ['instrument']

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
        return ['CCDBIN1', 'CCDBIN2']

    def pypeit_file_keys(self):
        """
        Define the list of keys to be output into a standard PypeIt file.

        Returns:
            :obj:`list`: The list of keywords in the relevant
            :class:`~pypeit.metadata.PypeItMetaData` instance to print to the
            :ref:`pypeit_file`.
        """
        return super().pypeit_file_keys() 

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
        if ftype in ['science', 'standard']:
            return good_exp & (fitstbl['idname'] == 'object')
        if ftype == 'bias':
            return good_exp & (fitstbl['idname'] == 'zero')
        if ftype == 'pixelflat': #Internal Flats
            return good_exp & (fitstbl['idname'] == 'flat') 
        if ftype in ['trace', 'illumflat']: 
            return good_exp & (fitstbl['idname'] == 'flat') 
        if ftype in ['pinhole', 'dark']:
            # Don't type pinhole or dark frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc','tilt']:
            return (good_exp & 
                    (fitstbl['mirror'] == 'Lamps') &
                    (fitstbl['lampstat02'] == '1') )
        log.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

    @property
    def norders(self):
        """
        Number of orders for this spectograph. Should only defined for
        echelle spectrographs, and it is undefined for the base class.
        """
        return 105

    @property
    def order_spat_pos(self):
        """
        Return the expected spatial position of each echelle order.
        """
        return np.array(
         [0.23185659, 0.24071052, 0.24951271, 0.25820542, 0.26682471, 0.27534693,
          0.28383979, 0.29222108, 0.30054427, 0.3087838 , 0.31695649, 0.32504069,
          0.33307599, 0.34101296, 0.34890994, 0.35671584, 0.36445721, 0.37212708,
          0.37972129, 0.38724958, 0.3947124 , 0.40210802, 0.40943847, 0.41669141,
          0.42387731, 0.43102231, 0.438077  , 0.44507518, 0.45202159, 0.45889284,
          0.46570195, 0.47244716, 0.47913928, 0.48576406, 0.49233616, 0.4988397 ,
          0.50528551, 0.51168016, 0.51800638, 0.52427919, 0.53049915, 0.53665685,
          0.5427665 , 0.5488149 , 0.55481995, 0.56076293, 0.56665304, 0.5724885 ,
          0.57827453, 0.58401064, 0.58969965, 0.59533909, 0.60092929, 0.60646368,
          0.61194805, 0.61738147, 0.62276776, 0.62810237, 0.63339203, 0.638642  ,
          0.64384236, 0.64900529, 0.65412626, 0.659205  , 0.66423787, 0.66922865,
          0.67418058, 0.6790963 , 0.68397457, 0.68880943, 0.69361243, 0.69838367,
          0.70311866, 0.70781911, 0.71248932, 0.71713185, 0.72174557, 0.72632912,
          0.73088776, 0.73542203, 0.73993791, 0.74442825, 0.7488972 , 0.75334996,
          0.7577895 , 0.76221721, 0.76662915, 0.77103302, 0.77543287, 0.7798291 ,
          0.78422839, 0.7886239 , 0.79302601, 0.79743542, 0.80185765, 0.80629519,
          0.81075452, 0.81523994, 0.81975688, 0.82430991, 0.82890433, 0.83354863,
          0.83825048, 0.84301529, 0.84784117])  #, 0.85275266, 0.85776007])

    @property
    def orders(self):
        """
        Return the order number for each echelle order.
        """
        return np.arange(160, 160-self.norders, -1, dtype=int)

    @property
    def spec_min_max(self):
        """
        Return the minimum and maximum spectral pixel expected for the
        spectral range of each order.
        """
        spec_max = np.zeros(self.norders) + 1800
        #spec_max[50:108] = 1500
        #spec_max[90:108] = 1350
        spec_min = np.zeros(self.norders) + 200
        #spec_max[50:108] = 750
        #spec_max[90:108] = 850
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
        # TODO: Figure out the order-dependence of the updated plate scale
        # From the X-Shooter P113 manual, average over all orders. No order-dependent values given.
        plate_scale = 0.245*np.ones_like(order_vec)
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
