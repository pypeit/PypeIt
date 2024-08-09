# Code to make the NIRSPEC DRP
"""
Module for Keck/NIRSPEC specific methods.

.. include:: ../include/links.rst
"""
import numpy as np

from pypeit import msgs
from pypeit import io
from pypeit.images import detector_container
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.core import parse
from pypeit.spectrographs import spectrograph


class KeckNIRSPECSpectrograph(spectrograph.Spectrograph):
    """
    Parent to handle the new Keck/NIRSPEC detector (i.e. post late 2018) specific code
    """
    ndet = 1
    telescope = telescopes.KeckTelescopePar()
    camera = 'NIRSPEC'
    comment = 'POST upgrade (~Dec 2018)'
    url = 'https://www2.keck.hawaii.edu/inst/nirspec/'
    header_name = 'NIRSPEC'

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
        pypeit_keys += ['comb_id', 'bkg_id', 'shift']
        return pypeit_keys


class KeckNIRSPECSpectrographOld(spectrograph.Spectrograph):
    """
    Parent to handle the old Keck/NIRSPEC detector (i.e. pre late 2018) specific code
    """
    ndet = 1
    telescope = telescopes.KeckTelescopePar()
    camera = 'NIRSPEC'
    comment = 'PRE upgrade (~Dec 2018)'
    url = 'https://www2.keck.hawaii.edu/inst/nirspec/'
    header_name = 'NIRSPEC'

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
        pypeit_keys += ['comb_id', 'bkg_id', 'shift']
        return pypeit_keys


class KeckNIRSPECHighSpectrograph(KeckNIRSPECSpectrograph):
    """
    Child to handle NIRSPEC high-dispersion specific code
    """
    name = 'keck_nirspec_high'
    supported = True
    comment = 'POST upgrade (~Dec 2018), High-dispersion grating, see :ref:`nirspec_high_howto`'
    pypeline = 'Echelle'
    ech_fixed_format = False
    lamps_list = []
    filter1 = ''
    filter2 = ''

    def get_detector_par(self, det, hdu=None):
        """
        Return metadata for the selected detector.
        see https://www2.keck.hawaii.edu/inst/nirspec/Specifications.html
        and https://www2.keck.hawaii.edu/inst/nirspec/images/spec_highres.png

        Args:
            det (:obj:`int`):
                1-indexed detector number.  This is not used because NIRSPEC
                only has one detector!
            hdu (`astropy.io.fits.HDUList`_, optional):
                The open fits file with the raw image of interest.  If not
                provided, frame-dependent parameters are set to a default.

        Returns:
            :class:`~pypeit.images.detector_container.DetectorContainer`:
            Object with the detector metadata.
        """

        detector_dict = dict(
            det=1,
            binning         ='1,1',  # No binning allowed
            dataext         = 0,
            specaxis        = 0,
            specflip        = True,
            spatflip        = False,
            platescale      = 0.13,
            darkcurr        = 2520.0,  # e-/pixel/hour  (=0.7 e-/pixel/s)
            saturation      = 1.0e5,
            nonlinear       = 0.9,  # docs say linear to 90,000 but our flats are usually higher
            numamplifiers   = 1,
            mincounts       = -1e10,
            gain            = np.atleast_1d(3.01),
            ronoise         = np.atleast_1d(11.56),
            datasec         = np.atleast_1d('[:,:]'),
            oscansec        = None, #np.atleast_1d('[:,:]')
            )

        return detector_container.DetectorContainer(**detector_dict)

    def order_platescale(self, order_vec, binning=None):
        """
        Return the platescale for each echelle order.

        Note that NIRSPEC has no binning.

        Args:
            order_vec (`numpy.ndarray`_):
                The vector providing the order numbers.
            binning (:obj:`str`, optional):
                The string defining the spectral and spatial binning. **This
                is always ignored.**

        Returns:
            `numpy.ndarray`_: An array with the platescale for each order
            provided by ``order``.
        """
        pscale = self.get_detector_par(1).platescale
        return np.full(order_vec.size, pscale)

    @classmethod
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
        return ['filter1', 'filter2', 'echangle', 'xdangle', 'decker']

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
        return ['SCIFILT1', 'SCIFILT2', 'ECHLPOS', 'DISPPOS', 'SLITNAME']

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
        par['calibrations']['wavelengths']['rms_thresh_frac_fwhm'] = 3.0 #0.20  # Might be grating dependent..
        par['calibrations']['wavelengths']['sigdetect']=5.0
        par['calibrations']['wavelengths']['fwhm']= 3.0
        par['calibrations']['wavelengths']['fwhm_fromlines']= False
        par['calibrations']['wavelengths']['n_final']= 4
        par['calibrations']['wavelengths']['lamps'] = ['NIRSPEC-ArNeKrXe']#, 'ThAr']
        #par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['method'] = 'echelle'
#        par['calibrations']['wavelengths']['ech_fix_format'] = True
        par['calibrations']['wavelengths']['cc_thresh'] = 0.5
        par['calibrations']['wavelengths']['cc_local_thresh'] = 0.50
        # Echelle parameters
        par['calibrations']['wavelengths']['echelle'] = True
        par['calibrations']['wavelengths']['ech_nspec_coeff'] = 4
        par['calibrations']['wavelengths']['ech_norder_coeff'] = 4
        par['calibrations']['wavelengths']['ech_sigrej'] = 3.0
        par['calibrations']['wavelengths']['ech_separate_2d'] = False

        
        # Trace ID parameters
        par['calibrations']['slitedges']['edge_thresh'] = 100.0
        par['calibrations']['slitedges']['fit_order'] = 4
        par['calibrations']['slitedges']['max_shift_adj'] = 0.5
        par['calibrations']['slitedges']['trace_thresh'] = 10.
        par['calibrations']['slitedges']['left_right_pca'] = True
        par['calibrations']['slitedges']['length_range'] = 0.3
        par['calibrations']['slitedges']['max_nudge'] = 10.
        par['calibrations']['slitedges']['overlap'] = True
        par['calibrations']['slitedges']['dlength_range'] = 0.25

        # Flats
        par['calibrations']['flatfield']['tweak_slits_thresh'] = 0.80
        par['calibrations']['flatfield']['slit_illum_relative'] = False#True

        # Extraction
        par['reduce']['skysub']['bspline_spacing'] = 0.8
        par['reduce']['extraction']['sn_gauss'] = 4.0
        par['reduce']['findobj']['fof_link'] = 0.2

        # Flexure
        par['flexure']['spec_method'] = 'skip'

        par['scienceframe']['process']['sigclip'] = 20.0
        par['scienceframe']['process']['satpix'] ='nothing'

        # Should be we be illumflattening?

        # Flats
        turn_off = dict(use_biasimage=False, use_overscan=False,
                        use_darkimage=False, use_specillum=False, use_illumflat = False, use_pixelflat=False) #use_illumflat=True, 
        par.reset_all_processimages_par(**turn_off)

        # Specify if cleaning cosmic ray hits/bad pixels

        # The settings below enable NIRSPEC dark subtraction from the
        # traceframe and pixelflatframe, but enforce that this bias won't be
        # subtracted from other images. It is a hack for now, because
        # eventually we want to perform this operation with the dark frame
        # class, and we want to attach individual sets of darks to specific
        # images.
        #par['calibrations']['biasframe']['useframe'] = 'bias'
        #par['calibrations']['traceframe']['process']['bias'] = 'force'
        #par['calibrations']['pixelflatframe']['process']['bias'] = 'force'
        #par['calibrations']['arcframe']['process']['bias'] = 'skip'
        #par['calibrations']['tiltframe']['process']['bias'] = 'skip'
        #par['calibrations']['standardframe']['process']['bias'] = 'skip'
        #par['scienceframe']['process']['bias'] = 'skip'

        # Set the default exposure time ranges for the frame typing
        par['calibrations']['standardframe']['exprng'] = [None, None]
        par['calibrations']['arcframe']['exprng'] = [20, None]
        par['calibrations']['darkframe']['exprng'] = [None, None]
        par['scienceframe']['exprng'] = [None, None]

        # Sensitivity function parameters
        par['sensfunc']['algorithm'] = 'IR'
        par['sensfunc']['polyorder'] = 8
        par['sensfunc']['IR']['telgridfile'] = 'TellPCA_9300_55100_R60000.fits'

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

        headarr = self.get_headarr(scifile)
        self.filter1 = self.get_meta_value(headarr, 'filter1')
        self.filter2 = self.get_meta_value(headarr, 'filter2')
        decker = self.get_meta_value(headarr, 'decker')
        self.lamps_list = par['calibrations']['wavelengths']['lamps']
        
        # wavelength calibration
        supported_filters = ['NIRSPEC-1', 'NIRSPEC-3', 'NIRSPEC-5', 'NIRSPEC-7', 'Kband-new', 'KL']
        if (self.filter1 not in supported_filters) and (self.filter2 not in supported_filters):
            msgs.warn(f'Filter {self.filter1} or {self.filter2} may not be supported!!')
        
        if self.filter1 == 'Kband-new' or self.filter2 == 'NIRSPEC-7':
            par['calibrations']['wavelengths']['n_final'] = 3
            par['calibrations']['wavelengths']['ech_nspec_coeff'] = 2
            par['calibrations']['wavelengths']['cc_thresh'] = 0.5
            par['calibrations']['wavelengths']['cc_local_thresh'] = 0.5
            par['calibrations']['wavelengths']['cc_offset_minmax'] = 0.25
            par['calibrations']['wavelengths']['cc_percent_ceil'] = 99.9
            par['calibrations']['wavelengths']['echelle_pad'] = 1
            par['calibrations']['wavelengths']['stretch_func'] = 'quadratic'

            par['calibrations']['slitedges']['overlap'] = False
            
        if self.filter1 == 'KL' or self.filter2 == 'KL':
            par['calibrations']['wavelengths']['n_final'] = 2
            par['calibrations']['wavelengths']['ech_nspec_coeff'] = 2
            par['calibrations']['wavelengths']['cc_thresh'] = 0.5
            par['calibrations']['wavelengths']['cc_local_thresh'] = 0.5
            par['calibrations']['wavelengths']['cc_offset_minmax'] = 0.25
            par['calibrations']['wavelengths']['cc_percent_ceil'] = 99.9
            par['calibrations']['wavelengths']['echelle_pad'] = 1
            par['calibrations']['wavelengths']['stretch_func'] = 'quadratic'

            par['calibrations']['slitedges']['overlap'] = False

        if self.filter2 == 'NIRSPEC-5':
            par['calibrations']['wavelengths']['n_final'] = 3
            par['calibrations']['wavelengths']['ech_nspec_coeff'] = 3
            par['calibrations']['wavelengths']['cc_thresh'] = 0.5
            par['calibrations']['wavelengths']['cc_local_thresh'] = 0.5
            par['calibrations']['wavelengths']['cc_offset_minmax'] = 0.25
            par['calibrations']['wavelengths']['cc_percent_ceil'] = 70.0
            par['calibrations']['wavelengths']['echelle_pad'] = 1
            par['calibrations']['wavelengths']['stretch_func'] = 'quadratic'
            if self.get_meta_value(headarr, 'xdangle') == 36.72:
                par['calibrations']['slitedges']['rm_slits'] = '1:1100:1925'

        if self.filter2 == 'NIRSPEC-3':
            par['calibrations']['wavelengths']['n_final'] = 3
            par['calibrations']['wavelengths']['ech_nspec_coeff'] = 3
            par['calibrations']['wavelengths']['cc_thresh'] = 0.5
            par['calibrations']['wavelengths']['cc_local_thresh'] = 0.5
            par['calibrations']['wavelengths']['cc_offset_minmax'] = 0.25
            par['calibrations']['wavelengths']['cc_percent_ceil'] = 99.9
            par['calibrations']['wavelengths']['echelle_pad'] = 0
            par['calibrations']['wavelengths']['stretch_func'] = 'quadratic'

        if self.filter2 == 'NIRSPEC-1':
            par['calibrations']['wavelengths']['n_final'] = 3
            par['calibrations']['wavelengths']['ech_nspec_coeff'] = 3
            par['calibrations']['wavelengths']['cc_thresh'] = 0.5
            par['calibrations']['wavelengths']['cc_local_thresh'] = 0.5
            par['calibrations']['wavelengths']['cc_offset_minmax'] = 0.25
            par['calibrations']['wavelengths']['cc_percent_ceil'] = 99.9
            par['calibrations']['wavelengths']['echelle_pad'] = 1
            par['calibrations']['wavelengths']['stretch_func'] = 'quadratic'
            #par['calibrations']['slitedges']['rm_slits'] = '1:'


        if decker == '0.144x12':
            par['calibrations']['wavelengths']['fwhm'] = 3.0
        if decker == '0.288x12' or decker == '0.288x24':
            par['calibrations']['wavelengths']['fwhm'] = 3.0
        if decker == '0.432x12' or decker == '0.432x24':
            par['calibrations']['wavelengths']['fwhm'] = 4.5
        if decker == '0.576x12':
            par['calibrations']['wavelengths']['fwhm'] = 6.0
        if decker == '0.720x12' or decker == '0.720x24':
            par['calibrations']['wavelengths']['fwhm'] = 7.5

        # Return
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
        self.meta['target'] = dict(ext=0, card='TARGNAME')
        self.meta['decker'] = dict(ext=0, card='SLITNAME')
        self.meta['binning'] = dict(ext=0, card=None, default='1,1')

        self.meta['mjd'] = dict(ext=0, card='MJD')
        self.meta['exptime'] = dict(ext=0, card='TRUITIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Extras for config and frametyping
        self.meta['dispname'] = dict(ext=0, card='OBSMODE')
        self.meta['hatch'] = dict(ext=0, card='CALMPOS')
        self.meta['frameno'] = dict(ext=0, card='FRAMENUM')
        self.meta['idname'] = dict(ext=0, card='IMTYPE')
        self.meta['instrument'] = dict(ext=0, card='INSTRUME')
        self.meta['filter1'] = dict(ext=0, card='SCIFILT1')
        self.meta['filter2'] = dict(ext=0, card='SCIFILT2')
        self.meta['echangle'] = dict(ext=0, card='ECHLPOS', rtol=1e-3)
        self.meta['xdangle'] = dict(ext=0, card='DISPPOS', rtol=1e-3)

        # Lamps
        lamp_names = ['NEON', 'ARGON', 'KRYPTON', 'XENON', 'ETALON']
        for kk,lamp_name in enumerate(lamp_names):
            self.meta['lampstat{:02d}'.format(kk+1)] = dict(ext=0, card=lamp_name)
        self.meta['lampstat06'] = dict(ext=0, card = 'HALOGEN')

    def get_echelle_angle_files(self):
        """ Pass back the files required
        to run the echelle method of wavecalib

        Returns:
            list: List of files
        """
        if self.filter1 == 'Kband-new' or self.filter1 == 'KL' or self.filter2 == '':
            band = self.filter1
        else:
            band = self.filter2
        lamps_list = np.copy(self.lamps_list)

        #msgs.info(lamps_list, 'Xe' in lamps_list[0])
        #msgs.info('filter1 = ', filter1)
        if 'Xe' in lamps_list[0]:
            if band == 'NIRSPEC-1':
                angle_fits_file = 'keck_nirspec_y_angle_fits.fits'
                composite_arc_file = 'keck_nirspec_y_composite_arc.fits'
            if band == 'NIRSPEC-3':
                angle_fits_file = 'keck_nirspec_j_angle_fits.fits'
                composite_arc_file = 'keck_nirspec_j_composite_arc.fits'
            if band == 'NIRSPEC-5':
                angle_fits_file = 'keck_nirspec_h_angle_fits.fits'
                composite_arc_file = 'keck_nirspec_h_composite_arc.fits'
            if band == 'Kband-new':
                angle_fits_file = 'keck_nirspec_k_angle_fits.fits'
                composite_arc_file = 'keck_nirspec_k_composite_arc.fits'
            if band == 'NIRSPEC-7':
                angle_fits_file = 'keck_nirspec_k_angle_fits.fits'
                composite_arc_file = 'keck_nirspec_k_composite_arc.fits'
            if band == 'KL':
                angle_fits_file = 'keck_nirspec_l_angle_fits.fits'
                composite_arc_file = 'keck_nirspec_l_composite_arc.fits'
        elif 'OH' in lamps_list[0]:
            msgs.info('Using OH Lines')
            if band == 'NIRSPEC-1':
                angle_fits_file = 'keck_nirspec_y_OH_angle_fits.fits'
                composite_arc_file = 'keck_nirspec_y_composite_OH.fits'
            if band == 'NIRSPEC-3':
                angle_fits_file = 'keck_nirspec_j_OH_angle_fits.fits'
                composite_arc_file = 'keck_nirspec_j_composite_OH.fits'
            if band == 'NIRSPEC-5':
                angle_fits_file = 'keck_nirspec_h_OH_angle_fits.fits'
                composite_arc_file = 'keck_nirspec_h_composite_OH.fits'

        return [angle_fits_file, composite_arc_file]
    
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
        hatch = np.copy(fitstbl['hatch'].data)#.data.astype(int)
        #define science filters so we know the band we're working in for the wavelength calss
        if ftype in ['science', 'standard']:
            return good_exp & self.lamps(fitstbl, 'off') & (hatch == 'Out') 
        if ftype in 'dark':
            return good_exp & self.lamps(fitstbl, 'off') & (hatch == 'In') 
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            # Flats and trace frames are typed together
            return good_exp & self.lamps(fitstbl, 'dome') & (hatch == 'In')
        if ftype == 'pinhole':
            # Don't type pinhole frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc', 'tilt']:
            # TODO: This is a kludge.  Allow science frames to also be
            # classified as arcs
            is_arc = self.lamps(fitstbl, 'arcs') & (hatch == 'In')
            good_exp[is_arc] = True 
            is_obj = self.lamps(fitstbl, 'off') & (hatch == 'Out') 
            good_exp[is_obj] = fitstbl['exptime'].data[is_obj] > 60.0
            return good_exp & (is_arc | is_obj)
        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

    def lamps(self, fitstbl, status):
        """
        Check the lamp status.

        Args:
            fitstbl (`astropy.table.Table`_):
                The table with the fits header meta data.
            status (:obj:`str`):
                The status to check. Can be ``'off'``, ``'arcs'``, or
                ``'dome'``.

        Returns:
            `numpy.ndarray`_: A boolean array selecting fits files that meet
            the selected lamp status.

        Raises:
            ValueError:
                Raised if the status is not one of the valid options.
        """

        if status == 'off':
            # Check if all are off
            lamp_stat = [k for k in fitstbl.keys() if 'lampstat' in k]
            retarr = np.zeros((len(lamp_stat), len(fitstbl)), dtype=bool)
            for kk, key in enumerate(lamp_stat):
                retarr[kk,:] = fitstbl[key] == 'Off'
            return np.all(retarr, axis=0)
        if status == 'arcs':
            # Check if any arc lamps are on
            lamp_stat = [ 'lampstat{0:02d}'.format(i) for i in range(1,6) ]
            retarr = np.zeros((len(lamp_stat), len(fitstbl)))
            for kk, key in enumerate(lamp_stat):
                retarr[kk,:] = fitstbl[key] == 'On'
            return np.any(retarr, axis=0)
        if status == 'dome':
            return fitstbl['lampstat06'] == 'On'

        raise ValueError('No implementation for status = {0}'.format(status))


class KeckNIRSPECHighSpectrographOld(KeckNIRSPECSpectrographOld):
    """
    Child to handle NIRSPEC high-dispersion pre-upgrade specific code
    """
    name = 'keck_nirspec_high_old'
    supported = True
    comment = 'PRE-upgrade (~Dec 2018),High-dispersion grating, see :ref:`nirspec_high_howto`'
    pypeline = 'Echelle'
    ech_fixed_format = False
    lamps_list = []
    filter = ''

    def get_detector_par(self, det, hdu=None):
        """
        Return metadata for the selected detector.
        see https://web.archive.org/web/20161107225756/https://www2.keck.hawaii.edu/inst/nirspec/Specifications.html
        and https://web.archive.org/web/20161102133624/https://www2.keck.hawaii.edu//inst/nirspec/nirspec-spec-highres.pdf

        Args:
            det (:obj:`int`):
                1-indexed detector number.  This is not used because NIRSPEC
                only has one detector!
            hdu (`astropy.io.fits.HDUList`_, optional):
                The open fits file with the raw image of interest.  If not
                provided, frame-dependent parameters are set to a default.

        Returns:
            :class:`~pypeit.images.detector_container.DetectorContainer`:
            Object with the detector metadata.
        """
        detector_dict = dict(
            det=1,
            binning         ='1,1',  # No binning allowed
            dataext         = 0,
            specaxis        = 1,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.193,
            darkcurr        = 0.8,  # this is e-/s/pixel
            saturation      = 1.0e6,
            nonlinear       = 0.9,  # docs say linear to 90,000 but our flats are usually higher
            numamplifiers   = 1,
            mincounts       = -1e10,
            gain            = np.atleast_1d(5.8),
            ronoise         = np.atleast_1d(23.0),
            datasec         = np.atleast_1d('[:,:]'),
            oscansec        = None, #np.atleast_1d('[:,:]')
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
        par['calibrations']['wavelengths']['rms_thresh_frac_fwhm'] = 3.0 #0.20  # Might be grating dependent..
        par['calibrations']['wavelengths']['sigdetect']=5.0
        par['calibrations']['wavelengths']['fwhm']= 3.0
        par['calibrations']['wavelengths']['fwhm_fromlines']= False
        par['calibrations']['wavelengths']['n_final']= 4
        par['calibrations']['wavelengths']['lamps'] = ['NIRSPEC-ArNeKrXe']#, 'ThAr']
        #par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['method'] = 'echelle'
#        par['calibrations']['wavelengths']['ech_fix_format'] = True
        par['calibrations']['wavelengths']['cc_thresh'] = 0.70
        par['calibrations']['wavelengths']['cc_local_thresh'] = 0.50
        # Echelle parameters
        par['calibrations']['wavelengths']['echelle'] = True
        par['calibrations']['wavelengths']['ech_nspec_coeff'] = 4
        par['calibrations']['wavelengths']['ech_norder_coeff'] = 4
        par['calibrations']['wavelengths']['ech_sigrej'] = 3.0
        par['calibrations']['wavelengths']['ech_separate_2d'] = False

        
        # Trace ID parameters
        par['calibrations']['slitedges']['edge_thresh'] = 50.0
        par['calibrations']['slitedges']['fit_order'] = 8
        par['calibrations']['slitedges']['max_shift_adj'] = 0.5
        par['calibrations']['slitedges']['trace_thresh'] = 10.
        par['calibrations']['slitedges']['left_right_pca'] = True
        par['calibrations']['slitedges']['length_range'] = 0.3
        par['calibrations']['slitedges']['max_nudge'] = 10.
        par['calibrations']['slitedges']['overlap'] = True
        par['calibrations']['slitedges']['dlength_range'] = 0.1

        # Flats
        par['calibrations']['flatfield']['tweak_slits_thresh'] = 0.80

        # Extraction
        par['reduce']['skysub']['bspline_spacing'] = 0.8
        par['reduce']['extraction']['sn_gauss'] = 4.0
        par['reduce']['findobj']['fof_link'] = 0.2

        # Flexure
        par['flexure']['spec_method'] = 'skip'

        par['scienceframe']['process']['sigclip'] = 20.0
        par['scienceframe']['process']['satpix'] ='nothing'

        # Should be we be illumflattening?

        # Flats
        turn_off = dict(use_illumflat=False, use_biasimage=False, use_overscan=False,
                        use_darkimage=False, use_pixelflat=False)
        par.reset_all_processimages_par(**turn_off)

        #turn_off = dict(use_biasimage=False, use_overscan=False)
        #par.reset_all_processimages_par(**turn_off)

        # Specify if cleaning cosmic ray hits/bad pixels

        # The settings below enable NIRSPEC dark subtraction from the
        # traceframe and pixelflatframe, but enforce that this bias won't be
        # subtracted from other images. It is a hack for now, because
        # eventually we want to perform this operation with the dark frame
        # class, and we want to attach individual sets of darks to specific
        # images.
        #par['calibrations']['biasframe']['useframe'] = 'bias'
        #par['calibrations']['traceframe']['process']['bias'] = 'force'
        #par['calibrations']['pixelflatframe']['process']['bias'] = 'force'
        #par['calibrations']['arcframe']['process']['bias'] = 'skip'
        #par['calibrations']['tiltframe']['process']['bias'] = 'skip'
        #par['calibrations']['standardframe']['process']['bias'] = 'skip'
        #par['scienceframe']['process']['bias'] = 'skip'

        # Set the default exposure time ranges for the frame typing
        par['calibrations']['standardframe']['exprng'] = [None, None]
        par['calibrations']['arcframe']['exprng'] = [20, None]
        par['calibrations']['darkframe']['exprng'] = [None, None]
        par['scienceframe']['exprng'] = [None, None]

        # Sensitivity function parameters
        par['sensfunc']['algorithm'] = 'IR'
        par['sensfunc']['polyorder'] = 8
        #par['sensfunc']['IR']['telgridfile'] = 'TelFit_MaunaKea_3100_26100_R20000.fits'
        par['sensfunc']['IR']['telgridfile'] = 'TellPCA_9300_55100_R60000.fits'
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

        headarr = self.get_headarr(scifile)
        self.filter1 = self.get_meta_value(headarr, 'filter1')
        self.filter2 = self.get_meta_value(headarr, 'filter2')
        self.lamps_list = par['calibrations']['wavelengths']['lamps']
        decker = self.get_meta_value(headarr, 'decker')


        # wavelength calibration
        supported_filters = ['NIRSPEC-1', 'NIRSPEC-3', 'NIRSPEC-5', 'NIRSPEC-7', 'KL']
        if (self.filter1 not in supported_filters) and (self.filter2 not in supported_filters):
            msgs.warn(f'Filter {self.filter1} or {self.filter2} may not be supported!!')
        
        if self.filter2 == 'NIRSPEC-7':
            par['calibrations']['wavelengths']['n_final'] = 3
            par['calibrations']['wavelengths']['ech_nspec_coeff'] = 3
            par['calibrations']['wavelengths']['cc_thresh'] = 0.5
            par['calibrations']['wavelengths']['cc_local_thresh'] = 0.5
            par['calibrations']['wavelengths']['cc_offset_minmax'] = 0.25
            par['calibrations']['wavelengths']['cc_percent_ceil'] = 99.9
            par['calibrations']['wavelengths']['echelle_pad'] = 1
            par['calibrations']['wavelengths']['stretch_func'] = 'quadratic'

            par['calibrations']['slitedges']['overlap'] = False

        if self.filter1 == 'KL' or self.filter2 == 'KL':
            par['calibrations']['wavelengths']['n_final'] = 2
            par['calibrations']['wavelengths']['ech_nspec_coeff'] = 2
            par['calibrations']['wavelengths']['cc_thresh'] = 0.5
            par['calibrations']['wavelengths']['cc_local_thresh'] = 0.5
            par['calibrations']['wavelengths']['cc_offset_minmax'] = 0.25
            par['calibrations']['wavelengths']['cc_percent_ceil'] = 99.9
            par['calibrations']['wavelengths']['echelle_pad'] = 1
            par['calibrations']['wavelengths']['stretch_func'] = 'quadratic'

            par['calibrations']['slitedges']['overlap'] = False

        if self.filter2 == 'NIRSPEC-5':
            par['calibrations']['wavelengths']['n_final'] = 3
            par['calibrations']['wavelengths']['ech_nspec_coeff'] = 3
            par['calibrations']['wavelengths']['cc_thresh'] = 0.5
            par['calibrations']['wavelengths']['cc_local_thresh'] = 0.5
            par['calibrations']['wavelengths']['cc_offset_minmax'] = 0.25
            par['calibrations']['wavelengths']['cc_percent_ceil'] = 70.0
            par['calibrations']['wavelengths']['echelle_pad'] = 1
            par['calibrations']['wavelengths']['stretch_func'] = 'quadratic'


        if self.filter2 == 'NIRSPEC-3':
            par['calibrations']['wavelengths']['n_final'] = 3
            par['calibrations']['wavelengths']['ech_nspec_coeff'] = 3
            par['calibrations']['wavelengths']['cc_thresh'] = 0.5
            par['calibrations']['wavelengths']['cc_local_thresh'] = 0.5
            par['calibrations']['wavelengths']['cc_offset_minmax'] = 0.25
            par['calibrations']['wavelengths']['cc_percent_ceil'] = 99.9
            par['calibrations']['wavelengths']['echelle_pad'] = 0
            par['calibrations']['wavelengths']['stretch_func'] = 'quadratic'

        if self.filter2 == 'NIRSPEC-1':
            par['calibrations']['wavelengths']['n_final'] = 3
            par['calibrations']['wavelengths']['ech_nspec_coeff'] = 3
            par['calibrations']['wavelengths']['cc_thresh'] = 0.5
            par['calibrations']['wavelengths']['cc_local_thresh'] = 0.5
            par['calibrations']['wavelengths']['cc_offset_minmax'] = 0.25
            par['calibrations']['wavelengths']['cc_percent_ceil'] = 99.9
            par['calibrations']['wavelengths']['echelle_pad'] = 1
            par['calibrations']['wavelengths']['stretch_func'] = 'quadratic'
        if decker == '0.144x12':
            par['calibrations']['wavelengths']['fwhm'] = 2.0
        if decker == '0.288x12' or decker == '0.288x24':
            par['calibrations']['wavelengths']['fwhm'] = 2.0
        if decker == '0.432x12' or decker == '0.432x24':
            par['calibrations']['wavelengths']['fwhm'] = 3.0
        if decker == '0.576x12':
            par['calibrations']['wavelengths']['fwhm'] = 4.0
        if decker == '0.720x12' or decker == '0.720x24':
            par['calibrations']['wavelengths']['fwhm'] = 5.0

        # Return
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
        self.meta['target'] = dict(ext=0, card='TARGNAME')
        self.meta['decker'] = dict(ext=0, card='SLITNAME')
        self.meta['binning'] = dict(ext=0, card=None, default='1,1')

        self.meta['mjd'] = dict(ext=0, card='MJD-OBS')
        self.meta['exptime'] = dict(ext=0, card='ITIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Extras for config and frametyping
        self.meta['dispname'] = dict(ext=0, card=None, default='Spectroscopy')
        self.meta['hatch'] = dict(ext=0, card='CALMPOS')
        self.meta['frameno'] = dict(ext=0, card='FILENUM')
        self.meta['idname'] = dict(ext=0, card=None, default = 'Spectrum')
        self.meta['instrument'] = dict(ext=0, card='CURRINST')
        self.meta['filter1'] = dict(ext=0, card='FILNAME')
        self.meta['filter2'] = dict(ext=0, card=None, default='')
        self.meta['echangle'] = dict(ext=0, card='ECHLPOS', rtol=1e-3)
        self.meta['xdangle'] = dict(ext=0, card='DISPPOS', rtol=1e-3)
        #self.meta['binning'] = dict(ext=0, card=None, default='1,1')
        #self.meta['imagrot'] = dict(ext=0, card='IROTPOS', rtol=1e-3)

        # Lamps
        lamp_names = ['NEON', 'ARGON', 'KRYPTON', 'XENON', 'ETALON']
        for kk,lamp_name in enumerate(lamp_names):
            self.meta['lampstat{:02d}'.format(kk+1)] = dict(ext=0, card=lamp_name)
        self.meta['lampstat06'] = dict(ext=0, card = 'FLAT')

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
        return ['filter1', 'filter2', 'echangle', 'xdangle']

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
        return ['SCIFILT1', 'SCIFILT2', 'ECHLPOS', 'DISPPOS', 'SLITNAME']

    def get_echelle_angle_files(self):
        """ Pass back the files required
        to run the echelle method of wavecalib

        Created for the pre-upgrade NIRSPEC

        Returns:
            list: List of files
        """
        band = self.filter1        
        lamps_list = np.copy(self.lamps_list)

        if 'Xe' in lamps_list[0]:
            if band == 'NIRSPEC-1':
                angle_fits_file = 'keck_nirspec_y_preupgrade_angle_fits.fits'
                composite_arc_file = 'keck_nirspec_y_preupgrade_composite_arc.fits'
            if band == 'NIRSPEC-3':
                angle_fits_file = 'keck_nirspec_j_preupgrade_angle_fits.fits'
                composite_arc_file = 'keck_nirspec_j_preupgrade_composite_arc.fits'
            if band == 'NIRSPEC-5':
                angle_fits_file = 'keck_nirspec_h_preupgrade_angle_fits.fits'
                composite_arc_file = 'keck_nirspec_h_preupgrade_composite_arc.fits'
            if band == 'NIRSPEC-7':
                angle_fits_file = 'keck_nirspec_k_preupgrade_angle_fits.fits'
                composite_arc_file = 'keck_nirspec_k_preupgrade_composite_arc.fits'
        elif 'OH' in lamps_list[0]:
            print('Using OH Lines')
            if band == 'NIRSPEC-1':
                angle_fits_file = 'keck_nirspec_y_preupgrade_OH_angle_fits.fits'
                composite_arc_file = 'keck_nirspec_y_preupgrade_composite_OH.fits'
            if band == 'NIRSPEC-3':
                angle_fits_file = 'keck_nirspec_j_preupgrade_OH_angle_fits.fits'
                composite_arc_file = 'keck_nirspec_j_preupgrade_composite_OH.fits'

        return [angle_fits_file, composite_arc_file]

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
        hatch = np.copy(fitstbl['hatch'].data)#.data.astype(int)

        if ftype in ['science']:
            return good_exp & self.lamps(fitstbl, 'off') & (hatch == '0') 
        if ftype in 'dark':
            return good_exp & self.lamps(fitstbl, 'off') & (hatch == '1') 
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            # Flats and trace frames are typed together
            return good_exp & self.lamps(fitstbl, 'dome') & (hatch == '1') 
        if ftype == 'pinhole':
            # Don't type pinhole frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc', 'tilt']:
            # TODO: This is a kludge.  Allow science frames to also be
            # classified as arcs
            is_arc = self.lamps(fitstbl, 'arcs') & (hatch == '1') 
            good_exp[is_arc] = True
            is_obj = self.lamps(fitstbl, 'off') & (hatch == '0') 
            good_exp[is_obj] = fitstbl['exptime'].data[is_obj] > 60.0
            return good_exp & (is_arc | is_obj)
        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

    def lamps(self, fitstbl, status):
        """
        Check the lamp status.

        Args:
            fitstbl (`astropy.table.Table`_):
                The table with the fits header meta data.
            status (:obj:`str`):
                The status to check. Can be ``'off'``, ``'arcs'``, or
                ``'dome'``.

        Returns:
            `numpy.ndarray`_: A boolean array selecting fits files that meet
            the selected lamp status.

        Raises:
            ValueError:
                Raised if the status is not one of the valid options.
        """
        if status == 'off':
            # Check if all are off
            lamp_stat = [k for k in fitstbl.keys() if 'lampstat' in k]
            retarr = np.zeros((len(lamp_stat), len(fitstbl)), dtype=bool)
            for kk, key in enumerate(lamp_stat):
                retarr[kk,:] = fitstbl[key] == '0'
            return np.all(retarr, axis=0)
        if status == 'arcs':
            # Check if any arc lamps are on
            lamp_stat = [ 'lampstat{0:02d}'.format(i) for i in range(1,6) ]
            retarr = np.zeros((len(lamp_stat), len(fitstbl)))
            for kk, key in enumerate(lamp_stat):
                retarr[kk,:] = fitstbl[key] == '1'
            return np.any(retarr, axis=0)
        if status == 'dome':
            return fitstbl['lampstat06'] == '1'

        raise ValueError('No implementation for status = {0}'.format(status))

    def order_platescale(self, order_vec, binning=None):
        """
        Return the platescale for each echelle order.

        Note that NIRSPEC has no binning.

        Args:
            order_vec (`numpy.ndarray`_):
                The vector providing the order numbers.
            binning (:obj:`str`, optional):
                The string defining the spectral and spatial binning. **This
                is always ignored.**

        Returns:
            `numpy.ndarray`_: An array with the platescale for each order
            provided by ``order``.
        """
        pscale = self.get_detector_par(1).platescale
        return np.full(order_vec.size, pscale)

    def get_rawimage(self, raw_file, det):
        """
        Read raw images and generate a few other bits and pieces that are key
        for image processing.

        .. warning::

            - When reading multiple detectors for a mosaic, this function
              expects all detector arrays to have exactly the same shape.

        Parameters
        ----------
        raw_file : :obj:`str`, `Path`_
            File to read
        det : :obj:`int`, :obj:`tuple`
            1-indexed detector(s) to read.  An image mosaic is selected using a
            :obj:`tuple` with the detectors in the mosaic, which must be one of
            the allowed mosaics returned by :func:`allowed_mosaics`.

        Returns
        -------
        detector_par : :class:`~pypeit.images.detector_container.DetectorContainer`, :class:`~pypeit.images.mosaic.Mosaic`
            Detector metadata parameters for one or more detectors.
        raw_img : `numpy.ndarray`_
            Raw image for this detector.  Shape is 2D if a single detector image
            is read and 3D if multiple detectors are read.  E.g., the image from
            the first detector in the tuple is accessed using ``raw_img[0]``.
        hdu : `astropy.io.fits.HDUList`_
            Opened fits file
        exptime : :obj:`float`
            Exposure time *in seconds*.
        rawdatasec_img : `numpy.ndarray`_
            Data (Science) section of the detector as provided by setting the
            (1-indexed) number of the amplifier used to read each detector
            pixel. Pixels unassociated with any amplifier are set to 0.  Shape
            is identical to ``raw_img``.
        oscansec_img : `numpy.ndarray`_
            Overscan section of the detector as provided by setting the
            (1-indexed) number of the amplifier used to read each detector
            pixel. Pixels unassociated with any amplifier are set to 0.  Shape
            is identical to ``raw_img``.
        """
        # Check extension and then open
        self._check_extensions(raw_file)
        hdu = io.fits_open(raw_file)

        # Validate the entered (list of) detector(s)
        nimg, _det = self.validate_det(det)

        # Grab the detector or mosaic parameters
        mosaic = None if nimg == 1 else self.get_mosaic_par(det, hdu=hdu)
        detectors = [self.get_detector_par(det, hdu=hdu)] if nimg == 1 else mosaic.detectors

        # Grab metadata from the header
        # NOTE: These metadata must be *identical* for all images when reading a
        # mosaic
        headarr = self.get_headarr(hdu)

        # save filter being used:
        self.filter1 = self.get_meta_value(headarr, 'filter1')

        # Exposure time (used by RawImage)
        # NOTE: This *must* be (converted to) seconds.
        exptime = self.get_meta_value(headarr, 'exptime')

        # Rawdatasec, oscansec images
        binning = self.get_meta_value(headarr, 'binning')
        # NOTE: This means that `specaxis` must be the same for all detectors in
        # a mosaic
        if detectors[0]['specaxis'] == 1:
            binning_raw = (',').join(binning.split(',')[::-1])
        else:
            binning_raw = binning

        raw_img = [None]*nimg
        rawdatasec_img = [None]*nimg
        oscansec_img = [None]*nimg
        for i in range(nimg):

            # Raw image
            raw_img[i] = hdu[detectors[i]['dataext']].data.astype(float)
            # Raw data from some spectrograph (i.e. FLAMINGOS2) have an addition
            # extention, so I add the following two lines. It's easier to change
            # here than writing another get_rawimage function in the
            # spectrograph file.
            # TODO: This feels dangerous, but not sure what to do about it...
            #if raw_img[i].ndim != 2:
            #    raw_img[i] = np.squeeze(raw_img[i])
            if raw_img[i].ndim != 2:
                msgs.error(f"Raw images must be 2D; check extension {detectors[i]['dataext']} "
                           f"of {raw_file}.")

            for section in ['datasec', 'oscansec']:

                # Get the data section
                # Try using the image sections as header keywords
                # TODO -- Deal with user windowing of the CCD (e.g. Kast red)
                #  Code like the following maybe useful
                #hdr = hdu[detector[det - 1]['dataext']].header
                #image_sections = [hdr[key] for key in detector[det - 1][section]]
                # Grab from Detector
                image_sections = detectors[i][section]
                #if not isinstance(image_sections, list):
                #    image_sections = [image_sections]
                # Always assume normal FITS header formatting
                one_indexed = True
                include_last = True

                # Initialize the image (0 means no amplifier)
                pix_img = np.zeros(raw_img[i].shape, dtype=int)
                for j in range(detectors[i]['numamplifiers']):

                    if image_sections is not None:  # and image_sections[i] is not None:
                        # Convert the data section from a string to a slice
                        datasec = parse.sec2slice(image_sections[j], one_indexed=one_indexed,
                                                  include_end=include_last, require_dim=2,
                                                  binning=binning_raw)
                        # Assign the amplifier
                        pix_img[datasec] = j+1

                # Finish
                if section == 'datasec':
                    rawdatasec_img[i] = pix_img.copy()
                else:
                    oscansec_img[i] = pix_img.copy()

        if nimg == 1:
            # Return single image
            return detectors[0], raw_img[0], hdu, exptime, rawdatasec_img[0], oscansec_img[0]

        if any([img.shape != raw_img[0].shape for img in raw_img[1:]]):
            msgs.error('All raw images in a mosaic must have the same shape.')
        # Return all images for mosaic
        return mosaic, np.array(raw_img), hdu, exptime, np.array(rawdatasec_img), np.array(oscansec_img)


class KeckNIRSPECLowSpectrograph(KeckNIRSPECSpectrograph):
    """
    Child to handle NIRSPEC low-dispersion specific code
    """
    name = 'keck_nirspec_low'
    supported = True
    comment = 'POST-upgrade (~Dec 2018), Low-dispersion grating'
    pypeline = 'MultiSlit'

    def get_detector_par(self, det, hdu=None):
        """
        Return metadata for the selected detector.
        see https://www2.keck.hawaii.edu/inst/nirspec/Specifications.html
        and https://www2.keck.hawaii.edu/inst/nirspec/images/spec_lowres.png

        Args:
            det (:obj:`int`):
                1-indexed detector number.  This is not used because NIRSPEC
                only has one detector!
            hdu (`astropy.io.fits.HDUList`_, optional):
                The open fits file with the raw image of interest.  If not
                provided, frame-dependent parameters are set to a default.

        Returns:
            :class:`~pypeit.images.detector_container.DetectorContainer`:
            Object with the detector metadata.
        """
        detector_dict = dict(
            det=1,
            binning         ='1,1',  # No binning allowed
            dataext         = 0,
            specaxis        = 0,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.098,
            darkcurr        = 2520.0,  # e-/pixel/hour  (=0.7 e-/pixel/s)
            saturation      = 100000.,
            nonlinear       = 0.9,  # docs say linear to 90,000 but our flats are usually higher
            numamplifiers   = 1,
            mincounts       = -1e10,
            gain            = np.atleast_1d(3.01),
            ronoise         = np.atleast_1d(11.56),
            datasec         = np.atleast_1d('[:,:]'),
            oscansec        = None, #np.atleast_1d('[:,:]')
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
        par['calibrations']['wavelengths']['rms_thresh_frac_fwhm'] = 0.06 #0.20  # Might be grating dependent..
        par['calibrations']['wavelengths']['sigdetect']=5.0
        par['calibrations']['wavelengths']['fwhm']= 3.5
        par['calibrations']['wavelengths']['n_final']= 4
        par['calibrations']['wavelengths']['lamps'] = ['OH_NIRES']
        #par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['method'] = 'holy-grail'
        # Reidentification parameters
        #par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_nires.fits'
        par['calibrations']['slitedges']['edge_thresh'] = 200.
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'

        # Flats
        par['calibrations']['flatfield']['tweak_slits_thresh'] = 0.80

        # Extraction
        par['reduce']['skysub']['bspline_spacing'] = 0.8
        par['reduce']['extraction']['sn_gauss'] = 4.0

        # Flexure
        par['flexure']['spec_method'] = 'skip'

        par['scienceframe']['process']['sigclip'] = 20.0
        par['scienceframe']['process']['satpix'] ='nothing'

        # Should be we be illumflattening?

        # Flats
        turn_off = dict(use_illumflat=False, use_biasimage=False, use_overscan=False,
                        use_darkimage=False)
        par.reset_all_processimages_par(**turn_off)

        #turn_off = dict(use_biasimage=False, use_overscan=False)
        #par.reset_all_processimages_par(**turn_off)


        # The settings below enable NIRSPEC dark subtraction from the
        # traceframe and pixelflatframe, but enforce that this bias won't be
        # subtracted from other images. It is a hack for now, because
        # eventually we want to perform this operation with the dark frame
        # class, and we want to attach individual sets of darks to specific
        # images.
        #par['calibrations']['biasframe']['useframe'] = 'bias'
        #par['calibrations']['traceframe']['process']['bias'] = 'force'
        #par['calibrations']['pixelflatframe']['process']['bias'] = 'force'
        #par['calibrations']['arcframe']['process']['bias'] = 'skip'
        #par['calibrations']['tiltframe']['process']['bias'] = 'skip'
        #par['calibrations']['standardframe']['process']['bias'] = 'skip'
        #par['scienceframe']['process']['bias'] = 'skip'

        # Set the default exposure time ranges for the frame typing
        par['calibrations']['standardframe']['exprng'] = [None, 20]
        par['calibrations']['arcframe']['exprng'] = [20, None]
        par['calibrations']['darkframe']['exprng'] = [20, None]
        par['scienceframe']['exprng'] = [20, None]

        # Sensitivity function parameters
        par['sensfunc']['algorithm'] = 'IR'
        par['sensfunc']['polyorder'] = 8
        par['sensfunc']['IR']['telgridfile'] = 'TellPCA_3000_26000_R25000.fits'
        par['sensfunc']['IR']['pix_shift_bounds'] = (-8.0,8.0)

        # Telluric parameters
        par['telluric']['pix_shift_bounds'] = (-8.0,8.0)

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
        self.meta['target'] = dict(ext=0, card='TARGNAME')
        self.meta['decker'] = dict(ext=0, card='SLITNAME')
        self.meta['binning'] = dict(ext=0, card=None, default='1,1')

        self.meta['mjd'] = dict(ext=0, card='MJD-OBS')
        self.meta['exptime'] = dict(ext=0, card='ELAPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Extras for config and frametyping
        self.meta['dispname'] = dict(ext=0, card='DISPERS')
        self.meta['hatch'] = dict(ext=0, card='CALMPOS')
        self.meta['idname'] = dict(ext=0, card='IMAGETYP')
        self.meta['instrument'] = dict(ext=0, card='INSTRUME')
        # Lamps
        lamp_names = ['NEON', 'ARGON', 'KRYPTON', 'XENON', 'ETALON', 'FLAT']
        for kk,lamp_name in enumerate(lamp_names):
            self.meta['lampstat{:02d}'.format(kk+1)] = dict(ext=0, card=lamp_name)

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
        return ['decker', 'dispname']

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
        return ['SLITNAME', 'DISPERS']

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
        hatch = fitstbl['hatch'].data.astype(int)
        if ftype in ['science', 'standard']:
            return good_exp & self.lamps(fitstbl, 'off') & (hatch == 0) \
                        & (fitstbl['idname'] == 'object')
        if ftype in ['bias', 'dark']:
            return good_exp & self.lamps(fitstbl, 'off') & (hatch == 0) \
                        & (fitstbl['idname'] == 'dark')
        if ftype in ['pixelflat', 'trace']:
            # Flats and trace frames are typed together
            return good_exp & self.lamps(fitstbl, 'dome') & (hatch == 1) \
                        & (fitstbl['idname'] == 'flatlamp')
        if ftype == 'pinhole':
            # Don't type pinhole frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc', 'tilt']:
            # TODO: This is a kludge.  Allow science frames to also be
            # classified as arcs
            is_arc = self.lamps(fitstbl, 'arcs') & (hatch == 1) \
                            & (fitstbl['idname'] == 'arclamp')
            is_obj = self.lamps(fitstbl, 'off') & (hatch == 0) \
                        & (fitstbl['idname'] == 'object')
            return good_exp & (is_arc | is_obj)
        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

    def lamps(self, fitstbl, status):
        """
        Check the lamp status.

        Args:
            fitstbl (`astropy.table.Table`_):
                The table with the fits header meta data.
            status (:obj:`str`):
                The status to check. Can be ``'off'``, ``'arcs'``, or
                ``'dome'``.

        Returns:
            `numpy.ndarray`_: A boolean array selecting fits files that meet
            the selected lamp status.

        Raises:
            ValueError:
                Raised if the status is not one of the valid options.
        """
        if status == 'off':
            # Check if all are off
            lamp_stat = [k for k in fitstbl.keys() if 'lampstat' in k]
            retarr = np.zeros((len(lamp_stat), len(fitstbl)), dtype=bool)
            for kk, key in enumerate(lamp_stat):
                retarr[kk,:] = fitstbl[key].data.astype(int) == 0
            return np.all(retarr, axis=0)
        if status == 'arcs':
            # Check if any arc lamps are on
            lamp_stat = [ 'lampstat{0:02d}'.format(i) for i in range(1,6) ]
            retarr = np.zeros((len(lamp_stat), len(fitstbl)))
            for kk, key in enumerate(lamp_stat):
                retarr[kk,:] = fitstbl[key].data.astype(int) == 1
            return np.any(retarr, axis=0)
        if status == 'dome':
            return fitstbl['lampstat06'].data.astype(int) == 1

        raise ValueError('No implementation for status = {0}'.format(status))

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

        # Edges of the detector are junk
        msgs.info("Custom bad pixel mask for NIRSPEC")
        bpm_img[:, :20] = 1.
        bpm_img[:, 1000:] = 1.

        return bpm_img


# NOTE: Including this unfinished class causes the docs to fault because the
# necessary meta data (using the init_meta method) is not defined.  I (KBW) am
# commenting this out until this spectrograph class is ready to be fully
# implemented.
#
#class KeckNIRSPECLowSpectrographOld(KeckNIRSPECSpectrographOld):
#    """
#    Child to handle NIRSPEC low-dispersion pre-upgrade specific code
#    """
#    name = 'keck_nirspec_low_old'
#    supported = False
#    comment = 'PRE-upgrade (~Dec 2018), Low-dispersion grating'
#    pypeline = 'MultiSlit'
#
#    # TODO: this is not tested, but added here for completeness
#
#    def get_detector_par(self, det, hdu=None):
#        """
#        Return metadata for the selected detector.
#        see https://web.archive.org/web/20161107225756/https://www2.keck.hawaii.edu/inst/nirspec/Specifications.html
#        and https://web.archive.org/web/20161102133622/https://www2.keck.hawaii.edu//inst/nirspec/nirspec-spec-lowres.pdf
#
#        Args:
#            det (:obj:`int`):
#                1-indexed detector number.  This is not used because NIRSPEC
#                only has one detector!
#            hdu (`astropy.io.fits.HDUList`_, optional):
#                The open fits file with the raw image of interest.  If not
#                provided, frame-dependent parameters are set to a default.
#
#        Returns:
#            :class:`~pypeit.images.detector_container.DetectorContainer`:
#            Object with the detector metadata.
#        """
#        detector_dict = dict(
#            det=1,
#            binning='1,1',  # No binning allowed
#            dataext=0,
#            specaxis=0,
#            specflip=False,
#            spatflip=False,
#            platescale=0.193,
#            darkcurr=2880.0,  # this is e-/pixel/hour... == 0.8 e-/s/pixel
#            saturation=1.0e6,
#            nonlinear=0.9,  # docs say linear to 90,000 but our flats are usually higher
#            numamplifiers=1,
#            mincounts=-1e10,
#            gain=np.atleast_1d(5.8),
#            ronoise=np.atleast_1d(23.0),
#            datasec=np.atleast_1d('[:,:]'),
#            oscansec=None,  # np.atleast_1d('[:,:]')
#        )
#
#        return detector_container.DetectorContainer(**detector_dict)
