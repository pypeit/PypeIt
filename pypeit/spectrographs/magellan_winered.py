"""
Module for Magellan/WINERED specific methods.

.. include:: ../include/links.rst
"""
import numpy as np

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container


class MagellanWINEREDSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Magellan/WINERED specific code

    """
    ndet = 1
    telescope = telescopes.MagellanTelescopePar()
    camera = 'WINERED'
    url = 'http://lihweb.kyoto-su.ac.jp/WINERED/index.html'
    header_name = 'WINERED'
    name = 'magellan_winered'
    pypeline = 'Echelle'
    ech_fixed_format = True
    supported = True
    comment = 'Magellan WINERED is a near-IR (0.91-1.35um) spectrograph at the Magellan Clay telescope. Currently only supports HIRES-Y mode.'

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
        self.meta['dichroic'] = dict(ext=0, card=None, default='default')
        self.meta['binning'] = dict(ext=0, card=None, default='1,1')

        self.meta['mjd'] = dict(ext=0, card='ACQTIME1')
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Extras for config and frametyping
        self.meta['slitwid'] = dict(ext=0, card='SLIT')
        self.meta['dispname'] = dict(ext=0, card='INSTMODE')
        self.meta['idname'] = dict(ext=0, card='DATA-TYP')
        self.meta['instrument'] = dict(ext=0, card='INSTRUME')
        self.meta['lampstat01'] = dict(ext=0, card='CMPLAMP')
        self.meta['lampstat02'] = dict(ext=0, card='INSFLAT')

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
        # TODO :: Come back and pull gain/RONOISE from the fits header
        # TODO :: update all parameter values here
        # Detector 1
        detector_dict = dict(
            binning         = '1,1',
            det             = 1,
            dataext         = 0,
            specaxis        = 1,
            specflip        = True,
            spatflip        = False,
            platescale      = 0.18,
            darkcurr        = 3.06,  # e-/pixel/hour  (=0.00085 e-/pixel/s)
            #saturation      = 20000., # high gain is 20000 ADU, low gain is 32000 ADU
            saturation      = 100000., # This is an arbitrary value.
            nonlinear       = 1.0, # high gain mode, low gain is 0.875
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d(1.2), # high gain mode, low gain is 3.8 e-/DN
            ronoise         = np.atleast_1d(5.0), # for high gain mode and SUTR read modes with exptime ~ 900s
            datasec         = np.atleast_1d('[5:2044,5:2044]'),
            oscansec        = np.atleast_1d('[5:2044,:5]')
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
        # 1D wavelength solution with OH lines
        # par['calibrations']['wavelengths']['rms_thresh_frac_fwhm'] = 0.25
        # par['calibrations']['wavelengths']['sigdetect']=[5,5,10,10,10,20,30,30,30,30,30,10,30,30,60,30,30,10,20,30,10]
        # par['calibrations']['wavelengths']['n_first']=2
        # par['calibrations']['wavelengths']['n_final']=[3,2,3,2,4,4,4,3,4,4,4,3,4,4,4,4,4,4,6,6,4]
        # par['calibrations']['wavelengths']['lamps'] = ['OH_FIRE_Echelle']
        # par['calibrations']['wavelengths']['method'] = 'reidentify'
        # par['calibrations']['wavelengths']['cc_thresh'] = 0.35
        # par['calibrations']['wavelengths']['reid_arxiv'] = 'magellan_fire_echelle.fits'

        # Echelle parameters
        # par['calibrations']['wavelengths']['echelle'] = True
        # par['calibrations']['wavelengths']['ech_nspec_coeff'] = 4
        # par['calibrations']['wavelengths']['ech_norder_coeff'] = 6
        # par['calibrations']['wavelengths']['ech_sigrej'] = 3.0

        # Set the parameters for the science frames
        par['scienceframe']['process']['sigclip'] = 20.0
        par['scienceframe']['process']['satpix'] ='nothing'

        # Set slits and tilts parameters
        par['calibrations']['tilts']['tracethresh'] = 5
        par['calibrations']['slitedges']['edge_thresh'] = 3.
        par['calibrations']['slitedges']['trace_thresh'] = 10.
        par['calibrations']['slitedges']['fit_order'] = 5
        par['calibrations']['slitedges']['max_shift_adj'] = 0.5
        par['calibrations']['slitedges']['fit_min_spec_length'] = 0.5
        par['calibrations']['slitedges']['left_right_pca'] = True
        par['calibrations']['slitedges']['pca_order'] = 3

        # Model entire slit
        par['reduce']['extraction']['model_full_slit'] = True  # local sky subtraction operates on entire slit
        par['reduce']['findobj']['maxnumber_sci'] = 2  # Slit is narrow so allow one object per order
        par['reduce']['findobj']['maxnumber_std'] = 1  # Slit is narrow so allow one object per order

        # Processing steps
        turn_off = dict(use_illumflat=False, use_biasimage=False, use_overscan=False,
                        use_darkimage=False)
        par.reset_all_processimages_par(**turn_off)
        # Do not correct for flexure
        par['flexure']['spec_method'] = 'skip'

        # Set the default exposure time ranges for the frame typing
        par['calibrations']['standardframe']['exprng'] = [None, 60]
        par['calibrations']['arcframe']['exprng'] = [1, 20]
        par['calibrations']['darkframe']['exprng'] = [1, None]
        par['scienceframe']['exprng'] = [10, None]

        # Sensitivity function parameters
        # Sensitivity function parameters
        par['sensfunc']['algorithm'] = 'IR'
        par['sensfunc']['polyorder'] = 5
        par['sensfunc']['IR']['maxiter'] = 2
        # place holder for telgrid file
        par['sensfunc']['IR']['telgridfile'] = 'TellPCA_3000_26000_R15000.fits'

        # Coadding
        par['coadd1d']['wave_method'] = 'log10'


        return par

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
            lamp_stat = ['lampstat01', 'lampstat02', 'target']
            # Check if all are off
            chk1 = np.array([(fitstbl[k] == 'off') for k in fitstbl.keys() if k in lamp_stat])
            chk2 = np.array([(fitstbl['idname'] == 'DOMEFLAT') & (fitstbl['target'] == 'off')])
            return np.all(chk1&chk2, axis=0)
        if status == 'arcs':
            # Check if any arc lamps are on
            arc_lamp_stat = [ 'lampstat01']
            return np.any(np.array([ fitstbl[k] == 'on' for k in fitstbl.keys()
                                            if k in arc_lamp_stat]), axis=0)
        if status == 'dome':
            return np.any(np.array([fitstbl['target'] == 'on']), axis=0)
        if status == 'internal':
            # Check if any internal lamps are on
            dome_lamp_stat = [ 'lampstat02' ]
            return np.any(np.array([ fitstbl[k] == 'on' for k in fitstbl.keys()
                                            if k in dome_lamp_stat]), axis=0)
        raise ValueError('No implementation for status = {0}'.format(status))

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
        if ftype in ['pinhole', 'bias']:
            # No pinhole or bias frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['dark', 'lampoffflats']:
            return good_exp & (fitstbl['idname'] == 'DOMEFLAT') & self.lamps(fitstbl, 'off')
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            return good_exp & (fitstbl['idname'] == 'DOMEFLAT') & \
                   (self.lamps(fitstbl, 'internal') | self.lamps(fitstbl, 'dome'))
        if ftype == 'standard':
            return good_exp & (fitstbl['idname'] == 'OBJECT')
        if ftype == 'science':
            return good_exp & (fitstbl['idname'] == 'OBJECT')
        if ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['idname'] == 'COMPARISON')
        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

    @property
    def norders(self):
        """
        Number of orders for this spectrograph. Should only defined for
        echelle spectrographs, and it is undefined for the base class.
        """
        self.check_disperser()
        if 'HIRES-Y' in self.dispname:
            return 31  # This includes orders that are only partially on the detector
        elif 'HIRES-J' in self.dispname:
            return 25  # TODO :: This is a guess based on the online documentation (probably doesn't include partial orders)
        elif 'WIDE' in self.dispname:
            return 20  # TODO :: This is a guess based on the online documentation (probably doesn't include partial orders)
        else:
            msgs.error('Unrecognized disperser')

    @property
    def order_spat_pos(self):
        """
        Return the expected spatial position of each echelle order.
        """
        ord_spat_pos = np.array([0.078125, 0.13769531, 0.19189453, 0.24414062, 0.29296875,
                                 0.34179688, 0.38330078, 0.42724609, 0.46582031, 0.50439453,
                                 0.54199219, 0.57763672, 0.61279297, 0.6484375 , 0.68457031,
                                 0.71875   , 0.75439453, 0.79443359, 0.83789062, 0.88671875,
                                 0.94091797])
        return ord_spat_pos

    @property
    def orders(self):
        """
        Return the order number for each echelle order.
        """
        return np.arange(31, 10, -1, dtype=int)

    @property
    def spec_min_max(self):
        """
        Return the minimum and maximum spectral pixel expected for the
        spectral range of each order.
        """
        spec_max = np.asarray([2048,2048,2048,2048,2048,2048,2048,2048,2048,2048,2048,2048,2048,2048,2048,2048,
                               2048,2048,2048,2048,2048])
        spec_min = np.asarray([ 500,   0,   0,   0,   0,   0,   0,    0,   0,   0,  0,   0,   0,   0,   0,   0,
                                  0,   0,   0,   0,   0])
        return np.vstack((spec_min, spec_max))

    def order_platescale(self, order_vec, binning=None):
        """
        Return the platescale for each echelle order.

        Note that FIRE has no binning.

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
        return np.full(order_vec.size, 0.27)

    # @property
    # def dloglam(self):
    #     """
    #     Return the logarithmic step in wavelength for output spectra.
    #     """
    #     # This number was determined using the resolution and sampling quoted on the FIRE website
    #     R = 6000.0 * 2.7
    #     dloglam = 1.0 / R / np.log(10.0)
    #     return dloglam
    #
    # @property
    # def loglam_minmax(self):
    #     """
    #     Return the base-10 logarithm of the first and last wavelength for
    #     output spectra.
    #     """
    #     return np.log10(8000.0), np.log10(25700)
