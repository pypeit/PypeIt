"""
Implements APF-specific functions

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import os

import numpy as np
from astropy.time import Time
from IPython import embed

from pypeit import msgs
from pypeit import telescopes
from pypeit import io
from pypeit.core import framematch
from pypeit.core import parse
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container



class APFLevySpectrograph(spectrograph.Spectrograph):
    """
    Child to handle APF specific code.

    This spectrograph is partially supported. The 
    data will be reduced but it will require some user
    intervention by editing the pypeit file.
    """
    ndet = 1
    telescope = telescopes.APFTelescopePar()
    pypeline = 'Echelle'
    name = 'apf_levy'
    camera = 'apf'
    header_name = 'apf'
    ech_fixed_format = True
    supported = True

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.
        
        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of ``PypeIt`` methods.
        """
        par = super().default_pypeit_par()

        par['calibrations']['slitedges']['edge_thresh'] = 1.
        par['calibrations']['slitedges']['fit_order'] = 4
        par['calibrations']['slitedges']['max_shift_adj'] = 0.5
        par['calibrations']['slitedges']['left_right_pca'] = True
        par['calibrations']['slitedges']['smash_range'] = [0.35,0.65]

        par['calibrations']['tilts']['tracethresh'] = 20

        # 1D wavelength solution
        par['calibrations']['wavelengths']['reid_arxiv'] = 'apf_levy_ech.fits'
        par['calibrations']['wavelengths']['lamps'] = ['ThAr_HARPS']
        par['calibrations']['wavelengths']['sigdetect'] = 5.0
        par['calibrations']['wavelengths']['fwhm'] = 2.5
        par['calibrations']['wavelengths']['rms_thresh_frac_fwhm'] = 0.2
        # Reidentification parameters
        par['calibrations']['wavelengths']['method'] = 'reidentify'
        par['calibrations']['wavelengths']['refframe'] = 'observed'
        # Echelle parameters
        par['calibrations']['wavelengths']['echelle'] = True
        par['calibrations']['wavelengths']['ech_nspec_coeff'] = 4
        par['calibrations']['wavelengths']['ech_norder_coeff'] = 4
        par['calibrations']['wavelengths']['ech_sigrej'] = 3.0

        par['calibrations']['flatfield']['slit_illum_finecorr'] = False
        par['calibrations']['flatfield']['tweak_slits'] = False
        par['calibrations']['flatfield']['spat_samp'] = 0.7
        par['calibrations']['flatfield']['slit_trim'] = 0


        # Processing steps
        turn_off = dict(use_biasimage=False,
                        use_darkimage=False,
                        use_illumflat=False)
        par.reset_all_processimages_par(**turn_off)
        # Do not correct for flexure
        par['flexure']['spec_method'] = 'skip'

        # no sky subtraction on standard stars
        par['reduce']['skysub']['global_sky_std'] = False

        # skip sky subtraction when searching for objects
        # this is because the sky subtraction is not very good with narrow
        # slits and the usual APF target is bright
        par['reduce']['findobj']['skip_skysub'] = True

        # there can be a 1" overlap between the orders, so we need to
        # trim the edges of the orders to avoid problems with the
        # extraction
        par['reduce']['findobj']['find_trim_edge'] = [2, 2]
        par['reduce']['findobj']['maxnumber_sci'] = 1
        par['reduce']['findobj']['maxnumber_std'] = 1
        par['reduce']['findobj']['find_fwhm'] = 4.0
        par['reduce']['findobj']['snr_thresh'] = 3.0

        par['reduce']['extraction']['boxcar_radius'] = 1.728 # 4 pixel radius
        par['reduce']['skysub']['mask_by_boxcar'] = True

        return par


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

        binning = "1,1"
        if hdu:
            # the CCD can only be binned 1x1 or 2x2
            # the square binning means both keywords will be the same
            # rbin means binning in the row, cbin is column (r is spatial)
            # finally, for reasons know only to god and Richard Stover,
            # 1 pixel binning is 0 and 2 pixel binning is 1 maybe this is
            # number of extra pixels on the detector?
            binning = self.get_meta_value(self.get_headarr(hdu), 'binning')

        detector_dict = dict(
            binning=binning,
            det=det,
            dataext=0,
            specaxis=0,
            specflip=True,
            spatflip=True,
            platescale=0.432, # SV made a very fast camera and the instrument takes a f/3 beam
            saturation=65535.,
            mincounts=-1e10,
            nonlinear=0.99, # the full well is like 300k and the gain is 1.031
            numamplifiers=1,
            gain=np.asarray([1.031]),
            ronoise=np.asarray([3.75]),
            xgap=0.,
            ygap=0.,
            ysize=1.,
            darkcurr=0.0008,
            # These are rows, columns on the raw frame, 1-indexed
            datasec=np.asarray(['[:, 1:2048]']),
            oscansec=np.asarray(['[:, 2049:2080]']),  # oscan is in the spatial direction
        )
        return detector_container.DetectorContainer(**detector_dict)

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
            time = headarr[0]['DATE-BEG']
            ttime = Time(time, format='isot')
            return ttime.mjd

        if meta_key == 'decker':
            decker_str = headarr[0]['DECKRNAM']
            if ":8" in decker_str:
                return '8.0'
            elif ":3" in decker_str:
                return '3.0'
            elif "Pinhole" in decker_str:
                return 'Pinhole'
            else:
                msgs.error(f"Unrecognized decker {decker_str}")

        if meta_key == 'binning':
            return f"{headarr[0]['RBIN']+1},{headarr[0]['CBIN']+1}"

        msgs.error("Not ready for this compound meta")

    def configuration_keys(self):
        """
        Return the default parameters to use for this instrument.
        
        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of ``PypeIt`` methods.
        """

        return ['binning', 'decker']

    def config_independent_frames(self):
        """
        Define frame types that are independent of the fully defined
        instrument configuration.

        This method returns a dictionary where the keys of the dictionary are
        the list of configuration-independent frame types. The value of each
        dictionary element can be set to one or more metadata keys that can
        be used to assign each frame type to a given configuration group. See
        :func:`~pypeit.metadata.PypeItMetaData.set_configurations` and how it
        interprets the dictionary values, which can be None.

        Returns:
            :obj:`dict`: Dictionary where the keys are the frame types that
            are configuration-independent and the values are the metadata
            keywords that can be used to assign the frames to a configuration
            group.
        """
        return {'bias': 'binning', 'dark': 'binning', 'pixelflat': 'binning',
                'pinhole': 'binning', 'illumflat': 'binning'}


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
        # TODO: Fit something
        # Current values are
        # Order Value
        # 58 0.43346
        # 66 0.43767
        # 77 0.43551
        # 93 0.42944
        # 108 0.42552
        # 124 0.43146

        if binning:
            _,bin_spat = binning.split(",")
            bin_spat = float(bin_spat)
        else:
            bin_spat = 1.0
        pscale = self.get_detector_par(1).platescale
        pscale = pscale * bin_spat
        return np.full(order_vec.size, pscale)

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the ``PypeIt``-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=0, card='RA')
        self.meta['dec'] = dict(ext=0, card='DEC')
        self.meta['target'] = dict(ext=0, card='TOBJECT')
        self.meta['decker'] = dict(ext=0, card=None, compound=True)
        self.meta['dispname'] = dict(ext=0, card=None, default='default')
        self.meta['mjd'] = dict(ext=0, card=None, compound=True)
        self.meta['binning'] = dict(ext=0, card=None, compound=True)

        self.meta['instrument'] = dict(ext=0, card='VERSION')
        self.meta['idname'] = dict(ext=0, card='OBJECT')
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
#        self.meta['dispname'] = dict(ext=0, card='ECHNAME')
        # Extras for config and frametyping

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

        # 'science' category
        if ftype == 'science':
            return good_exp & self.is_science(fitstbl)
        if ftype == 'bias':
            return good_exp & (fitstbl['idname'] == 'Bias')
        if ftype == 'dark':
            return good_exp & (fitstbl['idname'] == 'Dark')
        if ftype in ['pixelflat']:
            return good_exp & (fitstbl['idname'] == 'WideFlat')
        if ftype in ['trace']:
            return good_exp & (fitstbl['idname'] == 'NarrowFlat') & (fitstbl['decker'] == '3.0') \
                        | (fitstbl['idname'] == 'WideFlat') & (fitstbl['decker'] == '8.0')
        if ftype in ['tilt']:
            return good_exp & (fitstbl['idname'] == 'ThAr') & \
                ((fitstbl['decker'] == '8.0') | (fitstbl['decker'] == '3.0'))
        if ftype in ['arc']:
            return good_exp & (fitstbl['idname'] == 'ThAr')
        if ftype in ['pinhole']:
            return good_exp & (fitstbl['idname'] == 'NarrowFlat') & (fitstbl['decker'] == 'Pinhole')

        msgs.warn(f'Cannot determine if frames are of type {ftype}.')
        return np.zeros(len(fitstbl), dtype=bool)

    def is_science(self, fitstbl):
        """
        Return a boolean array selecting science frames.
        """

        return np.logical_not(np.isin(fitstbl['idname'], ['WideFlat', 'NarrowFlat', \
                                                          'ThAr', 'Dark', 'Bias', 'Iodine']))

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
        decker = self.get_meta_value(scifile, 'decker')

        if decker == '3.0':
            par['reduce']['trim_edge'] = [0, 0]
            par['reduce']['extraction']['model_full_slit'] = True
            par['reduce']['skysub']['no_local_sky'] = True
            par['reduce']['findobj']['find_trim_edge'] = [0, 0]
            par['calibrations']['slitedges']['pad'] = 5
            par['reduce']['extraction']['sn_gauss'] = 400
            # basically always use the Gaussian model for optimal extraction

        binning = self.get_meta_value(scifile, 'binning')

        if binning == "2,2":
            par['calibrations']['slitedges']['min_edge_side_sep'] = 2.0

        return par

    @property
    def norders(self):
        """
        Number of orders for this spectograph. Should only defined for
        echelle spectrographs, and it is undefined for the base class.
        """
        return 55

    @property
    def orders(self):
        """
        Return the order number for each echelle order.
        """
        return np.arange(122, 67, -1, dtype=int)

    @property
    def order_spat_pos(self):
        """
        Return the expected spatial position of each echelle order.
        """
        ord_spat_pos = np.array(
        #[0.06748562, 0.08327305, 0.09882858, 0.11422751, 0.1294643 ,
        [0.11422751, 0.1294643 ,
        0.14456395, 0.15949527, 0.17425991, 0.18887026, 0.20333023,
        0.21763396, 0.23178534, 0.2457807 , 0.25963618, 0.27333362,
        0.28688716, 0.30029984, 0.31356982, 0.32669448, 0.33967947,
        0.35252741, 0.36523925, 0.37781168, 0.3902439 , 0.40254701,
        0.41471288, 0.42674351, 0.43865832, 0.45044573, 0.46211223,
        0.47365838, 0.48506905, 0.4963403 , 0.50754044, 0.51861294,
        0.52956091, 0.54040772, 0.55114714, 0.56177537, 0.57228576,
        0.58268098, 0.59298235, 0.60318236, 0.61331737, 0.62333595,
        0.63325607, 0.64306738, 0.65280306, 0.66246134, 0.672025  ,
        0.68153262, 0.69093696, 0.70028011, 0.70953112, 0.71869543,
        0.72783612, 0.73692905, 0.74593691])
        #0.75489132, 0.76384089,
        #0.77268862, 0.781545  , 0.79036394, 0.79922398, 0.80807998,
        #0.81680164])

        return ord_spat_pos

    def get_rawimage(self, raw_file, det, sec_includes_binning=True):
        """ Read the image
        """
        # Check for file; allow for extra .gz, etc. suffix
        if not os.path.isfile(raw_file):
            msgs.error(f'{raw_file} not found!')
        hdu = io.fits_open(raw_file)

        head0 = hdu[0].header

        datasec = head0['DATASEC']
        datasec = datasec[1:-1] # trim [ ]
        xs , ys = datasec.split(",")
        yb, ye = ys.split(":")
        xb, xe = xs.split(":")
        xb = int(xb) - 1
        yb = int(yb) - 1
        xe = int(xe)
        ye = int(ye)


        # Grab the data
        full_image = hdu[0].data.astype(float)
        rawdatasec_img = np.zeros_like(full_image, dtype=int)
        oscansec_img = np.zeros_like(full_image, dtype=int)

        # Data
        rawdatasec_img[yb:ye, xb:xe] = 1
        # Overscan
        oscansec_img[yb:ye, xe:] = 1


        return self.get_detector_par(1, hdu=hdu), \
                full_image, hdu, head0['EXPTIME'], rawdatasec_img, oscansec_img
