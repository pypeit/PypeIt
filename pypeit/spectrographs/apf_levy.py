"""
Implements APF-specific functions

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import numpy as np
from astropy.time import Time

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container


class APFLevySpectrograph(spectrograph.Spectrograph):
    """
    Child to handle APF specific code.

    This spectrograph is not yet supported.
    """
    ndet = 1
    telescope = telescopes.APFTelescopePar()
    pypeline = 'Echelle'
    name = 'apf_levy'
    camera = 'apf'
    header_name = 'apf'
    ech_fixed_format = True

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.
        
        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of ``PypeIt`` methods.
        """
        par = super().default_pypeit_par()

        par['calibrations']['slitedges']['edge_thresh'] = 10.
        par['calibrations']['slitedges']['fit_order'] = 5
        par['calibrations']['slitedges']['max_shift_adj'] = 0.5
        par['calibrations']['slitedges']['left_right_pca'] = True
        par['calibrations']['slitedges']['smash_range'] = [0.4,0.6]

        par['calibrations']['tilts']['tracethresh'] = 20
        # Bias


        # 1D wavelength solution
        par['calibrations']['wavelengths']['lamps'] = ['ThAr']
        par['calibrations']['wavelengths']['rms_threshold'] = 0.25
        par['calibrations']['wavelengths']['sigdetect'] = 5.0
        par['calibrations']['wavelengths']['fwhm'] = 2.0
        # Reidentification parameters
        #par['calibrations']['wavelengths']['method'] = 'reidentify'
        #par['calibrations']['wavelengths']['ech_fix_format'] = True
        # Echelle parameters
        par['calibrations']['wavelengths']['echelle'] = True
        par['calibrations']['wavelengths']['ech_nspec_coeff'] = 4
        par['calibrations']['wavelengths']['ech_norder_coeff'] = 4
        par['calibrations']['wavelengths']['ech_sigrej'] = 3.0


        # Processing steps
        turn_off = dict(use_biasimage=False,
                        use_darkimage=False)
        par.reset_all_processimages_par(**turn_off)
        # Do not correct for flexure
        par['flexure']['spec_method'] = 'skip'

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
        detector_dict = dict(
            binning='1,1',
            det=1,
            dataext=0,
            specaxis=0,
            specflip=True,
            spatflip=True,
            platescale=0.39, # SV made a very fast camera and the instrument takes a f/3 beam
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
        # its a fixed format spectrometer
        # different deckers are used for different kinds of calibrations
        # we will treat deckers separately
        return ['decker']

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
        plate_scale = np.zeros_like(order_vec)
        plate_scale += (0.43346 + 0.43767 + 0.43551 + 0.42944 + 0.42552 + 0.43146)/6.0

        return plate_scale

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
        self.meta['binning'] = dict(ext=0, card=None, default='1,1')
        self.meta['dispname'] = dict(ext=0, card=None, default='default')
        self.meta['mjd'] = dict(ext=0, card=None, compound=True)

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
            # Flats and trace frames are typed together
            return good_exp & (fitstbl['idname'] == 'WideFlat')
        if ftype in ['trace', 'illumflat']:
            return good_exp & ((fitstbl['idname'] == 'WideFlat') |
                                   (fitstbl['idname'] == 'NarrowFlat'))
        if ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['idname'] == 'ThAr')

        msgs.warn(f'Cannot determine if frames are of type {ftype}.')
        return np.zeros(len(fitstbl), dtype=bool)

    def is_science(self, fitstbl):
        """
        Return a boolean array selecting science frames.
        """

        rv = fitstbl['idname'] != 'WideFlat'

        for filetype in ['NarrowFlat','ThAr','Dark','Bias','Iodine']:
            rv = rv & (fitstbl['idname'] != filetype)

        return rv

    @property
    def norders(self):
        """
        Number of orders for this spectograph. Should only defined for
        echelle spectrographs, and it is undefined for the base class.
        """
        return 65

    @property
    def orders(self):
        """
        Return the order number for each echelle order.
        """
        return np.arange(125, 60, dtype=int)

    @property
    def order_spat_pos(self):
        """
        Return the expected spatial position of each echelle order.
        """
        ord_spat_pos = np.array([0.0983199235264478, 0.11375312819654312,
                                0.12897539482511905, 0.1440953815504298,
                                0.15900762855967568, 0.17377299180935876,
                                0.18838232501574215, 0.20284463310496173,
                                0.21714890931272507, 0.23129979938915052,
                                0.24529189806355084, 0.25914692702252085,
                                0.2728407669050262, 0.28640472166394293,
                                0.2998128953546085, 0.31308216800464767,
                                0.32620673184319937, 0.33919142310076345,
                                0.352043519170773, 0.3647491532650974,
                                0.377320973331427, 0.38975546681481527,
                                0.4020557752605035, 0.41422063996111946,
                                0.4262569359010734, 0.43817113193728346,
                                0.44995809937279335, 0.46162488123623735,
                                0.47317467844888994, 0.4845794119440806,
                                0.49585022440661947, 0.5070556525631347,
                                0.5181287361328383, 0.5290728748588746,
                                0.5399230762066759, 0.5506557239902183,
                                0.5612878806046232, 0.5717958659670593,
                                0.5821946454493663, 0.5924898391530614,
                                0.602694572800143, 0.6128232035249103,
                                0.6228467082726102, 0.6327686913140586,
                                0.6425759357758352, 0.6523113583754457,
                                0.661967513378642, 0.6715374133541441,
                                0.681038638586849, 0.6904473777247495,
                                0.6997878996249236, 0.709041703663812,
                                0.7182069577174326, 0.7273478239726879,
                                0.7364428418264531, 0.745445883042169,
                                0.7544037473189316, 0.7633508225256459,
                                0.7722047243309123, 0.7810504844199905,
                                0.7898802294237828, 0.7987155388555277,
                                0.8076016765076058, 0.8161601047976773,
                                0.8248325393949485])

        return ord_spat_pos


# def apf_read_chip(hdu):
#     """ Read the APF detector

#     Parameters
#     ----------
#     hdu : HDUList

#     Returns
#     -------
#     data : ndarray
#     oscan : ndarray
#     """

#     # Extract datasec from header
#     datsec = hdu[0].header['DATASEC']
#     postpix = hdu[0].header['COVER']


#     x1_dat, x2_dat, y1_dat, y2_dat = np.array(parse.load_sections(datsec)).flatten()
#     x1_det, x2_det, y1_det, y2_det = np.array(parse.load_sections(detsec)).flatten()

#     # This rotates the image to be increasing wavelength to the top
#     #data = np.rot90((hdu[0].data).T, k=2)
#     #nx=data.shape[0]
#     #ny=data.shape[1]

#     # Science data
#     fullimage = hdu[0].data
#     data = fullimage[x1_dat:x2_dat,y1_dat:y2_dat]

#     # Overscan
#     oscan = fullimage[:,y2_dat:]

#     # Return
#     return data, oscan
