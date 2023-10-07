"""
Implements APF-specific functions

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import glob

from IPython import embed

import numpy as np
from astropy.time import Time

from pypeit import msgs
from pypeit import telescopes
from pypeit import io
from pypeit.core import parse
from pypeit.core import framematch
from pypeit.par import pypeitpar
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

        par['calibrations']['tilts']['tracethresh'] = 20
        # Bias


        # 1D wavelength solution
        par['calibrations']['wavelengths']['lamps'] = ['ThAr']
        par['calibrations']['wavelengths']['rms_threshold'] = 0.25
        par['calibrations']['wavelengths']['sigdetect'] = 5.0
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
        self.meta['decker'] = dict(ext=0, card='DECKRNAM')
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
        if ftype in ['pixelflat','illumflat']:
            # Flats and trace frames are typed together
            return good_exp & (fitstbl['idname'] == 'WideFlat') 
        if ftype in ['trace']:
            return good_exp & ((fitstbl['idname'] == 'WideFlat') |
                                   (fitstbl['idname'] == 'Iodine'))
        if ftype in ['pinhole']:
            return good_exp & (fitstbl['idname'] == 'NarrowFlat') 
        if ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['idname'] == 'ThAr')

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

    def is_science(self, fitstbl):
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
        return np.arange(60, 125, dtype=int)

    @property
    def order_spat_pos(self):
        """
        Return the expected spatial position of each echelle order.
        """
        ord_spat_pos = np.array([0.18358168, 0.19230169, 0.20103488, 0.20985631, 0.21863824,
                                 0.22747433, 0.23634767, 0.24524769, 0.25419597, 0.26322385,
                                 0.27228656, 0.28142214, 0.29059795, 0.29985087, 0.30918132,
                                 0.31858318, 0.32809034, 0.33764746, 0.34730427, 0.35703641,
                                 0.36684215, 0.37676406, 0.38676816, 0.39689161, 0.40708169,
                                 0.41738262, 0.42777566, 0.43828581, 0.44890403, 0.45963836,
                                 0.47047902, 0.48142425, 0.49248485, 0.50366581, 0.51494656,
                                 0.52635061, 0.53787141, 0.54955084, 0.56133274, 0.57323834,
                                 0.58526665, 0.59743135, 0.60972577, 0.62215249, 0.63471998,
                                 0.64742603, 0.66026556, 0.67323877, 0.68636125, 0.69962374,
                                 0.71302498, 0.7265794 , 0.7402736 , 0.75412483, 0.76811608,
                                 0.78228211, 0.79658665, 0.81104268, 0.8256534 , 0.84044735,
                                 0.85533756, 0.87046297, 0.88571258, 0.90111377, 0.91669895])

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

