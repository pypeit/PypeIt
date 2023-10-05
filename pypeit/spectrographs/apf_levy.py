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
        return []

    
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
            return good_exp & (fitstbl['idname'] == 'WideFlat') 
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
        ord_spat_pos = np.array([0.078125, 0.13769531, 0.19189453, 0.24414062, 0.29296875,
                                 0.34179688, 0.38330078, 0.42724609, 0.46582031, 0.50439453,
                                 0.54199219, 0.57763672, 0.61279297, 0.6484375 , 0.68457031,
                                 0.71875   , 0.75439453, 0.79443359, 0.83789062, 0.88671875,
                                 0.94091797])
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

