"""
Implements HIRES-specific functions, including reading in slitmask design
files.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import glob

from IPython import embed

import numpy as np

from scipy import interpolate

from pypeit import msgs
from pypeit import telescopes
from pypeit import io
from pypeit.core import parse
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph


class APFLevySpectrograph(spectrograph.Spectrograph):
    """
    Child to handle APF/Levy specific code.

    This spectrograph is not yet supported.
    """
    ndet = 1
    telescope = telescopes.APFTelescopePar()
    pypeline = 'Echelle'


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
        return ['dispname', 'decker', 'binning']

 

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
        self.meta['target'] = dict(ext=0, card='OBJECT')
        self.meta['decker'] = dict(ext=0, card='DECKRNAM')
        self.meta['binning'] = dict(card=None, compound=True)

        self.meta['mjd'] = dict(card=None, compound=True)
        self.meta['idname'] = dict(ext=0, card='OBJECT')
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        # self.meta['dispname'] = dict(ext=0, card='ECHNAME')
        # Extras for config and frametyping
        
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
            binspec = headarr[0]['RBIN']
            binspatial = headarr[0]['CBIN']
            return parse.binning2string(binspec, binspatial)
        elif meta_key == 'mjd':
            ttime = Time(headarr[1]['DATE-OBS'], format='isot')
            return ttime.mjd
        else:
            msgs.error("Not ready for this compound meta")

    def get_detector_par(self, det, hdu=None):
        """
        Return metadata for the selected detector.

        Args:
            det (:obj:`int`):
                1-indexed detector number.  This is not used because NIRES only
                has one detector!
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
            dataext         = 0,
            specaxis        = 0,
            specflip        = False,
            spatflip=False,
            platescale      = 0.39,
            darkcurr        = 0.0008,
            saturation      = 65535, 
            nonlinear       = 0.,
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d(1.03),
            ronoise         = np.atleast_1d(3.75),
            datasec         = np.atleast_1d('[:,0:2048]'),
            oscansec        = np.atleast_1d('[:,2052:2080]')
            )
        return detector_container.DetectorContainer(**detector_dict)

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
        if ftype == 'bias':
            return good_exp & (fitstbl['idname'] == 'Bias')
        if ftype == 'dark':
            return good_exp & (fitstbl['idname'] == 'Dark')
        if ftype in ['pixelflat']:
            return good_exp & (fitstbl['idname'] == 'WideFlat') 
        if ftype in ['trace']:
            return good_exp & (fitstbl['idname'] == 'NarrowFlat')
        if ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['idname'] == 'ThAr')
        if ftype == 'science':
            return good_exp 

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)


    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.
        
        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of ``PypeIt`` methods.
        """
        par = super().default_pypeit_par()
        
        # Adjustments to parameters from Xshooter-VIS
        # we have pixel flats normally, but using them will be a trick I think
        # Processing steps
         
        turn_on = dict(use_biasimage=False, use_overscan=True, overscan_method='median',
                       use_darkimage=False, use_illumflat=False, use_pixelflat=False,
                       use_specillum=False)
        par.reset_all_processimages_par(**turn_on)

        par['calibrations']['traceframe']['process']['overscan_method'] = 'median'
        par['scienceframe']['process']['use_biasimage']=True
        par['scienceframe']['process']['use_illumflat']=True
        par['scienceframe']['process']['use_pixelflat']=True
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
        par['calibrations']['wavelengths']['lamps'] = ['ThAr_lines.dat']

        par['calibrations']['wavelengths']['rms_threshold'] = 0.50
        par['calibrations']['wavelengths']['sigdetect'] = 5.0
        par['calibrations']['wavelengths']['n_final'] = [3] + 13*[4] + [3]

        par['calibrations']['wavelengths']['fwhm'] = 2.0

        par['calibrations']['wavelengths']['method'] = 'reidentify'
        
        # par['calibrations']['wavelengths']['reid_arxiv'] = 'vlt_xshooter_vis1x1.fits'
        par['calibrations']['wavelengths']['cc_thresh'] = 0.50
        par['calibrations']['wavelengths']['cc_local_thresh'] = 0.50
        par['calibrations']['wavelengths']['ech_fix_format'] = True
        
        # Echelle parameters
        par['calibrations']['wavelengths']['echelle'] = True
        par['calibrations']['wavelengths']['ech_nspec_coeff'] = 4
        par['calibrations']['wavelengths']['ech_norder_coeff'] = 4
        par['calibrations']['wavelengths']['ech_sigrej'] = 3.0
        
        # Flats
        par['calibrations']['flatfield']['tweak_slits_thresh'] = 0.90
        par['calibrations']['flatfield']['tweak_slits_maxfrac'] = 0.10

        # local sky subtraction operates on entire slit
        # par['reduce']['extraction']['model_full_slit'] = True
        # Mask 3 edges pixels since the slit is short, insted of default (5,5)
        par['reduce']['findobj']['find_trim_edge'] = [3,3]
        # Continnum order for determining thresholds
        par['reduce']['findobj']['find_npoly_cont'] = 0
        # Don't attempt to fit a continuum to the trace rectified image
        par['reduce']['findobj']['find_cont_fit'] = False

        par['reduce']['extraction']['boxcar_radius'] = 2.0  # arcsec

        # Set the default exposure time ranges for the frame typing
        par['calibrations']['standardframe']['exprng'] = [None, 60]
        par['calibrations']['arcframe']['exprng'] = [10, None]
        par['calibrations']['tiltframe']['exprng'] = [10, None]
        par['calibrations']['darkframe']['exprng'] = [10, None]
        par['scienceframe']['exprng'] = [5, None]


        return par

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

        msgs.info("Using hard-coded BPM for APF/Levy")
        bpm_img[:, 0] = 1

        return bpm_img

    

