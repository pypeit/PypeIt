"""
Module for the SOAR/Goodman instrument

.. include:: ../include/links.rst
"""
from IPython.terminal.embed import embed
from pkg_resources import resource_filename

import numpy as np

from astropy.time import Time

from pypeit import msgs
from pypeit import telescopes
from pypeit import io
from pypeit.core import framematch
from pypeit.spectrographs import spectrograph
from pypeit.core import parse
from pypeit.images import detector_container


class SOARGoodmanSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Goodman specific code for each camera
    """
    ndet = 1
    telescope = telescopes.SOARTelescopePar()

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
        return ['dispname', 'decker', 'binning', 'dispangle']

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the ``PypeIt``-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=1, card='RA')
        self.meta['dec'] = dict(ext=1, card='DEC')
        self.meta['target'] = dict(ext=1, card='OBJECT')
        self.meta['decker'] = dict(ext=1, card='SLIT')
        self.meta['binning'] = dict(card=None, compound=True)
        self.meta['exptime'] = dict(ext=1, card='EXPTIME')
        self.meta['mjd'] = dict(card=None, compound=True)
        self.meta['airmass'] = dict(ext=1, card='AIRMASS')
        # Extras for config and frametyping
        self.meta['dispname'] = dict(ext=1, card='GRATING')
        self.meta['dispangle'] = dict(ext=1, card='GRT_ANG', rtol=1e-3)
        self.meta['idname'] = dict(ext=1, card='OBSTYPE')
        # used for arc and continuum lamps
        self.meta['lampstat01'] = dict(ext=1, card='LAMP_HGA')
        self.meta['lampstat02'] = dict(ext=1, card='LAMP_NE')
        self.meta['lampstat03'] = dict(ext=1, card='LAMP_AR')
        self.meta['lampstat04'] = dict(ext=1, card='LAMP_FE')
        self.meta['lampstat05'] = dict(ext=1, card='LAMP_CU')
        self.meta['lampstat06'] = dict(ext=1, card='LAMP_QUA')
        self.meta['lampstat07'] = dict(ext=1, card='LAMP_BUL')
        self.meta['lampstat08'] = dict(ext=1, card='LAMP_DOM')

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
            binspec, binspatial = [int(item) for item in headarr[1]['CCDSUM'].split(' ')]
            return parse.binning2string(binspec, binspatial)
        elif meta_key == 'mjd':
            ttime = Time(headarr[1]['DATE-OBS'], format='isot')
            return ttime.mjd
        else:
            msgs.error("Not ready for this compound meta")

#    def pypeit_file_keys(self):
#        """
#        Define the list of keys to be output into a standard ``PypeIt`` file.
#
#        Returns:
#            :obj:`list`: The list of keywords in the relevant
#            :class:`~pypeit.metadata.PypeItMetaData` instance to print to the
#            :ref:`pypeit_file`.
#        """
#        return super().pypeit_file_keys() + ['slitwid']

class SOARGoodmanRedSpectrograph(SOARGoodmanSpectrograph):
    name = 'soar_goodman_red'
    camera = 'red'
    comment = 'Supported gratings: M1, M2 and 2x2 binning'
    supported = True

    def get_detector_par(self, hdu, det):
        """
        Return metadata for the selected detector.

        Args:
            hdu (`astropy.io.fits.HDUList`_):
                The open fits file with the raw image of interest.
            det (:obj:`int`):
                1-indexed detector number.

        Returns:
            :class:`~pypeit.images.detector_container.DetectorContainer`:
            Object with the detector metadata.
        """
        header = hdu[1].header
        # Binning
        binning = self.get_meta_value(self.get_headarr(hdu), 'binning')  # Could this be detector dependent??

        # Detector 1
        detector_dict = dict(
            binning         = binning,
            det             = 1,
            dataext         = 1,
            specaxis        = 1,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.15,
            darkcurr        = 0.00008,  # e-/s/pix
            saturation      = 65535.,
            nonlinear       = 1.0,
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d(header['GAIN']),
            ronoise         = np.atleast_1d(header['RDNOISE']),
            )
        
        # Only tested for 2x2
        if binning == '2,2':
            # parse TRIMSEC
            col0 = int(header['TRIMSEC'][1:].split(':')[0])
            dsec = f"[:,{col0*2}:]"  # rows, columns on the raw frame
            detector_dict['datasec'] = np.atleast_1d(dsec)
            # Overscan
            osec = f"[:,1:{int(col0*2)-2}:]"
            detector_dict['oscansec'] = np.atleast_1d(osec)
        else:
            msgs.error("Ask the developers to add your binning.  Or add it yourself.")
        
        # Return
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

        # Turn off bias and turn on overscan
        turn_off_on = dict(use_biasimage=False, use_darkimage=False, use_overscan=True)
        par.reset_all_processimages_par(**turn_off_on)

        # Ignore PCA
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'

        # Set pixel flat combination method
        #par['calibrations']['pixelflatframe']['process']['combine'] = 'median'
        # Change the wavelength calibration method
        par['calibrations']['wavelengths']['method'] = 'holy-grail'
        #par['calibrations']['wavelengths']['method'] = 'reidentify'
        par['calibrations']['wavelengths']['lamps'] = ['NeI', 'ArI', 'HgI']
        # Wavelengths
        # 1D wavelength solution
        par['calibrations']['wavelengths']['rms_threshold'] = 0.5
        par['calibrations']['wavelengths']['sigdetect'] = 5.
        par['calibrations']['wavelengths']['fwhm']= 5.0

        #par['calibrations']['wavelengths']['n_first'] = 3
        #par['calibrations']['wavelengths']['n_final'] = 5
        #par['calibrations']['wavelengths']['sigdetect'] = 10.0
        #par['calibrations']['wavelengths']['wv_cen'] = 4859.0
        #par['calibrations']['wavelengths']['disp'] = 0.2

        # Set the default exposure time ranges for the frame typing
        #par['calibrations']['biasframe']['exprng'] = [None, 1]
        #par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        #par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['arcframe']['exprng'] = [None, 30]
        par['calibrations']['standardframe']['exprng'] = [None, 120]
        par['scienceframe']['exprng'] = [90, None]

        # Extraction
        #par['reduce']['skysub']['bspline_spacing'] = 0.8
        #par['reduce']['skysub']['no_poly'] = True
        #par['reduce']['skysub']['bspline_spacing'] = 0.6
        #par['reduce']['skysub']['joint_fit'] = False
        #par['reduce']['skysub']['global_sky_std']  = False
#
#        par['reduce']['extraction']['sn_gauss'] = 4.0
#        par['reduce']['findobj']['sig_thresh'] = 5.0
#        par['reduce']['skysub']['sky_sigrej'] = 5.0
#        par['reduce']['findobj']['find_trim_edge'] = [5,5]

        # Sensitivity function parameters
        #par['sensfunc']['polyorder'] = 7

        # Do not correct for flexure
#        par['flexure']['spec_method'] = 'skip'

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

        msgs.info("Using hard-coded BPM for SOAR/Goodman")
        bpm_img[:, 0] = 1

        return bpm_img

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
        if ftype in ['science']:
            return good_exp & (fitstbl['idname'] == 'SPECTRUM')
        if ftype in ['standard']:
            # Don't type pinhole or dark frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype == 'bias':
            # Don't type bias
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            return good_exp & (fitstbl['idname'] == 'LAMPFLAT')
        if ftype in ['pinhole', 'dark']:
            # Don't type pinhole or dark frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['idname'] == 'ARC')
        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)



