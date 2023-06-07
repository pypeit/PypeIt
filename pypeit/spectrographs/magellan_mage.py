"""
Module for Magellan/MAGE specific methods.

.. include:: ../include/links.rst
"""

from IPython import embed

import numpy as np

from astropy.time import Time

from pypeit import msgs
from pypeit import telescopes
from pypeit import io
from pypeit.core import framematch
from pypeit.core import parse
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container

class MagellanMAGESpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Magellan/MAGE specific code
    """
    ndet = 1
    name = 'magellan_mage'
    camera = 'MagE'
    url = 'https://www.lco.cl/?epkb_post_type_1=mage'
    header_name = 'MagE'
    telescope = telescopes.MagellanTelescopePar()
    pypeline = 'Echelle'
    ech_fixed_format = True
    supported = True
    comment = 'See :doc:`mage`'

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
        # TODO: Could this be detector dependent??
        binning = '1,1' if hdu is None else self.get_meta_value(self.get_headarr(hdu), 'binning')

        # Detector 1
        detector_dict = dict(
            binning         = binning,
            det             = 1,
            dataext         = 0,
            specaxis        = 1,
            specflip        = True,
            spatflip        = False,
            # plate scale in arcsec/pixel
            platescale      = 0.3,
            # electrons/pixel/hour. From: http://www.lco.cl/telescopes-information/magellan/instruments/mage/the-mage-spectrograph-user-manual
            darkcurr        = 1.00,
            saturation      = 65535.,
            # CCD is linear to better than 0.5 per cent up to digital saturation (65,536 DN including bias) in the Fast readout mode.
            nonlinear       = 0.99,
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d(1.02), # depends on the readout
            ronoise         = np.atleast_1d(2.9), # depends on the readout
            datasec         = np.atleast_1d('[1:1024, 1:2048]'),
            oscansec        = np.atleast_1d('[1:1024, 2049:2176]'),
            )
        # Taken from the MASE paper: https://arxiv.org/pdf/0910.1834.pdf
        #self.norders = 15
        # 20-6
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

        # Bias
        #par['calibrations']['biasframe']['useframe'] = 'overscan'
        # Wavelengths
        # 1D wavelength solution
        par['calibrations']['wavelengths']['rms_threshold'] = 0.20  # Might be grating dependent..
        par['calibrations']['wavelengths']['sigdetect'] = 5.0
        par['calibrations']['wavelengths']['lamps'] = ['ThAr_MagE']

        par['calibrations']['wavelengths']['method'] = 'reidentify'
        par['calibrations']['wavelengths']['cc_thresh'] = 0.50
        par['calibrations']['wavelengths']['cc_local_thresh'] = 0.50

        # Reidentification parameters
        par['calibrations']['wavelengths']['reid_arxiv'] = 'magellan_mage.fits'
#        par['calibrations']['wavelengths']['ech_fix_format'] = True
        # Echelle parameters
        par['calibrations']['wavelengths']['echelle'] = True
        par['calibrations']['wavelengths']['ech_nspec_coeff'] = 4
        par['calibrations']['wavelengths']['ech_norder_coeff'] = 4
        par['calibrations']['wavelengths']['ech_sigrej'] = 3.0

        par['scienceframe']['process']['sigclip'] = 20.0
        par['scienceframe']['process']['satpix'] = 'nothing'

        # Set slits and tilts parameters
        par['calibrations']['tilts']['tracethresh'] = 10. #[10]*self.norders
        par['calibrations']['slitedges']['fit_order'] = 5
        par['calibrations']['slitedges']['max_shift_adj'] = 3.
        par['calibrations']['slitedges']['edge_thresh'] = 10.  # Tough to get the bluest orders
        par['calibrations']['slitedges']['left_right_pca'] = True
        par['calibrations']['slitedges']['fit_min_spec_length'] = 0.3  # Allow for a short detected blue order
        # Find object parameters
        par['reduce']['findobj']['find_trim_edge'] = [4,4]    # Slit is too short to trim 5,5 especially with 2x binning
        par['reduce']['findobj']['maxnumber_sci'] = 2  # Slit is narrow so allow one object per order
        par['reduce']['findobj']['maxnumber_std'] = 1  # Slit is narrow so allow one object per order
        par['reduce']['extraction']['model_full_slit'] = True  # local sky subtraction operates on entire slit


        # Always flux calibrate, starting with default parameters
        # Do not correct for flexure
        par['flexure']['spec_method'] = 'skip'
        # Set the default exposure time ranges for the frame typing
        par['calibrations']['standardframe']['exprng'] = [None, 20]
        par['calibrations']['arcframe']['exprng'] = [20, None]
        par['calibrations']['darkframe']['exprng'] = [20, None]
        par['scienceframe']['exprng'] = [20, None]
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
        self.meta['target'] = dict(ext=0, card='OBJECT')
        #TODO: Check decker is correct
        self.meta['decker'] = dict(ext=0, card='SLITNAME')
        self.meta['binning'] = dict(card=None, compound=True)
#        self.meta['binning'] = dict(ext=0, card='BINNING')
        self.meta['mjd'] = dict(ext=0, card=None, compound=True)
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Extras for config and frametyping
        self.meta['dispname'] = dict(ext=0, card='INSTRUME')
        self.meta['idname'] = dict(ext=0, card='EXPTYPE')
        self.meta['instrument'] = dict(ext=0, card='INSTRUME')

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
            binspatial, binspec = parse.parse_binning(headarr[0]['BINNING'])
            return parse.binning2string(binspec, binspatial)
        elif meta_key == 'mjd':
            time = '{:s}T{:s}'.format(headarr[0]['UT-DATE'], headarr[0]['UT-TIME'])
            ttime = Time(time, format='isot')
            return ttime.mjd
        else:
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
        return []

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
        if ftype in ['pinhole', 'dark']:
            # No pinhole or bias or dark frames
            return np.zeros(len(fitstbl), dtype=bool)
        elif ftype in ['bias']:
            return fitstbl['idname'] == 'Bias'
        elif ftype in ['pixelflat', 'trace']:
            return fitstbl['idname'] == 'Flat'
        elif ftype in ['arc']:
            return fitstbl['idname'] == 'ThAr-Lamp'
        else:
            return (fitstbl['idname'] == 'Object') \
                        & framematch.check_frame_exptime(fitstbl['exptime'], exprng)

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

        # Get the binning
        msgs.info("Custom bad pixel mask for MAGE")
        hdu = io.fits_open(filename)
        binspatial, binspec = parse.parse_binning(hdu[0].header['BINNING'])
        hdu.close()
        # Do it
        bpm_img[:, :10//binspatial] = 1.  # Setting BPM on the edge of the detector often leads to false edges
        bpm_img[:, 1020//binspatial:] = 1.
        # Return
        return bpm_img

    @property
    def norders(self):
        """
        Number of orders for this spectograph. Should only defined for
        echelle spectrographs, and it is undefined for the base class.
        """
        return 12   # 20-6

    @property
    def order_spat_pos(self):
        """
        Return the expected spatial position of each echelle order.
        """
        ord_spat_pos =  np.array([0.316, 0.399, 0.475, 0.545, 0.609, 0.669, 0.723, 0.774, 0.823,
                                  0.869, 0.915, 0.965])
        return ord_spat_pos

    @property
    def orders(self):
        """
        Return the order number for each echelle order.
        """
        return  np.arange(17, 5, -1, dtype=int)

    @property
    def spec_min_max(self):
        """
        Return the minimum and maximum spectral pixel expected for the
        spectral range of each order.
        """
        spec_max = np.full(self.norders, np.inf)
        spec_min = np.full(self.norders, -np.inf)
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
        norders = len(order_vec)
        binspatial, binspec = parse.parse_binning(binning)
        return np.full(norders, 0.30*binspatial)


