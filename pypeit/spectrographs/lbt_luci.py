"""
Module for LBT/LUCI specific methods.

.. include:: ../include/links.rst
"""

from IPython import embed

import numpy as np

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container


class LBTLUCISpectrograph(spectrograph.Spectrograph):
    """
    Class to handle LBT/LUCI specific code
    """
    ndet = 1
    telescope = telescopes.LBTTelescopePar()
    url = 'https://scienceops.lbto.org/luci/'

#    def __init__(self):
#        super().__init__()
#        self.timeunit = 'isot'

#    @classmethod
#    def default_pypeit_par(cls):
#        """
#        Return the default parameters to use for this instrument.
#        
#        Returns:
#            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
#            all of ``PypeIt`` methods.
#        """
#        par = super().default_pypeit_par()
#
#        # Processing steps
#        turn_off = dict(use_illumflat=False, use_biasimage=False, use_overscan=False,
#                        use_darkimage=False)
#        par.reset_all_processimages_par(**turn_off)
#
#        par['calibrations']['biasframe']['exprng'] = [None, 1]
#        par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
#        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
#        par['calibrations']['pixelflatframe']['exprng'] = [0, None]
#        par['calibrations']['traceframe']['exprng'] = [0, None]
#        par['calibrations']['arcframe']['exprng'] = [None, 60]
#        par['calibrations']['standardframe']['exprng'] = [1, 200]
#        par['scienceframe']['exprng'] = [200, None]
#        return par

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the ``PypeIt``-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=0, card='OBJRA')
        self.meta['dec'] = dict(ext=0, card='OBJDEC')
        self.meta['target'] = dict(ext=0, card='OBJECT')
        self.meta['decker'] = dict(ext=0, card='MASKID')
        self.meta['binning'] = dict(ext=0, card=None, default='1,1')
        self.meta['filter1'] = dict(ext=0, card='FILTERS')
        self.meta['idname'] = dict(card=None, compound=True)
        self.meta['mjd'] = dict(ext=0, card='MJD-OBS')
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        self.meta['dispname'] = dict(ext=0, card='GRATNAME')
        self.meta['instrument'] = dict(ext=0, card='INSTRUME')

    # TODO: Deal with isot time here.
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
        # Populate the idname based on the header information of LUCI
        # This is an implicit way of pre-typing without adding too many
        # variables to the self.meta.
        if meta_key == 'idname':
            targetname = (headarr[0].get('OBJECT'))
            dispname = (headarr[0].get('GRATNAME'))
            calib_unit = (headarr[0].get('CALIB'))
            filter1 = (headarr[0].get('FILTER1'))
            filter2 = (headarr[0].get('FILTER2'))
            lamp1 = (headarr[0].get('LAMP1'))
            lamp2 = (headarr[0].get('LAMP2'))
            lamp3 = (headarr[0].get('LAMP3'))
            lamp4 = (headarr[0].get('LAMP4'))
            lamp5 = (headarr[0].get('LAMP5'))
            lamp6 = (headarr[0].get('LAMP6'))

            # object frame -> will be typed as science
            # This currently includes sky flats, science and standard images
            # We will guess standards using the beginning of their names.
            if ((dispname != 'Mirror') and
                (calib_unit == False) and
                (lamp1 == False) and
                (lamp2 == False) and
                (lamp3 == False) and
                (lamp4 == False) and
                (lamp5 == False) and
                (lamp6 == False)):
                if (targetname[:3] == 'HIP' or
                    targetname[:2] == 'HD' or
                    targetname[:5] == 'Feige'):
                    return 'standard'
                else:
                    return 'object'
            # flat frames -> will be typed as pixelflat, trace
            elif ((calib_unit == True) and
                  ((lamp4 == True) or
                   (lamp5 == True) or
                   (lamp6 == True))):
                return 'flat'
            # arcs -> will be typed as arc, tilt
            elif ((dispname != 'Mirror') and
                  (calib_unit == True) and
                  ((lamp1 == True) or
                   (lamp2 == True) or
                   (lamp3 == True))):
                return 'arc'
            # pixelflat off -> will be typed as bias
            elif ((dispname != 'Mirror') and
                (calib_unit == True) and
                (lamp1 == False) and
                (lamp2 == False) and
                (lamp3 == False) and
                (lamp4 == False) and
                (lamp5 == False) and
                (lamp6 == False) and
                (filter1 != 'blind') and
                (filter2 != 'blind')):
                return 'flat_off'
            # darks -> will not be typed currently
            elif ((filter1 == 'blind') or
                  (filter2 == 'blind')):
                return 'dark'
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
        return ['decker', 'dispname']

    def pypeit_file_keys(self):
        """
        Define the list of keys to be output into a standard ``PypeIt`` file.

        Returns:
            :obj:`list`: The list of keywords in the relevant
            :class:`~pypeit.metadata.PypeItMetaData` instance to print to the
            :ref:`pypeit_file`.
        """
        pypeit_keys = super().pypeit_file_keys()
        # TODO: Why are these added here? See
        # pypeit.metadata.PypeItMetaData.set_pypeit_cols
        # TODO: This should only add idname
        pypeit_keys += ['calib', 'comb_id', 'bkg_id', 'idname']
        return pypeit_keys

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
        # ATTENTION: Standards have to be added manually for LUCI because
        # there is not unique flag that allows to distinguish between targets
        # and standards
        if ftype in ['science']:
            return good_exp & (fitstbl['idname'] == 'object')
        if ftype in ['standard']:
            return good_exp & (fitstbl['idname'] == 'standard')
        if ftype == 'bias':
            # for NIR data we type off lamp flats as biases
            return good_exp & (fitstbl['idname'] == 'flat_off')
        if ftype in ['pixelflat', 'trace']:
            # Flats and trace frames are typed together
            return good_exp & (fitstbl['idname'] == 'flat')
        if ftype in ['dark']:
            # NOT Typing dark frames
            # return np.zeros(len(fitstbl), dtype=bool)
            # for testing dark typing uncommen the following line and comment
            # out the previous line
            return good_exp & (fitstbl['idname'] == 'dark')
        if ftype in ['arc', 'tilt']:
            return (good_exp & ((fitstbl['idname'] == 'object') |
                    (fitstbl['idname'] == 'arc')))

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)


class LBTLUCI1Spectrograph(LBTLUCISpectrograph):
    """
    Child to handle LBT/LUCI1 specific code
    """
    name = 'lbt_luci1'
    camera = 'LUCI1'
    header_name = 'LUCI1'
    supported = True

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
            binning         = '1,1',
            det             = 1,
            dataext         = 0,
            specaxis        = 1,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.25,
            # Dark current nominally is < 360 electrons per hours
            # but the dark subtraction will effectively bring this to 0
            darkcurr        = 0.03,
            # Saturation is 55000, but will be set to dummy value for
            # now
            saturation      = 1e+8,
            # NIR detectors are non-linear even in lower percentages
            # of the full well, thus for precision measurements one
            # should take into account a general non-linearity
            # correction.
            nonlinear       = 0.80,
            mincounts       = -1e10,
            # In fact there are 32 amplifiers, which gain and ronoise
            # are extremely similar to each other, thus it will be
            # mimicked as 1
            numamplifiers   = 1,
            # The readout noise for LUCI are different for
            # different readout modes. MER noise given here.
            gain= np.atleast_1d(2.21),
            ronoise         = np.atleast_1d(5.1),
            datasec=np.atleast_1d('[5:2044,5:2044]'),
            # For Luci the first 4 pixels on each side can
            # technically be used for as a biassec. This is not
            # included here.
            oscansec= np.atleast_1d('[5:2044,1:4]'),
            )
            
        detector = detector_container.DetectorContainer(**detector_dict)
        
        if hdu is None:
            return detector
        
        camera = hdu[0].header['CAMERA']
        if camera == 'N3.75 Camera':
            detector.platescale = 0.1178
        
        readmode = hdu[0].header['READMODE']
        if readmode == 'LIR':
            detector.ronoise = np.atleast_1d(9.6)
            
        return detector

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.
        
        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of ``PypeIt`` methods.
        """
        par = super().default_pypeit_par()

        # Wavelengths
        # 1D wavelength solution
        par['calibrations']['wavelengths'][
            'rms_threshold'] = 0.30  # 0.20  # Might be grating dependent..
        par['calibrations']['wavelengths']['sigdetect'] = 5.0
        par['calibrations']['wavelengths']['fwhm'] = 4.0
        par['calibrations']['wavelengths']['n_final'] = 4
        par['calibrations']['wavelengths']['lamps'] = ['OH_NIRES']
        #par['calibrations']['wavelengths']['nonlinear_counts'] = \
        #self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['method'] = 'holy-grail'
        # Reidentification parameters
        par['calibrations']['slitedges']['minimum_slit_length'] = 10.
        par['calibrations']['slitedges']['edge_thresh'] = 30.
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'
        
        # Large chunk of slit is lost with default tweak threshold.
        par['calibrations']['flatfield']['tweak_slits_thresh'] = 0.85

        # Extraction
        # Model full slit currently turned off
        par['reduce']['extraction']['model_full_slit'] = False
        # Tailored profile nsigma parameter for the standard, trying 100 (30
        # was standard
        par['reduce']['extraction']['std_prof_nsigma'] = 100.
        # Perform global sky subtraction for standard stars
        par['reduce']['skysub']['global_sky_std'] = True
        par['reduce']['skysub']['bspline_spacing'] = 0.8
        par['reduce']['extraction']['sn_gauss'] = 4.0

        # Processing steps
        turn_off = dict(use_illumflat=False, use_biasimage=False, use_overscan=False,
                        use_darkimage=False)
        par.reset_all_processimages_par(**turn_off)

        # Flexure
        # Parameters should work for long-slit N1.8 camera exposures
        # N3.75 camera and/or multi-slit may require careful adjustment
        par['scienceframe']['process']['spat_flexure_correct'] = True
        par['scienceframe']['process']['spat_flexure_maxshift'] = 100
        par['scienceframe']['process']['spat_flexure_cont_samp'] = 2
        par['scienceframe']['process']['spat_flexure_sigdetect'] = 2.
        par['calibrations']['tiltframe']['process']['spat_flexure_correct'] = True
        par['calibrations']['tiltframe']['process']['spat_flexure_maxshift'] = 100
        par['calibrations']['tiltframe']['process']['spat_flexure_cont_samp'] = 2
        par['calibrations']['tiltframe']['process']['spat_flexure_sigdetect'] = 2.


        par['scienceframe']['process']['sigclip'] = 20.0
        par['scienceframe']['process']['satpix'] = 'nothing'
        # par['scienceframe']['process']['satpix'] = 'reject'

        return par
        
    def get_rawimage(self, raw_file, det):
        """
        Read raw images and generate a few other bits and pieces that are key
        for image processing.

        .. warning::

            - When reading multiple detectors for a mosaic, this function
              expects all detector arrays to have exactly the same shape.

        Parameters
        ----------
        raw_file : :obj:`str`
            File to read
        det : :obj:`int`, :obj:`tuple`
            1-indexed detector(s) to read.  An image mosaic is selected using a
            :obj:`tuple` with the detectors in the mosaic, which must be one of
            the allowed mosaics returned by :func:`allowed_mosaics`.

        Returns
        -------
        detector : :class:`~pypeit.images.detector_container.DetectorContainer`, :class:`~pypeit.images.mosaic.Mosaic`
            Detector metadata parameters for one or more detectors.
        raw : `numpy.ndarray`_
            Raw image for this detector.  Shape is 2D if a single detector image
            is read and 3D if multiple detectors are read.  E.g., the image from
            the first detector in the tuple is accessed using ``raw_img[0]``.
        hdu : `astropy.io.fits.HDUList`_
            Opened fits file
        texp : :obj:`float`
            Exposure time *in seconds*.
        rawdatasec : `numpy.ndarray`_
            Data (Science) section of the detector as provided by setting the
            (1-indexed) number of the amplifier used to read each detector
            pixel. Pixels unassociated with any amplifier are set to 0.  Shape
            is identical to ``raw_img``.
        oscansec : `numpy.ndarray`_
            Overscan section of the detector as provided by setting the
            (1-indexed) number of the amplifier used to read each detector
            pixel. Pixels unassociated with any amplifier are set to 0.  Shape
            is identical to ``raw_img``.
        """
        
        detector, raw, hdu, texp, datasec, oscansec = super().get_rawimage(raw_file, det)
        
        # Non-linearity correction
        # See: https://scienceops.lbto.org/luci/instrument-characteristics/detector/
        # I assume that the correction applies to each DIT, not the full exposure.
        ndit = hdu[0].header['NDIT']
        raw = ndit*(raw/ndit + 2.767e-6*((raw/ndit)**2.0))
        
        return detector, raw, hdu, texp, datasec, oscansec

# TODO: OUT OF DATE
#    def check_headers(self, headers):
#        """
#        Check headers match expectations for an LBT LUCI1 exposure.
#
#        See also
#        :func:`pypeit.spectrographs.spectrograph.Spectrograph.check_headers`.
#
#        Args:
#            headers (list):
#                A list of headers read from a fits file
#        """
#        expected_values = { '0.INSTRUME': 'LUCI1',
#                            '0.NAXIS': 2 }
#        super(LBTLUCI1Spectrograph, self).check_headers(headers,
#                                                              expected_values=expected_values)

class LBTLUCI2Spectrograph(LBTLUCISpectrograph):
    """
    Child to handle LBT/LUCI2 specific code
    """
    name = 'lbt_luci2'
    camera = 'LUCI2'
    header_name = 'LUCI2'
    supported = True

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
            binning         = '1,1',
            det             = 1,
            dataext         = 0,
            specaxis        = 1,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.25,
            darkcurr        = 0.006,
            # Saturation is 55000, but will be set to dummy value for
            # now
            saturation=1e+8,
            nonlinear       = 0.80,
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d(2.15),
            # Read noise here assumes MER readout mode
            ronoise         = np.atleast_1d(4.5),
            datasec= np.atleast_1d('[5:2044,5:2044]'),
            oscansec= np.atleast_1d('[5:2044,1:4]'),
            )
            
            
        detector = detector_container.DetectorContainer(**detector_dict)
        
        if hdu is None:
            return detector
        
        camera = hdu[0].header['CAMERA']
        if camera == 'N3.75 Camera':
            detector.platescale = 0.119
        
        readmode = hdu[0].header['READMODE']
        if readmode == 'LIR':
            detector.ronoise = np.atleast_1d(9.2)
            
        return detector

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.
        
        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of ``PypeIt`` methods.
        """
        par = super().default_pypeit_par()

        # Wavelengths
        # 1D wavelength solution
        par['calibrations']['wavelengths'][
            'rms_threshold'] = 0.30  # 0.20  # Might be grating dependent..
        par['calibrations']['wavelengths']['sigdetect'] = 5.0
        par['calibrations']['wavelengths']['fwhm'] = 4.0
        par['calibrations']['wavelengths']['n_final'] = 4
        par['calibrations']['wavelengths']['lamps'] = ['OH_NIRES']
        #par['calibrations']['wavelengths']['nonlinear_counts'] = \
        #    self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['method'] = 'holy-grail'
        
        # Large chunk of slit is lost with default tweak threshold.
        par['calibrations']['flatfield']['tweak_slits_thresh'] = 0.85

        par['calibrations']['slitedges']['minimum_slit_length'] = 10.
        par['calibrations']['slitedges']['edge_thresh'] = 30.
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'
        #par['calibrations']['slitedges']['fit_order'] = 8

        # Extraction
        # Model full slit currently turned on
        par['reduce']['extraction']['model_full_slit'] = False
        # Tailored profile nsigma parameter for the standard
        par['reduce']['extraction']['std_prof_nsigma'] = 100.
        # Perform global sky subtraction for standard stars
        par['reduce']['skysub']['global_sky_std'] = True
        par['reduce']['skysub']['bspline_spacing'] = 0.8
        par['reduce']['extraction']['sn_gauss'] = 4.0

        # Processing steps
        turn_off = dict(use_illumflat=False, use_biasimage=False, use_overscan=False,
                        use_darkimage=False)
        par.reset_all_processimages_par(**turn_off)

        # Flexure
        par['flexure']['spec_method'] = 'skip'
        par['scienceframe']['process']['spat_flexure_correct'] = True
        par['scienceframe']['process']['spat_flexure_maxshift'] = 100
        par['scienceframe']['process']['spat_flexure_cont_samp'] = 2
        par['scienceframe']['process']['spat_flexure_sigdetect'] = 2.
        par['calibrations']['tiltframe']['process']['spat_flexure_correct'] = True
        par['calibrations']['tiltframe']['process']['spat_flexure_maxshift'] = 100
        par['calibrations']['tiltframe']['process']['spat_flexure_cont_samp'] = 2
        par['calibrations']['tiltframe']['process']['spat_flexure_sigdetect'] = 2.

        par['scienceframe']['process']['sigclip'] = 20.0
        par['scienceframe']['process']['satpix'] = 'nothing'
        # par['scienceframe']['process']['satpix'] = 'reject'

        return par
        
    def get_rawimage(self, raw_file, det):
        """
        Read raw images and generate a few other bits and pieces that are key
        for image processing.

        .. warning::

            - When reading multiple detectors for a mosaic, this function
              expects all detector arrays to have exactly the same shape.

        Parameters
        ----------
        raw_file : :obj:`str`
            File to read
        det : :obj:`int`, :obj:`tuple`
            1-indexed detector(s) to read.  An image mosaic is selected using a
            :obj:`tuple` with the detectors in the mosaic, which must be one of
            the allowed mosaics returned by :func:`allowed_mosaics`.

        Returns
        -------
        detector : :class:`~pypeit.images.detector_container.DetectorContainer`, :class:`~pypeit.images.mosaic.Mosaic`
            Detector metadata parameters for one or more detectors.
        raw : `numpy.ndarray`_
            Raw image for this detector.  Shape is 2D if a single detector image
            is read and 3D if multiple detectors are read.  E.g., the image from
            the first detector in the tuple is accessed using ``raw_img[0]``.
        hdu : `astropy.io.fits.HDUList`_
            Opened fits file
        texp : :obj:`float`
            Exposure time *in seconds*.
        rawdatasec : `numpy.ndarray`_
            Data (Science) section of the detector as provided by setting the
            (1-indexed) number of the amplifier used to read each detector
            pixel. Pixels unassociated with any amplifier are set to 0.  Shape
            is identical to ``raw_img``.
        oscansec : `numpy.ndarray`_
            Overscan section of the detector as provided by setting the
            (1-indexed) number of the amplifier used to read each detector
            pixel. Pixels unassociated with any amplifier are set to 0.  Shape
            is identical to ``raw_img``.
        """
        
        detector, raw, hdu, texp, datasec, oscansec = super().get_rawimage(raw_file, det)
        
        # Non-linearity correction
        # See: https://scienceops.lbto.org/luci/instrument-characteristics/detector/
        ndit = hdu[0].header['NDIT']
        raw = ndit*(raw/ndit+2.898e-6*((raw/ndit)**2.0))
        
        return detector, raw, hdu, texp, datasec, oscansec

# TODO: OUT OF DATE
#    def check_headers(self, headers):
#        """
#        Check headers match expectations for an LBT LUCI1 exposure.
#
#        See also
#        :func:`pypeit.spectrographs.spectrograph.Spectrograph.check_headers`.
#
#        Args:
#            headers (list):
#                A list of headers read from a fits file
#        """
#        expected_values = { '0.INSTRUME': 'LUCI2',
#                            '0.NAXIS': 2 }
#        super(LBTLUCI1Spectrograph, self).check_headers(headers,
#                                                              expected_values=expected_values)


