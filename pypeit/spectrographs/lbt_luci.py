"""
Module for LBT/LUCI specific methods.

.. include:: ../include/links.rst
"""
from IPython import embed

from pathlib import Path

import numpy as np

from astropy.io import fits
from astropy.table import Table

from pypeit import log
from pypeit import PypeItError
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container
from pypeit.par import parset


class LBTLUCISpectrograph(spectrograph.Spectrograph):
    """
    Class to handle LBT/LUCI specific code

    The provided default reduction parameters have been tailored and tested
    for LUCI1 and LUCI2 using the G200 low resolution grating and the zJspec
    and HKspec filters.

    """
    ndet = 1
    telescope = telescopes.LBTTelescopePar()
    url = 'https://scienceops.lbto.org/luci/'

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)

        # Target meta
        # CORE
        self.meta['ra'] = dict(ext=0, card='OBJRA')
        self.meta['dec'] = dict(ext=0, card='OBJDEC')
        self.meta['target'] = dict(ext=0, card='OBJECT')

        # Instrument meta
        # CORE
        self.meta['dispname'] = dict(ext=0, card='GRATNAME')
        self.meta['decker'] = dict(ext=0, card='MASKID')
        self.meta['binning'] = dict(ext=0, card=None, default='1,1')
        # additional
        self.meta['instrument'] = dict(ext=0, card='INSTRUME')
        self.meta['filter1'] = dict(ext=0, card='FILTER2')
        self.meta['camera'] = dict(ext=0, card='CAMERA')
        # populating dispangle with the grating wavelength
        self.meta['dispangle'] = dict(ext=0, card='GRATWLEN', rtol=0.1)

        # Observation meta
        # CORE
        self.meta['mjd'] = dict(ext=0, card='MJD-OBS')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        # additional
        self.meta['readmode'] = dict(ext=0, card='READMODE')
        self.meta['savemode'] = dict(ext=0, card='SAVEMODE')
        self.meta['dit'] = dict(ext=0, card='DIT')
        self.meta['ndit'] = dict(ext=0, card='NDIT')

        # Define compound meta
        self.meta['idname'] = dict(card=None, compound=True)

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
      
            print ("%s %s %s %s %s %s %s %s %s %s %s" % (targetname,dispname,calib_unit,filter1,filter2,lamp1,lamp2,lamp3,lamp4,lamp5,lamp6))

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
            # pixelflat off -> will be typed as lampoffflats
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

        else:
            raise PypeItError("Not ready for this compound meta")

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
        return ['decker', 'dispname', 'dispangle', 'filter1', 'camera']

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
        return ['MASKID', 'GRATNAME']

    def pypeit_file_keys(self):
        """
        Define the list of keys to be output into a standard PypeIt file.

        Returns:
            :obj:`list`: The list of keywords in the relevant
            :class:`~pypeit.metadata.PypeItMetaData` instance to print to the
            :ref:`pypeit_file`.
        """
        pypeit_keys = super().pypeit_file_keys()
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
        if ftype == 'lampoffflats':
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

        log.debug('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

# Detector information from official LBT LUCI website
# https://sites.google.com/a/lbto.org/luci/instrument-characteristics/detector
 # Parameter 	 LUCI1 	 LUCI2
 # Gain (e-/ADU) 	 2.21 	 2.15
 # Dark Current (e-/sec/pix) 	  0.03 	 0.006
 # Read Noise (e-) 	 LIR 	  9.6 	9.2
 # MER 	  5.1 	 4.5
 # Full Well Capacity (e-) 	  TBD 	 122,000
 # Linearity (% at 40k ADU) 	  11.6% 	 11.1%
 # Crosstalk at Saturation 	  TBD 	 ≤0.2%
 # Persistence (% after 5min) 	  TBD 	 0.03%
 # Minimum DIT (sec) 	 LIR 	  2.503 	 2.503
 # MER 	  6.31 	 6.31


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

        # Get detector parameters from hdu header if available
        if hdu is None:
            # needed for automated construction of the detector table provided
            # in the docs
            platescale = None
            ronoise = None
        else:
            readmode = self.get_meta_value(self.get_headarr(hdu), 'readmode')

            if readmode == 'LIR':
                ronoise = np.atleast_1d(9.6)
            elif readmode == 'MER':
                ronoise = np.atleast_1d(5.1)
            else:
                raise PypeItError("Read mode not recognized (options: LIR, MER)")

            camera = self.get_meta_value(self.get_headarr(hdu), 'camera')
            if camera == 'N1.8 Camera':
                platescale = 0.2500
            elif camera == 'N3.75 Camera':
                platescale = 0.1190
            elif camera == 'N30 Camera': # currently untested but should work in principle
                platescale = 0.0150
            else:
                raise PypeItError("Camera not recognized (options: N1.8, N3.75, N30)")

        # Detector 1
        detector_dict = dict(
            binning='1,1',
            det=1,
            dataext=0,
            specaxis=1,
            specflip=False,
            spatflip=False,
            platescale=platescale,
            darkcurr=108.0, #e-/pix/hour
            # Saturation is 55000, but will be set to dummy value for
            # now as integrated exposures over multiple detector integrations
            # will provide higher counts.
            saturation=1e+8,
            # NIR detectors are non-linear even in lower percentages
            # of the full well, thus for precision measurements one
            # should take into account a general non-linearity
            # correction.
            nonlinear=0.80,
            mincounts=-1e10,
            # In fact there are 32 amplifiers, which gain and ronoise
            # are extremely similar to each other, thus it will be
            # mimicked as 1
            numamplifiers=1,
            gain=np.atleast_1d(2.21),
            ronoise=ronoise,  # variable populated from readmode meta above
            datasec=np.atleast_1d('[5:2044,5:2044]'),
            # For Luci the first 4 pixels on each side can
            # technically be used for as a biassec. This is not
            # included here.
            oscansec=np.atleast_1d('[5:2044,1:4]'),
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
        par['calibrations']['wavelengths']['rms_thresh_frac_fwhm'] = 0.07
        par['calibrations']['wavelengths']['sigdetect'] = 5.0
        par['calibrations']['wavelengths']['fwhm'] = 4.0
        par['calibrations']['wavelengths']['n_final'] = 4
        par['calibrations']['wavelengths']['lamps'] = ['OH_NIRES']
        par['calibrations']['wavelengths']['method'] = 'holy-grail'

        # Reidentification parameters
        par['calibrations']['slitedges']['minimum_slit_length'] = 10.
        par['calibrations']['slitedges']['edge_thresh'] = 30.
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'

        # Large chunk of long slit is lost with default tweak threshold.
        par['calibrations']['flatfield']['tweak_slits_thresh'] = 0.85

        # Extraction
        # Model full slit currently turned off
        par['reduce']['extraction']['model_full_slit'] = False
        # Tailored profile nsigma parameter for the standard, trying 100 (30
        # was standard
        par['reduce']['extraction']['std_prof_nsigma'] = 100.
        # Perform global sky subtraction for standard stars
        par['reduce']['skysub']['global_sky_std'] = False
        par['reduce']['findobj']['maxnumber_std'] = 1
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
        par['scienceframe']['process']['spat_flexure_sigdetect'] = 100
        par['calibrations']['tiltframe']['process']['spat_flexure_correct'] = True
        par['calibrations']['tiltframe']['process']['spat_flexure_sigdetect'] = 100

        par['scienceframe']['process']['sigclip'] = 20.0
        par['scienceframe']['process']['satpix'] = 'nothing'
        
        # Sensitivity function parameters
        par['sensfunc']['algorithm'] = 'IR'
        par['sensfunc']['IR']['telgridfile'] = 'TellPCA_3000_26000_R10000.fits'


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


    def config_specific_par(
            self,
            inp:str|list|Path|fits.Header|Table,
            inp_par:parset.ParSet|None=None
        ) -> parset.ParSet:
        """
        Modify the PypeIt parameters to hard-wired values used for
        specific instrument configurations.

        Args:
            inp (:obj:`str`, :obj:`list`, `Path`_, `astropy.io.fits.Header`_, `astropy.table.Table`_):
                Input filename, an `astropy.io.fits.Header`_ object, or a list
                of `astropy.io.fits.Header`_ objects.  Or a row from the
                metadata table.
            inp_par (:class:`~pypeit.par.parset.ParSet`, optional):
                Parameter set used for the full run of PypeIt.  If None,
                use :func:`default_pypeit_par`.

        Returns:
            :class:`~pypeit.par.parset.ParSet`: The PypeIt parameter set
            adjusted for configuration specific parameter values.
        """
        # Start with instrument-wide parameters
        par = super().config_specific_par(inp, inp_par=inp_par)

        # Adjust parameters based on filters, decker, and disperser angle used
        dispname = self.get_meta_value(inp, 'dispname')
        filter = self.get_meta_value(inp, 'filter1')
        cenwave = self.get_meta_value(inp, 'dispangle')

        if dispname == 'G200 LoRes' and filter == 'zJspec' and cenwave > \
                1.165 and cenwave < 1.175:
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = \
                'lbt_luci1_g200_zJ.fits'

        elif dispname == 'G200 LoRes' and filter == 'HKspec' and cenwave > \
                1.925 and cenwave < 1.935:
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = \
                'lbt_luci1_g200_HK.fits'

        return par

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
        # Get detector parameters from hdu header if available
        if hdu is None:
            # needed for automated construction of the detector table provided
            # in the docs
            platescale = None
            ronoise = None
        else:
            readmode = self.get_meta_value(self.get_headarr(hdu), 'readmode')

            if readmode == 'LIR':
                ronoise = np.atleast_1d(9.2)
            elif readmode == 'MER':
                ronoise = np.atleast_1d(4.5)
            else:
                raise PypeItError("Read mode not recognized (options: LIR, MER)")
                raise ValueError()

            camera = self.get_meta_value(self.get_headarr(hdu), 'camera')
            if camera == 'N1.8 Camera':
                platescale = 0.2500
            elif camera == 'N3.75 Camera':
                platescale = 0.1178
            elif camera == 'N30 Camera': # currently untested but should work in principle
                platescale = 0.0150
            else:
                raise PypeItError("Camera not recognized (options: N1.8, N3.75, N30)")

        # Detector 1
        detector_dict = dict(
            binning='1,1',
            det=1,
            dataext=0,
            specaxis=1,
            specflip=False,
            spatflip=False,
            platescale=platescale,
            darkcurr=21.6, #e-/pix/hour
            # Saturation is 55000, but will be set to dummy value for
            # now as integrated exposures over multiple detector integrations
            # will provide higher counts.
            saturation=1e+8,
            nonlinear=0.80,
            mincounts=-1e10,
            numamplifiers=1,
            gain=np.atleast_1d(2.15),
            ronoise=ronoise,  # variable populated from readmode meta above
            datasec= np.atleast_1d('[5:2044,5:2044]'),
            oscansec= np.atleast_1d('[5:2044,1:4]'),
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
        par['calibrations']['wavelengths']['rms_thresh_frac_fwhm'] = 0.07
        par['calibrations']['wavelengths']['sigdetect'] = 5.0
        par['calibrations']['wavelengths']['fwhm'] = 4.0
        par['calibrations']['wavelengths']['n_final'] = 4
        par['calibrations']['wavelengths']['lamps'] = ['OH_NIRES']
        par['calibrations']['wavelengths']['method'] = 'holy-grail'

        # Large chunk of slit is lost with default tweak threshold.
        par['calibrations']['flatfield']['tweak_slits_thresh'] = 0.85

        par['calibrations']['slitedges']['minimum_slit_length'] = 10.
        par['calibrations']['slitedges']['edge_thresh'] = 30.
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'

        # Extraction
        # Model full slit currently turned on
        par['reduce']['extraction']['model_full_slit'] = False
        # Tailored profile nsigma parameter for the standard
        par['reduce']['extraction']['std_prof_nsigma'] = 100.
        # Perform global sky subtraction for standard stars
        par['reduce']['skysub']['global_sky_std'] = True
        par['reduce']['findobj']['maxnumber_std'] = 1
        par['reduce']['skysub']['bspline_spacing'] = 0.8
        par['reduce']['extraction']['sn_gauss'] = 4.0


        # Processing steps
        turn_off = dict(use_illumflat=False, use_biasimage=False, use_overscan=False,
                        use_darkimage=False)
        par.reset_all_processimages_par(**turn_off)

        # Flexure
        par['flexure']['spec_method'] = 'skip'
        par['scienceframe']['process']['spat_flexure_correct'] = True
        par['scienceframe']['process']['spat_flexure_sigdetect'] = 100
        par['calibrations']['tiltframe']['process']['spat_flexure_correct'] = True
        par['calibrations']['tiltframe']['process']['spat_flexure_sigdetect'] = 100

        par['scienceframe']['process']['sigclip'] = 20.0
        par['scienceframe']['process']['satpix'] = 'nothing'
        # par['scienceframe']['process']['satpix'] = 'reject'
        
        # Sensitivity function parameters
        par['sensfunc']['algorithm'] = 'IR'
        par['sensfunc']['IR']['telgridfile'] = 'TellPCA_3000_26000_R10000.fits'

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

    def config_specific_par(
            self,
            inp:str|list|Path|fits.Header|Table,
            inp_par:parset.ParSet|None=None
        ) -> parset.ParSet:
        """
        Modify the PypeIt parameters to hard-wired values used for
        specific instrument configurations.

        Args:
            inp (:obj:`str`, :obj:`list`, `Path`_, `astropy.io.fits.Header`_, `astropy.table.Table`_):
                Input filename, an `astropy.io.fits.Header`_ object, or a list
                of `astropy.io.fits.Header`_ objects.  Or a row from the
                metadata table.
            inp_par (:class:`~pypeit.par.parset.ParSet`, optional):
                Parameter set used for the full run of PypeIt.  If None,
                use :func:`default_pypeit_par`.

        Returns:
            :class:`~pypeit.par.parset.ParSet`: The PypeIt parameter set
            adjusted for configuration specific parameter values.
        """
        # Start with instrument-wide parameters
        par = super().config_specific_par(inp, inp_par=inp_par)

        # Adjust parameters based on filters and decker used
        dispname = self.get_meta_value(inp, 'dispname')
        filter = self.get_meta_value(inp, 'filter1')

        if dispname == 'G200 LoRes' and filter == 'zJspec':
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = \
                'lbt_luci2_g200_zJ.fits'

        elif dispname == 'G200 LoRes' and filter == 'HKspec':
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = \
                'lbt_luci2_g200_HK.fits'

        return par
