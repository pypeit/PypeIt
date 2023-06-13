"""
Module for the SOAR/Goodman instrument

.. include:: ../include/links.rst
"""
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
    url = 'https://noirlab.edu/science/programs/ctio/instruments/goodman-high-throughput-spectrograph'
    allowed_extensions = [".fz"]

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
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
        return ['GRATING', 'SLIT', 'CCDSUM', 'GRT_ANG']

#    def pypeit_file_keys(self):
#        """
#        Define the list of keys to be output into a standard PypeIt file.
#
#        Returns:
#            :obj:`list`: The list of keywords in the relevant
#            :class:`~pypeit.metadata.PypeItMetaData` instance to print to the
#            :ref:`pypeit_file`.
#        """
#        return super().pypeit_file_keys() + ['slitwid']

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
            # Check if all are off
            return np.all(np.array([ (np.char.lower(fitstbl[k]) == 'false') | (np.char.lower(fitstbl[k]) == 'none')
                                        for k in fitstbl.keys() if 'lampstat' in k]), axis=0)
        if status == 'arc':
            # Check if any arc lamps are on
            arc_lamp_stat = [ 'lampstat{0:02d}'.format(i) for i in range(1,5) ]
            return np.any(np.array([ np.char.lower(fitstbl[k]) == 'true' for k in fitstbl.keys()
                                            if k in arc_lamp_stat]), axis=0)
        if status == 'dome':
            # Check if any dome lamps are on
            dome_lamp_stat = [ 'lampstat{0:02d}'.format(i) for i in [6, 8] ]
            return np.any(np.array([ np.char.lower(fitstbl[k]) == 'true' for k in fitstbl.keys()
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
        if ftype in ['science']:
            return good_exp & (fitstbl['idname'] == 'SPECTRUM') & self.lamps(fitstbl, 'off')
        if ftype in ['standard']:
            # Don't type pinhole or dark frames
            return np.zeros(len(fitstbl), dtype=bool) & self.lamps(fitstbl, 'off')
        if ftype == 'bias':
            # Don't type bias
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            return good_exp & self.lamps(fitstbl, 'dome')
        if ftype in ['pinhole', 'dark']:
            # Don't type pinhole or dark frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc', 'tilt']:
            return good_exp & self.lamps(fitstbl, 'arc')
        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)


class SOARGoodmanRedSpectrograph(SOARGoodmanSpectrograph):
    name = 'soar_goodman_red'
    camera = 'red'
    comment = 'Supported gratings: 400_SYZY at M1 and M2 tilts'
    supported = True

    def get_detector_par(self, det, hdu=None):
        """
        Return metadata for the selected detector.

        .. warning::

            Many of the necessary detector parameters are read from the file
            header, meaning the ``hdu`` argument is effectively **required** for
            SOAR/Goodman-Red.  The optional use of ``hdu`` is only viable for
            automatically generated documentation.

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
        if hdu is None:
            binning = '2,2'
            gain = None
            ronoise = None
            datasec = None
            oscansec = None
        else:
            # TODO: Could this be detector dependent??
            binning = self.get_meta_value(self.get_headarr(hdu), 'binning')
            gain = np.atleast_1d(hdu[1].header['GAIN'])
            ronoise = np.atleast_1d(hdu[1].header['RDNOISE'])
            datasec = None
            oscansec = None

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
            gain            = gain,
            ronoise         = ronoise,
            datasec         = datasec,
            oscansec        = oscansec
            )

        if hdu is None:
            return detector_container.DetectorContainer(**detector_dict)

        # Only tested for 2x2
        if binning == '2,2':
            # parse TRIMSEC
            col0 = int(hdu[1].header['TRIMSEC'][1:].split(':')[0])
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
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        # Turn off bias and turn on overscan
        turn_off_on = dict(use_biasimage=False, use_darkimage=False, use_overscan=True)
        par.reset_all_processimages_par(**turn_off_on)

        # Ignore PCA
        par['calibrations']['slitedges']['bound_detector'] = True
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'

        # Always correct for flexure, starting with default parameters
        par['flexure']['spec_method'] = 'boxcar'

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

        #par['sensfunc']['algorithm'] = 'IR'
        par['sensfunc']['IR']['telgridfile'] = 'TelFit_LasCampanas_3100_26100_R20000.fits'

        # TODO: Temporary fix for failure mode.  Remove once Ryan provides a
        # fix.
        par['calibrations']['flatfield']['slit_illum_finecorr'] = False

        return par


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
        # Start with instrument wide
        par = super().config_specific_par(scifile, inp_par=inp_par)

        # Wavelength calibrations
        # Here is a useful website with an arc atlas
        # http://soartelescope.org/soar/content/goodman-comparison-lamps
        if self.get_meta_value(scifile, 'dispname') == '400_SYZY':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'soar_goodman_red_400_SYZY.fits'
            par['calibrations']['wavelengths']['method'] = 'full_template'
        elif self.get_meta_value(scifile, 'dispname') == '600_SYZY_OLD':
            par['calibrations']['wavelengths']['lamps'] = ['NeI', 'ArI', 'HgI']
            par['calibrations']['wavelengths']['reid_arxiv'] = 'soar_goodman_red_600_SYZY_OLD.fits'
            par['calibrations']['wavelengths']['method'] = 'full_template'

        # Return
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
                Processed bias frame used to identify bad pixels

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


class SOARGoodmanBlueSpectrograph(SOARGoodmanSpectrograph):
    name = 'soar_goodman_blue'
    camera = 'blue'
    comment = 'Supported gratings: 400_SYZY at M1 tilt'
    supported = True

    def get_detector_par(self, det, hdu=None):
        """
        Return metadata for the selected detector.

        .. warning::

            Many of the necessary detector parameters are read from the file
            header, meaning the ``hdu`` argument is effectively **required** for
            SOAR/Goodman-Blue.  The optional use of ``hdu`` is only viable for
            automatically generated documentation.

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
        if hdu is None:
            binning = '2,2'
            gain = None
            ronoise = None
            datasec = None
            oscansec = None
        else:
            # TODO: Could this be detector dependent??
            binning = self.get_meta_value(self.get_headarr(hdu), 'binning')
            gain = np.atleast_1d(hdu[1].header['GAIN'])
            ronoise = np.atleast_1d(hdu[1].header['RDNOISE'])
            datasec = None
            oscansec = None

        # Detector 1
        detector_dict = dict(
            binning=binning,
            det=1,
            dataext=1,
            specaxis=1,
            specflip=False,
            spatflip=False,
            platescale=0.15,
            darkcurr=0.00008,  # e-/s/pix
            saturation=65535.,
            nonlinear=1.0,
            mincounts=-1e10,
            numamplifiers=1,
            gain=gain,
            ronoise=ronoise,
            datasec=np.asarray(['[:,20:4112]']),
            oscansec=np.asarray(['[:,2:16]'])
        )

        if hdu is None:
            return detector_container.DetectorContainer(**detector_dict)

        # Return
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

        # Turn off bias and turn on overscan
        turn_off_on = dict(use_biasimage=False, use_darkimage=False, use_overscan=True)
        par.reset_all_processimages_par(**turn_off_on)

        # Ignore PCA
        par['calibrations']['slitedges']['bound_detector'] = True
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'

        # Always correct for flexure, starting with default parameters
        par['flexure']['spec_method'] = 'boxcar'

        # Set pixel flat combination method
        # par['calibrations']['pixelflatframe']['process']['combine'] = 'median'
        # Change the wavelength calibration method
        par['calibrations']['wavelengths']['method'] = 'holy-grail'
        # par['calibrations']['wavelengths']['method'] = 'reidentify'
        par['calibrations']['wavelengths']['lamps'] = ['NeI', 'ArI', 'HgI']
        # Wavelengths
        # 1D wavelength solution
        par['calibrations']['wavelengths']['rms_threshold'] = 0.5
        par['calibrations']['wavelengths']['sigdetect'] = 5.
        par['calibrations']['wavelengths']['fwhm'] = 5.0

        # par['calibrations']['wavelengths']['n_first'] = 3
        # par['calibrations']['wavelengths']['n_final'] = 5
        # par['calibrations']['wavelengths']['sigdetect'] = 10.0
        # par['calibrations']['wavelengths']['wv_cen'] = 4859.0
        # par['calibrations']['wavelengths']['disp'] = 0.2

        # Set the default exposure time ranges for the frame typing
        # par['calibrations']['biasframe']['exprng'] = [None, 1]
        # par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        # par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['arcframe']['exprng'] = [None, 30]
        par['calibrations']['standardframe']['exprng'] = [None, 120]
        par['scienceframe']['exprng'] = [90, None]

        # par['sensfunc']['algorithm'] = 'IR'
        par['sensfunc']['IR']['telgridfile'] = 'TelFit_LasCampanas_3100_26100_R20000.fits'

        return par

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
        # Start with instrument wide
        par = super().config_specific_par(scifile, inp_par=inp_par)

        # Wavelength calibrations
        # Here is a useful website with an arc atlas
        # http://soartelescope.org/soar/content/goodman-comparison-lamps
        if self.get_meta_value(scifile, 'dispname') == '400_SYZY':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'soar_goodman_blue_400_SYZY.fits'
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['flatfield']['slit_illum_finecorr'] = False  # Turn this off due to junk in the unilluminated part of the detector

        # Return
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
                Processed calibration frame used to identify bad pixels

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
