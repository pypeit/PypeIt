"""
Module for MMT/Blue Channel specific methods.

.. include:: ../include/links.rst
"""
import glob
import numpy as np
from astropy.io import fits
from astropy.time import Time

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph
from pypeit.core import parse
from pypeit.images import detector_container


class MMTBlueChannelSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle MMT/Blue Channel specific code
    """
    ndet = 1
    name = 'mmt_bluechannel'
    header_name = 'mmtbluechan'
    telescope = telescopes.MMTTelescopePar()
    camera = 'Blue_Channel'
    url = 'http://www.mmto.org/instrument-suite/blue-red-channel-spectrographs/blue-channel-details/'
    supported = True

    def get_detector_par(self, det, hdu=None):
        """
        Return metadata for the selected detector.

        .. warning::

            Many of the necessary detector parameters are read from the file
            header, meaning the ``hdu`` argument is effectively **required** for
            MMT/BlueChannel.  The optional use of ``hdu`` is only viable for
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
            binning = '1,1'
            gain = None
            ronoise = None
            darkcurr = None
            datasec = None
            oscansec = None
        else:
            binning = self.get_meta_value(self.get_headarr(hdu), 'binning')
            gain = np.atleast_1d(hdu[0].header['GAIN'])
            ronoise = np.atleast_1d(hdu[0].header['RDNOISE'])
            darkcurr = hdu[0].header['DARKCUR']
            datasec = np.atleast_1d(hdu[0].header['DATASEC'])
            oscansec = np.atleast_1d(hdu[0].header['BIASSEC'])

        # Detector 1
        detector_dict = dict(
            binning         = binning,
            det             = 1,
            dataext         = 0,
            specaxis        = 0,
            specflip        = False,
            spatflip        = False,
            xgap            = 0.,
            ygap            = 0.,
            ysize           = 1.,
            platescale      = 0.3,
            darkcurr        = darkcurr, #header['DARKCUR'],
            saturation      = 65535.,
            nonlinear       = 0.95,  # need to look up and update
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = gain, #np.atleast_1d(header['GAIN']),
            ronoise         = ronoise, #np.atleast_1d(header['RDNOISE']),
            # note that the header entries use the binned sizes
# TODO: These aren't needed because the read_rawimage sets these directly,
# right?
            datasec         = datasec, #np.atleast_1d(header['DATASEC']),
            oscansec        = oscansec #np.atleast_1d(header['BIASSEC'])
        )

        return detector_container.DetectorContainer(**detector_dict)

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
        self.meta['decker'] = dict(ext=0, card='APERTURE')
        self.meta['dichroic'] = dict(ext=0, card='INSFILTE')
        self.meta['binning'] = dict(ext=0, card=None, compound=True)
        self.meta['mjd'] = dict(ext=0, card=None, compound=True)
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')

        # Extras for config and frametyping
        self.meta['dispname'] = dict(ext=0, card='DISPERSE')
        self.meta['idname'] = dict(ext=0, card='IMAGETYP')
        self.meta['dispangle'] = dict(ext=0, card=None, compound=True, rtol=0.002)
        self.meta['cenwave'] = dict(ext=0, card=None, compound=True, rtol=5.0)
        self.meta['filter1'] = dict(ext=0, card='INSFILTE')

        # used for arc and continuum lamps
        self.meta['lampstat01'] = dict(ext=0, card=None, compound=True)
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
            """
            Binning in blue channel headers is space-separated rather than comma-separated.
            """
            binspec, binspatial = headarr[0]['CCDSUM'].split()
            binning = parse.binning2string(binspec, binspatial)
            return binning
        elif meta_key == 'cenwave':
            if headarr[0]['CENWAVE'] == 'moving':
                cenwave = None
            else:
                cenwave = int(headarr[0]['CENWAVE'])
            return cenwave
        elif meta_key == 'dispangle':
            if headarr[0]['TILTPOS'] == 'moving':
                dispangle = None
            else:
                dispangle = float(headarr[0]['TILTPOS'])
            return dispangle
        elif meta_key == 'mjd':
            """
            Need to combine 'DATE-OBS' and 'UT' headers and then use astropy to make an mjd.
            """
            date = headarr[0]['DATE-OBS']
            ut = headarr[0]['UT']
            ttime = Time(f"{date}T{ut}", format='isot')
            return ttime.mjd
        elif meta_key == 'lampstat01':
            """
            If the comparison mirror is in, there will be a 'COMPLAMP' header entry containing the lamps
            that are turned on. However, if the comparison mirror is out, then this header entry doesn't exist.
            So need to test for it and set to 'Off' if it's not there.
            """
            if 'COMPLAMP' in headarr[0]:
                return headarr[0]['COMPLAMP']
            else:
                return 'off'

        msgs.error(f"Not ready for compound meta, {meta_key}, for MMT Blue Channel.")

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
        return ['dispname', 'dispangle', 'filter1']

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
        return ['DISPERSE', 'TILTPOS', 'INSFILTE']

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
        par['calibrations']['wavelengths']['rms_threshold'] = 0.5
        par['calibrations']['wavelengths']['sigdetect'] = 5.
        par['calibrations']['wavelengths']['fwhm_fromlines'] = True

        # Parse the lamps used from the image header.
        par['calibrations']['wavelengths']['lamps'] = ['use_header']
        par['calibrations']['wavelengths']['method'] = 'holy-grail'

        # Processing steps
        turn_off = dict(use_illumflat=False, use_biasimage=False, use_darkimage=False)
        par.reset_all_processimages_par(**turn_off)

        # Extraction
        par['reduce']['skysub']['bspline_spacing'] = 0.8
        par['reduce']['extraction']['sn_gauss'] = 4.0

        ## Do not perform global sky subtraction for standard stars
        par['reduce']['skysub']['global_sky_std']  = False

        # cosmic ray rejection parameters for science frames
        par['scienceframe']['process']['sigclip'] = 5.0
        par['scienceframe']['process']['objlim'] = 2.0

        # Set the default exposure time ranges for the frame typing

        # Appropriate exposure times for Blue Channel can vary a lot depending
        # on grating and wavelength. E.g. 300 and 500 line gratings need very
        # short exposures for flats to avoid saturation, but the 1200 and 832
        # can use much longer exposures due to the higher resolution and the
        # continuum lamp not being very bright in the blue/near-UV.
        par['calibrations']['pixelflatframe']['exprng'] = [None, 600]
        par['calibrations']['traceframe']['exprng'] = [None, 600]
        par['calibrations']['standardframe']['exprng'] = [None, 600]
        par['calibrations']['arcframe']['exprng'] = [1, None]
        par['calibrations']['darkframe']['exprng'] = [300, None]
        par['calibrations']['illumflatframe']['exprng'] = [None, 3600]

        # less than 30 sec implies conditions are bright enough for scattered
        # light to be significant which affects the illumination of the slit.
        par['calibrations']['illumflatframe']['exprng'] = [1, None]

        # Need to specify this for long-slit data
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'
        par['calibrations']['slitedges']['bound_detector'] = True

        # Sensitivity function parameters
        par['sensfunc']['polyorder'] = 7

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
        par = super().config_specific_par(scifile, inp_par=inp_par)

        grating = self.get_meta_value(scifile, 'dispname')
        cenwave = self.get_meta_value(scifile, 'cenwave')

        if grating in ['300GPM', '500GPM', '800GPM', '1200GPM']:
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = f"mmt_bluechannel_{grating}.fits"

        # the 832 GPM grating can be used in 1st or 2nd order and therefore needs two templates.
        # the blue, 2nd order setting has a usable range from 3200-5500 A while the red, 1st order
        # setting is usable from 6400-10000 A. of course, why one would use a "blue" channel that far
        # into the red is a valid question...
        if grating == '832GPM' and cenwave < 6000:
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = f"mmt_bluechannel_{grating}_order2.fits"

        if grating == '832GPM' and cenwave >= 6000:
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = f"mmt_bluechannel_{grating}_order1.fits"

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

        if det == 1:
            msgs.info("Using hard-coded BPM for  Blue Channel")

            bpm_img[-1, :] = 1

        else:
            msgs.error(f"Invalid detector number, {det}, for MMT Blue Channel (only one detector).")

        return bpm_img

    def get_lamps(self, fitstbl):
        """
        Extract the list of arc lamps used from header

        .. note::

            Blue channel uses a variety of lamps depending on grating and wavelength range. HeNeAr covers
            the vast majority of cases, but ThAr, HgCd, and CuAr have important use cases.

        Args:
            fitstbl (`astropy.table.Table`_):
                The table with the metadata for one or more arc frames.
        Returns:
            lamps (:obj:`list`) : List used arc lamps
        """
        lampspecs = fitstbl['lampstat01']
        lamps = []

        for lampstr in lampspecs:
            if 'Ne' in lampstr:
                lamps += ['NeI']
            if 'HeAr' in lampstr:
                lamps += ['HeI', 'ArI', 'ArII']
            if 'ThAr' in lampstr:
                # this is a hack to work around non-functional ThAr lamp at time test data was taken
                lamps += ['ArI', 'ArII']
            if 'HgCd' in lampstr:
                lamps += ['HgI', 'CdI']

        return list(set(lamps))

    def pypeit_file_keys(self):
        """
        Define the list of keys to be output into a standard PypeIt file.

        Returns:
            :obj:`list`: The list of keywords in the relevant
            :class:`~pypeit.metadata.PypeItMetaData` instance to print to the
            :ref:`pypeit_file`.
        """
        return super().pypeit_file_keys() + ['cenwave','lampstat01']

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
        if ftype == 'bias':
            return fitstbl['idname'] == 'zero'
        if ftype == 'dark':
            return fitstbl['idname'] == 'dark'
        if ftype == 'science':
            return good_exp & (fitstbl['lampstat01'] == 'off') & (fitstbl['idname'] == 'object') & (fitstbl['target'] != 'skyflat')
        if ftype == 'standard':
            return (
                good_exp
                & (fitstbl['lampstat01'] == 'off')
                & (fitstbl['idname'] == 'object')
                & (fitstbl['target'] != 'skyflat')
                & (fitstbl['decker'] == '5.0x180')
            )
        if ftype == 'arc':
            # should flesh this out to include all valid arc lamp combos
            return (
                good_exp
                & (fitstbl['lampstat01'] != 'off')
                & (fitstbl['lampstat01'] != 'BC')
                & (fitstbl['idname'] == 'comp')
                & (fitstbl['decker'] != '5.0x180')
                & (fitstbl['target'] != 'focus')
            )
        if ftype == 'tilt':
            # should flesh this out to include all valid arc lamp combos
            return (
                good_exp
                & (fitstbl['lampstat01'] != 'off')
                & (fitstbl['lampstat01'] != 'BC')
                & (fitstbl['idname'] == 'comp')
                & (fitstbl['decker'] != '5.0x180')
            )
        if ftype in ['trace', 'pixelflat']:
            # i think the bright lamp, BC, is the only one ever used for this. imagetyp should always be set to flat, but sometimes not.
            return good_exp & (fitstbl['lampstat01'] == 'BC')
        if ftype == 'illumflat':
            # i think the bright lamp, BC, is the only one ever used for this. imagetyp should always be set to flat.
            return good_exp & (fitstbl['lampstat01'] == 'off') & (fitstbl['target'] == 'skyflat')

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

    def get_rawimage(self, raw_file, det):
        """
        Read raw images and generate a few other bits and pieces
        that are key for image processing.

        Parameters
        ----------
        raw_file : :obj:`str`
            File to read
        det : :obj:`int`
            1-indexed detector to read

        Returns
        -------
        detector_par : :class:`pypeit.images.detector_container.DetectorContainer`
            Detector metadata parameters.
        raw_img : `numpy.ndarray`_
            Raw image for this detector.
        hdu : `astropy.io.fits.HDUList`_
            Opened fits file
        exptime : :obj:`float`
            Exposure time read from the file header
        rawdatasec_img : `numpy.ndarray`_
            Data (Science) section of the detector as provided by setting the
            (1-indexed) number of the amplifier used to read each detector
            pixel. Pixels unassociated with any amplifier are set to 0.
        oscansec_img : `numpy.ndarray`_
            Overscan section of the detector as provided by setting the
            (1-indexed) number of the amplifier used to read each detector
            pixel. Pixels unassociated with any amplifier are set to 0.
        """
        # Check for file; allow for extra .gz, etc. suffix
        fil = glob.glob(raw_file + '*')
        if len(fil) != 1:
            msgs.error("Found {:d} files matching {:s}".format(len(fil)))

        # Read FITS image
        msgs.info("Reading MMT Blue Channel file: {:s}".format(fil[0]))
        hdu = fits.open(fil[0])
        hdr = hdu[0].header

        # we're flipping FITS x/y to pypeit y/x here. pypeit wants blue on the
        # bottom, slit bottom on the right...
        rawdata = np.fliplr(hdu[0].data.astype(float).transpose())

        exptime = hdr['EXPTIME']

        # TODO Store these parameters in the DetectorPar.
        # Number of amplifiers
        detector_par = self.get_detector_par(det if det is not None else 1, hdu=hdu)
        numamp = detector_par['numamplifiers']

        # First read over the header info to determine the size of the output array...
        datasec = hdr['DATASEC']
        xdata1, xdata2, ydata1, ydata2 = np.array(parse.load_sections(datasec, fmt_iraf=False)).flatten()

        # Get the overscan section
        biassec = hdr['BIASSEC']
        xbias1, xbias2, ybias1, ybias2 = np.array(parse.load_sections(biassec, fmt_iraf=False)).flatten()

        # allocate output arrays and fill in with mask values
        rawdatasec_img = np.zeros_like(rawdata, dtype=int)
        oscansec_img = np.zeros_like(rawdata, dtype=int)

        # trim bad sections at beginning of data and bias sections
        rawdatasec_img[xdata1+2:xdata2, ydata1:ydata2-1] = 1
        oscansec_img[xbias1+2:xbias2, ybias1:ybias2-1] = 1

        return detector_par, rawdata, hdu, exptime, rawdatasec_img, oscansec_img


