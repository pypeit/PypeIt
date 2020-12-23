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

from pkg_resources import resource_filename

class MMTBlueChannelSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle MMT/Blue Channel specific code
    """
    ndet = 1
    name = 'mmt_bluechannel'
    telescope = telescopes.MMTTelescopePar()
    camera = 'Blue_Channel'
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
        header = hdu[0].header

        # Binning
        binning = self.get_meta_value(self.get_headarr(hdu), 'binning')

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
            darkcurr        = header['DARKCUR'],
            saturation      = 65535.,
            nonlinear       = 0.95,  # need to look up and update
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d(header['GAIN']),
            ronoise         = np.atleast_1d(header['RDNOISE']),
            # note that the header entries use the binned sizes
            datasec         = np.atleast_1d(header['DATASEC']),
            oscansec        = np.atleast_1d(header['BIASSEC'])
        )

        return detector_container.DetectorContainer(**detector_dict)

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
        self.meta['decker'] = dict(ext=0, card='APERTURE')
        self.meta['dichroic'] = dict(ext=0, card='INSFILTE')
        self.meta['binning'] = dict(ext=0, card=None, compound=True)
        self.meta['mjd'] = dict(ext=0, card=None, compound=True)
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')

        # Extras for config and frametyping
        self.meta['dispname'] = dict(ext=0, card='DISPERSE')
        self.meta['idname'] = dict(ext=0, card='IMAGETYP')

        # used for arc and continuum lamps
        self.meta['lampstat01'] = dict(ext=0, card=None, compound=True)

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
        par['calibrations']['wavelengths']['rms_threshold'] = 0.5
        par['calibrations']['wavelengths']['sigdetect'] = 5.
        par['calibrations']['wavelengths']['fwhm']= 5.0
        # HeNeAr is by far most commonly used, though ThAr is used for some situations.
        par['calibrations']['wavelengths']['lamps'] = ['ArI', 'ArII', 'HeI', 'NeI']
        par['calibrations']['wavelengths']['method'] = 'holy-grail'

        # Processing steps
        turn_off = dict(use_biasimage=False, use_darkimage=False)
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
        par['calibrations']['pixelflatframe']['exprng'] = [None, 100]
        par['calibrations']['traceframe']['exprng'] = [None, 100]
        par['calibrations']['standardframe']['exprng'] = [None, 600]
        par['calibrations']['arcframe']['exprng'] = [10, None]
        par['calibrations']['darkframe']['exprng'] = [300, None]

        # less than 30 sec implies conditions are bright enough for scattered
        # light to be significant which affects the illumination of the slit.
        par['calibrations']['illumflatframe']['exprng'] = [30, None]

        # Need to specify this for long-slit data
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'

        # Sensitivity function parameters
        par['sensfunc']['polyorder'] = 7

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

        if det == 1:
            msgs.info("Using hard-coded BPM for  Blue Channel")

            bpm_img[-1, :] = 1

        else:
            msgs.error(f"Invalid detector number, {det}, for MMT Blue Channel (only one detector).")

        return bpm_img

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
        return ['dispname']

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
        if ftype in ['science', 'standard']:
            return good_exp & (fitstbl['lampstat01'] == 'off') & (fitstbl['idname'] == 'object')
        if ftype in ['arc', 'tilt']:
            # should flesh this out to include all valid arc lamp combos
            return good_exp & (fitstbl['lampstat01'] != 'off') & (fitstbl['idname'] == 'comp') & (fitstbl['decker'] != '5.0x180')
        if ftype in ['trace', 'pixelflat']:
            # i think the bright lamp, BC, is the only one ever used for this. imagetyp should always be set to flat.
            return good_exp & (fitstbl['lampstat01'] == 'BC') & (fitstbl['idname'] == 'flat')
        if ftype in ['illumflat']:
            # these can be set to flat or object depending if they're twilight or dark sky
            return good_exp & (fitstbl['idname'] in ['flat', 'object']) & (fitstbl['lampstat01'] == 'off')

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
        detector_par = self.get_detector_par(hdu, det if det is None else 1)
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


