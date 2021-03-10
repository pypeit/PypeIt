""" Module for Magellan/IMACS specific codes
"""
import glob
import numpy as np

from astropy.io import fits
from astropy.time import Time

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.spectrographs import spectrograph
from pypeit.core import parse
from pypeit.images import detector_container

from pkg_resources import resource_filename

class MagellanIMACSSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Magellan/IMACS specific code
    """
    ndet = 8
    telescope = telescopes.MMTTelescopePar()

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
        self.meta['dispname'] = dict(ext=0, card='DISPERSR')
        self.meta['decker'] = dict(ext=0, card='SLITMASK')
        self.meta['binning'] = dict(ext=0, card='BINNING', default='1x1')

        self.meta['mjd'] = dict(ext=0, card=None, compound=True)
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Extras for config and frametyping
        self.meta['idname'] = dict(ext=0, card='EXPTYPE')
        self.meta['dispangle'] = dict(ext=0, card='G-ANGLE', rtol=1e-5)


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
        # TODO: This should be how we always deal with timeunit = 'isot'. Are
        # we doing that for all the relevant spectrographs?
        if meta_key == 'mjd':
            try:
                time = headarr[1]['DATE-OBS']
            except:
                time = headarr[0]['DATE-OBS']
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
        return ['dispname', 'decker', 'binning', 'dispangle']


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

        if ftype == 'science':
            return good_exp & (fitstbl['idname'] == 'Object')
        if ftype == 'standard':
            return good_exp & (fitstbl['idname'] == 'Object')
        if ftype in ['arc', 'tilt']:
            #arc1 = good_exp & (fitstbl['target'] == 'Arc')
            #arc2 = (fitstbl['dispname'] == 'Gri-300-26.7') & \
            #       (fitstbl['idname'] == 'Object') &  \
            #       (fitstbl['exptime'] >1200)
            return good_exp & (fitstbl['target'] == 'Arc')
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            return good_exp & (fitstbl['target'] == 'Flat')

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
        raw_file = raw_file.replace('c1.fits','c{:01d}.fits'.format(det))
        fil = glob.glob(raw_file + '*')
        if len(fil) != 1:
            msgs.error("Found {:d} files matching {:s}".format(len(fil)))

        # Read
        msgs.info("Reading IMACS F2 file: {:s}".format(fil[0]))
        hdu = fits.open(fil[0])
        head1 = fits.getheader(fil[0], 0)
        rawdata = fits.getdata(fil[0], 0)*1.0

        # get the x and y binning factors...
        detector_par = self.get_detector_par(hdu, det if det is not None else 1)
        msgs.work('Need to tweak with binned data.')
        xbin, ybin = parse.parse_binning(detector_par['binning'])


        # First read over the header info to determine the size of the output array...
        datasec = head1['DATASEC']
        xdata1, xdata2, ydata1, ydata2 = np.array(parse.load_sections(datasec, fmt_iraf=True)).flatten()

        # Get the overscan section
        biassec = head1['BIASSEC']
        xbias1, xbias2, ybias1, ybias2 = np.array(parse.load_sections(biassec, fmt_iraf=True)).flatten()

        # allocate output arrays and fill in with mask values
        rawdatasec_img = np.zeros_like(rawdata, dtype=int)
        oscansec_img = np.zeros_like(rawdata, dtype=int)
        rawdatasec_img[xdata1:xdata2, ydata1:ydata2] = 1
        oscansec_img[xbias1:xbias2, ybias1:ybias2] = 1

        if det>4:
            rawdata = np.flipud(np.fliplr(rawdata))
            rawdatasec_img = np.flipud(np.fliplr(rawdatasec_img))
            oscansec_img = np.flipud(np.fliplr(oscansec_img))

        # Need the exposure time
        try:
            exptime = hdu[self.meta['exptime']['ext']].header[self.meta['exptime']['card']]
        except:
            exptime = head1[self.meta['exptime']['card']]
        # Return, transposing array back to orient the overscan properly
        return detector_par,rawdata, hdu, exptime, rawdatasec_img, oscansec_img

class MagellanIMACSF2Spectrograph(MagellanIMACSSpectrograph):
    """
    Child to handle IMACS/F2 specific code
    """
    name = 'magellan_imacsf2'
    camera = 'IMACS'
    supported = True
    comment = 'IMACS f2 camera'

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
        # Binning
        # TODO: Could this be detector dependent?
        binning = self.get_meta_value(self.get_headarr(hdu), 'binning')

        # Detector 1
        detector_dict1 = dict(
            binning         = binning,
            det             = 1,
            dataext         = 0,
            specaxis        = 0,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.2,
            darkcurr        = 2.28,
            saturation      = 65535., # ADU
            nonlinear       = 0.95,
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d(1.56),
            ronoise         = np.atleast_1d(5.4),
            )
        # Detector 2
        detector_dict2 = detector_dict1.copy()
        detector_dict2.update(dict(
            det=2,
            darkcurr=1.,
            gain            = np.atleast_1d(1.56),
            ronoise         = np.atleast_1d(5.6),
        ))
        # Detector 3
        detector_dict3 = detector_dict1.copy()
        detector_dict3.update(dict(
            det=3,
            darkcurr=1.,
            gain            = np.atleast_1d(1.68),
            ronoise         = np.atleast_1d(5.4),
        ))
        # Detector 4
        detector_dict4 = detector_dict1.copy()
        detector_dict4.update(dict(
            det=4,
            darkcurr=1.,
            gain            = np.atleast_1d(1.59),
            ronoise         = np.atleast_1d(6.8),
        ))
        # Detector 5
        detector_dict5 = detector_dict1.copy()
        detector_dict5.update(dict(
            det=5,
            darkcurr=1.,
            gain            = np.atleast_1d(1.58),
            ronoise         = np.atleast_1d(5.6),
        ))
        # Detector 6
        detector_dict6 = detector_dict1.copy()
        detector_dict6.update(dict(
            det=6,
            darkcurr=1.,
            gain            = np.atleast_1d(1.61),
            ronoise         = np.atleast_1d(5.9),
        ))
        # Detector 7
        detector_dict7 = detector_dict1.copy()
        detector_dict7.update(dict(
            det=7,
            darkcurr=1.,
            gain            = np.atleast_1d(1.59),
            ronoise         = np.atleast_1d(6.3),
        ))
        # Detector 8
        detector_dict8 = detector_dict1.copy()
        detector_dict8.update(dict(
            det=8,
            darkcurr=1.,
            gain            = np.atleast_1d(1.65),
            ronoise         = np.atleast_1d(6.7),
        ))
        detectors = [detector_dict1, detector_dict2, detector_dict3, detector_dict4,
                     detector_dict5, detector_dict6, detector_dict7, detector_dict8]
        # Return
        return detector_container.DetectorContainer(**detectors[det-1])

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.

        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of ``PypeIt`` methods.
        """
        par = super().default_pypeit_par()

        # Spectral flexure correction
        par['flexure']['spec_method'] = 'boxcar'
        # Set wave tilts order
        par['calibrations']['slitedges']['edge_thresh'] = 3.
        par['calibrations']['slitedges']['fit_order'] = 3
        par['calibrations']['slitedges']['minimum_slit_gap'] = 1.0
        par['calibrations']['slitedges']['minimum_slit_length_sci'] = 50.
        par['calibrations']['slitedges']['det_buffer'] = 40

        # 1D wavelength solution
        par['calibrations']['wavelengths']['rms_threshold'] = 0.5
        par['calibrations']['wavelengths']['sigdetect'] = 10.
        par['calibrations']['wavelengths']['fwhm']= 5.0
        par['calibrations']['wavelengths']['method'] = 'holy-grail'

        # Set the default exposure time ranges for the frame typing
        par['calibrations']['standardframe']['exprng'] = [10, 150]
        par['calibrations']['arcframe']['exprng'] = [0, 10]
        par['calibrations']['tiltframe']['exprng'] = [0, 10]
        par['calibrations']['darkframe']['exprng'] = [0, None]
        par['scienceframe']['exprng'] = [150, None]


        # Do not require bias frames
        turn_off = dict(use_biasimage=False)
        par.reset_all_processimages_par(**turn_off)

        # Sensitivity function parameters
        par['sensfunc']['algorithm'] = 'IR'
        par['sensfunc']['polyorder'] = 7
        par['sensfunc']['IR']['telgridfile'] = resource_filename('pypeit',
                                                '/data/telluric/TelFit_LasCampanas_3100_26100_R20000.fits')


        return par


    def config_specific_par(self, scifile, inp_par=None):
        """
        Modify the ``PypeIt`` parameters to hard-wired values used for
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

        headarr = self.get_headarr(scifile)

        # Turn PCA off for long slits
        if ('Long' in self.get_meta_value(headarr, 'decker')) :
            par['calibrations']['slitedges']['sync_predict'] = 'nearest'
        # Turn on the use of mask design
        else:
            par['calibrations']['slitedges']['use_maskdesign'] = True
            # Since we use the slitmask info to find the alignment boxes, I don't need `minimum_slit_length_sci`
            par['calibrations']['slitedges']['minimum_slit_length_sci'] = None
            # Sometime the added missing slits at the edge of the detector are to small to be useful.
            par['calibrations']['slitedges']['minimum_slit_length'] = 2.
            # Since we use the slitmask info to add and remove traces, 'minimum_slit_gap' may undo the matching effort.
            par['calibrations']['slitedges']['minimum_slit_gap'] = 0.

        # Templates
        #if self.get_meta_value(headarr, 'dispname') == 'Gri-300-26.7':
        #    par['calibrations']['wavelengths']['lamps'] = ['OH_MODS']
        #else:
        par['calibrations']['wavelengths']['lamps'] = ['HeI', 'NeI', 'ArI', 'ArII']

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
        # Get the empty bpm: force is always True
        bpm_img = self.empty_bpm(filename, det, shape=shape)

        # Fill in bad pixels if a master bias frame is provided
        if msbias is not None:
            return self.bpm_frombias(msbias, det, bpm_img)

        return bpm_img
