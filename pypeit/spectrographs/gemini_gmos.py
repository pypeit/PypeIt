"""
Module for Gemini GMOS specific methods.

.. include:: ../include/links.rst
"""
import glob
from pkg_resources import resource_filename

from IPython import embed

import numpy as np

from pypeit import msgs
from pypeit.spectrographs import spectrograph
from pypeit import telescopes
from pypeit import io
from pypeit.core import framematch
from pypeit.core import parse
from pypeit.images import detector_container
from pypeit.par import pypeitpar


class GeminiGMOSSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Gemini/GMOS specific code. This is a base class that
    should not be instantiated.
    """
    ndet = 3

    def __init__(self):
        super().__init__()
#        self.timeunit = 'isot'  # Synthesizes date+time
        self.nod_shuffle_pix = None # Nod & Shuffle

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
        self.meta['decker'] = dict(ext=0, card='MASKNAME')
        self.meta['binning'] = dict(card=None, compound=True)  # Uses CCDSUM

        self.meta['mjd'] = dict(ext=0, card='OBSEPOCH')
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Extras for config and frametyping
        self.meta['dispname'] = dict(ext=0, card='GRATING')
        self.meta['dispangle'] = dict(ext=0, card='CENTWAVE', rtol=1e-5)
        self.meta['dichroic'] = dict(ext=0, card='FILTER1')

        self.meta['datasec'] = dict(ext=1, card='DATASEC')

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
            binspatial, binspec = parse.parse_binning(headarr[1]['CCDSUM'])
            binning = parse.binning2string(binspec, binspatial)
            return binning

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
            return good_exp & (fitstbl['target'] != 'CuAr') & (fitstbl['target'] != 'GCALflat') \
                    & (fitstbl['target'] != 'Bias')
            #& (fitstbl['idname'] == 'OBJECT')
        if ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['target'] == 'CuAr')#& (fitstbl['idname'] == 'ARC')
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            return good_exp & (fitstbl['target'] == 'GCALflat')#& (fitstbl['idname'] == 'FLAT')
        if ftype == 'bias':
            return good_exp & (fitstbl['target'] == 'Bias')#& (fitstbl['idname'] == 'BIAS')

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

        par['calibrations']['slitedges']['edge_thresh'] = 20.
        par['calibrations']['slitedges']['fit_order'] = 3

        # 1D wavelength solution
        par['calibrations']['wavelengths']['rms_threshold'] = 0.40  # Might be grating dependent..
        par['calibrations']['wavelengths']['sigdetect'] = 5.  # Doesn't work for reddest chip
        par['calibrations']['wavelengths']['lamps'] = ['CuI', 'ArI', 'ArII']
        par['calibrations']['wavelengths']['method'] = 'full_template'
        par['calibrations']['wavelengths']['nsnippet'] = 1  # 3 detectors splitting is already a lot

        par['calibrations']['tilts']['tracethresh'] = 10.  # Deals with faint CuAr lines

        #   IF YOU CHANGE THIS, YOU WILL NEED TO DEAL WITH THE OVERSCAN GOING ALONG ROWS
        #for key in par['calibrations'].keys():
        #    if 'frame' in key:
        #        par['calibrations'][key]['process']['overscan'] = 'median'

        # Overscan subtract the images
        #par['calibrations']['biasframe']['useframe'] = 'overscan'

        # Alter the method used to combine pixel flats
        par['calibrations']['pixelflatframe']['process']['combine'] = 'median'

        # Always correct for flexure
        par['flexure']['spec_method'] = 'boxcar'
        # Splice detectors 1,2,3 when creating sensitivity function
        par['sensfunc']['multi_spec_det'] = [1,2,3]

        # Set the default exposure time ranges for the frame typing
        #par['scienceframe']['exprng'] = [30, None]

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
        if 'arcsec' in self.get_meta_value(headarr, 'decker'):
            par['calibrations']['slitedges']['sync_predict'] = 'nearest'

        # Allow for various binning
        binning = parse.parse_binning(self.get_meta_value(headarr, 'binning'))
        par['calibrations']['wavelengths']['fwhm'] = 8.0 / binning[1]

        return par

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
        return super().configuration_keys() + ['dispangle', 'datasec']

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
        # TODO: I don't remember what we decided about this use of glob...
        # Check for file; allow for extra .gz, etc. suffix
        fil = glob.glob(raw_file + '*')
        if len(fil) != 1:
            msgs.error("Found {:d} files matching {:s}".format(len(fil)))

        # Read
        msgs.info("Reading GMOS file: {:s}".format(fil[0]))
        hdu = io.fits_open(fil[0])
        head0 = hdu[0].header
        head1 = hdu[1].header

        # TODO: I don't understand this comment, why not use self.ndet?

        # Number of amplifiers (could pull from DetectorPar but this avoids
        # needing the spectrograph, e.g. view_fits)
        numamp = (len(hdu) - 1) // 3

        # get the x and y binning factors...
        binning = head1['CCDSUM']
        xbin, ybin = [int(ibin) for ibin in binning.split(' ')]

        # First read over the header info to determine the size of the output array...
        datasec = head1['DATASEC']
        x1, x2, y1, y2 = np.array(parse.load_sections(datasec, fmt_iraf=False)).flatten()
        biassec = head1['BIASSEC']
        b1, b2, b3, b4 = np.array(parse.load_sections(biassec, fmt_iraf=False)).flatten()
        nxb = b2 - b1 + 1

        # determine the output array size...
        nx = (x2 - x1 + 1) * numamp + nxb * numamp
        ny = y2 - y1 + 1

        # allocate output array...
        array = np.zeros((nx, ny))
        rawdatasec_img = np.zeros_like(array, dtype=int)
        oscansec_img = np.zeros_like(array, dtype=int)

        # TODO: Why is this stuff here and not in the relevant subclass?
        if numamp == 2:  # E2V
            if det == 1:  # BLUEST DETECTOR
                order = range(6, 4, -1)
            elif det == 2:  # NEXT
                order = range(3, 5)
            elif det == 3:  # REDDEST DETECTOR
                order = range(1, 3)
        elif numamp == 4:  # Hamamatsu
            if det == 1:  # BLUEST DETECTOR
                order = range(12, 8, -1)
            elif det == 2:  # BLUEST DETECTOR
                order = range(8, 4, -1)
            elif det == 3:  # BLUEST DETECTOR
                order = range(4, 0, -1)
        else:
            embed()

        # insert extensions into master image...
        for kk, jj in enumerate(order):
            # grab complete extension...
            data, overscan, datasec, biassec, x1, x2 = gemini_read_amp(hdu, jj)
            # insert components into output array...
            inx = data.shape[0]
            xs = inx * kk
            xe = xs + inx

            # insert data...
            # Data section
            #section = '[{:d}:{:d},:]'.format(xs * xbin, xe * xbin)  # Eliminate lines
            #dsec.append(section)
            array[xs:xe, :] = np.flipud(data)
            rawdatasec_img[xs:xe, :] = kk+1

            # ; insert postdata...
            xs = nx - numamp * nxb + kk * nxb
            xe = xs + nxb

            #osection = '[{:d}:{:d},:]'.format(xs * xbin, xe * xbin)  # TRANSPOSED FOR WHAT COMES
            #osec.append(osection)
            array[xs:xe, :] = overscan
            oscansec_img[xs:xe, :] = kk+1

        # Need the exposure time
        exptime = hdu[self.meta['exptime']['ext']].header[self.meta['exptime']['card']]

        # Transpose now (helps with debuggin)
        array = array.T
        rawdatasec_img = rawdatasec_img.T
        oscansec_img = oscansec_img.T

        # TODO: Move to the relevant subclass.
        # Hack me
        if self.name == 'gemini_gmos_north_ham_ns' \
                and head0['object'] in ['GCALflat', 'CuAr', 'Bias'] \
                and self.nod_shuffle_pix is not None:
            # TODO -- Should double check NOD&SHUFFLE was not on
            row1, row2 = 1456, 2812 # NEED TO FIGURE OUT HOW TO GENERALIZE THIS
            nodpix = self.nod_shuffle_pix
            # Shuffle me
            array[row1-nodpix:row2-nodpix,:] = array[row1:row2,:]

        # Return, transposing array back to orient the overscan properly
        return self.get_detector_par(hdu, det if det is None else 1), \
                array, hdu, exptime, rawdatasec_img, oscansec_img


class GeminiGMOSSHamSpectrograph(GeminiGMOSSpectrograph):
    """
    Child to handle Gemini/GMOS-S instrument with Hamamatsu detector
    """
    name = 'gemini_gmos_south_ham'
    camera = 'GMOS-S'
    telescope = telescopes.GeminiSTelescopePar()
    supported = True
    comment = 'Hamamatsu detector (R400, B600, R831); see :doc:`gemini_gmos`'

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
        # TODO: Could this be detector dependent??
        binning = self.get_meta_value(self.get_headarr(hdu), 'binning')

        # Detector 1
        detector_dict1 = dict(
            binning         = binning,
            det             = 1,
            dataext         = 1,  # Not sure this is used
            specaxis        = 1,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.080,
            darkcurr        = 0.0,
            saturation      = 129000.,
            nonlinear       = 0.95,
            mincounts       = -1e10,
            numamplifiers   = 4,
            gain            = np.atleast_1d([1.83]*4),
            ronoise         = np.atleast_1d([3.98]*4),
            )
        # Detector 2
        detector_dict2 = dict(
            binning         = binning,
            det             = 2,
            dataext         = 2,  # Not sure this is used
            specaxis        = 1,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.080,
            darkcurr        = 0.0,
            saturation      = 123000.,
            nonlinear       = 0.95,
            mincounts       = -1e10,
            numamplifiers   = 4,
            gain            = np.atleast_1d([1.83]*4),
            ronoise         = np.atleast_1d([3.98]*4),
            )
        # Detector 3
        detector_dict3 = dict(
            binning         = binning,
            det             = 3,
            dataext         = 3,  # Not sure this is used
            specaxis        = 1,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.080,
            darkcurr        = 0.0,
            saturation      = 125000.,
            nonlinear       = 0.95,
            mincounts       = -1e10,
            numamplifiers   = 4,
            gain            = np.atleast_1d([1.83]*4),
            ronoise         = np.atleast_1d([3.98]*4),
            )
        detectors = [detector_dict1, detector_dict2, detector_dict3]
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
        par['sensfunc']['IR']['telgridfile'] \
                = resource_filename('pypeit',
                                    '/data/telluric/TelFit_LasCampanas_3100_26100_R20000.fits')
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

        # Add to it
        if det == 1:
            msgs.info("Using hard-coded BPM for det=1 on GMOSs")

            # TODO: Fix this
            # Get the binning
            hdu = io.fits_open(filename)
            binning = hdu[1].header['CCDSUM']
            hdu.close()

            # Apply the mask
            xbin = int(binning.split(' ')[0])
            badc = 616//xbin
            bpm_img[badc,:] = 1
        elif det == 2:
            msgs.info("Using hard-coded BPM for det=2 on GMOSs")

            # Get the binning
            hdu = io.fits_open(filename)
            binning = hdu[1].header['CCDSUM']
            hdu.close()

            # Apply the mask
            xbin = int(binning.split(' ')[0])
            if xbin != 2:
                msgs.error("Not prepared for GMOS data wihout 2x binning!")
            # Up high
            badr = (902*2)//xbin # Transposed
            bpm_img[badr:badr+(3*2)//xbin,:] = 1
            # Down low
            badr = (161*2)//xbin # Transposed
            bpm_img[badr,:] = 1
        elif det == 3:
            msgs.info("Using hard-coded BPM for det=3 on GMOSs")

            # Get the binning
            hdu = io.fits_open(filename)
            binning = hdu[1].header['CCDSUM']
            hdu.close()

            # Apply the mask
            xbin = int(binning.split(' ')[0])
            if xbin != 2:
                embed()
            badr = (281*2)//xbin # Transposed
            bpm_img[badr:badr+(2*2)//xbin,:] = 1

        return bpm_img

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
        # Start with instrument wide
        par = super().config_specific_par(scifile, inp_par=inp_par)

        if self.get_meta_value(scifile, 'dispname')[0:4] == 'R400':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gemini_gmos_r400_ham.fits'
        elif self.get_meta_value(scifile, 'dispname')[0:4] == 'B600':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gemini_gmos_b600_ham.fits'
        #
        return par



class GeminiGMOSNSpectrograph(GeminiGMOSSpectrograph):
    """
    Child to handle Gemini/GMOS-N instrument
    """
    telescope = telescopes.GeminiNTelescopePar()
    camera = 'GMOS-N'


class GeminiGMOSNHamSpectrograph(GeminiGMOSNSpectrograph):
    """
    Child to handle Gemini/GMOS-N instrument with Hamamatsu detector
    Used since February 2017
    """
    name = 'gemini_gmos_north_ham'
    supported = True
    comment = 'Hamamatsu detector (R400, B600, R831); Used since Feb 2017; see :doc:`gemini_gmos`'

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
        # TODO: Could this be detector dependent?
        # Binning
        binning = self.get_meta_value(self.get_headarr(hdu), 'binning')

        # Detector 1
        detector_dict1 = dict(
            binning         = binning,
            det             = 1,
            dataext         = 1,  # Not sure this is used
            specaxis        = 1,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.0807,
            darkcurr        = 0.0,
            saturation      = 129000.,
            nonlinear       = 0.95,
            mincounts       = -1e10,
            numamplifiers   = 4,
            gain            = np.atleast_1d([1.63]*4),
            ronoise         = np.atleast_1d([4.14]*4),
            )
        # Detector 2
        detector_dict2 = dict(
            binning         = binning,
            det             = 2,
            dataext         = 2,  # Not sure this is used
            specaxis        = 1,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.0807,
            darkcurr        = 0.0,
            saturation      = 123000.,
            nonlinear       = 0.95,
            mincounts       = -1e10,
            numamplifiers   = 4,
            gain            = np.atleast_1d([1.63]*4),
            ronoise         = np.atleast_1d([4.14]*4),
            )
        # Detector 3
        detector_dict3 = dict(
            binning         = binning,
            det             = 3,
            dataext         = 3,  # Not sure this is used
            specaxis        = 1,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.0807,
            darkcurr        = 0.0,
            saturation      = 125000.,
            nonlinear       = 0.95,
            mincounts       = -1e10,
            numamplifiers   = 4,
            gain            = np.atleast_1d([1.63]*4),
            ronoise         = np.atleast_1d([4.14]*4),
            )
        detectors = [detector_dict1, detector_dict2, detector_dict3]
        # Return
        return detector_container.DetectorContainer(**detectors[det-1])

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
        # Start with instrument wide
        par = super().config_specific_par(scifile, inp_par=inp_par)

        if self.get_meta_value(scifile, 'dispname')[0:4] == 'R400':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gemini_gmos_r400_ham.fits'
        elif self.get_meta_value(scifile, 'dispname')[0:4] == 'B600':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gemini_gmos_b600_ham.fits'
        elif self.get_meta_value(scifile, 'dispname')[0:4] == 'R831':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gemini_gmos_r831_ham.fits'
        return par


class GeminiGMOSNHamNSSpectrograph(GeminiGMOSNHamSpectrograph):
    """
    Child to handle Gemini/GMOS-N instrument with Hamamatsu detector
    and Nod+Shuffle in an not-really NS manner (for now)
    """
    name = 'gemini_gmos_north_ham_ns'
    supported = True
    comment = 'Same as gemini_gmos_north_ham when used in nod-and-shuffle mode; ' \
              'see :doc:`gemini_gmos`'

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
        # Slurp the NOD&Shuffle
        headarr = self.get_headarr(scifile)
        self.nod_shuffle_pix = headarr[0]['NODPIX']
        #
        return par


class GeminiGMOSNE2VSpectrograph(GeminiGMOSNSpectrograph):
    """
    Child to handle Gemini/GMOS-N instrument with E2V detector
    Used until February 2017
    """
    name = 'gemini_gmos_north_e2v'
    supported = True
    comment = 'E2V detector; see :doc:`gemini_gmos`'

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
        # TODO: Could this be detector dependent?
        # Binning
        binning = self.get_meta_value(self.get_headarr(hdu), 'binning')

        # Detector 1
        detector_dict1 = dict(
            binning         = binning,
            det             = 1,
            dataext         = 1,  # Not sure this is used
            specaxis        = 1,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.0728,  # arcsec per pixel
            darkcurr        = 0.0,
            saturation      = 110900.,
            nonlinear       = 0.95,
            mincounts       = -1e10,
            numamplifiers   = 2,
            gain            = np.atleast_1d([2.27]*2),
            ronoise         = np.atleast_1d([3.32]*2),
            )
        # Detector 2
        detector_dict2 = dict(
            binning         = binning,
            det             = 2,
            dataext         = 2,  # Not sure this is used
            specaxis        = 1,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.0728,
            darkcurr        = 0.0,
            saturation      = 115500.,
            nonlinear       = 0.95,
            mincounts       = -1e10,
            numamplifiers   = 2,
            gain            = np.atleast_1d([2.27]*2),
            ronoise         = np.atleast_1d([3.32]*2),
            )
        # Detector 3
        detector_dict3 = dict(
            binning         = binning,
            det             = 3,
            dataext         = 3,  # Not sure this is used
            specaxis        = 1,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.0728,
            darkcurr        = 0.0,
            saturation      = 116700.,
            nonlinear       = 0.95,
            mincounts       = -1e10,
            numamplifiers   = 2,
            gain            = np.atleast_1d([2.27]*2),
            ronoise         = np.atleast_1d([3.32]*2),
            )
        detectors = [detector_dict1, detector_dict2, detector_dict3]
        # Return
        return detector_container.DetectorContainer(**detectors[det-1])

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

        if self.get_meta_value(scifile, 'dispname')[0:4] == 'R400':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gemini_gmos_r400_e2v.fits'
        #
        return par

# TODO: Someone please check the docstring
def gemini_read_amp(inp, ext):
    """
    Read one amplifier of an Gemini GMOS multi-extension FITS image

    Parameters
    ----------
    inp: :obj:`tuple`
        A two-tuple with either the filename and extension ``(str,int)`` with
        the data to read or the already opened `astropy.io.fits.HDUList`_
        object and extension ``(hdu,int)``.

    Returns
    -------
    data : `numpy.ndarray`_
        2D array with the science region of the raw image.
    overscan : `numpy.ndarray`_
        2D array with the overscan region of the raw image.
    datasec : :obj:`str`
        String representation of the section in the raw image with the
        science data.
    baissec : :obj:`str`
        String representation of the section in the raw image with the
        overscan.
    x1 : :obj:`int`
        Starting pixel along the first axis with the science data in the raw
        image.
    y1 : :obj:`int`
        Starting pixel along the second axis with the science data in the raw
        image.
    """
    # Parse input
    hdu = io.fits_open(inp) if isinstance(inp, str) else inp

    # get entire extension...
    temp = hdu[ext].data.transpose()
    tsize = temp.shape
    nxt = tsize[0]

    # parse the DETSEC keyword to determine the size of the array.
    header = hdu[ext].header
    detsec = header['DETSEC']
    x1, x2, y1, y2 = np.array(parse.load_sections(detsec, fmt_iraf=False)).flatten()

    # parse the DATASEC keyword to determine the size of the science region (unbinned)
    datasec = header['DATASEC']
    xdata1, xdata2, ydata1, ydata2 \
            = np.array(parse.load_sections(datasec, fmt_iraf=False)).flatten()

    # grab the components...
    data = temp[xdata1-1:xdata2,:]

    # Overscan
    biassec = header['BIASSEC']
    xdata1, xdata2, ydata1, ydata2 \
            = np.array(parse.load_sections(biassec, fmt_iraf=False)).flatten()
    overscan = temp[xdata1-1:xdata2,:]

    # Return
    return data, overscan, datasec, biassec, x1, x2



