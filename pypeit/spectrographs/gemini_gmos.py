"""
Module for Gemini GMOS specific methods.

.. include:: ../include/links.rst
"""
import numpy as np

from pypeit import msgs
from pypeit.spectrographs import spectrograph
from pypeit import telescopes
from pypeit import io
from pypeit.core import framematch
from pypeit.core import parse
from pypeit.images import detector_container
from pypeit.images.mosaic import Mosaic
from pypeit.core.mosaic import build_image_mosaic_transform

from IPython import embed

class GeminiGMOSMosaicLookUp:
    """
    Provides the geometry required to mosaic Gemini GMOS data.

    This is purposely a direct copy of the data provided by v3.0.0 of the Gemini
    `DRAGONS
    <https://github.com/GeminiDRSoftware/DRAGONS/blob/a59e6ff5c8ca79bc64ecd690ac50e4a91278530b/geminidr/gmos/lookups/geometry_conf.py#L26>`__
    software package.  Specifically, if you have DRAGONS installed,
    :attr:`pypeit.spectrographs.gemini_gmos.GeminiGMOSMosaicLookUp.geometry`
    should be identical to:

    .. code-block:: python

        from geminidr.gmos.lookups.geometry_conf import geometry
    
    Updating to any changes made to the DRAGONS version requires by-hand editing
    of the ``PypeIt`` code.
    """
    geometry = {
        # GMOS-N
        'EEV9273-16-03EEV9273-20-04EEV9273-20-03': {'default_shape': (2048, 4608),
                                                    (0, 0): {'shift': (-2087.50,-1.58),
                                                            'rotation': -0.004},
                                                    (2048, 0): {},
                                                    (4096, 0): {'shift': (2088.8723, -1.86),
                                                                'rotation': -0.046}
                                                            },
        'e2v 10031-23-05,10031-01-03,10031-18-04': {'default_shape': (2048, 4608),
                                                    (0, 0): {'shift': (-2087.7,-0.749),
                                                            'rotation': -0.009},
                                                    (2048, 0): {},
                                                    (4096, 0): {'shift': (2087.8014, 2.05),
                                                                'rotation': -0.003}
                                                    },
        'BI13-20-4k-1,BI12-09-4k-2,BI13-18-4k-2': {'default_shape': (2048, 4224),
                                                   (0, 0): {'shift': (-2115.95, -0.21739),
                                                            'rotation': -0.004},
                                                   (2048, 0): {},
                                                   (4096, 0): {'shift': (2115.48, 0.1727),
                                                               'rotation': -0.00537}
                                                   },

    # GMOS-S
        'EEV8056-20-03EEV8194-19-04EEV8261-07-04': {'default_shape': (2048, 4608),
                                                    (0, 0): {'shift': (-2086.44, 5.46),
                                                            'rotation': -0.01},
                                                    (2048, 0): {},
                                                    (4096, 0): {'shift': (2092.53, 9.57),
                                                                'rotation': 0.02}
                                                    },
        'EEV2037-06-03EEV8194-19-04EEV8261-07-04': {'default_shape': (2048, 4608),
                                                    (0, 0): {'shift': (-2086.49, -0.22),
                                                            'rotation': 0.011},
                                                    (2048, 0): {},
                                                    (4096, 0): {'shift': (2089.31, 2.04),
                                                                'rotation': 0.012}
                                                    },
        'BI5-36-4k-2,BI11-33-4k-1,BI12-34-4k-1': {'default_shape': (2048, 4224),
                                                  (0, 0): {'shift': (-2110.2, 0.71),
                                                           'rotation': 0.},
                                                  (2048, 0): {},
                                                  (4096, 0): {'shift': (2109., -0.73),
                                                              'rotation': 0.}
                                                  },
    }


class GeminiGMOSSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Gemini/GMOS specific code. This is a base class that
    should not be instantiated.
    """
    ndet = 3
    detid = None

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
            # binning in the raw frames
            ccdsum = headarr[1].get('CCDSUM')
            if ccdsum is not None:
                binspatial, binspec = parse.parse_binning(ccdsum)
                binning = parse.binning2string(binspec, binspatial)
            else:
                # binning in the spec2d file
                binning = headarr[0].get('BINNING')
            if binning is None:
                msgs.error('Binning not found')
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

        # Always default to reducing as a mosaic
        par['rdx']['detnum'] = [(1,2,3)]

        par['calibrations']['slitedges']['follow_span'] = 80
        par['calibrations']['slitedges']['edge_thresh'] = 100.
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

        # TODO: Note the default is now to mosaic the detectors.  This means the
        # user will need to set this if they ever reduce single detectors at a
        # time.
#        # Splice detectors 1,2,3 when creating sensitivity function
#        par['sensfunc']['multi_spec_det'] = [1,2,3]
        # NOTE: New syntax uses the detector names
#        par['sensfunc']['multi_spec_det'] = ['DET01','DET02','DET03']

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
        par['calibrations']['wavelengths']['fwhm_fromlines'] = True

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

    def hdu_read_order(self):
        """
        Return the order in which to read HDU extensions for this instrument.
        """
        pass

    def get_rawimage(self, raw_file, det):
        """
        Read raw images and generate a few other bits and pieces
        that are key for image processing.

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
        detector_par : :class:`~pypeit.images.detector_container.DetectorContainer`, :class:`~pypeit.images.mosaic.Mosaic`
            Detector metadata parameters for one or more detectors.
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
        # Read
        msgs.info(f'Attempting to read GMOS file: {raw_file}')
        # NOTE: io.fits_open checks that the file exists
        hdu = io.fits_open(raw_file)
        head0 = hdu[0].header
        head1 = hdu[1].header

        # Validate the entered (list of) detector(s)
        nimg, _det = self.validate_det(det)

        # Grab the detector or mosaic parameters
        mosaic = None if nimg == 1 else self.get_mosaic_par(det, hdu=hdu)
        detectors = [self.get_detector_par(det, hdu=hdu)] if nimg == 1 else mosaic.detectors

        # Number of amplifiers is hard-coded as follows
        numamp = (len(hdu) - 1) // self.ndet
        if numamp != detectors[0].numamplifiers:
            msgs.error(f'Unexpected number of amplifiers for {self.name} based on number of '
                       f'extensions in {raw_file}.')

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
        array = np.zeros((nimg, nx, ny))
        rawdatasec_img = np.zeros_like(array, dtype=int)
        oscansec_img = np.zeros_like(array, dtype=int)

        # Get the HDU read order for this instrument
        order = self.hdu_read_order()
        for ii in range(nimg):

            # insert extensions into master image...
            for kk, jj in enumerate(order[_det[ii]-1]):
                # grab complete extension...
                data, overscan, datasec, biassec, x1, x2 = gemini_read_amp(hdu, jj)
                # insert components into output array...
                inx = data.shape[0]
                xs = inx * kk
                xe = xs + inx

                # insert data...
                array[ii,xs:xe,:] = np.flipud(data)
                rawdatasec_img[ii,xs:xe, :] = kk+1

                # ; insert postdata...
                xs = nx - numamp * nxb + kk * nxb
                xe = xs + nxb

                array[ii,xs:xe,:] = overscan
                oscansec_img[ii,xs:xe,:] = kk+1

        # Need the exposure time
        exptime = self.get_meta_value(self.get_headarr(hdu), 'exptime')

        # Transpose now (helps with debuggin)
        array = np.transpose(array, axes=(0,2,1))
        rawdatasec_img = np.transpose(rawdatasec_img, axes=(0,2,1))
        oscansec_img = np.transpose(oscansec_img, axes=(0,2,1))

        # Handle returning both single and multiple images
        if nimg == 1:
            return detectors[0], array[0], hdu, exptime, rawdatasec_img[0], oscansec_img[0]
        return mosaic, array, hdu, exptime, rawdatasec_img, oscansec_img

    def get_mosaic_par(self, mosaic, hdu=None, msc_order=0):
        """
        Return the hard-coded parameters needed to construct detector mosaics
        from unbinned images.

        The parameters expect the images to be trimmed and oriented to follow
        the ``PypeIt`` shape convention of ``(nspec,nspat)``.  For returned
        lists, the length of the list is the same as the number of detectors in
        the mosaic, and they are ordered by the detector number.

        Args:
            mosaic (:obj:`tuple`):
                Tuple of detector numbers used to construct the mosaic.  Must be
                one among the list of possible mosaics as hard-coded by the
                :func:`allowed_mosaics` function.
            hdu (`astropy.io.fits.HDUList`_, optional):
                The open fits file with the raw image of interest.  If not
                provided, frame-dependent detector parameters are set to a
                default.  BEWARE: If ``hdu`` is not provided, the binning is
                assumed to be `1,1`, which will cause faults if applied to
                binned images!
            msc_order (:obj:`int`, optional):
                Order of the interpolation used to construct the mosaic.

        Returns:
            :class:`~pypeit.images.mosaic.Mosaic`: Object with the mosaic *and*
            detector parameters.
        """
        if self.detid is None:
            return None, None, None

        # Validate the entered (list of) detector(s)
        nimg, _ = self.validate_det(mosaic)

        # Index of mosaic in list of allowed detector combinations
        mosaic_id = self.allowed_mosaics.index(mosaic)+1

        # Get the detectors
        detectors = np.array([self.get_detector_par(det, hdu=hdu) for det in mosaic])
        # Binning *must* be consistent for all detectors
        if any(d.binning != detectors[0].binning for d in detectors[1:]):
            msgs.error('Binning is somehow inconsistent between detectors in the mosaic!')

        # Collect the offsets and rotations for *all unbinned* detectors in the
        # full instrument, ordered by the number of the detector.  Detector
        # numbers must be sequential and 1-indexed.
        # NOTE: These lines use the directly copied metadata from the Gemini
        # DRAGONS software and then adjusts them so that they are in "PypeIt
        # format".  See the mosaic documentattion.
        expected_shape = GeminiGMOSMosaicLookUp.geometry[self.detid]['default_shape']
        shift = np.array([( GeminiGMOSMosaicLookUp.geometry[self.detid][(4096,0)]['shift'][1],
                           -GeminiGMOSMosaicLookUp.geometry[self.detid][(4096,0)]['shift'][0]),
                          (0.,0.),
                          ( GeminiGMOSMosaicLookUp.geometry[self.detid][(0,0)]['shift'][1],
                           -GeminiGMOSMosaicLookUp.geometry[self.detid][(0,0)]['shift'][0])])
        rotation = np.array([GeminiGMOSMosaicLookUp.geometry[self.detid][(4096,0)]['rotation'],
                             0.,
                             GeminiGMOSMosaicLookUp.geometry[self.detid][(0,0)]['rotation']])

        # The binning and process image shape must be the same for all images in
        # the mosaic
        binning = tuple(int(b) for b in detectors[0].binning.split(','))
        shape = tuple(n // b for n, b in zip(expected_shape, binning))

        msc_sft = [None]*nimg
        msc_rot = [None]*nimg
        msc_tfm = [None]*nimg
        for i, d in enumerate([det-1 for det in mosaic]):
            msc_sft[i] = shift[d]
            msc_rot[i] = rotation[d]
            msc_tfm[i] = build_image_mosaic_transform(shape, msc_sft[i], msc_rot[i], binning)

        return Mosaic(mosaic_id, detectors, shape, np.array(msc_sft), np.array(msc_rot),
                      np.array(msc_tfm), msc_order)

    @property
    def allowed_mosaics(self):
        """
        Return the list of allowed detector mosaics.

        Gemini GMOS only allows for mosaicing all three detectors.

        Returns:
            :obj:`list`: List of tuples, where each tuple provides the 1-indexed
            detector numbers that can be combined into a mosaic and processed by
            ``PypeIt``.
        """
        return [(1,2,3)]

    @property
    def default_mosaic(self):
        return self.allowed_mosaics[0]


class GeminiGMOSSHamSpectrograph(GeminiGMOSSpectrograph):
    """
    Child to handle Gemini/GMOS-S instrument with Hamamatsu detector
    """
    name = 'gemini_gmos_south_ham'
    camera = 'GMOS-S'
    header_name = 'GMOS-S'
    telescope = telescopes.GeminiSTelescopePar()
    supported = True
    comment = 'Hamamatsu detector (R400, B600, R831); see :doc:`gemini_gmos`'
    detid = 'BI5-36-4k-2,BI11-33-4k-1,BI12-34-4k-1'

    def hdu_read_order(self):
        """
        Return the order in which to read HDU extensions for this instrument.

        Returns:
            :obj:`tuple`: A tuple of three iterables that provide the order of
            the HDU extensions for reading the data from this instrument.
        """
        # Order expected is "blue" detector, "green" detector, "red" detector
        return (range(12, 8, -1), range(8, 4, -1), range(4, 0, -1))

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
        par['sensfunc']['algorithm'] = 'IR'
        par['sensfunc']['IR']['telgridfile'] = 'TelFit_LasCampanas_3100_26100_R20000.fits'
        # Bound the detector with slit edges if no edges are found. These data are often trimmed
        # so we implement this here as the default.
        par['calibrations']['slitedges']['bound_detector'] = True
        return par

    def bpm(self, filename, det, shape=None, msbias=None):
        """
        Generate a default bad-pixel mask.

        Even though they are both optional, either the precise shape for
        the image (``shape``) or an example file that can be read to get
        the shape (``filename`` using :func:`get_image_shape`) *must* be
        provided.

        Args:
            filename (:obj:`str`):
                An example file to use to get the image shape.  **Cannot** be
                None.
            det (:obj:`int`, :obj:`tuple`):
                1-indexed detector(s) to read.  An image mosaic is selected
                using a :obj:`tuple` with the detectors in the mosaic, which
                must be one of the allowed mosaics returned by
                :func:`allowed_mosaics`.
            shape (:obj:`tuple`, optional):
                Processed image shape.  If ``filename`` is None, this *must* be
                provided; otherwise, this is ignored.
            msbias (:class:`~pypeit.images.pypeitimage.PypeItImage`, optional):
                Master bias frame.  If provided, it is used by
                :func:`~pypeit.spectrographs.spectrograph.Spectrograph.bpm_frombias`
                to identify bad pixels.

        Returns:
            `numpy.ndarray`_: An integer array with a masked value set to 1 and
            an unmasked value set to 0.
        """
        # Validate the entered (list of) detector(s)
        nimg, _det = self.validate_det(det)
        _det = list(_det)
        # Call the base-class method to generate the empty bpm
        bpm_img = super().bpm(filename, det, shape=shape, msbias=msbias)
        # NOTE: expand_dims does *not* copy the array.  We can edit it directly
        # because we've created it inside this function.
        _bpm_img = np.expand_dims(bpm_img, 0) if nimg == 1 else bpm_img

        # Get the binning
        # TODO: By definition this line means that filename cannot be None.
        # TODO: We're opening the file too many times...
        hdu = io.fits_open(filename)
        # TODO: Why aren't we usig get_meta_value for this?
        binning = hdu[1].header['CCDSUM']
        xbin = int(binning.split(' ')[0])
        hdu.close()

        # Add the detector-specific, hard-coded bad columns
        if 1 in _det:
            msgs.info("Using hard-coded BPM for det=1 on GMOSs")
            i = _det.index(1)
            # Apply the mask
            badc = 616//xbin
            _bpm_img[i,badc,:] = 1
        if 2 in _det:
            msgs.info("Using hard-coded BPM for det=2 on GMOSs")
            i = _det.index(2)
            # Apply the mask
            # Up high
            badr = (902*2)//xbin # Transposed
            _bpm_img[i,badr:badr+(3*2)//xbin,:] = 1
            # Down low
            badr = (161*2)//xbin # Transposed
            _bpm_img[i,badr,:] = 1
        if 3 in _det:
            msgs.info("Using hard-coded BPM for det=3 on GMOSs")
            i = _det.index(3)
            # Apply the mask
            badr = (281*2)//xbin # Transposed
            _bpm_img[i,badr:badr+(2*2)//xbin,:] = 1
        # Done
        return _bpm_img[0] if nimg == 1 else _bpm_img

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
    header_name = 'GMOS-N'


class GeminiGMOSNHamSpectrograph(GeminiGMOSNSpectrograph):
    """
    Child to handle Gemini/GMOS-N instrument with Hamamatsu detector
    Used since February 2017
    """
    name = 'gemini_gmos_north_ham'
    supported = True
    comment = 'Hamamatsu detector (R400, B600, R831); Used since Feb 2017; see :doc:`gemini_gmos`'
    detid = 'BI13-20-4k-1,BI12-09-4k-2,BI13-18-4k-2'

    def hdu_read_order(self):
        """
        Return the order in which to read HDU extensions for this instrument.

        Returns:
            :obj:`tuple`: A tuple of three iterables that provide the order of
            the HDU extensions for reading the data from this instrument.
        """
        # Order expected is "blue" detector, "green" detector, "red" detector
        return (range(12, 8, -1), range(8, 4, -1), range(4, 0, -1))

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
        # TODO: Could this be detector dependent?
        # Binning
        binning = '1,1' if hdu is None else self.get_meta_value(self.get_headarr(hdu), 'binning')

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

    def __init__(self):
        super().__init__()
        self.nod_shuffle_pix = None # Nod & Shuffle

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

    def get_rawimage(self, raw_file, det):
        """
        Read raw images and generate a few other bits and pieces
        that are key for image processing.

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
        detector_par : :class:`~pypeit.images.detector_container.DetectorContainer`, :class:`~pypeit.images.mosaic.Mosaic`
            Detector metadata parameters for one or more detectors.
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
        detpar, array, hdu, exptime, rawdatasec_img, oscansec_img \
                = super().get_rawimage(raw_file, det)
        # TODO: Actually assign as follows here?
        # self.nod_shuffle_pix = hdu[0].header['NODPIX']

        if self.nod_shuffle_pix is None \
                or hdu[0].header['object'] not in ['GCALflat', 'CuAr', 'Bias']:
            # No need to adjust output
            return detpar, array, hdu, exptime, rawdatasec_img, oscansec_img

        nimg = 1 if array.ndim == 2 else array.shape[0]
        if nimg == 1:
            _detpar = [detpar]
            # NOTE: expand_dims does *not* copy the array.  We can edit it
            # directly because we've created it inside this function.
            _array = np.expand_dims(array, 0)
        elif isinstance(detpar, Mosaic):
            _detpar = detpar.detectors
            _array = array
        else:
            _detpar = detpar
            _array = array
        for i in range(nimg):
            xbin, ybin = parse.parse_binning(_detpar[i].binning)
            # TODO: Should double check NOD&SHUFFLE was not on
            nodpix = int(self.nod_shuffle_pix/xbin)
            #48 is a solid value for the unusful rows in GMOS data
            row1, row2 = nodpix + int(48/xbin), 2*nodpix+int(48/xbin)
            # Shuffle me
            _array[i,row1-nodpix:row2-nodpix,:] = _array[i,row1:row2,:]
        if nimg == 1:
            array = _array[0]

        return detpar, array, hdu, exptime, rawdatasec_img, oscansec_img

class GeminiGMOSNE2VSpectrograph(GeminiGMOSNSpectrograph):
    """
    Child to handle Gemini/GMOS-N instrument with E2V detector
    Used until February 2017
    """
    name = 'gemini_gmos_north_e2v'
    supported = True
    comment = 'E2V detector; see :doc:`gemini_gmos`'
    # TODO: Check this is correct
    detid = 'e2v 10031-23-05,10031-01-03,10031-18-04'

    def hdu_read_order(self):
        """
        Return the order in which to read HDU extension for each detector.

        Returns:
            :obj:`tuple`: A tuple of three iterables that provide the order of
            the HDU extensions for reading the data from this instrument.
        """
        # Order expected is "blue" detector, "green" detector, "red" detector
        return (range(6, 4, -1), range(3, 5), range(1, 3))

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
        # TODO: Could this be detector dependent?
        # Binning
        binning = '1,1' if hdu is None else self.get_meta_value(self.get_headarr(hdu), 'binning')

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
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gemini_gmos_r400_e2v_mosaic.fits'
            # The blue wavelengths are *faint*
            #   But redder observations may prefer something closer to the default
            par['calibrations']['wavelengths']['sigdetect'] = 1.  
        # Return
        return par

# TODO: Someone please check the docstring
def gemini_read_amp(inp, ext):
    """
    Read one amplifier of an Gemini GMOS multi-extension FITS image

    Parameters
    ----------
    inp : :obj:`str`, `astropy.io.fits.HDUList`_
        The file name of a fits file with the data or the already opened
        `astropy.io.fits.HDUList`_ object.
    ext : :obj:`str`, :obj:`int`
        The name or index of the extension in the list of HDUs with the relevant
        data.

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



