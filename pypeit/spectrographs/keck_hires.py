"""
Module for Keck/HIRES

.. include:: ../include/links.rst
"""
import os

import numpy as np
from scipy.io import readsav

from pkg_resources import resource_filename

from astropy.table import Table

from pypeit import msgs
from pypeit import telescopes
from pypeit import io
from pypeit.core import parse
from pypeit.core import framematch
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container
from pypeit.par import pypeitpar
from pypeit.images.mosaic import Mosaic
from pypeit.core.mosaic import build_image_mosaic_transform


from IPython import embed



class HIRESMosaicLookUp:
    """
    Provides the geometry required to mosaic Keck DEIMOS data.
    Similar to :class:`~pypeit.spectrographs.gemini_gmos.GeminiGMOSMosaicLookUp`

    """
    geometry = {
        'MSC01': {'default_shape': (6168, 3990),
                  'blue_det': {'shift': (-2048.0 -12.0, 0.0), 'rotation': 0.},
                  'green_det': {'shift': (0., 0.), 'rotation': 0.},
                  'red_det': {'shift': (2048.0 + 12.0, 0.), 'rotation': 0.}},
    }



class KECKHIRESSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle KECK/HIRES specific code.

    This spectrograph is not yet supported.
    """

    ndet = 3
    name = 'keck_hires'
    telescope = telescopes.KeckTelescopePar()
    camera = 'HIRES'
    header_name = 'HIRES'
    pypeline = 'Echelle'
    ech_fixed_format = False
    supported = True


    # Place holder taken from X-shooter VIS
    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.

        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of ``PypeIt`` methods.
        """
        par = super().default_pypeit_par()

        # Adjustments to parameters for VIS
        turn_on = dict(use_biasimage=False, use_overscan=True, overscan_method='median',
                       use_darkimage=False, use_illumflat=False, use_pixelflat=False,
                       use_specillum=False)
        par.reset_all_processimages_par(**turn_on)
        # X-SHOOTER arcs/tilts are also have different binning with bias
        # frames, so don't use bias frames. Don't use the biases for any
        # calibrations since it appears to be a different amplifier readout
        par['calibrations']['traceframe']['process']['overscan_method'] = 'median'

        # Right now we are using the overscan and not biases becuase the
        # standards are read with a different read mode and we don't yet have
        # the option to use different sets of biases for different standards,
        # or use the overscan for standards but not for science frames
        # TODO testing
        par['scienceframe']['process']['use_biasimage'] = False
        par['scienceframe']['process']['use_illumflat'] = False
        par['scienceframe']['process']['use_pixelflat'] = False
        par['calibrations']['standardframe']['process']['use_illumflat'] = False
        par['calibrations']['standardframe']['process']['use_pixelflat'] = False
        # par['scienceframe']['useframe'] ='overscan'

        par['calibrations']['slitedges']['edge_thresh'] = 8.0
        par['calibrations']['slitedges']['fit_order'] = 8
        par['calibrations']['slitedges']['max_shift_adj'] = 0.5
        par['calibrations']['slitedges']['trace_thresh'] = 10.
        par['calibrations']['slitedges']['left_right_pca'] = True
        par['calibrations']['slitedges']['length_range'] = 0.3

        # These are the defaults
        par['calibrations']['tilts']['tracethresh'] = 15
        par['calibrations']['tilts']['spat_order'] = 3
        par['calibrations']['tilts']['spec_order'] = 5  # [5, 5, 5] + 12*[7] # + [5]

        # 1D wavelength solution
        par['calibrations']['wavelengths']['lamps'] = ['ThAr_XSHOOTER_VIS']
        # This is for 1x1 binning. TODO GET BINNING SORTED OUT!!
        par['calibrations']['wavelengths']['rms_threshold'] = 0.50
        par['calibrations']['wavelengths']['sigdetect'] = 5.0
        par['calibrations']['wavelengths']['n_final'] = [3] + 13 * [4] + [3]
        # This is for 1x1 binning. Needs to be divided by binning for binned data!!
        par['calibrations']['wavelengths']['fwhm'] = 11.0
        # Reidentification parameters
        par['calibrations']['wavelengths']['method'] = 'full_template'
        # TODO: the arxived solution is for 1x1 binning. It needs to be
        # generalized for different binning!
        par['calibrations']['wavelengths']['reid_arxiv'] = 'vlt_xshooter_vis1x1.fits'
        par['calibrations']['wavelengths']['cc_thresh'] = 0.50
        par['calibrations']['wavelengths']['cc_local_thresh'] = 0.50
        par['calibrations']['wavelengths']['ech_fix_format'] = True
        # Echelle parameters
        par['calibrations']['wavelengths']['echelle'] = True
        par['calibrations']['wavelengths']['ech_nspec_coeff'] = 4
        par['calibrations']['wavelengths']['ech_norder_coeff'] = 4
        par['calibrations']['wavelengths']['ech_sigrej'] = 3.0

        # Flats
        par['calibrations']['flatfield']['tweak_slits_thresh'] = 0.90
        par['calibrations']['flatfield']['tweak_slits_maxfrac'] = 0.10

        # Extraction
        par['reduce']['skysub']['bspline_spacing'] = 0.5
        par['reduce']['skysub']['global_sky_std'] = False
        # local sky subtraction operates on entire slit
        par['reduce']['extraction']['model_full_slit'] = True
        # Mask 3 edges pixels since the slit is short, insted of default (5,5)
        par['reduce']['findobj']['find_trim_edge'] = [3, 3]
        # Continnum order for determining thresholds

        # Sensitivity function parameters
        par['sensfunc']['algorithm'] = 'IR'
        par['sensfunc']['polyorder'] = [9, 11, 11, 9, 9, 8, 8, 7, 7, 7, 7, 7, 7, 7, 7]
        par['sensfunc']['IR']['telgridfile'] = 'TelFit_MaunaKea_3100_26100_R20000.fits'
        return par

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the ``PypeIt``-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=0, card='RA', required_ftypes=['science', 'standard'])
        self.meta['dec'] = dict(ext=0, card='DEC', required_ftypes=['science', 'standard'])
        self.meta['target'] = dict(ext=0, card='OBJECT')
        self.meta['decker'] = dict(ext=0, card='DECKNAME')
        self.meta['binning'] = dict(card=None, compound=True)
        self.meta['mjd'] = dict(ext=0, card='MJD')
        self.meta['exptime'] = dict(ext=0, card='ELAPTIME') # This may depend on the old/new detector
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        #self.meta['dispname'] = dict(ext=0, card='ECHNAME')
        # Extras for config and frametyping
        self.meta['hatch'] = dict(ext=0, card='HATOPEN')
        self.meta['dispname'] = dict(ext=0, card='XDISPERS')
        self.meta['filter1'] = dict(ext=0, card='FIL1NAME')
        self.meta['echangle'] = dict(ext=0, card='ECHANGL', rtol=1e-3)
        self.meta['xdangle'] = dict(ext=0, card='XDANGL', rtol=1e-3)
        self.meta['idname'] = dict(ext=0, card='IMAGETYP')
        self.meta['frameno'] = dict(ext=0, card='FRAMENO')
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
            # TODO JFH Is this correct or should it be flipped?
            binspatial, binspec = parse.parse_binning(headarr[0]['BINNING'])
            binning = parse.binning2string(binspec, binspatial)
            return binning
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
        return ['filter1', 'echangle', 'xdangle']


    def pypeit_file_keys(self):
        """
        Define the list of keys to be output into a standard ``PypeIt`` file.

        Returns:
            :obj:`list`: The list of keywords in the relevant
            :class:`~pypeit.metadata.PypeItMetaData` instance to print to the
            :ref:`pypeit_file`.
        """
        return super().pypeit_file_keys() + ['frameno']




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
        # TODO: Allow for 'sky' frame type, for now include sky in
        # 'science' category
        if ftype == 'science':
            return good_exp & (fitstbl['idname'] == 'Object')
        if ftype == 'standard':
            return good_exp & ((fitstbl['idname'] == 'Std') | (fitstbl['idname'] == 'Object'))
        if ftype == 'bias':
            return good_exp & (fitstbl['idname'] == 'Bias')
        if ftype == 'dark':
            return good_exp & (fitstbl['idname'] == 'Dark')
        if ftype in ['pixelflat', 'trace']:
            # Flats and trace frames are typed together
            return good_exp & ((fitstbl['idname'] == 'Flat') | (fitstbl['idname'] == 'IntFlat'))
        if ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['idname'] == 'Line')

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)


    def get_rawimage(self, raw_file, det, spectrim=20):
        """
        Read raw images and generate a few other bits and pieces
        that are key for image processing.

        Based on readmhdufits.pro

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
        # TODO -- Put a check in here to avoid data using the
        #  original CCD (1 chip)


        # Check for file; allow for extra .gz, etc. suffix
        if not os.path.isfile(raw_file):
            msgs.error(f'{raw_file} not found!')
        hdu = io.fits_open(raw_file)

        head0 = hdu[0].header

        # Get post, pre-pix values
        precol = head0['PRECOL']
        postpix = head0['POSTPIX']
        preline = head0['PRELINE']
        postline = head0['POSTLINE']
        detlsize = head0['DETLSIZE']
        x0, x_npix, y0, y_npix = np.array(parse.load_sections(detlsize)).flatten()

        # get the x and y binning factors...
        #binning = head0['BINNING']

        binning = self.get_meta_value(self.get_headarr(hdu), 'binning')
        # TODO: JFH I think this works fine
        if binning != '3,1':
            msgs.warn("This binning for HIRES might not work.  But it might..")

        # We are flipping this because HIRES stores the binning oppostire of the (binspec, binspat) pypeit convention.
        binspatial, binspec = parse.parse_binning(head0['BINNING'])
        # Validate the entered (list of) detector(s)
        nimg, _det = self.validate_det(det)

        # Grab the detector or mosaic parameters
        mosaic = None if nimg == 1 else self.get_mosaic_par(det, hdu=hdu)
        detectors = [self.get_detector_par(det, hdu=hdu)] if nimg == 1 else mosaic.detectors

        # get the chips to read in
        # DP: I don't know if this needs to still exist. I believe det is never None
        if det is None:
            chips = range(self.ndet)
        else:
            chips = [d-1 for d in _det]  # Indexing starts at 0 here

        # get final datasec and oscan size (it's the same for every chip so
        # it's safe to determine it outsize the loop)

        # Create final image
        if det is None:
            # JFH: TODO is this a good idea?
            image = np.zeros((x_npix, y_npix + 4 * postpix))
            rawdatasec_img = np.zeros_like(image, dtype=int)
            oscansec_img = np.zeros_like(image, dtype=int)
        else:
            data, oscan = hires_read_1chip(hdu, chips[0] + 1)
            image = np.zeros((nimg, data.shape[0], data.shape[1] + oscan.shape[1]))
            rawdatasec_img = np.zeros_like(image, dtype=int)
            oscansec_img = np.zeros_like(image, dtype=int)


        # Loop over the chips
        for ii, tt in enumerate(chips):
            image_ii, oscan_ii = hires_read_1chip(hdu, tt + 1)

            # Indexing
            x1, x2, y1, y2, o_x1, o_x2, o_y1, o_y2 = indexing(tt, postpix, det=det, xbin=binspatial, ybin=binspec)

            # Fill
            image[ii, y1:y2, x1:x2] = image_ii
            image[ii, o_y1:o_y2, o_x1:o_x2] = oscan_ii
            rawdatasec_img[ii, y1:y2-spectrim//binspec, x1:x2] = 1  # Amp
            oscansec_img[ii, o_y1:o_y2-spectrim//binspec, o_x1:o_x2] = 1  # Amp

        exptime = hdu[self.meta['exptime']['ext']].header[self.meta['exptime']['card']]

        # Return
        # Handle returning both single and multiple images
        if nimg == 1:
            return detectors[0], image[0], hdu, exptime, rawdatasec_img[0], oscansec_img[0]
        return mosaic, image, hdu, exptime, rawdatasec_img, oscansec_img


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

        # Validate the entered (list of) detector(s)
        nimg, _ = self.validate_det(mosaic)

        # Index of mosaic in list of allowed detector combinations
        mosaic_id = self.allowed_mosaics.index(mosaic)+1
        detid = f'MSC0{mosaic_id}'

        # Get the detectors
        detectors = np.array([self.get_detector_par(det, hdu=hdu) for det in mosaic])
        # Binning *must* be consistent for all detectors
        if any(d.binning != detectors[0].binning for d in detectors[1:]):
            msgs.error('Binning is somehow inconsistent between detectors in the mosaic!')

        # Collect the offsets and rotations for *all unbinned* detectors in the
        # full instrument, ordered by the number of the detector.  Detector
        # numbers must be sequential and 1-indexed.
        # See the mosaic documentattion.
        msc_geometry = HIRESMosaicLookUp.geometry
        expected_shape = msc_geometry[detid]['default_shape']
        shift = np.array([(msc_geometry[detid]['blue_det']['shift'][0], msc_geometry[detid]['blue_det']['shift'][1]),
                          (msc_geometry[detid]['green_det']['shift'][0], msc_geometry[detid]['green_det']['shift'][1]),
                          (msc_geometry[detid]['red_det']['shift'][0], msc_geometry[detid]['red_det']['shift'][1])])

        rotation = np.array([msc_geometry[detid]['blue_det']['rotation'], msc_geometry[detid]['green_det']['rotation'],
                             msc_geometry[detid]['red_det']['rotation']])

        # The binning and process image shape must be the same for all images in
        # the mosaic
        binning = tuple(int(b) for b in detectors[0].binning.split(','))
        shape = tuple(n // b for n, b in zip(expected_shape, binning))

        msc_sft = [None]*nimg
        msc_rot = [None]*nimg
        msc_tfm = [None]*nimg

        for i in range(nimg):
            msc_sft[i] = shift[i]
            msc_rot[i] = rotation[i]
            # binning is here in the PypeIt convention of (binspec, binspat), but the mosaic tranformations
            # occur in the raw data frame, which flips spectral and spatial
            msc_tfm[i] = build_image_mosaic_transform(shape, msc_sft[i], msc_rot[i], tuple(reversed(binning)))

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
        binning = '1,1' if hdu is None else self.get_meta_value(self.get_headarr(hdu), 'binning')

        # Detector 1
        # TODO -- Fill these in properly for HIRES
        detector_dict1 = dict(
            binning         = binning,
            det             = 1,
            dataext         = 1,
            specaxis        = 0,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.135,
            darkcurr        = 0.0,
            saturation      = 65535.,
            nonlinear       = 0.86,
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d([1.]),
            ronoise         = np.atleast_1d([3.]),
            )

        # Detector 2. HACK
        detector_dict2 = detector_dict1.copy()
        detector_dict2.update(dict(
            det=2,
            dataext=2,
            gain=np.atleast_1d([1.]),
            ronoise=np.atleast_1d([3.])
        ))


        # Detector 3,. HACK
        detector_dict3 = detector_dict1.copy()
        detector_dict3.update(dict(
            det=3,
            dataext=3,
            gain=np.atleast_1d([1.]),
            ronoise=np.atleast_1d([3.])
        ))

        # Instantiate
        detector_dicts = [detector_dict1, detector_dict2, detector_dict3]
        detector = detector_container.DetectorContainer(
            **detector_dicts[det-1])

        if hdu is None:
            return detector

        # Return
        return detector

# TODO: This spectrograph is way out of date
#    def load_raw_img_head(self, raw_file, det=None, **null_kwargs):
#        """
#        Wrapper to the raw image reader for HIRES
#
#        Args:
#            raw_file (:obj:`str`):
#                filename
#            det (:obj:`int`, optional):
#              Desired detector.  Despite default value, cannot be
#              ``None`` (todo: set a sensible default).
#            **null_kwargs:
#              Captured and never used
#
#        Returns:
#            tuple: Raw image and header
#
#        """
#        raw_img, head0, _ = read_hires(raw_file, det=det)
#
#        return raw_img, head0
#
#    def get_image_section(self, inp=None, det=1, section='datasec'):
#        """
#        Return a string representation of a slice defining a section of
#        the detector image.
#
#        Overwrites base class function to use :func:`read_hires` to get
#        the image sections.
#
#        .. todo ::
#            - It is really ineffiecient.  Can we parse
#              :func:`read_hires` into something that can give you the
#              image section directly?
#
#        This is done separately for the data section and the overscan
#        section in case one is defined as a header keyword and the other
#        is defined directly.
#
#        Args:
#            inp (:obj:`str`, `astropy.io.fits.Header`_, optional):
#                String providing the file name to read, or the relevant
#                header object.  Default is None, meaning that the
#                detector attribute must provide the image section
#                itself, not the header keyword.
#            det (:obj:`int`, optional):
#                1-indexed detector number.
#            section (:obj:`str`, optional):
#                The section to return.  Should be either 'datasec' or
#                'oscansec', according to the
#                :class:`pypeitpar.DetectorPar` keywords.
#
#        Returns:
#            tuple: Returns three objects: (1) A list of string
#            representations for the image sections, one string per
#            amplifier.  The sections are *always* returned in PypeIt
#            order: spectral then spatial.  (2) Boolean indicating if the
#            slices are one indexed.  (3) Boolean indicating if the
#            slices should include the last pixel.  The latter two are
#            always returned as True following the FITS convention.
#        """
#        # Read the file
#        if inp is None:
#            msgs.error('Must provide Keck HIRES file to get image section.')
#        elif not os.path.isfile(inp):
#            msgs.error('File {0} does not exist!'.format(inp))
#        temp, head0, secs = read_hires(inp, det)
#        if section == 'datasec':
#            return secs[0], False, False
#        elif section == 'oscansec':
#            return secs[1], False, False
#        else:
#            raise ValueError('Unrecognized keyword: {0}'.format(section))
#
#     def get_datasec_img(self, filename, det=1, force=True):
#         """
#         Create an image identifying the amplifier used to read each pixel.
#
#         Args:
#             filename (str):
#                 Name of the file from which to read the image size.
#             det (:obj:`int`, optional):
#                 Detector number (1-indexed)
#             force (:obj:`bool`, optional):
#                 Force the image to be remade
#
#         Returns:
#             `numpy.ndarray`: Integer array identifying the amplifier
#             used to read each pixel.
#         """
#         if self.datasec_img is None or force:
#             # Check the detector is defined
#             self._check_detector()
#             # Get the image shape
#             raw_naxis = self.get_raw_image_shape(filename, det=det)
#
#             # Binning is not required because read_hires accounts for it
# #            binning = self.get_meta_value(filename, 'binning')
#
#             data_sections, one_indexed, include_end, transpose \
#                     = self.get_image_section(filename, det, section='datasec')
#
#             # Initialize the image (0 means no amplifier)
#             self.datasec_img = np.zeros(raw_naxis, dtype=int)
#             for i in range(self.detector[det-1]['numamplifiers']):
#                 # Convert the data section from a string to a slice
#                 datasec = parse.sec2slice(data_sections[i], one_indexed=one_indexed,
#                                           include_end=include_end, require_dim=2,
#                                           transpose=transpose) #, binning=binning)
#                 # Assign the amplifier
#                 self.datasec_img[datasec] = i+1
#         return self.datasec_img

# This is deprecated.
# class KECKHIRESRSpectrograph(KECKHIRESSpectrograph):
#     """
#     Child to handle KECK/HIRES-R specific code
#
#     .. warning::
#         Spectrograph not yet supported
#     """
#     name = 'keck_hires_red'
#     camera = 'HIRES_R'
#
#     @classmethod
#     def default_pypeit_par(cls):
#         """
#         Return the default parameters to use for this instrument.
#
#         Returns:
#             :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
#             all of ``PypeIt`` methods.
#         """
#         par = super().default_pypeit_par()
#
#         # Adjustments to slit and tilts for NIR
#         par['calibrations']['slitedges']['edge_thresh'] = 600.
#         par['calibrations']['slitedges']['fit_order'] = 5
#         par['calibrations']['slitedges']['max_shift_adj'] = 0.5
#         par['calibrations']['slitedges']['left_right_pca'] = True
#
#         par['calibrations']['tilts']['tracethresh'] = 20
#         # Bias
#         par['calibrations']['biasframe']['useframe'] = 'bias'
#
#         # 1D wavelength solution
#         par['calibrations']['wavelengths']['lamps'] = ['ThAr']
#         #par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
#         par['calibrations']['wavelengths']['rms_threshold'] = 0.25
#         par['calibrations']['wavelengths']['sigdetect'] = 5.0
#         # Reidentification parameters
#         #par['calibrations']['wavelengths']['method'] = 'reidentify'
#         #par['calibrations']['wavelengths']['reid_arxiv'] = 'vlt_xshooter_nir.json'
#         par['calibrations']['wavelengths']['ech_fix_format'] = True
#         # Echelle parameters
#         par['calibrations']['wavelengths']['echelle'] = True
#         par['calibrations']['wavelengths']['ech_nspec_coeff'] = 4
#         par['calibrations']['wavelengths']['ech_norder_coeff'] = 4
#         par['calibrations']['wavelengths']['ech_sigrej'] = 3.0
#
#         # Always correct for flexure, starting with default parameters
#         par['flexure'] = pypeitpar.FlexurePar()
#         par['scienceframe']['process']['sigclip'] = 20.0
#         par['scienceframe']['process']['satpix'] ='nothing'
#         par['calibrations']['standardframe']['exprng'] = [None, 600]
#         par['scienceframe']['exprng'] = [600, None]
#
#         return par



def indexing(itt, postpix, det=None,xbin=1,ybin=1):
    """
    Some annoying book-keeping for instrument placement.

    Parameters
    ----------
    itt : int
    postpix : int
    det : int, optional

    Returns
    -------

    """
    # Deal with single chip
    if det is not None:
        tt = 0
    else:
        tt = itt
    ii = int(np.round(2048/xbin))
    jj = int(np.round(4096/ybin))
    # y indices
    y1, y2 = 0, jj
    o_y1, o_y2 = y1, y2

    # x
    x1, x2 = (tt%4)*ii, (tt%4 + 1)*ii
    if det is None:
        o_x1 = 4*ii + (tt%4)*postpix
    else:
        o_x1 = ii + (tt%4)*postpix
    o_x2 = o_x1 + postpix

    # Return
    return x1, x2, y1, y2, o_x1, o_x2, o_y1, o_y2

def hires_read_1chip(hdu,chipno):
    """ Read one of the HIRES detectors

    Parameters
    ----------
    hdu : HDUList
    chipno : int

    Returns
    -------
    data : ndarray
    oscan : ndarray
    """

    # Extract datasec from header
    datsec = hdu[chipno].header['DATASEC']
    detsec = hdu[chipno].header['DETSEC']
    postpix = hdu[0].header['POSTPIX']
    precol = hdu[0].header['PRECOL']

    x1_dat, x2_dat, y1_dat, y2_dat = np.array(parse.load_sections(datsec)).flatten()
    x1_det, x2_det, y1_det, y2_det = np.array(parse.load_sections(detsec)).flatten()

    # This rotates the image to be increasing wavelength to the top
    #data = np.rot90((hdu[chipno].data).T, k=2)
    #nx=data.shape[0]
    #ny=data.shape[1]

    # Science data
    fullimage = hdu[chipno].data
    data = fullimage[x1_dat:x2_dat,y1_dat:y2_dat]

    # Overscan
    oscan = fullimage[:,y2_dat:]

    # Flip as needed
    if x1_det > x2_det:
        data = np.flipud(data)
        oscan = np.flipud(oscan)
    if y1_det > y2_det:
        data = np.fliplr(data)
        oscan = np.fliplr(oscan)

    # Return
    return data, oscan

# THIS is deprecated
def read_hires(hdu, det=None, spectrim=20):
    """
    Read a raw HIRES data frame (one or more detectors).

    Data are unpacked from the multi-extension HDU.  Function is
    based :func:`pypeit.spectrographs.keck_lris.read_lris`, which
    was based on the IDL procedure ``readmhdufits.pro``.
    
    Parameters
    ----------
    raw_file : str
        Filename
    spectrim: int, optional
       Number of unbinned pixels to trim off spectral direction because top side of detector is not illuiminated.
       Default is 20.



    Returns
    -------
    image : ndarray
        Raw image including overscan
    header : FITS header
    rawdatasec_img : `numpy.ndarray`_
        Data (Science) section of the detector as provided by setting the
        (1-indexed) number of the amplifier used to read each detector pixel.
        Pixels unassociated with any amplifier are set to 0.
    oscansec_img : `numpy.ndarray`_
        Overscan section of the detector as provided by setting the
        (1-indexed) number of the amplifier used to read each detector pixel.
        Pixels unassociated with any amplifier are set to 0.
    hdu : List of HUDs

    """

    head0 = hdu[0].header

    # Get post, pre-pix values
    precol = head0['PRECOL']
    postpix = head0['POSTPIX']
    preline = head0['PRELINE']
    postline = head0['POSTLINE']
    detlsize = head0['DETLSIZE']
    x0, x_npix, y0, y_npix = np.array(parse.load_sections(detlsize)).flatten()

    # Create final image
    if det is None:
        image = np.zeros((x_npix,y_npix+4*postpix))

    # Setup for datasec, oscansec
    dsec = []
    osec = []

    # get the x and y binning factors...
    binning = head0['BINNING']
    if binning != '3,1':
        msgs.warn("This binning for HIRES might not work.  But it might..")


    binspatial, binspec = parse.parse_binning(head0['BINNING'])

    # HIRES detectors
    nchip = 3


    if det is None:
        chips = range(nchip)
    else:
        chips = [det-1] # Indexing starts at 0 here
    # Loop
    for tt in chips:
        data, oscan = hires_read_1chip(hdu, tt+1)

        # One detector??
        if det is not None:
            image = np.zeros((data.shape[0],data.shape[1]+oscan.shape[1]))

        # Indexing
        x1, x2, y1, y2, o_x1, o_x2, o_y1, o_y2 = indexing(tt, postpix, det=det,xbin=binspatial, ybin=binspec)

        # Fill
        image[y1:y2, x1:x2] = data
        image[o_y1:o_y2, o_x1:o_x2] = oscan

        # Rawdatasec
        rawdatasec_img = np.zeros_like(image, dtype=int)
        rawdatasec_img[y1:y2-spectrim//binspec, x1:x2] = 1
        oscansec_img = np.zeros_like(image, dtype=int)
        oscansec_img[o_y1:o_y2-spectrim//binspec, o_x1:o_x2] = 1

        # Sections
        #idsec = '[{:d}:{:d},{:d}:{:d}]'.format(y1, y2, x1, x2)
        #iosec = '[{:d}:{:d},{:d}:{:d}]'.format(o_y1, o_y2, o_x1, o_x2)
        #dsec.append(idsec)
        #osec.append(iosec)

    # Return
    return image, head0, rawdatasec_img, oscansec_img, hdu


def grab_arctempl_dict(arc_meta:dict, det, ORDRS=None):
    # Read template file
    templ_table_file = os.path.join(
        resource_filename('pypeit', 'data'), 'arc_lines',
        'hires', 'hires_templ.dat')
    tbl = Table.read(templ_table_file, format='ascii')
    
    gd = tbl['XDISP'] == arc_meta['dispname']
    cut_tbl = tbl[gd]

    # Unpack for convenience
    chip = cut_tbl['Chip']
    echa = cut_tbl['ECH']
    xda = cut_tbl['XDAng']

    # Match
    # ;; Closest XDANGL irrespective of binning
    if ORDRS is None:
        # Tolerances on ECHA and XDA
        tols = [[0.05, 0.2], [0.05, 0.4], [0.1, 0.4]]
        for tol in tols:
            #; Best = EDANGL LT 0.05 and XDANGL LT 0.2
            mtch = np.where((np.abs(arc_meta['echangle']-echa) < tol[0]) &
                   (np.abs(arc_meta['xdangle']-xda) < tol[1]) &
                   (det == chip) )[0]
            if len(mtch) > 0:
                break
        if len(mtch) == 0:
            mtch = np.where(det == chip)[0]
    else: 
        msgs.error("NEED TO PORT THE IDL CODE BELOW")
    '''
      ;; Specified orders
      gdo = where(min(ordrs) GE ordi AND $
                  max(ordrs) LE ordf, ngdo)
      if ngdo EQ 0 then begin
         mtch = where(hires.cross EQ xdisp, nmt)
         gdo = where(min(ordrs) GE ordi[mtch] AND $
                     (max(ordrs) < maxo) LE ordf[mtch], ngdo)
         if ngdo EQ 0 then begin
            print, 'hires_arctempl:  No archived wavelengths fitting your' + $
                   'setup.  Contact JXP ASAP (xavier@ucolick.org)!'
            stop
         endif
         mtch = mtch[gdo]
      endif else mtch = gdo 
    '''

    #; Closet ECHANGL
    imn = np.argmin(np.abs(arc_meta['echangle']-echa[mtch]))
    allx = np.where(np.abs(echa[mtch]-echa[mtch[imn]]) < 0.001)[0]

    if len(allx) != 1: #;; Close XDANGL
        imne = np.argmin(np.abs(arc_meta['xdangle']-xda[mtch[allx]]))
        idx = mtch[allx[imne]]
    else:
        idx = mtch[allx[0]]

    # Return the filename (without path)
    # Return
    return dict(cut_tbl[idx])

def load_hires_template(template_file:str):
    dat_dict = readsav(template_file)

    # Chop down to good orders
    n_orders = len(dat_dict['guess_ordr'])
    # Return
    return dat_dict['guess_ordr'][::-1], dat_dict['sv_aspec'][:n_orders,:][::-1, :].T