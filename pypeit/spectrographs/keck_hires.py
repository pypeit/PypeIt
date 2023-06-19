"""
Module for Keck/HIRES

.. include:: ../include/links.rst
"""
import os

from IPython import embed



import numpy as np
from scipy.io import readsav

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


class HIRESMosaicLookUp:
    """
    Provides the geometry required to mosaic Keck HIRES data.
    Similar to :class:`~pypeit.spectrographs.gemini_gmos.GeminiGMOSMosaicLookUp`

    """
    # Original
    geometry = {
        'MSC01': {'default_shape': (6168, 3990),
                  'blue_det': {'shift': (-2048.0 - 41.0, 0.0), 'rotation': 0.},
                  'green_det': {'shift': (0., 0.), 'rotation': 0.},
                  'red_det': {'shift': (2048.0 + 53.0, 0.), 'rotation': 0.}},
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
    url = 'https://www2.keck.hawaii.edu/inst/hires/'
    header_name = 'HIRES'
    url = 'https://www2.keck.hawaii.edu/inst/hires/'
    pypeline = 'Echelle'
    ech_fixed_format = False
    supported = False
    # TODO before support = True
    # 1. Implement flat fielding
    # 2. Test on several different setups
    # 3. Implement PCA extrapolation into the blue 


    # TODO: Place holder parameter set taken from X-shooter VIS for now.
    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.

        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        par['rdx']['detnum'] = [(1,2,3)]

        # Adjustments to parameters for Keck HIRES
        turn_on = dict(use_biasimage=False, use_overscan=True, overscan_method='median',
                       use_darkimage=False, use_illumflat=False, use_pixelflat=False,
                       use_specillum=False)
        par.reset_all_processimages_par(**turn_on)
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
        par['calibrations']['slitedges']['max_nudge'] = 0.
        par['calibrations']['slitedges']['overlap'] = True
        par['calibrations']['slitedges']['dlength_range'] = 0.25

        # These are the defaults
        par['calibrations']['tilts']['tracethresh'] = 15
        par['calibrations']['tilts']['spat_order'] = 3
        par['calibrations']['tilts']['spec_order'] = 5  # [5, 5, 5] + 12*[7] # + [5]

        # 1D wavelength solution
        par['calibrations']['wavelengths']['lamps'] = ['ThAr']
        # This is for 1x1 binning. TODO GET BINNING SORTED OUT!!
        par['calibrations']['wavelengths']['rms_threshold'] = 0.50
        par['calibrations']['wavelengths']['sigdetect'] = 5.0
        par['calibrations']['wavelengths']['n_final'] = 4 #[3] + 13 * [4] + [3]
        # This is for 1x1 binning. Needs to be divided by binning for binned data!!
        par['calibrations']['wavelengths']['fwhm'] = 8.0
        # Reidentification parameters
        par['calibrations']['wavelengths']['method'] = 'echelle'
        # TODO: the arxived solution is for 1x1 binning. It needs to be
        # generalized for different binning!
        #par['calibrations']['wavelengths']['reid_arxiv'] = 'vlt_xshooter_vis1x1.fits'
        par['calibrations']['wavelengths']['cc_thresh'] = 0.50
        par['calibrations']['wavelengths']['cc_local_thresh'] = 0.50
#        par['calibrations']['wavelengths']['ech_fix_format'] = True
        # Echelle parameters
        par['calibrations']['wavelengths']['echelle'] = True
        par['calibrations']['wavelengths']['ech_nspec_coeff'] = 4
        par['calibrations']['wavelengths']['ech_norder_coeff'] = 4
        par['calibrations']['wavelengths']['ech_sigrej'] = 3.0
        par['calibrations']['wavelengths']['ech_separate_2d'] = True

        # Flats
        par['calibrations']['flatfield']['tweak_slits_thresh'] = 0.90
        par['calibrations']['flatfield']['tweak_slits_maxfrac'] = 0.10

        # Extraction
        par['reduce']['skysub']['bspline_spacing'] = 0.6
        par['reduce']['skysub']['global_sky_std'] = False
        # local sky subtraction operates on entire slit
        par['reduce']['extraction']['model_full_slit'] = True
        # Mask 3 edges pixels since the slit is short, insted of default (5,5)
        par['reduce']['findobj']['find_trim_edge'] = [3, 3]
        # Continnum order for determining thresholds

        # Sensitivity function parameters
        par['sensfunc']['algorithm'] = 'IR'
        par['sensfunc']['polyorder'] = 11 #[9, 11, 11, 9, 9, 8, 8, 7, 7, 7, 7, 7, 7, 7, 7]
        par['sensfunc']['IR']['telgridfile'] = 'TelFit_MaunaKea_3100_26100_R20000.fits'
        return par

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
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
        # This may depend on the old/new detector
        self.meta['exptime'] = dict(ext=0, card='ELAPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        #self.meta['dispname'] = dict(ext=0, card='ECHNAME')
        # Extras for config and frametyping
        self.meta['hatch'] = dict(ext=0, card='HATOPEN')
        self.meta['dispname'] = dict(ext=0, card='XDISPERS')
        self.meta['filter1'] = dict(ext=0, card='FIL1NAME')
        self.meta['echangle'] = dict(ext=0, card='ECHANGL', rtol=1e-3)
        self.meta['xdangle'] = dict(ext=0, card='XDANGL', rtol=1e-3)
#        self.meta['idname'] = dict(ext=0, card='IMAGETYP')
        # NOTE: This is the native keyword.  IMAGETYP is from KOA.
        self.meta['idname'] = dict(ext=0, card='OBSTYPE')
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
        return ['FIL1NAME', 'ECHANGL', 'XDANGL']

    def pypeit_file_keys(self):
        """
        Define the list of keys to be output into a standard PypeIt file.

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
            return good_exp & (fitstbl['idname'] == 'Object')
        if ftype == 'bias':
            return good_exp & (fitstbl['idname'] == 'Bias')
        if ftype == 'dark':
            return good_exp & (fitstbl['idname'] == 'Dark')
        if ftype in ['pixelflat', 'trace']:
            # Flats and trace frames are typed together
            return good_exp & (fitstbl['idname'] == 'IntFlat')
        if ftype in ['arc', 'tilt']:
            # Arc and tilt frames are typed together
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
#        # TODO: JFH I think this works fine
#        if binning != '3,1':
#            msgs.warn("This binning for HIRES might not work.  But it might..")

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
        the PypeIt shape convention of ``(nspec,nspat)``.  For returned
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

        Keck/HIRES only allows for mosaicing all three detectors.

        Returns:
            :obj:`list`: List of tuples, where each tuple provides the 1-indexed
            detector numbers that can be combined into a mosaic and processed by
            PypeIt.
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
            nonlinear       = 0.7, # Website says 0.6, but we'll push it a bit
            mincounts       = -1e10,
            numamplifiers   = 1,
            ronoise         = np.atleast_1d([2.8]),
            )

        # Detector 2. 
        detector_dict2 = detector_dict1.copy()
        detector_dict2.update(dict(
            det=2,
            dataext=2,
            ronoise=np.atleast_1d([3.1])
        ))


        # Detector 3,. 
        detector_dict3 = detector_dict1.copy()
        detector_dict3.update(dict(
            det=3,
            dataext=3,
            ronoise=np.atleast_1d([3.1])
        ))

        # Set gain 
        # https://www2.keck.hawaii.edu/inst/hires/instrument_specifications.html
        if hdu is None or hdu[0].header['CCDGAIN'].strip() == 'low':
            detector_dict1['gain'] = np.atleast_1d([1.9])
            detector_dict2['gain'] = np.atleast_1d([2.1])
            detector_dict3['gain'] = np.atleast_1d([2.1])
        elif hdu[0].header['CCDGAIN'].strip() == 'high':
            detector_dict1['gain'] = np.atleast_1d([0.78])
            detector_dict2['gain'] = np.atleast_1d([0.86])
            detector_dict3['gain'] = np.atleast_1d([0.84])
        else:
            msgs.error("Bad CCDGAIN mode for HIRES")
            
        # Instantiate
        detector_dicts = [detector_dict1, detector_dict2, detector_dict3]
        return detector_container.DetectorContainer( **detector_dicts[det-1])

    def get_echelle_angle_files(self):
        """ Pass back the files required
        to run the echelle method of wavecalib

        Returns:
            list: List of files
        """
        angle_fits_file = 'keck_hires_angle_fits.fits'
        composite_arc_file = 'keck_hires_composite_arc.fits'

        return [angle_fits_file, composite_arc_file]
        

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
        det = self.get_detector_par(1)
        binspectral, binspatial = parse.parse_binning(binning)

        # Assume no significant variation (which is likely true)
        return np.ones_like(order_vec)*det.platescale*binspatial

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
