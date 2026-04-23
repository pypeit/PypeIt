"""
Module for JWST NIRSpec specific methods.

.. include:: ../include/links.rst
"""
import copy

import numpy as np
from astropy.io import fits

from pypeit import log
from pypeit import PypeItError
from pypeit import telescopes
from pypeit import utils
from pypeit import io
from pypeit.core import framematch
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container
from pypeit.images.mosaic import Mosaic
from pypeit.core.mosaic import build_image_mosaic_transform
from IPython import embed


class JWSTNIRSpecSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle JWST NIRSpec specific code
    """
    ndet = 2
    name = 'jwst_nirspec'
    header_name = 'NIRSPEC'
    telescope = telescopes.JWSTTelescopePar()
    camera = 'NIRSPEC'
    url = 'https://jwst-docs.stsci.edu/jwst-near-infrared-spectrograph'
    pypeline = 'NIRSpecSlit'
    #supported = True
    allowed_extensions = ['rate.fits','rate.fits.gz' , 'uncal.fits.gz', 'uncal.fits', '.fits', '.fits.gz']

    def rawfile_basename(self, filename, targname=None, slitname=None):
        """
        Return the basename of a raw file, which is used for naming output
        files by the function :func:`~pypeit.metadata.construct_basename`.
        This can be spectrograph-dependent if specific changes need to be made.

        Here we strip ``_assign_wcs``, ``_nrs1``, and ``_nrs2`` suffixes that
        are typical of JWST NIRSpec data products.

        Args:
            filename (:obj:`str`):
                Input raw fits filename.
            slitname (:obj:`str`, optional):
                Slit name to be added in the basename for per-slit outputs.
                If None, no slit name will be added.

        Returns:
            :obj:`str`:
            The basename of the input file.

        """

        _filename = filename.split('.fits')[0]

        for tag in ['_assign_wcs', '_nrs1', '_nrs2']:
            _filename = _filename.replace(tag, '')

        if targname is not None:
            _filename = _filename + '-' + targname

        # Embed slit name in basename (for per-slit outputs)
        if slitname is not None:
            _filename = _filename + f'_{slitname}'

        return _filename

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

        # Detector 1, i.e. NRS1 from
        # https://jwst-docs.stsci.edu/jwst-near-infrared-spectrograph/nirspec-instrumentation/nirspec-detectors/nirspec-detector-performance
        detector_dict1 = dict(
            binning='1,1',
            det=1,
            dataext=1,
            specaxis=1,
            specflip=False,
            spatflip=False,
            xgap=180.,
            ygap=0.,
            ysize=1.,
            platescale=0.1,
            darkcurr=33.12,  # e-/pixel/hour  (=0.0092 e-/pixel/s)
            saturation=55100.,
            nonlinear=0.95,  # need to look up and update
            mincounts=-1e10,
            numamplifiers=1,
            gain=np.atleast_1d(0.996),
            ronoise=np.atleast_1d(5.17),
            datasec=np.atleast_1d('[:,:]'),
            oscansec=None
        )

        # Detector 2
        detector_dict2 = detector_dict1.copy()
        detector_dict2.update(dict(
            det=2,
            dataext=1,
            darkcurr=20.52,  # e-/pixel/hour,  (=0.0057 e-/pixel/s)
            saturation=60400.,
            gain=np.atleast_1d(1.137),
            ronoise=np.atleast_1d(6.60),
        ))
        detector_dicts = [detector_dict1, detector_dict2]
        return detector_container.DetectorContainer(**detector_dicts[det-1])

    def get_mosaic_par(self, mosaic, hdu=None, msc_ord=0):
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
            msc_ord (:obj:`int`, optional):
                Order of the interpolation used to construct the mosaic.

        Returns:
            :class:`~pypeit.images.mosaic.Mosaic`: Object with the mosaic *and*
            detector parameters.
        """

        # Validate the entered (list of) detector(s)
        nimg, _ = self.validate_det(mosaic)

        # Index of mosaic in list of allowed detector combinations
        mosaic_id = self.allowed_mosaics.index(mosaic)+1

        # Get the detectors
        detectors = np.array([self.get_detector_par(det, hdu=hdu) for det in mosaic])
        # TODO: Implement proper mosaic geometry for NIRSpec NRS1+NRS2.

        # expected_shape = (2048, 4608)
        # shift = np.array(shift = [(0., 0.), (0, self.image[0].shape[0]+ int(self.get_detector_par(1).xgap))])
        expected_shape = (1031, 25)
        shift = [(0.0, 0.0), (0, 605)]
        rotation = np.array([0., 0.])

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
            # binning is here in the PypeIt convention of (binspec, binspat), but the mosaic tranformations
            # occur in the raw data frame, which has spatial in y and spectral in x
            msc_tfm[i] = build_image_mosaic_transform(shape, msc_sft[i], msc_rot[i], tuple(reversed(binning)))


        # For now, use placeholder values; the actual mosaic is built by
        # make_mosaic() using slit_slices-based spatial offsets.
        return Mosaic(mosaic_id, detectors, shape, np.array(msc_sft), np.array(msc_rot),
                      np.array(msc_tfm), msc_ord)

    # def validate_det(self, det):
    #     """
    #     Validate the detector specification and return the number of images
    #     and a standardized detector tuple.
    #
    #     Args:
    #         det (:obj:`int`, :obj:`tuple`):
    #             Detector number or tuple of detector numbers.
    #
    #     Returns:
    #         :obj:`tuple`: Number of images and the validated detector
    #         specification.
    #     """
    #     if isinstance(det, (int, np.integer)):
    #         return 1, (det,)
    #     elif isinstance(det, (tuple, list)):
    #         return len(det), tuple(det)
    #     else:
    #         raise PypeItError(f'Invalid detector specification: {det}')

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.

        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        # turn_off = dict(use_biasimage=False, use_overscan=False, use_darkimage=False, use_illumflat=False)
        turn_off = dict(use_illumflat=False, use_biasimage=False, use_overscan=False, use_darkimage=False,
                        subtract_scattlight=False, spat_flexure_correct=False, use_pixelflat=True,
                        use_specillum=False, apply_gain=False, trim=True)
        par.reset_all_processimages_par(**turn_off)

        # Reduce
        par['reduce']['trim_edge'] = [0,0]

        # Object finding
        par['reduce']['findobj']['find_trim_edge'] = [0,0]
        par['reduce']['findobj']['maxnumber_sci'] = 2
        par['reduce']['findobj']['snr_thresh'] = 10.0
        par['reduce']['findobj']['trace_npoly'] = 5
        par['reduce']['findobj']['find_fwhm'] = 2.0

        # Sky-subtraction
        par['reduce']['skysub']['bspline_spacing'] = 5.0 # JWST sky is smooth
        par['reduce']['skysub']['max_mask_frac'] = 0.95
        par['reduce']['skysub']['mask_by_boxcar'] = True
        par['reduce']['skysub']['sky_sigrej'] = 4.0

        # Extraction
        par['reduce']['extraction']['model_full_slit'] = True
        par['reduce']['extraction']['sn_gauss'] = 5.0
        par['reduce']['extraction']['boxcar_radius'] = 0.2 # extent in calwebb is 0.55" source and on NIRSpec website
        par['reduce']['extraction']['use_2dmodel_mask'] = False # Don't use 2d mask in local skysub

        # Cosmic ray rejection parameters for science frames
        par['scienceframe']['process']['sigclip'] = 5.0
        par['scienceframe']['process']['objlim'] = 2.0
        par['scienceframe']['process']['mask_cr'] = False # Turn off for now since we coadd.

        # identify science frames
        par['scienceframe']['exprng'] = [0.1, None]

        # Skip reference frame correction for now.
        par['calibrations']['wavelengths']['refframe'] = 'observed'

        return par

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the ``PypeIt``-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=0, card='TARG_RA')
        self.meta['dec'] = dict(ext=0, card='TARG_DEC')
        self.meta['target'] = dict(ext=0, card=None, compound=True)
        self.meta['mode'] = dict(ext=0, card='EXP_TYPE')
        self.meta['decker'] = dict(ext=0, card='APERNAME')

        self.meta['binning'] = dict(ext=0, card=None, default='1,1')
        self.meta['mjd'] = dict(ext=0, card='EXPMID')
        self.meta['exptime'] = dict(ext=0, card='EFFEXPTM')
        self.meta['airmass'] = dict(ext=0, card=None, compound=True)

        # Extras for config and frametyping
        self.meta['dispname'] = dict(ext=0, card='GRATING')
        self.meta['filter1'] = dict(ext=0, card='FILTER')
        self.meta['idname'] = dict(ext=0, card=None, compound=True)
        self.meta['dithpat'] = dict(ext=0, card=None, compound=True)
        self.meta['dithpos'] = dict(ext=0, card=None, compound=True)
        self.meta['dithoff'] = dict(ext=0, card=None, compound=True)

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
        # get the filename if available, it is used in several cases
        fname = headarr[0].get('FILENAME')

        if meta_key == 'airmass':
            return 0.

        if meta_key == 'dithpos':
            # the dither information is stored in all frame type, but we want it only for science frames
            if fname is not None and '_assign_wcs' not in fname:
                return None
            return headarr[0].get('PATT_NUM')

        elif meta_key == 'dithoff':
            # the dither information is stored in all frame type, but we want it only for science frames
            if fname is not None and '_assign_wcs' not in fname:
                return 0
            return round(headarr[0].get('YOFFSET'), 3) if headarr[0].get('YOFFSET') is not None else 0.0

        if meta_key == 'dithpat':
            # the dither information is stored in all frame type, but we want it only for science frames
            if fname is not None and '_assign_wcs' not in fname:
                return None
            exp_type = headarr[0].get('EXP_TYPE')
            if exp_type == 'NRS_MSASPEC':
                return headarr[0].get('NOD_TYPE')
            elif exp_type == 'NRS_FIXEDSLIT':
                return headarr[0].get('PATTTYPE')
            else:
                log.warning(f'Cannot determine dithering pattern for EXP_TYPE={exp_type}.')
                return None
        elif meta_key == 'target':
            # The target name is stored in the TARGPROP header card
            # for all NIRSpec modes.
            return headarr[0].get('TARGNAME') if headarr[0].get('TARGNAME') is not None and \
                                                 len(headarr[0].get('TARGNAME').strip()) > 0 \
                                             else headarr[0].get('TARGPROP')
        elif meta_key == 'idname':
            if fname is not None:
                if '_assign_wcs' in fname:
                    return 'science'
                elif '_interpolatedflat' in fname:
                    return 'interpolatedflat'
                elif '_cal' in fname:
                    return 'calib'
                else:
                    return None
            else:
                log.warning("Cannot determine idname from header. Setting to None.")
                return None
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
        return ['dispname', 'filter1', 'decker', 'target']

    def pypeit_file_keys(self):
        """
        Define the list of keys to be output into a standard ``PypeIt`` file.

        Returns:
            :obj:`list`: The list of keywords in the relevant
            :class:`~pypeit.metadata.PypeItMetaData` instance to print to the
            :ref:`pypeit_file`.
        """
        pypeit_keys = super().pypeit_file_keys()
        pypeit_keys.remove('airmass')
        pypeit_keys.remove('binning')
        return pypeit_keys + ['dithpat', 'dithpos', 'dithoff']


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
            return good_exp & (fitstbl['idname'] == 'science')
        if ftype in ['pixelflat']:
            return good_exp & (fitstbl['idname'] == 'interpolatedflat')
        if ftype in ['arc', 'tilt', 'trace']:
            return good_exp & (fitstbl['idname'] == 'calib')
        log.debug('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)



    @property
    def allowed_mosaics(self):
        """
        Return the list of allowed detector mosaics.

        JWST/NIRSpec only allows for mosaicing the NRS1 and NRS2 detectors.

        Returns:
            :obj:`list`: List of tuples, where each tuple provides the 1-indexed
            detector numbers that can be combined into a mosaic and processed by
            ``PypeIt``.
        """
        return [(1,2)]

    @property
    def default_mosaic(self):
        return self.allowed_mosaics[0]
    
    def get_rawimage(self, raw_file, det, extname='SCI'):
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

        raw_files = np.atleast_1d(raw_file)

        file_list = [str(utils.find_single_file(f'{rfile}*', required=True)) for rfile in raw_files]
        for fil in file_list:
             # Check extension
            self._check_extensions(fil)
        # Read
        log.info(f'Reading JWST/NIRSpec file(s): \n{"\n".join(file_list)}')

        # get hdul, headarr, and exptime (we use one of the frames to get this info,
        # since they should all be the same for the frames in a mosaic)
        hdul = fits.open(file_list[0])
        headarr = self.get_headarr(hdul)
        # Exposure time (used by RawImage)
        exptime = self.get_meta_value(headarr, 'exptime')


        # Validate the entered (list of) detector(s)
        nimg, _det = self.validate_det(det)

        # Grab the detector or mosaic parameters
        mosaic = None if nimg == 1 else self.get_mosaic_par(det, hdu=None)
        detectors = [self.get_detector_par(det, hdu=None)] if nimg == 1 else mosaic.detectors

        # Read the image(s)
        # check that the number of files found matches the number of detectors in the mosaic
        if nimg==2 and len(file_list) != nimg:
            raise PypeItError(f'Expected {nimg} files for mosaic with detectors {det}, but found {len(file_list)} files: {file_list}')
        raw_img = [None]*nimg
        rawdatasec_img = [None]*nimg
        oscansec_img = [None]*nimg
        for i in range(nimg):
            indx = np.where([f'nrs{detectors[i].det}' in n for n in file_list])[0]
            if len(indx) == 0:
                raise PypeItError(f'Could not find file for detector nrs{detectors[i].det} in mosaic with detectors {det}.')
            _hdu = io.fits_open(file_list[indx[0]], ignore_missing_end=True, output_verify='ignore', ignore_blank=True)
            # Raw image
            raw_img[i] = _hdu[extname].data.astype(float)
            rawdatasec_img[i] = np.zeros_like(raw_img[i], dtype=int) + int(detectors[i].det)
            oscansec_img = np.zeros_like(raw_img[i], dtype=int)

        if nimg == 1:
            return detectors[0], raw_img[0], hdul, exptime, rawdatasec_img[0], oscansec_img[0]
        return mosaic, np.array(raw_img), hdul, exptime, np.array(rawdatasec_img), np.array(oscansec_img)

    def group_rawfiles(self, raw_files, det=1):
        """
        Group raw files. This can be useful for mosaic reductions, when
        the detectors that need to be mosaiced are in different files.
        This is spectrograph-specific, and it is not defined for all
        spectrographs. The use of this method is not restricted to
        mosaics, but it is expected to be most useful for that.

        Args:
            raw_files (:obj:`list`):
                List of raw files to group.
            det (:obj:`int`, :obj:`tuple`):
                The single detector or set of detectors in a mosaic to process.
        Returns:
            :obj:`list`: List or List of lists of raw files.
                For a non-mosaic reduction, the list would simply be
                a list of individual raw files that match the provided detector.
                For a mosaic reduction, the list would be a list of sublists,
                where the files in each sublist are grouped together for processing.
        """

        if det in [1,2]:
            detname = f'nrs{det}'
            # select only the files that have detname in their name
            return [f for f in raw_files if detname in f]
        if det in self.allowed_mosaics:
            grouped_files = {}
            # group files with the same basename and different detectors
            for _file in raw_files:
                key = str(_file).split('nrs')[0]
                if key not in grouped_files:
                    grouped_files[key] = []
                grouped_files[key].append(_file)
            return list(grouped_files.values())
        else:
            log.warning(f'Detector {det} not supported.')
            return raw_files

    # def make_mosaic(self, img_list, det, slit_slices):
    #     """
    #     Create a mosaic image from the provided list of images.
    #     The images are assumed to be trimmed and oriented to follow
    #     the PypeIt shape convention of ``(nspec,nspat)``.
    #
    #     Args:
    #         img_list (:obj:`list` or `numpy.ndarray`_):
    #             List of images to be combined into a mosaic.  The images must
    #             be trimmed and oriented to follow the PypeIt shape
    #             convention of ``(nspec,nspat)``.
    #         det (:obj:`tuple`):
    #             Tuple of detector numbers used to construct the mosaic.  Must be
    #             one among the list of possible mosaics as hard-coded by the
    #             :func:`allowed_mosaics` function.
    #         slit_slices (:obj:`list`):
    #             List of slices for the slit in the form
    #             ``[(spec_lo, spec_hi), (spat_lo, spat_hi)]`` for each detector
    #             in the mosaic.  The slices are used to determine the spatial
    #             offset between the two detectors in the mosaic.
    #     Returns:
    #         `numpy.ndarray`_: The mosaic image constructed from the provided
    #         list of images.  The image is trimmed and oriented to
    #         follow the PypeIt shape convention of ``(nspec,nspat)``.
    #     """
    #
    #
    #     nimg, _ = self.validate_det(det)
    #
    #     if nimg == 1:
    #         raise PypeItError('Mosaic cannot be made with only one detector!')
    #     else:
    #         if len(img_list) != 2:
    #             raise PypeItError('Mosaic can only be made with two images!')
    #         if det not in self.allowed_mosaics:
    #             raise PypeItError(f'Mosaic with detectors {det} is not allowed! '
    #                        f'Allowed mosaics are: {self.allowed_mosaics}')
    #
    #     detector_gap = int(self.get_detector_par(1).xgap)
    #     spat_offset = (slit_slices[1][1].start - slit_slices[0][1].start)
    #     spec_lo1, spec_hi1 = 0, img_list[0].shape[0]
    #     spec_lo2, spec_hi2 = img_list[0].shape[0] + detector_gap, \
    #                          img_list[0].shape[0] + detector_gap + img_list[1].shape[0]
    #     shape = (img_list[0].shape[0] + img_list[1].shape[0] + detector_gap,
    #              np.max([img_list[0].shape[1], img_list[1].shape[1]]) + np.abs(spat_offset))
    #
    #     if spat_offset >= 0:
    #         # Detector nrs2 starts at a larger spat value than detector nrs1
    #         spat_lo1, spat_hi1 = 0, img_list[0].shape[1]  # nrs1 not shifted spatially
    #         spat_lo2, spat_hi2 = spat_offset, spat_offset + img_list[1].shape[1] # nrs2 shifted spatially
    #     else:
    #         # Detector nrs1 starts at larger spat value than detector nrs1
    #         spat_lo1, spat_hi1 = np.abs(spat_offset), img_list[0].shape[1] + np.abs(spat_offset)
    #         spat_lo2, spat_hi2 = 0, img_list[1].shape[1]
    #
    #
    #     nrs1_slice = np.s_[spec_lo1: spec_hi1, spat_lo1: spat_hi1]
    #     nrs2_slice = np.s_[spec_lo2: spec_hi2, spat_lo2: spat_hi2]
    #
    #     # Create the mosaic
    #     mosaic = np.zeros(shape, dtype=img_list[0].dtype)
    #     mosaic[nrs1_slice] = img_list[0]
    #     mosaic[nrs2_slice] = img_list[1]
    #
    #     return mosaic

    @staticmethod
    def get_slit_slice(cal_data, slit_name):
        """
        Get the pixel slice for a named slit from JWST calibration data.

        Args:
            cal_data: JWST MultiSlitModel calibration data
            slit_name (:obj:`str`): The slit name, e.g. 'S200A1'.

        Returns:
            :obj:`tuple`: A tuple of two slices ``(spec_slice, spat_slice)``.
        """
        slit_names = np.array([slit.name for slit in cal_data.slits])
        if slit_name not in slit_names:
            raise PypeItError(f'Slit name {slit_name} not found in '
                              f'calibration data {cal_data.meta.filename}')
        #TODO: THIS VALUES DO NOT SEEM ACCURATE. FIX IT
        indx = np.where(slit_names == slit_name)[0][0]
        this_slit = cal_data.slits[indx]

        # --- Get slit geometry ---
        slit_xstart = int(this_slit.xstart)
        slit_xsize = int(this_slit.xsize)
        slit_ystart = int(this_slit.ystart)
        slit_ysize = int(this_slit.ysize)

        # Convert slit position to 0-indexed full-frame coordinates
        # xstart = slit_xstart - 1 + subarray_xstart - 1
        # ystart = slit_ystart - 1 + subarray_ystart - 1
        xstart = slit_xstart - 1
        ystart = slit_ystart - 1
        xstop = xstart + slit_xsize
        ystop = ystart + slit_ysize


        # spec_lo = this_slit.xstart - 1
        # spec_hi = spec_lo + this_slit.xsize
        # spat_lo = this_slit.ystart - 1
        # spat_hi = spat_lo + this_slit.ysize

        # spat_start, spat_end, spec_start, spec_end = ystart, ystop, xstart, xstop
        return ystart, ystop, xstart, xstop
