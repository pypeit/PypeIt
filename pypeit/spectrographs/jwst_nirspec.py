"""
Module for JWST NIRSpec specific methods.

.. include:: ../include/links.rst
"""
import copy
import functools
import importlib.util
from io import BytesIO
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.time import Time

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
    allowed_extensions = [
        '_assign_wcs.fits', '_interpolatedflat.fits', '_interpolatedflat_fs.fits', '_cal.fits',
        'msa.fits', 'rate.fits', 'rate.fits.gz', 'uncal.fits', 'uncal.fits.gz', '.fits', '.fits.gz'
    ]

    def rawfile_basename(self, filename, targname=None, slitname=None, mjd=None):
        """
        Return the basename of a raw file, which is used for naming output
        files by the function :func:`~pypeit.metadata.construct_basename`.

        The ``_nrs1`` and ``_nrs2`` detector tags are removed so that NRS1 and
        NRS2 files from the same exposure map to the same basename.  The
        result is then passed to :func:`~pypeit.outputfiles.construct_basename`,
        which strips the JWST data-product suffix (e.g., ``_cal.fits``,
        ``_assign_wcs.fits``; see :attr:`allowed_extensions`) and appends the
        target, camera, observation time, and slit name.

        Args:
            filename (:obj:`str`, `Path`_):
                Input raw fits filename.
            targname (:obj:`str`, optional):
                Target name to be added in the basename.
                If None, no target name will be added.
            slitname (:obj:`str`, optional):
                Slit name to be added in the basename for per-slit outputs.
                If None, no slit name will be added.
            mjd (:obj:`float`, optional):
                The MJD of the observation.  If None, no observation time will
                be added.

        Returns:
            :obj:`str`:
            The basename of the input file.

        """
        _filename = Path(filename).name.replace('_nrs1', '').replace('_nrs2', '')
        return super().rawfile_basename(
            _filename, targname=targname, slitname=slitname, mjd=mjd
        )

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
            # NOTE: For NIRSpec, xgap is used as the gap (in pixels, along the
            # dispersion direction) between NRS1 and NRS2 when building the
            # mosaic.  146 px is the separation between the full 2048-pixel
            # arrays (NRS1 column 2048 to NRS2 column 1) given by the CRDS FPA
            # reference file (e.g., jwst_nirspec_fpa_0011.asdf).  Excluding
            # the 4 reference columns on each side, this is the ~2.8 mm (~156
            # px) gap between responsive columns quoted by Jakobsen et al.
            # (2022, A&A, 661, A80).
            xgap=146.,
            ygap=0.,
            ysize=1.,
            # Spatial-direction (along-slit) platescale: 105 mas/pixel
            # (Jakobsen et al. 2022). The platescale in the dispersion
            # direction is 103 mas/pixel.
            platescale=0.105,
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
        Return a placeholder :class:`~pypeit.images.mosaic.Mosaic` object for
        the NRS1+NRS2 detector pair.

        Because the correct mosaic geometry depends on the per-slit bounding
        boxes (available only after reading the data), this method returns a
        :class:`~pypeit.images.mosaic.Mosaic` with **placeholder** shifts and
        transforms computed from a fixed ``(2048, 2048)`` detector shape.
        The placeholder is sufficient to bootstrap the ``RawImage`` and
        calibration pipelines; the actual geometry is overwritten in-place
        once the slit extents are known:

        * :func:`~pypeit.images.rawimage.NIRSpecRawImage.build_mosaic`
          recomputes the transforms from :attr:`~pypeit.images.rawimage.NIRSpecRawImage.slit_slices`
          before calling :func:`~pypeit.core.mosaic.build_image_mosaic`.
        * :func:`~pypeit.calibrations.NIRSpecSlitCalibrations._build_calibimage`
          likewise derives the geometry from ``self.slit_slices`` for each
          calibration image.

        Args:
            mosaic (:obj:`tuple`):
                Tuple of detector numbers used to construct the mosaic.  Must be
                one among the list of possible mosaics as hard-coded by the
                :func:`allowed_mosaics` function.
            hdu (`astropy.io.fits.HDUList`_, optional):
                The open fits file with the raw image of interest.  If not
                provided, frame-dependent detector parameters are set to a
                default.  BEWARE: If ``hdu`` is not provided, the binning is
                assumed to be ``1,1``, which will cause faults if applied to
                binned images!
            msc_ord (:obj:`int`, optional):
                Order of the interpolation used to construct the mosaic.

        Returns:
            :class:`~pypeit.images.mosaic.Mosaic`: Placeholder mosaic object
            whose ``shape``, ``shift``, and ``tform`` attributes will be
            overwritten with the true geometry at processing time.
        """

        # Validate the entered (list of) detector(s)
        nimg, _ = self.validate_det(mosaic)

        # Index of mosaic in list of allowed detector combinations
        mosaic_id = self.allowed_mosaics.index(mosaic)+1

        # Get the detectors
        detectors = np.array([self.get_detector_par(det, hdu=hdu) for det in mosaic])

        expected_shape = (2048, 2048)
        shift = [(0.0, 0.0), (0, 2048.)]
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


        # For now, use placeholder values; the mosaic is updated later using slit_slices-based spatial offsets.
        return Mosaic(mosaic_id, detectors, shape, np.array(msc_sft), np.array(msc_rot), np.array(msc_tfm), msc_ord)

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

        # this is just so the code does not complain if this is not set. This
        # parameter is not really used for NIRSpec. detnum is determined by where
        # the slit that is being reduced is located.
        par['rdx']['detnum'] = [(1,2)]

        # Reduce
        par['reduce']['trim_edge'] = [1,1]

        # Object finding
        par['reduce']['findobj']['find_trim_edge'] = [1,1]
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
        # Skip flexure correction.
        par['flexure']['spec_method'] = 'skip'

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

        # Extras for pypeit file
        self.meta['dateobs'] = dict(card=None, compound=True)
        self.meta['pid'] = dict(ext=0, card='PROGRAM')
        self.meta['obsnum'] = dict(ext=0, card='OBSERVTN')
        self.meta['visitnum'] = dict(ext=0, card='VISIT')
        self.meta['groupnum'] = dict(card=None, compound=True)
        self.meta['expnum'] = dict(card=None, compound=True)

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
            # The PypeIt output files (e.g., spec2d) record the dither offset
            # of the reduced slit, measured from the center of the slit (see
            # get_slit_dithoff); use it if available.
            if headarr[0].get('DITHOFF') is not None:
                return headarr[0]['DITHOFF']
            # YOFFSET is the offset (arcsec) along the slit in the SI ideal
            # frame, measured from the reference point of the aperture.  A
            # positive YOFFSET moves the trace to *larger* spatial pixels in
            # the PypeIt image, which is opposite to the sign convention
            # adopted by PypeIt for ``dithoff`` (e.g., used by coadd2d when
            # offsets='header'), so we flip the sign here.
            yoffset = headarr[0].get('YOFFSET')
            return 0.0 if yoffset is None else -round(yoffset, 3) + 0.0

        elif meta_key == 'dithpat':
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
        elif meta_key == 'dateobs':
            return Time(headarr[0].get('EXPMID'), format='mjd').isot if headarr[0].get('EXPMID') is not None else None
        elif meta_key == 'groupnum':
            return f"{headarr[0].get('VISITGRP')}{headarr[0].get('SEQ_ID')}{headarr[0].get('ACT_ID')}" if headarr[0].get('VISITGRP') is not None else None
        elif meta_key == 'expnum':
           return f"{int(headarr[0].get('EXPOSURE')):05d}" if headarr[0].get('EXPOSURE') is not None else None
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

    def get_comb_group(self, fitstbl):
        """
        Automatically assign combination groups by grouping files that share the
        same ``rawfile_basename``.

        This method is used in
        :func:`~pypeit.metadata.PypeItMetaData.set_combination_groups`, and
        directly modifies the ``comb_id`` and ``calib`` columns in the provided
        table.

        Specifically, for JWST NIRSpec each physical exposure produces two raw
        files — one per detector (``_nrs1`` / ``_nrs2``).  The
        :func:`rawfile_basename` method strips those suffixes so that both files
        share the same basename.  This method exploits that property to:

        1. Assign the same ``calib`` value to every file that shares a
           ``rawfile_basename`` (i.e., NRS1 and NRS2 of the same exposure are
           placed in the same calibration group).  A sequential integer counter
           (1-based, maximum 64) is used so that
           different exposures get different ``calib`` values.  Frames whose
           ``frametype`` is ``'None'`` are excluded from this assignment and
           their ``calib`` value is left at the default ``'0'``.  If the number
           of distinct basenames with a valid frametype exceeds 64, a warning is
           issued and the counter is capped at 64, meaning that different exposures
           may end up in the same calib group.

        2. Assign the same ``comb_id`` to science frames that share a
           ``rawfile_basename``, and — independently — to standard frames that
           share a ``rawfile_basename``.  This ensures both detector files of
           the same science (or standard) exposure are combined together.

        Args:
            fitstbl (`astropy.table.Table`_):
                The table with the metadata for all the frames.

        Returns:
            `astropy.table.Table`_: modified fitstbl.
        """
        # Compute rawfile_basename for every row (no targname/slitname/mjd so that
        # _nrs1 and _nrs2 files from the same exposure map to the same value).
        basenames = np.array([self.rawfile_basename(f) for f in fitstbl['filename']])
        unique_basenames = np.unique(basenames)

        # Flag rows whose frametype is 'None' (i.e. unclassified frames).
        # Their calib value is left at '0' and they are excluded from the
        # counter-based grouping below.
        has_frametype = np.array([ft != 'None' for ft in fitstbl['frametype']])

        # 1. Files with the same rawfile_basename get the same calib value.
        #    Frames with frametype 'None' are skipped (their calib stays '0').
        #    The calib counter is 1-based and is capped at 64.
        if 'calib' in fitstbl.keys():
            calib_counter = 1
            for basename in unique_basenames:
                same_basename = (basenames == basename) & has_frametype
                if not np.any(same_basename):
                    continue
                fitstbl['calib'][same_basename] = str(calib_counter)
                if calib_counter > 64:
                    log.warning(f'Number of unique calib groups exceeded the maximum of 64. '
                                f'Therefore the calib group counter was capped at 64. '
                                f'Consider splitting your dataset into smaller subsets.')
                    calib_counter = 64
                else:
                    calib_counter += 1

        # 2. Science frames and standard frames are handled separately.
        #    Within each frame type, files with the same rawfile_basename are
        #    assigned the same comb_id (the minimum value already set for that
        #    group, which keeps the IDs globally unique across groups).
        for frame_type in ['science', 'standard']:
            same_ftype = np.array([frame_type in _tab for _tab in fitstbl['frametype']])
            if not np.any(same_ftype):
                continue
            for basename in unique_basenames:
                in_group = same_ftype & (basenames == basename)
                if not np.any(in_group):
                    continue
                # Use the smallest existing comb_id in the group as the
                # representative value for all files sharing this basename.
                rep_comb_id = np.min(fitstbl['comb_id'][in_group])
                fitstbl['comb_id'][in_group] = rep_comb_id


        return fitstbl

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
        return pypeit_keys + ['pid', 'obsnum', 'visitnum', 'groupnum', 'expnum',
            'dateobs', 'dithpat', 'dithpos', 'dithoff']


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
        if ftype in ['pixelflat', 'trace']:
            return good_exp & (fitstbl['idname'] == 'interpolatedflat')
        if ftype in ['arc', 'tilt']:
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
    
    def get_rawimage(self, raw_file, det, keys='SCI'):
        """
        Read raw images and generate a few other bits and pieces
        that are key for image processing.

        This reads either the full-frame JWST images (e.g., the
        ``_assign_wcs.fits`` science frames) or, when ``keys`` selects the
        ``SCI`` extension of a slit (see
        :func:`~pypeit.calibrations.NIRSpecSlitCalibrations._build_calibimage`),
        the cutout of a single slit in the JWST ``_cal.fits`` or
        ``_interpolatedflat.fits`` files.  In the latter case, the slit is
        identified by the ``SLTNAME`` keyword of the HDU header, and the
        nominal platescale of the returned detector(s) (and mosaic) is
        replaced by the along-slit platescale of that slit, computed from
        the slit GWCS (see :func:`get_slit_platescale`).  For full-frame
        images, the slit is not known here, and the nominal platescale is
        returned; the platescale of the slit is then set by
        :func:`~pypeit.images.rawimage.NIRSpecRawImage.process`.

        Parameters
        ----------
        raw_file : :obj:`str` or :obj:`list` of :obj:`str`
            File(s) to read
        det : :obj:`int`, :obj:`tuple`
            1-indexed detector or tuple of detectors (for a mosaic) to read
        keys : :obj:`list` of or single int, :obj:`str` or :obj:`tuple` of (string, int)
            A list of keys or a single key identifying the HDU to read.  If a key is a tuple, it is of the
            form `(name, ver)` where `ver` is an `EXTVER`` value that must
            match the HDU being searched for. If keys is a list, it has to be of the same length as the number of
            raw files.

        Returns
        -------
        detector_par : :class:`~pypeit.images.detector_container.DetectorContainer`, :class:`~pypeit.images.mosaic.Mosaic`
            Detector (or mosaic) metadata parameters.
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

        # check the `keys` values
        if keys is None:
            keys = ['SCI']*nimg
        elif isinstance(keys, (str, int, tuple)):
            keys = [keys]*nimg
        elif isinstance(keys, list):
            if len(keys) != nimg:
                raise PypeItError(f'Number of keys provided ({len(keys)}) does not match the number of detectors ({nimg}).')
        else:
            raise PypeItError(f'Invalid keys specification: {keys}')

        # check that the number of files found matches the number of detectors in the mosaic
        if nimg==2 and len(file_list) != nimg:
            raise PypeItError(f'Expected {nimg} files for mosaic with detectors {det}, but found {len(file_list)} files: {file_list}')

        # Grab the detector or mosaic parameters
        mosaic = None if nimg == 1 else self.get_mosaic_par(det, hdu=None)
        detectors = [self.get_detector_par(det, hdu=None)] if nimg == 1 else mosaic.detectors

        # When reading the cutout of a single slit from the JWST _cal.fits
        # files, use the along-slit platescale of that slit.  For the
        # full-frame science images, the slit is not known here, and the
        # platescale is updated by NIRSpecRawImage.process.
        slit_name = hdul[keys[0]].header.get('SLTNAME')
        if slit_name is not None:
            platescale = self.get_slit_platescale(file_list, slit_name)
            if platescale is not None:
                for d in (detectors if mosaic is None else [*detectors, mosaic]):
                    d.platescale = platescale

        # Read the image(s)
        raw_img = [None]*nimg
        rawdatasec_img = [None]*nimg
        oscansec_img = [None]*nimg
        for i in range(nimg):
            indx = np.where([f'nrs{detectors[i].det}' in n for n in file_list])[0]
            if len(indx) == 0:
                raise PypeItError(f'Could not find file for detector nrs{detectors[i].det} in mosaic with detectors {det}.')
            _hdu = io.fits_open(file_list[indx[0]], ignore_missing_end=True, output_verify='ignore', ignore_blank=True)
            # Raw image
            _img = _hdu[keys[i]].data.astype(float)
            nofinite = np.logical_not(np.isfinite(_img))
            _img[nofinite] = 0.0
            raw_img[i] = _img
            rawdatasec_img[i] = np.zeros_like(raw_img[i], dtype=int) + int(detectors[i].det)
            oscansec_img[i] = np.zeros_like(raw_img[i], dtype=int)

        if nimg == 1:
            return detectors[0], raw_img[0], hdul, exptime, rawdatasec_img[0], oscansec_img[0]

        # For mosaics, stack into a uniform 3-D array when all shapes match,
        # or into a 1-D object array when NRS1 and NRS2 have different shapes.
        shapes = [img.shape for img in raw_img]
        if len(set(shapes)) == 1:
            return mosaic, np.array(raw_img), hdul, exptime, \
                   np.array(rawdatasec_img), np.array(oscansec_img)

        # Different shapes: pack into 1-D object arrays so RawImage.__init__
        # can still determine nimg from rawimage.shape[0].
        raw_img_arr = np.empty(nimg, dtype=object)
        rawdatasec_arr = np.empty(nimg, dtype=object)
        oscansec_arr = np.empty(nimg, dtype=object)
        for i in range(nimg):
            raw_img_arr[i] = raw_img[i]
            rawdatasec_arr[i] = rawdatasec_img[i]
            oscansec_arr[i] = oscansec_img[i]
        return mosaic, raw_img_arr, hdul, exptime, rawdatasec_arr, oscansec_arr

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
        

    def get_nirspec_slits(self, fitstbl, calib_ID, par):
        """
        Determine the list of NIRSpec slits to reduce for a given calibration group.

        Reads the JWST calibration FITS extensions  directly — without any ``jwst.datamodels``
        dependency — to discover the available slits, map each slit to its detector
        (NRS1/NRS2 or mosaic), and optionally filter by ``par['rdx']['maskIDs']``.

        The detector assignment is inferred from the filename: a file whose name
        contains ``'nrs1'`` is classified as NRS1, one containing ``'nrs2'`` as
        NRS2, and the primary–header ``DETECTOR`` keyword is used as a fallback.
        Slit metadata (name, source name, source identifier, source alias) is
        read from the per-extension FITS header keywords ``SLTNAME``,
        ``SRCNAME``, ``SRCID``, and ``SRCALIAS``, respectively.

        Args:
            fitstbl (:class:`~pypeit.metadata.PypeItMetaData`):
                The metadata table for this reduction run.
            calib_ID (:obj:`str`):
                The calibration group identifier.
            par (:class:`~pypeit.par.pypeitpar.PypeItPar`):
                The full parameter set for the reduction.

        Returns:
            :obj:`tuple`: A 2-tuple ``(gd_slits_sources, slit_det_map)``:

            - **gd_slits_sources** (:obj:`list` or ``None``): List of
              ``(slit_name, source_name, source_id, source_alias)`` tuples
              for each slit to be reduced.  ``None`` if no slits are found.
            - **slit_det_map** (:obj:`dict` or ``None``): Dictionary mapping
              each slit name to its detector — ``1`` for NRS1-only slits,
              ``2`` for NRS2-only slits, and ``(1, 2)`` for slits that span
              both detectors (mosaic).  ``None`` if no slits are found.
        """

        # get cal_files (i.e., _cal.fits), which have the `arc` frametype
        cal_files = fitstbl.find_frame_files('arc', calib_ID=calib_ID)

        if not cal_files:
            log.warning(f'No _cal.fits files found for calib_ID {calib_ID}.')
            return None, None

        # Collect per-slit info by reading FITS headers directly.
        slit_names_nrs1 = []
        slit_names_nrs2 = []
        slit_info = {}  # slit_name -> (source_name, source_id, source_alias)

        for cal_file in cal_files:
            with fits.open(cal_file) as hdul:
                # Determine detector from filename first, fall back to primary header.
                if 'nrs1' in str(cal_file).lower():
                    detector = 'NRS1'
                elif 'nrs2' in str(cal_file).lower():
                    detector = 'NRS2'
                else:
                    detector = hdul[0].header.get('DETECTOR', '').upper()

                for hdu in hdul:
                    slt_name = hdu.header.get('SLTNAME')
                    if slt_name is None:
                        continue

                    # Use two independent checks so a slit appearing in both NRS1
                    # and NRS2 files is correctly recorded in both lists (mosaic case).
                    if detector == 'NRS1' and slt_name not in slit_names_nrs1:
                        slit_names_nrs1.append(slt_name)
                    if detector == 'NRS2' and slt_name not in slit_names_nrs2:
                        slit_names_nrs2.append(slt_name)

                    # First occurrence wins for source metadata.
                    if slt_name not in slit_info:
                        src_name = hdu.header.get('SRCNAME', '')
                        src_id = hdu.header.get('SOURCEID', -1)
                        src_alias = hdu.header.get('SRCALIAS', '')
                        slit_info[slt_name] = (src_name, src_id, src_alias)

        slit_names_nrs1 = np.array(slit_names_nrs1) if slit_names_nrs1 else np.array([], dtype=str)
        slit_names_nrs2 = np.array(slit_names_nrs2) if slit_names_nrs2 else np.array([], dtype=str)

        if not len(slit_names_nrs1) and not len(slit_names_nrs2):
            log.warning(f'No slits with SLTNAME header keyword found in calibration files for calib_ID '
                        f'{calib_ID}. Skipping.')
            return None, None

        all_slit_names = np.unique(np.hstack([slit_names_nrs1, slit_names_nrs2]))
        slit_sources_uni = [(slt, *slit_info[slt]) for slt in all_slit_names]

        # Filter by maskIDs if specified.
        if par['rdx']['maskIDs'] is not None:
            maskIDs = [str(mid).strip() for mid in par['rdx']['maskIDs']]
            gd_slits_sources = [
                (slt, src_n, src_id, src_alias)
                for slt, src_n, src_id, src_alias in slit_sources_uni
                if str(slt).strip() in maskIDs
                or str(src_id) in maskIDs
                or str(src_alias) in maskIDs
            ]
            if not gd_slits_sources:
                log.warning(f'No slits or sources matched '
                            f'maskIDs={par["rdx"]["maskIDs"]}. Skipping.')
                return None, None
            if len(gd_slits_sources) < len(maskIDs):
                missing = [
                    mid for mid in maskIDs
                    if mid not in [s for s, _, _, _ in gd_slits_sources]
                    and mid not in [str(si) for _, _, si, _ in gd_slits_sources]
                    and mid not in [sa for _, _, _, sa in gd_slits_sources]
                ]
                log.warning(f'The following maskIDs were not found: '
                            f'{", ".join(missing)}. Processing will continue '
                            'on the available slits and sources.')
        else:
            gd_slits_sources = slit_sources_uni

        slit_strs = ', '.join([f'({s}, {n})' for s, n, _, _ in gd_slits_sources])
        log.info(f'Found {len(gd_slits_sources)} slit(s) to process '
                 f'(slit_name, src_name): {slit_strs}')

        # Build a per-slit detector map.
        slit_det_map = {}
        for slit_name, _, _, _ in gd_slits_sources:
            in_nrs1 = slit_name in slit_names_nrs1
            in_nrs2 = slit_name in slit_names_nrs2
            if in_nrs1 and not in_nrs2:
                slit_det_map[slit_name] = 1
            elif in_nrs2 and not in_nrs1:
                slit_det_map[slit_name] = 2
            else:
                slit_det_map[slit_name] = (1, 2)  # Slit spans both detectors (mosaic).

        return gd_slits_sources, slit_det_map

    def get_slit_platescale(self, cal_files, slit_name):
        """
        Compute the along-slit (spatial) platescale of a single NIRSpec slit
        from the GWCS stored in the JWST calibrated (``_cal.fits``) files.

        The NIRSpec platescale is not uniform: it ranges from ~0.103 to ~0.108
        arcsec/pixel across the MSA (depending on quadrant and shutter row)
        and is ~0.1053 arcsec/pixel for the fixed slits.  This method uses
        :func:`slit_platescale_from_gwcs` to measure it for the requested slit.
        When the slit falls on both detectors, the mean of the values
        measured on each detector is returned.

        The platescale is used to replace the nominal one for the
        calibration images of the slit (see :func:`get_rawimage`) and for the
        science frames (see
        :func:`~pypeit.calibrations.NIRSpecSlitCalibrations.get_slit_slices`
        and :func:`~pypeit.images.rawimage.NIRSpecRawImage.process`).

        .. warning::

            Reading the GWCS requires the ``stdatamodels`` package, which
            provides the ``asdf`` converters for the JWST-specific transforms.
            If ``stdatamodels`` is not installed (or the GWCS cannot be
            evaluated), this method returns None and the nominal
            platescale from :func:`get_detector_par` is used instead, which
            can be off by up to ~3% for MSA slits.

        Args:
            cal_files (:obj:`str`, `Path`_, :obj:`list`):
                One or more JWST ``_cal.fits`` files (e.g., one per detector
                for a mosaic) that include the slit.
            slit_name (:obj:`str`):
                Name of the slit (the ``SLTNAME`` header keyword).

        Returns:
            :obj:`float`: The along-slit platescale in arcsec/pixel, or None
            if it could not be computed.
        """
        nominal = self.get_detector_par(1).platescale
        if not stdatamodels_installed():
            if str(slit_name) in _STDATAMODELS_WARNED:
                return None
            if len(_STDATAMODELS_WARNED) == 0:
                log.warning('The stdatamodels package is NOT installed!  PypeIt requires it to '
                            'read the slit GWCS from the JWST _cal.fits files and compute the '
                            'along-slit platescale of each NIRSpec slit.  Falling back to the '
                            f'nominal platescale ({nominal} arcsec/pixel), which can be off by '
                            'up to ~3% for MSA slits and propagates to, e.g., the offsets '
                            'computed by pypeit_coadd_2dspec from the header dither offsets.  '
                            'We strongly recommend installing it with: pip install stdatamodels')
            else:
                log.warning(f'stdatamodels is not installed; using the nominal platescale '
                            f'({nominal} arcsec/pixel) for slit {slit_name}.')
            _STDATAMODELS_WARNED.add(str(slit_name))
            return None

        _cal_files = [cal_files] if isinstance(cal_files, (str, Path)) else list(cal_files)
        platescales = []
        for cal_file in _cal_files:
            try:
                platescales.append(slit_platescale_from_gwcs(cal_file, slit_name))
            except Exception as e:
                log.warning(f'Could not compute the platescale of slit {slit_name} from the GWCS '
                            f'in {Path(cal_file).name}: {type(e).__name__}: {e}')
        if len(platescales) == 0:
            log.warning(f'Using the nominal platescale ({nominal} arcsec/pixel) for slit '
                        f'{slit_name}.')
            return None

        return float(np.mean(platescales))

    def get_slit_dithoff(self, fitstbl, frame, slit_name):
        """
        Return the dither offset of a frame for a specific slit, measured
        from the center of the slit.

        The ``dithoff`` metadata of the NIRSpec frames is set by
        :func:`compound_meta` to the negative of the ``YOFFSET`` header
        keyword, which is measured from the reference point of the aperture.
        For MSA observations, this reference point is not the center of the
        slit, and it changes with the pointing, such that the offsets of
        frames taken in different visits or with different MSA
        configurations cannot be compared.  Instead, this method uses the
        position of the source in the slit (see
        :func:`slit_source_offset_from_gwcs`), which is measured from the
        center of the slit in every frame, and returns its negative, to
        follow the ``PypeIt`` sign convention for ``dithoff``.

        The position is measured using the JWST ``_cal.fits`` file(s) of the
        same exposure as the frame (i.e., with the same
        :func:`rawfile_basename`).  When the slit falls on both detectors,
        the mean of the values measured on each detector is used.

        .. warning::

            This requires the ``stdatamodels`` package to be installed.  If
            the offset cannot be computed, None is returned and the
            ``dithoff`` value in the metadata table (i.e., ``-YOFFSET``) is
            used, which is only reliable for frames taken in the same visit.

        Args:
            fitstbl (:class:`~pypeit.metadata.PypeItMetaData`):
                The metadata table for the frames being reduced.
            frame (:obj:`int`):
                0-indexed row in the metadata table with the frame.
            slit_name (:obj:`str`):
                Name of the slit (the ``SLTNAME`` header keyword).

        Returns:
            :obj:`float`: The dither offset in arcsec, or None if it could
            not be computed.
        """
        fallback = ('Using the dither offset from the YOFFSET header keyword, which is only '
                    'reliable for frames taken in the same visit.')
        if not stdatamodels_installed():
            log.warning(f'stdatamodels is not installed; cannot compute the position of the '
                        f'source in slit {slit_name}.  {fallback}')
            return None

        # The calibrated files of the same exposure
        basename = self.rawfile_basename(fitstbl.frame_paths(frame))
        cal_files = [f for f in fitstbl.find_frame_files('arc')
                     if self.rawfile_basename(f) == basename]

        offsets = []
        for cal_file in cal_files:
            with fits.open(cal_file) as hdul:
                slit_hdr = [h.header for h in hdul if h.name == 'SCI'
                            and h.header.get('SLTNAME', '').strip() == str(slit_name).strip()]
                if len(slit_hdr) == 0:
                    # Slit not on this detector
                    continue
                if not slit_has_source(hdul[0].header, slit_hdr[0]):
                    log.info(f'Slit {slit_name} does not include a planned source.  {fallback}')
                    return None
            try:
                offset = slit_source_offset_from_gwcs(cal_file, slit_name)
            except Exception as e:
                log.warning(f'Could not compute the position of the source in slit {slit_name} '
                            f'from the GWCS in {Path(cal_file).name}: {type(e).__name__}: {e}')
                continue
            if offset is not None:
                offsets.append(offset)

        if len(offsets) == 0:
            log.warning(f'Could not determine the position of the source in slit {slit_name} '
                        f'for {Path(fitstbl.frame_paths(frame)).name}.  {fallback}')
            return None

        dithoff = -round(float(np.mean(offsets)), 4) + 0.0
        log.info(f'Dither offset of slit {slit_name} (from the slit center): {dithoff:.4f} arcsec')
        return dithoff


_STDATAMODELS_WARNED = set()
"""
Slits for which the warning about the missing ``stdatamodels`` package has
been issued.  The full warning is issued only once per process, and a short
one once per slit.
"""


def slit_has_source(primary_hdr, slit_hdr):
    """
    Check if a NIRSpec slit includes the planned source, such that the source
    position recorded by the JWST pipeline (``SRCYPOS``) follows the dither
    pattern.

    For the slits that do not include a planned source (i.e., the fixed slits
    other than the one targeted by the observation, and the MSA background
    slitlets), the JWST pipeline sets ``SRCYPOS`` to 0 in all the exposures,
    regardless of the dither pattern.

    Args:
        primary_hdr (`astropy.io.fits.Header`_):
            Primary header of the JWST ``_cal.fits`` file.
        slit_hdr (`astropy.io.fits.Header`_):
            Header of the ``SCI`` extension of the slit.

    Returns:
        :obj:`bool`: True if the slit includes the planned source.
    """
    if primary_hdr.get('EXP_TYPE') == 'NRS_FIXEDSLIT':
        return str(slit_hdr.get('SLTNAME', '')).strip() \
                    == str(primary_hdr.get('FXD_SLIT', '')).strip()
    alias = str(slit_hdr.get('SRCALIAS') or '').strip().upper()
    ra, dec = slit_hdr.get('SRCRA'), slit_hdr.get('SRCDEC')
    return not alias.startswith('BKG') and ra not in (None, 0) and dec not in (None, 0)


def stdatamodels_installed():
    """
    Check if the ``stdatamodels`` package is installed.

    ``stdatamodels`` is needed to read the GWCS of NIRSpec slits (see
    :func:`slit_platescale_from_gwcs`).  It is not imported here: once
    installed, it registers its ``asdf`` converters for the JWST-specific
    transforms via ``asdf`` entry points.

    Returns:
        :obj:`bool`: True if ``stdatamodels`` is installed.
    """
    return importlib.util.find_spec('stdatamodels') is not None


def slit_platescale_from_gwcs(cal_file, slit_name, nsamp=15):
    """
    Measure the along-slit (spatial) platescale of a NIRSpec slit from the
    GWCS stored in a JWST calibrated (``_cal.fits``) file.

    The GWCS of each slit is read directly from the ``ASDF`` extension of the
    file using ``asdf``, without opening the file as a JWST datamodel.  The
    WCS is evaluated at ``nsamp`` positions along the spectral direction at
    the center of the slit and one pixel above it (along the slit, i.e.,
    detector +y in the JWST frame), and the median sky separation between
    the two sets of positions is returned.

    .. warning::

        This requires the ``stdatamodels`` package to be installed.  It
        provides the ``asdf`` converters for the JWST-specific transforms used
        by the NIRSpec GWCS (``tag:stsci.edu:jwst_pipeline/*`` tags) and
        registers them automatically with ``asdf``, but it does not need to
        be imported.

    Args:
        cal_file (:obj:`str`, `Path`_):
            JWST ``_cal.fits`` file.
        slit_name (:obj:`str`):
            Name of the slit (the ``SLTNAME`` header keyword).
        nsamp (:obj:`int`, optional):
            Number of positions sampled along the spectral direction.

    Returns:
        :obj:`float`: The along-slit platescale in arcsec/pixel.
    """
    return _slit_gwcs_measurements(str(Path(cal_file).resolve()), str(slit_name).strip(),
                                   int(nsamp))[0]


def slit_source_offset_from_gwcs(cal_file, slit_name):
    """
    Measure the position of the source along a NIRSpec slit, relative to the
    center of the slit, from the GWCS stored in a JWST calibrated
    (``_cal.fits``) file.

    The JWST pipeline records the position of the source in the slit
    (``source_ypos``, also in the ``SRCYPOS`` header keyword of the slit) in
    the ``slit_frame`` coordinate frame of the slit GWCS: in units of the slit
    length relative to the slit center for the fixed slits, and in units of
    the shutter pitch relative to the center of the shutter that contains the
    source for the MSA.  The detector -> ``slit_frame`` transform of the GWCS
    is used to find the detector rows of the source and of the center of the
    slit (the midpoint of ``slit_ymin`` and ``slit_ymax``) at the center of
    the slit in the spectral direction.  Their difference, converted to
    arcsec with the along-slit platescale (see
    :func:`slit_platescale_from_gwcs`), is returned.

    Contrary to the ``YOFFSET`` header keyword, which is measured from the
    reference point of the aperture (e.g., of the whole MSA for MSA
    observations) and therefore depends on the pointing, this offset is
    measured from the center of the slit in every exposure.  The offsets of
    different exposures can then be compared even if they were taken in
    different visits, with different pointings or MSA configurations.

    .. warning::

        This requires the ``stdatamodels`` package to be installed; see
        :func:`slit_platescale_from_gwcs`.

    Args:
        cal_file (:obj:`str`, `Path`_):
            JWST ``_cal.fits`` file.
        slit_name (:obj:`str`):
            Name of the slit (the ``SLTNAME`` header keyword).

    Returns:
        :obj:`float`: The offset in arcsec of the source from the center of
        the slit, along the slit, which is positive when the source is at
        larger spatial pixels in the ``PypeIt`` frame.  None if the source
        position is not defined.
    """
    platescale, offset_pix = _slit_gwcs_measurements(str(Path(cal_file).resolve()),
                                                     str(slit_name).strip(), 15)
    return None if offset_pix is None else offset_pix * platescale


def _slit_rows(wcs, slit_y, xpix):
    """
    Find the detector rows where the ``slit_frame`` y coordinate of a NIRSpec
    slit GWCS has the provided values.

    Args:
        wcs (`gwcs.wcs.WCS`):
            Slit GWCS.
        slit_y (array-like):
            Values of the ``slit_frame`` y coordinate.
        xpix (:obj:`float`):
            Detector column (spectral direction) at which to evaluate the
            transform.

    Returns:
        `numpy.ndarray`_: Detector rows, NaN for the values outside the slit.
    """
    (_, _), (y0, y1) = wcs.bounding_box.bounding_box(order='F')
    y = np.arange(y0, y1, 0.01)
    sy = wcs.get_transform('detector', 'slit_frame')(np.full_like(y, xpix), y)[1]
    indx = np.isfinite(sy)
    y, sy = y[indx], sy[indx]
    if y.size < 2:
        return np.full(np.size(slit_y), np.nan)
    srt = np.argsort(sy)
    rows = np.interp(slit_y, sy[srt], y[srt])
    outside = (np.asarray(slit_y) < sy.min()) | (np.asarray(slit_y) > sy.max())
    return np.where(outside, np.nan, rows)


@functools.lru_cache(maxsize=1024)
def _slit_gwcs_measurements(cal_file, _slit_name, nsamp):
    """
    Cached worker for :func:`slit_platescale_from_gwcs` and
    :func:`slit_source_offset_from_gwcs`.

    The GWCS of the slit is read only once to measure both quantities.  The
    result is cached because the same slit is processed for multiple
    calibration steps and frames during a reduction.

    Args:
        cal_file (:obj:`str`):
            Full path to the JWST ``_cal.fits`` file.
        _slit_name (:obj:`str`):
            Name of the slit, with no leading or trailing spaces.
        nsamp (:obj:`int`):
            Number of positions sampled along the spectral direction to
            measure the platescale.

    Returns:
        :obj:`tuple`: The along-slit platescale in arcsec/pixel and the
        offset of the source from the center of the slit in pixels (None if
        the source position is not defined).
    """
    # asdf is a dependency of stdatamodels
    import asdf

    with fits.open(cal_file) as hdul:
        # The slits in the ASDF tree are in the same order as the SCI extensions
        names = [h.header.get('SLTNAME', '').strip() for h in hdul if h.name == 'SCI']
        if _slit_name not in names:
            raise PypeItError(f'Slit {_slit_name} not found in {cal_file}.')
        try:
            buf = BytesIO(hdul['ASDF'].data.tobytes())
        except KeyError:
            raise PypeItError(f'No ASDF extension found in {cal_file}.')
    islit = names.index(_slit_name)

    # NOTE: Lazy loading is required because the data arrays referenced by the
    # ASDF tree are stored in the other FITS extensions, not in the ASDF one.
    with asdf.open(buf, lazy_load=True, lazy_tree=True) as af:
        slits = af.tree['slits']
        if str(slits[islit]['name']).strip() != _slit_name:
            # Unexpected order; fall back to (slower) searching by name
            islit = [str(s['name']).strip() for s in slits].index(_slit_name)
        slit = slits[islit]
        wcs = slit['meta']['wcs']

        # Platescale: sample the center of the slit, avoiding the edges of the
        # bounding box
        (x0, x1), (y0, y1) = wcs.bounding_box.bounding_box(order='F')
        x = np.linspace(x0, x1, nsamp + 2)[1:-1]
        y = np.full_like(x, 0.5 * (y0 + y1))
        ra0, dec0 = np.radians(wcs(x, y)[:2])
        ra1, dec1 = np.radians(wcs(x, y + 1)[:2])

        # Rows of the source and of the center of the slit
        src_ypos = slit.get('source_ypos')
        slit_ymin, slit_ymax = slit.get('slit_ymin'), slit.get('slit_ymax')
        offset_pix = None
        if src_ypos is not None and slit_ymin is not None and slit_ymax is not None:
            src_row, cen_row = _slit_rows(wcs, [src_ypos, 0.5 * (slit_ymin + slit_ymax)],
                                          0.5 * (x0 + x1))
            if np.isfinite(src_row) and np.isfinite(cen_row):
                offset_pix = float(src_row - cen_row)

    # Haversine angular separation
    hav = np.sin((dec1 - dec0) / 2)**2 + np.cos(dec0) * np.cos(dec1) * np.sin((ra1 - ra0) / 2)**2
    sep = np.degrees(2 * np.arcsin(np.sqrt(hav))) * 3600
    if not np.any(np.isfinite(sep)):
        raise PypeItError(f'GWCS of slit {_slit_name} in {cal_file} could not be evaluated.')
    platescale = float(np.nanmedian(sep))
    log.info(f'Platescale of slit {_slit_name} from the GWCS in {Path(cal_file).name}: '
             f'{platescale:.5f} arcsec/pixel')
    return platescale, offset_pix
