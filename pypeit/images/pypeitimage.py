""" Module for the PypeItImage include its Mask

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import numpy as np
import os
import inspect

from pypeit import msgs
from pypeit.images import imagebitmask
from pypeit.images.detector_container import DetectorContainer
from pypeit.images.mosaic import Mosaic
from pypeit.core import procimg
from pypeit.display import display
from pypeit import datamodel
from pypeit import utils
from pypeit import masterframe

from IPython import embed


class PypeItImage(datamodel.DataContainer):
    r"""
    Container class for processed ``PypeIt`` images and associated data.

    Image data can be either from a single detector, multiple detectors arranged
    in a 3D array, or multiple detectors resampled to a single mosaic image.
    The orientation of the images is always spectral pixels along the first axis
    and spatial pixels along the second axis; i.e., the shape is :math:`(N_{\rm
    spec},N_{\rm spat})`.

    Note that all arguments for instantiation are optional, meaning that an
    empty object can be created.

    Args:
        image (`numpy.ndarray`_, optional):
            Primary image data
        ivar (`numpy.ndarray`_, optional):
            Inverse variance image
        nimg (`numpy.ndarray`_, optional):
            If a combination of multiple images, this is the number of images
            that contributed to each pixel
        det_img (`numpy.ndarray`_, optional):
            If a detector mosaic, this image provides the detector that
            contributed to each pixel
        rn2img (`numpy.ndarray`_, optional):
            Read-noise-squared image
        base_var (`numpy.ndarray`_, optional):
            Base-level image variance, excluding count shot-noise
        img_scale (`numpy.ndarray`_, optional):
            Image count scaling applied (e.g., 1/flat-field)
        bpm (`numpy.ndarray`_, optional):
            Bad pixel mask
        crmask (`numpy.ndarray`_, optional):
            CR mask image
        fullmask (`numpy.ndarray`_, optional):
            Full image bitmask
        detector (:class:`~pypeit.images.detector_container.DetectorContainer`, :class:`~pypeit.images.mosaic.Mosaic`, optional):
            The detector or mosaic parameters
        PYP_SPEC (:obj:`str`, optional):
            PypeIt spectrograph name
        units (:obj:`str`, optional):
            (Unscaled) Pixel units (e- or ADU)
        exptime (:obj:`int`, :obj:`float`, optional):
            Effective exposure time (s)
        noise_floor (:obj:`float`, optional):
            Noise floor included in variance
        shot_noise (:obj:`bool`, optional):
            Shot-noise included in variance
        spat_flexure (:obj:`float`, optional):
            Shift, in spatial pixels, between this image and SlitTrace
        imgbitm (:obj:`str`, optional):
            List of BITMASK keys from :class:`~pypeit.images.imagebitmask.ImageBitMask`

    Attributes:
        files (:obj:`list`):
            List of source files
        rawheadlst (:obj:`list`):
            List containing headers of the raw image file.  If the source data
            is a stack of multiple files, this is the set of headers for the
            *last* file in the list.
        process_steps (:obj:`list`):
            List of steps executed during processing of the image data; see
            :class:`~pypeit.images.rawimage.RawImage.process`.
        master_key (:obj:`str`):
            Master key, only for Master frames
        master_dir (:obj:`str`):
            Master directory, only for Master frames
    """
#        head0 (`astropy.io.fits.Header`_):
#            Primary header of the source fits file.  If the image is a
#            combination of multiple files, this is None.

    version = '1.2.0'
    """Datamodel version number"""

    datamodel = {'image': dict(otype=np.ndarray, atype=np.floating, descr='Primary image data'),
                 'ivar': dict(otype=np.ndarray, atype=np.floating,
                              descr='Inverse variance image'),
                 'nimg': dict(otype=np.ndarray, atype=np.integer,
                              descr='If a combination of multiple images, this is the number of '
                                    'images that contributed to each pixel'),
                 'amp_img': dict(otype=np.ndarray, atype=np.integer,
                                 descr='Provides the amplifier that contributed to each pixel.  '
                                       'If this is a detector mosaic, this must be used in '
                                       'combination with ``det_img`` to select pixels for a '
                                       'given detector amplifier.'),
                 'det_img': dict(otype=np.ndarray, atype=np.integer,
                                 descr='If a detector mosaic, this image provides the detector '
                                       'that contributed to each pixel.'),
                 'rn2img': dict(otype=np.ndarray, atype=np.floating,
                                descr='Read noise squared image'),
                 'base_var': dict(otype=np.ndarray, atype=np.floating,
                                  descr='Base-level image variance, excluding count shot-noise'),
                 'img_scale': dict(otype=np.ndarray, atype=np.floating,
                                   descr='Image count scaling applied (e.g., 1/flat-field)'),
                 'bpm': dict(otype=np.ndarray, atype=np.integer, descr='Bad pixel mask'),
                 'crmask': dict(otype=np.ndarray, atype=np.bool_, descr='CR mask image'),
                 'fullmask': dict(otype=np.ndarray, atype=np.integer, descr='Full image bitmask'),
                 'detector': dict(otype=(DetectorContainer, Mosaic),
                                  descr='The detector or mosaic parameters'),
                 'PYP_SPEC': dict(otype=str, descr='PypeIt spectrograph name'),
                 'units': dict(otype=str, descr='(Unscaled) Pixel units (e- or ADU)'),
                 # TODO: Consider forcing exptime to be a float.
                 'exptime': dict(otype=(int, float), descr='Effective exposure time (s)'),
                 'noise_floor': dict(otype=float, descr='Noise floor included in variance'),
                 'shot_noise': dict(otype=bool, descr='Shot-noise included in variance'),
                 'spat_flexure': dict(otype=float,
                                      descr='Shift, in spatial pixels, between this image '
                                            'and SlitTrace'),
                 'imgbitm': dict(otype=str,
                                 descr='List of BITMASK keys from '
                                       ':class:`~pypeit.images.imagebitmask.ImageBitMask`')}
    """Data model components."""

    bitmask = imagebitmask.ImageBitMask()
    """Class mask attribute"""

    @classmethod
    def from_pypeitimage(cls, pypeitImage):
        """
        Generate an instance
        This enables building the Child from the Parent, e.g. a MasterFrame Image

        This is *not* a deepcopy

        Args:
            pypeitImage (:class:`PypeItImage`):

        Returns:
            pypeitImage (:class:`PypeItImage`):

        """
        _d = {}
        for key in pypeitImage.datamodel.keys():
            if pypeitImage[key] is None:
                continue
            _d[key] = pypeitImage[key]
        # Instantiate
        slf = cls(**_d)
        # Internals
        slf.process_steps = pypeitImage.process_steps
        slf.files = pypeitImage.files
        slf.rawheadlist = pypeitImage.rawheadlist
        slf.master_dir = pypeitImage.master_dir
        slf.master_key = pypeitImage.master_key
        # Return
        return slf

    def __init__(self, image=None, ivar=None, nimg=None, amp_img=None, det_img=None, rn2img=None,
                 base_var=None, img_scale=None, bpm=None, crmask=None, fullmask=None,
                 detector=None, spat_flexure=None, PYP_SPEC=None, units=None, exptime=None,
                 noise_floor=None, shot_noise=None, imgbitm=None):

        # Setup the DataContainer. Dictionary elements include
        # everything but self in the instantiation call.
        # NOTE: Instantiating this way means that all of the arguments must also
        # be members of the datamodel.
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        super(PypeItImage, self).__init__(d={k: values[k] for k in args[1:]})

    def _init_internals(self):
        """
        Initialize attributes that are not part of the datamodel.
        """
#        self.head0 = None
        self.process_steps = None
        self.files = None
        self.rawheadlist = None
        # Master stuff
        self.master_key = None
        self.master_dir = None

    def _validate(self):
        """
        Validate the object.

        Validation steps are:

            #. Check that the bitmasks used to create the image are the same as
               those defined by this version of
               :class:`~pypeit.images.imagebitmask.ImageBitMask`.

            #. Check that all the image data arrays (image, ivar, etc) have the
               same shape.

            #. Initialize the bit mask and incorporate any BPM (:attr:`bpm`) or
               cosmic-ray flags (:attr:`crmask`).

            #. Set the units to counts if none are provided.
        """
        if self.imgbitm is None:
            self.imgbitm = ','.join(list(self.bitmask.keys()))
        elif self.imgbitm != ','.join(list(self.bitmask.keys())):
            msgs.error("Input BITMASK keys differ from current data model!")

        if self.image is not None:
            for attr in ['ivar', 'nimg', 'amp_img', 'det_img', 'rn2img', 'base_var', 'img_scale',
                         'bpm', 'crmask', 'fullmask']:
                _arr = getattr(self, attr)
                if _arr is not None and _arr.shape != self.shape:
                    msgs.error(f'Attribute {attr} does not match image shape.')

        # Make sure the fullmask, bpm, and crmask are all consistent
        # TODO: These steps ignore any existing flags in fullmask and resets
        # them to current arrays.  The masking needs to be handled better so
        # that consistency is straight-forward.
        if self.fullmask is None:
            self.reinit_mask()
        if self.bpm is not None:
            self.update_mask('BPM', action='turn_off')
            self.update_mask('BPM', indx=self.bpm.astype(bool))
        if self.crmask is not None:
            self.update_mask('CR', action='turn_off')
            self.update_mask('CR', indx=self.crmask.astype(bool))

        # Make sure the units are defined
        if self.units is None:
            self.units = 'e-'

    def _bundle(self):
        """
        Package the datamodel for writing.

        Returns:
            :obj:`list`: A list of dictionaries, each list element is written to
            its own fits extension. See
            :class:`~pypeit.datamodel.DataContainer`.
        """
        # TODO: Add `files` and `process_steps`?
        d = []
        # Primary image
        d.append(dict(image=self.image))

        # Rest of the datamodel
        for key in self.keys():
            # TODO: Why are these skipped?
            if key in ['image', 'crmask', 'bpm']:
                continue
            # Skip None
            if self[key] is None:
                continue
            # Array?
            if self.datamodel[key]['otype'] == np.ndarray:
                tmp = {}
                tmp[key] = self[key]
                d.append(tmp)
            elif key == 'detector':
                d.append(dict(detector=self.detector))
            else: # Add to header of the primary image
                d[0][key] = self[key]
        # Return
        return d

    @property
    def shape(self):
        """
        Return the primary image shape.
        """
        return () if self.image is None else self.image.shape

    @property
    def is_multidetector(self):
        """
        Flag if the "image" has more than one detector image that *has not* been
        resampled into a mosaic image.
        """
        return isinstance(self.detector, Mosaic) and self.image.ndim == 3

    @property
    def is_mosaic(self):
        """
        Flag if the image is a mosaic of multiple detectors.
        """
        return isinstance(self.detector, Mosaic) and self.image.ndim == 2

    def build_crmask(self, par, subtract_img=None):
        """
        Identify and flag cosmic rays in the image.

        This is mainly a wrapper to :func:`pypeit.core.procimg.lacosmic`, with
        the boolean cosmic-ray mask is saved to :attr:`crmask`. 

        Args:
            par (:class:`~pypeit.par.pypeitpar.ProcessImagesPar`):
                Parameters that dictate the processing of the images.  See
                :class:`~pypeit.par.pypeitpar.ProcessImagesPar` for the
                defaults.
            subtract_img (`numpy.ndarray`_, optional):
                An image to subtract from the primary image *before* executing
                the cosmic-ray detection algorithm.  If provided, it must have
                the correct shape (see :func:`shape`).

        Returns:
            `numpy.ndarray`_: Copy of :attr:`crmask`.
        """
        if subtract_img is not None and subtract_img.shape != self.shape:
            msgs.error('In cosmic-ray detection, image to subtract has incorrect shape.')

        # Image to flag
        use_img = self.image if subtract_img is None else self.image - subtract_img
        # TODO: Is there error in `subtract_img`?  Should this be added to the
        # variance used by LACosmic?
        var = utils.inverse(self.ivar)

        # TODO: What flags should and should not be included in the "bpm" passed
        # to L.A.Cosmic?  For now, I'm doing the simple thing of just using the
        # bad pixel mask, but we should be constructing this from some or all of
        # the bits flagged in `fullmask`.
        _bpm = self.bpm
        saturation = self.map_detector_value('saturation')
        nonlinear = self.map_detector_value('nonlinear')

        # If the object has multiple images, need to flag each image individually
        if self.is_multidetector:
            self.crmask = [None] * self.shape[0]
            for i in range(self.shape[0]):
                self.crmask[i] = procimg.lacosmic(use_img[i], saturation=saturation[i],
                                                  nonlinear=nonlinear[i], bpm=_bpm[i],
                                                  varframe=var, maxiter=par['lamaxiter'],
                                                  grow=par['grow'],
                                                  remove_compact_obj=par['rmcompact'],
                                                  sigclip=par['sigclip'], sigfrac=par['sigfrac'],
                                                  objlim=par['objlim'])
            self.crmask = np.array(self.crmask)
            return self.crmask.copy()

        # Run LA Cosmic to get the cosmic ray mask and return a copy of the
        # result
        self.crmask = procimg.lacosmic(use_img, saturation=saturation, nonlinear=nonlinear,
                                       bpm=_bpm, varframe=var, maxiter=par['lamaxiter'],
                                       grow=par['grow'], remove_compact_obj=par['rmcompact'],
                                       sigclip=par['sigclip'], sigfrac=par['sigfrac'],
                                       objlim=par['objlim'])
        return self.crmask.copy()

    def map_detector_value(self, attr):
        """
        Provided a detector specific value, remap it as necessary for numpy
        operations with :attr:`image`.

        Args:
            attr (:obj:`str`):
                Attribute of
                :class:`~pypeit.images.detector_container.DetectorContainer` to
                map to an image, vector, or scalar.

        Returns:
            scalar, `numpy.ndarray`_: An array that can be mapped for use with
            :attr:`image` containing the detector values.
        """
        if isinstance(self.detector, DetectorContainer):
            data = self.detector[attr]
            if np.isscalar(data):
                return data

            # Must be defining the per-amplifier value
            if self.amp_img is None:
                msgs.error(f'To remap detector {attr}, object must have amp_img defined.')
            out = np.zeros(self.shape, dtype=type(data[0]))
            for j in range(len(data)):
                out[self.amp_img == j+1] = data[j]
            return out

        # The object holds data from multiple detectors.  Collate the detector
        # data.
        data = [d[attr] for d in self.detector.detectors]

        # Object holds data for multiple detectors in separate images
        if self.is_multidetector:
            if np.isscalar(data[0]):
                # Single value per detector, so just repeat it.
                return np.repeat(data, np.prod(self.shape[1:])).reshape(self.shape)
            # Must be defining the per-amplifier value
            if self.amp_img is None:
                msgs.error(f'To remap detector {attr}, object must have amp_img defined.')
            out = np.zeros(self.shape, dtype=type(data[0][0]))
            for i in range(self.detector.ndet):
                for j in range(len(data[i])):
                    out[i,self.amp_img[i] == j+1] = data[i][j]
            return out

        # Object holds data for multiple detectors that have been mosaiced into
        # a single image
        if self.is_mosaic:
            # Check for amplifier dependent output before entering loop
            if not np.isscalar(data[0]) and self.amp_img is None:
                # Must be defining the per-amplifier value
                msgs.error(f'To remap detector {attr}, object must have amp_img defined.')
            # Get the output type
            otype = type(data[0]) if np.isscalar(data[0]) else type(data[0][0])
            # Fill the array
            out = np.zeros(self.shape, dtype=otype)
            for i in range(self.detector.ndet):
                indx = self.det_img == self.detector.detectors[i].det
                if np.isscalar(data[i]):
                    out[indx] = data[i]
                    continue
                # Must be defining the per-amplifier value
                for j in range(len(data[i])):
                    out[indx & (self.amp_img == j+1)] = data[i][j]
            return out

        # Should not get here
        msgs.error('CODING ERROR: Bad logic in map_detector_value.')

    def build_mask(self, saturation=None, mincounts=None, slitmask=None, from_scratch=True):
        """
        Construct the bit value mask used during extraction.

        .. warning::

            By default, this **erases** any existing mask and starts from
            scratch!  See ``from_scratch``.

        The mask bit keys are defined by
        :class:`~pypeit.images.imagebitmask.ImageBitMask`.  Assuming an instance
        of :class:`~pypeit.images.pypeitimage.PypeItImage` called ``img``, any
        pixel with ``img.fullmask == 0`` is valid, otherwise the pixel has been
        masked.  To determine why a given pixel has been masked (see
        also :func:`~pypeit.images.pypeitimage.PypeItImage.select_flag`):

        .. code-block:: python

            reasons = img.bitmask.flagged_bits(img.fullmask[0,0])

        To get all the pixel masked for a specific set of reasons, e.g.:

        .. code-block:: python

            has_cr = img.select_flag(flag='CR')
            is_saturated = img.select_flag(flag='SATURATION')

        Args:
            saturation (:obj:`float`, :obj:`str`, `numpy.ndarray`_, optional):
                Saturation limit (i.e., maximum allowed value).  Any value above
                this is flagged with the SATURATION bit.  Units must match the
                image.  Can be a single float to use for all pixels or a
                `numpy.ndarray`_ for a pixel-dependent value (e.g., for an image
                mosaic); if the latter, the shape must match :attr:`shape`.  Can
                also be a string, but, if so, the string must be ``'default'``,
                and that means the saturation limit is set by :attr:`detector`.
                The tabulated saturation limit is assumed to be in ADU/DN, which
                will be converted to counts if :attr:`units` is `'e-'`, and the
                value used for the saturation includes the detector
                non-linearity threshold; see
                :func:`~pypeit.images.detector_container.DetectorContainer.nonlinear_counts`.
                If None, pixels are not flagged for saturation.
            mincounts (:obj:`float`, :obj:`str`, `numpy.ndarray`_, optional):
                The minimum valid value; units must match the image.  Any value
                below this is flagged with the MINCOUNTS bit.  Can be a single
                float to use for all pixels or a `numpy.ndarray`_ for a
                pixel-dependent value (e.g., for an image mosaic); if the
                latter, the shape must match :attr:`shape`.  Can also be a
                string, but, if so, the string must be ``'default'``, and that
                means the minimum count threshold is set by :attr:`detector`.
                The tabulated minimum counts is assumed to be in e-, which will
                be converted to ADU/DN if :attr:`units` is `'ADU'`.  If None,
                pixels are not flagged for a minimum value.
            slitmask (`numpy.ndarray`_, optional):
                Slit mask image.  Pixels not in a slit are flagged with the
                OFFSLITS bit; see :func:`update_mask_slitmask`.
            from_scratch (:obj:`bool`, optional):
                Build the mask from scratch.  That is, if :attr:`fullmask`
                already exists, reset it to 0 and reflag all the bits.  Bad
                pixels and cosmic rays must be provided by :attr:`bpm` and
                :attr:`crmask`.
        """
        # Check input
        if saturation is not None and isinstance(saturation, np.ndarray) \
                and saturation.shape != self.shape:
            msgs.error('Saturation array must have the same shape as the image.')
        if mincounts is not None and isinstance(mincounts, np.ndarray) \
                and mincounts.shape != self.shape:
            msgs.error('Minimum counts array must have the same shape as the image.')

        # Setup the saturation level 
        if isinstance(saturation, str):
            if saturation != 'default':
                msgs.error(f'Unknown saturation string: {saturation}')
            _saturation = self.map_detector_value('saturation') \
                            * self.map_detector_value('nonlinear')
            if self.units == 'e-':
                _saturation *= self.map_detector_value('gain')
        else:
            _saturation = saturation

        # Setup the minimum counts level 
        if isinstance(mincounts, str):
            if mincounts != 'default':
                msgs.error(f'Unknown mincounts string: {mincounts}')
            _mincounts = self.map_detector_value('mincounts')
            if self.units == 'ADU':
                _mincounts /= self.map_detector_value('gain')
        else:
            _mincounts = mincounts

        if from_scratch:
            # Instatiate the mask
            self.reinit_mask()

        # Add the bad pixel mask
        if self.bpm is not None:
            self.update_mask('BPM', indx=self.bpm.astype(bool))

        # Add the cosmic rays
        if self.crmask is not None:
            self.update_mask('CR', indx=self.crmask.astype(bool))

        # Saturated pixels
        if _saturation is not None:
            self.update_mask('SATURATION', indx=self.image>=_saturation)

        # Minimum counts
        if _mincounts is not None:
            self.update_mask('MINCOUNTS', indx=self.image<=_mincounts)

        # Undefined counts
        self.update_mask('IS_NAN', indx=np.logical_not(np.isfinite(self.image)))

        if self.ivar is not None:
            # Bad inverse variance values
            self.update_mask('IVAR0', indx=np.logical_not(self.ivar > 0.0))
            # Undefined inverse variances
            self.update_mask('IVAR_NAN', indx=np.logical_not(np.isfinite(self.ivar)))

        if slitmask is not None:
            self.update_mask_slitmask(slitmask)

    def update_mask_slitmask(self, slitmask):
        """
        Update a mask using the slitmask

        Args:
            slitmask (`numpy.ndarray`_):
                Slitmask with -1 values pixels *not* in a slit

        """
        if slitmask.shape != self.shape:
            msgs.error('Slit mask image must have the same shape as data image.')
        # Pixels excluded from any slit.
        # TODO: Check me on this: Shouldn't we be turning off the old mask
        # before setting the new one?
        self.update_mask('OFFSLITS', action='turn_off')
        self.update_mask('OFFSLITS', indx=slitmask==-1)

    def update_mask_cr(self, crmask_new):
        """
        Update the mask bits for cosmic rays

        The original are turned off and the new
        ones are turned on.

        Args:
            crmask_new (`numpy.ndarray`_):
                New CR mask
        """
        self.update_mask('CR', action='turn_off')
        self.update_mask('CR', indx=crmask_new.astype(bool))

    def update_mask(self, flag, indx=None, action='turn_on'):
        """
        Update :attr:`fullmask` by operating on the bits for the provided (list
        of) flags.

        This method alters :attr:`fullmask` in-place.

        Args:
            flag (:obj:`str`, array-like):
                One or more flags to turn on for the selected pixels.
            indx (`numpy.ndarray`_, optional):
                Boolean array selecting the pixels in :attr:`fullmask` to alter.
                Must be the same shape as :attr:`fullmask`.  If None, *all*
                pixels are altered.
            action (:obj:`str`, optional):
                The action to perform.  Must be ``'turn_on'`` or ``'turn_off'``.
        """
        if action not in ['turn_on', 'turn_off']:
            msgs.error(f'{action} is not a known bit action!')
        if indx is None:
            self.fullmask = getattr(self.bitmask, action)(self.fullmask, flag)
            return
        if indx.shape != self.fullmask.shape:
            msgs.error('Array selecting pixels to update must be the same shape as fullmask.')
        self.fullmask[indx] = getattr(self.bitmask, action)(self.fullmask[indx], flag)

    def reinit_mask(self):
        """
        Reset the mask; :attr:`fullmask` is reset to be 0 everywhere.
        """
        if self.image is not None:
            self.fullmask = np.zeros(self.image.shape,
                                     dtype=self.bitmask.minimum_dtype(asuint=True))

    def select_flag(self, flag=None, invert=False):
        """
        Return a boolean array that selects pixels masked with the specified
        bits in :attr:`fullmask`.

        For example, to create a bad-pixel mask based on which pixels have
        cosmic-ray detections, run:

        .. code-block:: python

            cr_bpm = self.select_flag(flag='CR')

        Or, to create a good-pixel mask for all pixels that are not flagged for
        any reason, run:

        .. code-block:: python

            gpm = self.select_flag(invert=True)

        Args:
            flag (:obj:`str`, array-like, optional):
                One or more flags to select when returning the boolean mask.  If
                None, pixels flagged for any reason are returned as True.
            invert (:obj:`bool`, optional):
                If False, the return mask is True for masked pixels, False for
                good pixels (i.e., a bad-pixel mask).  If True, invert the sense
                of the mask (i.e., create a good-pixel mask, True for good
                pixels, False for bad pixels).
    
        Returns:
            `numpy.ndarray`_: Boolean array where pixels with the selected bits
            flagged are returned as True (if ``invert`` is False); i.e., this is
            a boolean bad-pixel mask.  If ``flag`` is not provided, pixels
            flagged for any reason are returned as True.
        """
        mask = np.zeros(self.image.shape, dtype=bool) if self.fullmask is None \
                    else self.bitmask.flagged(self.fullmask, flag=flag)
        return np.logical_not(mask) if invert else mask

    def sub(self, other, par):
        """
        Subtract one PypeItImage from another
        Extras (e.g. ivar, masks) are included if they are present

        Args:
            other (:class:`PypeItImage`):
            par (:class:`pypeit.par.pypeitpar.ProcessImagesPar`):
                Parameters that dictate the processing of the images.  See
                :class:`pypeit.par.pypeitpar.ProcessImagesPar` for the defaults
        Returns:
            PypeItImage:
        """
        if not isinstance(other, PypeItImage):
            msgs.error("Misuse of the subtract method")
        # Images
        newimg = self.image - other.image

        # Mask time
        outmask_comb = self.select_flag(invert=True) & other.select_flag(invert=True)

        # Variance
        if self.ivar is not None:
            new_ivar = utils.inverse(utils.inverse(self.ivar) + utils.inverse(other.ivar))
            new_ivar[np.logical_not(outmask_comb)] = 0
        else:
            new_ivar = None

        # RN2
        if self.rn2img is not None and other.rn2img is not None:
            new_rn2 = self.rn2img + other.rn2img
        else:
            new_rn2 = None

        # Instantiate
        new_sciImg = PypeItImage(image=newimg, ivar=new_ivar, bpm=self.bpm, rn2img=new_rn2,
                                 detector=self.detector)
        # Files
        new_sciImg.files = self.files + other.files

        #TODO: KW properly handle adding the bits
        #crmask_diff = new_sciImg.build_crmask(par) if par['mask_cr'] else np.zeros_like(other.image, dtype=bool)
        # crmask_eff assumes evertything masked in the outmask_comb is a CR in the individual images
        # JFH changed to below because this was not respecting the desire not to mask_crs
        new_sciImg.crmask = (new_sciImg.build_crmask(par) | np.logical_not(outmask_comb)) if par['mask_cr'] else np.logical_not(outmask_comb)
        #new_sciImg.crmask = crmask_diff | np.logical_not(outmask_comb)
        # Note that the following uses the saturation and mincounts held in self.detector
        new_sciImg.build_mask()

        return new_sciImg

    def show(self):
        """
        Show the image in a ginga viewer.
        """
        if self.image is None:
            # TODO: This should fault.
            msgs.warn("No image to show!")
            return
        display.show_image(self.image, chname='image')

    def __repr__(self):
        repr = '<{:s}: '.format(self.__class__.__name__)
        # Image
        rdict = {}
        for attr in self.datamodel.keys():
            if hasattr(self, attr) and getattr(self, attr) is not None:
                rdict[attr] = True
            else:
                rdict[attr] = False
        repr += ' images={}'.format(rdict)
        repr = repr + '>'
        return repr


