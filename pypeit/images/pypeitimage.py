""" Module for the PypeItImage include its Mask

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import os
import inspect

from IPython import embed

import numpy as np

from pypeit import msgs
from pypeit.images.imagebitmask import ImageBitMaskArray
from pypeit.images.detector_container import DetectorContainer
from pypeit.images.mosaic import Mosaic
from pypeit.core import procimg
from pypeit.display import display
from pypeit import datamodel
from pypeit import utils
from pypeit import masterframe


class PypeItImage(datamodel.DataContainer):
    r"""
    Container class for processed PypeIt images and associated data.

    Image data can be either from a single detector, multiple detectors arranged
    in a 3D array, or multiple detectors resampled to a single mosaic image.
    The orientation of the images is always spectral pixels along the first axis
    and spatial pixels along the second axis; i.e., the shape is :math:`(N_{\rm
    spec},N_{\rm spat})`.

    Class instantiation only requires the ``image`` data; all other elements of
    the datamodel are completely optional.  If not defined, the datamodel
    elements are set to None, with the exception of the image mask (see below).
    Any code that uses a ``PypeItImage`` and accesses the datamodel components
    *must* be able to handle when any of the values are None.

    The datamodel components are:

    .. include:: ../include/class_datamodel_pypeitimage.rst

    Regardless of whether or not it is provided directly (see ``fullmask``), the
    image mask is always instantiated.  If ``fullmask`` is not provided, this
    instantiation begins with all values being unmasked.  Although the class
    allows direct access/manipulation of ``fullmask`` (see
    :class:`~pypeit.images.imagebitmask.ImageBitMaskArray`), convenience
    functions are provided to interface with the underlying object;  see
    :func:`update_mask` and :func:`select_flag`, as well as
    :class:`~pypeit.images.imagebitmask.ImageBitMask` for the valid flag names.
    
    Additionally, pixels can be flagged on instantiation as being part of a
    bad-pixel mask (``bpm``), a cosmic ray hit (``crmask``), or a user-level
    bad-pixel mask (``usermask``).  Importantly, if both ``fullmask`` and one of
    these masks are provided, the default behavior is to *add*, e.g., the
    ``crmask`` to the list of already flagged CRs provided by ``fullmask`` (if
    there are any).  To remove any of the relevant flags before including the
    new flags in ``fullmask`` use ``clean_mask``.  See the class parameters
    below.

    The following lists only those parameters that are *not* part of the class
    datamodel (see above).
    
    Args:
        bpm (`numpy.ndarray`_, optional):
            The image bad-pixel mask, which typically selects image values that
            contain detector artifacts.  This is a boolean array that must have
            the same shape as ``image``.
        crmask (`numpy.ndarray`_, optional):
            A boolean array selecting pixels that contain cosmic-ray hits.
            Shape must match ``image``.
        usermask (`numpy.ndarray`_, optional):
            A boolean array selecting pixels that the user wishes to ignore.
            Shape must match ``image``.
        clean_mask (:obj:`bool`, optional):
            If true, first remove any existing BPM, CR, or USER mask from
            ``fullmask`` before including the masks provided by ``bpm``,
            ``crmask``, or ``usermask``.

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

    version = '1.3.0'
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
                 'fullmask': dict(otype=ImageBitMaskArray, descr='Image mask'),
                 'detector': dict(otype=(DetectorContainer, Mosaic),
                                  descr='The detector (see :class:`~pypeit.images.detector_container.DetectorContainer`) '
                                        'or mosaic (see :class:`~pypeit.images.mosaic.Mosaic`) '
                                        'parameters'),
                 'PYP_SPEC': dict(otype=str, descr='PypeIt spectrograph name'),
                 'units': dict(otype=str, descr='(Unscaled) Pixel units (e- or ADU)'),
                 # TODO: Consider forcing exptime to be a float.
                 'exptime': dict(otype=(int, float), descr='Effective exposure time (s)'),
                 'noise_floor': dict(otype=float, descr='Noise floor included in variance'),
                 'shot_noise': dict(otype=bool, descr='Shot-noise included in variance'),
                 'spat_flexure': dict(otype=float,
                                      descr='Shift, in spatial pixels, between this image '
                                            'and SlitTrace')}
    """Data model components."""

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

    def __init__(self, image, ivar=None, nimg=None, amp_img=None, det_img=None, rn2img=None,
                 base_var=None, img_scale=None, fullmask=None, detector=None, spat_flexure=None,
                 PYP_SPEC=None, units=None, exptime=None, noise_floor=None, shot_noise=None,
                 bpm=None, crmask=None, usermask=None, clean_mask=False):

        if image is None:
            msgs.error('Must provide an image when instantiating PypeItImage.')

        # Instantiate as an empty DataContainer
        super().__init__()

        # Get the list of arguments and values provided to the instantiation call.
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        # Loop through them all and assign them directly to attributes of the
        # class if they're part of the datamodel and *not* None.
        for key in args[1:]:
            if key in self.datamodel and values[key] is not None:
                setattr(self, key, values[key])

        # Make sure all the relevant arrays have the same shape
        for attr in ['ivar', 'nimg', 'amp_img', 'det_img', 'rn2img', 'base_var', 'img_scale',
                     'fullmask']:
            _arr = getattr(self, attr)
            if _arr is not None and _arr.shape != self.shape:
                msgs.error(f'Attribute {attr} does not match image shape.')

        # Make sure the units are defined
        if self.units is None:
            self.units = 'e-'

        # Reset the masks if they're provided.
        if fullmask is None:
            self.reinit_mask()

        if bpm is not None:
            if not np.issubdtype(bpm.dtype, np.bool_) and not np.issubdtype(bpm.dtype, bool):
                msgs.error('CODING ERROR: bpm entry in PypeItImage must have boolean type')
            if clean_mask:
                self.update_mask('BPM', action='turn_off')
            self.update_mask('BPM', indx=bpm)
        if crmask is not None:
            if not np.issubdtype(crmask.dtype, np.bool_) and not np.issubdtype(crmask.dtype, bool):
                msgs.error('CODING ERROR: crmask entry in PypeItImage must have boolean type')
            if clean_mask:
                self.update_mask('CR', action='turn_off')
            self.update_mask('CR', indx=crmask)
        if usermask is not None:
            if not np.issubdtype(usermask.dtype, np.bool_) and not np.issubdtype(usermask.dtype, bool):
                msgs.error('CODING ERROR: usermask entry in PypeItImage must have boolean type')
            if clean_mask:
                self.update_mask('USER', action='turn_off')
            self.update_mask('USER', indx=crmask)

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

    def _bundle(self):
        """
        Package the datamodel for writing.

        Returns:
            :obj:`list`: A list of dictionaries, each list element is written to
            its own fits extension. See
            :class:`~pypeit.datamodel.DataContainer`.
        """
        keys = list(self.keys())
        # Primary image
        d = [dict(image=self.image)]
        keys.remove('image')

        # Rest of the datamodel
        for key in keys:
            # Skip None
            if self[key] is None:
                continue
            # Array?
            if self.datamodel[key]['otype'] == np.ndarray \
                    or isinstance(self[key], datamodel.DataContainer):
                d.append({key : self[key]})
            else: # Add to header of the primary image
                d[0][key] = self[key]

        # TODO: Add `files` and `process_steps`?

        # Return
        return d

    @classmethod
    def from_hdu(cls, hdu, hdu_prefix=None, chk_version=True):
        """
        Instantiate the object from an HDU extension.

        This overrides the base-class method. Overriding this method
        is preferrable to overriding the ``_parse`` method because it
        makes it easier to deal with the multiple
        :class:`~pypeit.datamodel.DataContainer` objects contained by
        :class:`PypeItImage`.

        Args:
            hdu (`astropy.io.fits.HDUList`_, `astropy.io.fits.ImageHDU`_, `astropy.io.fits.BinTableHDU`_):
                The HDU(s) with the data to use for instantiation.
            hdu_prefix (:obj:`str`, optional):
                Maintained for consistency with the base class but is
                not used by this method.
            chk_version (:obj:`bool`, optional):
                If True, raise an error if the datamodel version or
                type check failed. If False, throw a warning only.
        """
        # Set the default hdu prefix if none is provided
        if hdu_prefix is None:
            hdu_prefix = cls.hdu_prefix

        # Need to separately parse the mask because it is not called 'MASK'.
        # Set the mask extension name.
        mask_ext = 'FULLMASK' if hdu_prefix is None else f'{hdu_prefix}FULLMASK'

        # Get all the extensions
        ext = [h.name for h in hdu] if hasattr(hdu, '__len__') else [hdu.name]

        # Remove the mask extension if it's there
        if mask_ext in ext:
            ext.remove(mask_ext)

        # Parse everything *but* the mask extension
        d, version_passed, type_passed, parsed_hdus \
                = super()._parse(hdu, ext=ext, hdu_prefix=hdu_prefix)
        if not type_passed:
            msgs.error(f'The HDU(s) cannot be parsed by a {cls.__name__} object!')
        if not version_passed:
            _f = msgs.error if chk_version else msgs.warn
            _f(f'Current version of {cls.__name__} object in code (v{cls.version})'
               ' does not match version used to write your HDU(s)!')

        if mask_ext in hdu:
            # If the mask extension exists, parse it
            d['fullmask'] = ImageBitMaskArray.from_hdu(hdu[mask_ext], ext_pseudo='MASK',
                                                       chk_version=chk_version)

        # Instantiate
        return cls.from_dict(d=d)

    @property
    def shape(self):
        """
        Return the primary image shape.
        """
        return self.image.shape

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
            `numpy.ndarray`_: Boolean array with flagging pixels with identified
            cosmic rays; True mean a CR was flagged.
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
        # bad pixel mask, but what other flags from ``fullmask`` should be
        # included?
        bpm = self.fullmask.bpm

        # TODO: These saturation and non-linear values are typically for the raw
        # pixel value.  E.g., a saturation of 65535 is because the digitization
        # is for a 16-bit integer.  But this function is called on an image that
        # has likely already been bias-subtracted, etc.  Are we doing this
        # right?
        saturation = self.map_detector_value('saturation')
        nonlinear = self.map_detector_value('nonlinear')

        if self.is_multidetector:
            # If the object has multiple images, need to flag each image individually
            crmask = np.empty(self.shape, dtype=bool)
            for i in range(self.shape[0]):
                crmask[i] = procimg.lacosmic(use_img[i], saturation=saturation[i],
                                             nonlinear=nonlinear[i], bpm=bpm[i], varframe=var[i],
                                             maxiter=par['lamaxiter'], grow=par['grow'],
                                             remove_compact_obj=par['rmcompact'],
                                             sigclip=par['sigclip'], sigfrac=par['sigfrac'],
                                             objlim=par['objlim'])
        else:
            # Otherwise, just run LA Cosmic once
            crmask = procimg.lacosmic(use_img, saturation=saturation, nonlinear=nonlinear,
                                      bpm=bpm, varframe=var, maxiter=par['lamaxiter'],
                                      grow=par['grow'], remove_compact_obj=par['rmcompact'],
                                      sigclip=par['sigclip'], sigfrac=par['sigfrac'],
                                      objlim=par['objlim'])
        # Update the mask (this erases any existing CR mask!)
        self.update_mask_cr(crmask)
        # Return the result
        return crmask

    def map_detector_value(self, attr):
        """
        Provided a detector-specific value, remap it as necessary for numpy
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
        pixel with ``img.fullmask.mask == 0`` is valid, otherwise the pixel has been
        masked.  To determine why a given pixel has been masked (see
        also :func:`~pypeit.images.pypeitimage.PypeItImage.select_flag`):

        .. code-block:: python

            reasons = img.fullmask.flagged_bits([0,0])

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
                already exists, all bits are turned off except for BPM and CR.
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
            # Save the existing BPM and CR masks
            bpm = self.fullmask.bpm
            cr = self.fullmask.cr
            # Re-initialize the fullmask (erases all existing masks)
            self.reinit_mask()
            # Recover the BPM and CR masks
            self.update_mask('BPM', indx=bpm)
            self.update_mask('CR', indx=cr)

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
            indx (:obj:`tuple`, :obj:`slice`, `numpy.ndarray`_, optional):
                Object used to select elements of the mask array to at which to
                toggle the provided bit flags.  I.e., for the internal
                :attr:`fullmask`, ``fullmask.mask[indx]`` must be a valid (fancy
                indexing) operation.  If None, the action is performed for the
                full mask!
            action (:obj:`str`, optional):
                The action to perform.  Must be ``'turn_on'`` or ``'turn_off'``.
        """
        if action not in ['turn_on', 'turn_off']:
            msgs.error(f'{action} is not a known bit action!')
        if indx is None:
            getattr(self.fullmask, action)(flag)
        getattr(self.fullmask, action)(flag, select=indx)

    def reinit_mask(self):
        """
        Reset the mask; :attr:`fullmask` is reset to be 0 everywhere.
        """
        self.fullmask = ImageBitMaskArray(self.image.shape)

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
                None, pixels flagged for *any* reason are returned as True.
            invert (:obj:`bool`, optional):
                If False, the return mask is True for masked pixels, False for
                good pixels (i.e., a bad-pixel mask).  If True, invert the sense
                of the mask (i.e., create a good-pixel mask, True for good
                pixels, False for bad pixels).
    
        Returns:
            `numpy.ndarray`_: Boolean array where pixels with the selected bits
            flagged are returned as True (if ``invert`` is False); i.e., this is
            a boolean bad-pixel mask (or a good-pixel mask when ``invert`` is
            True).  If ``flag`` is not provided, pixels flagged for any reason
            are returned as True.
        """
        return self.fullmask.flagged(flag=flag, invert=invert)

    def __sub__(self, other):
        return self.sub(other)

    def sub(self, other):
        """
        Subtract this PypeItImage from another.
        
        The following operations are performed:

            - the image data is subtracted (images must have the same shape)
            - the inverse variance (:attr:`ivar`) is propagated
            - the number of images is combined (:attr:`nimg`)
            - the RN2 (:attr:`rn2img`) is propagated
            - the base variance (:attr:`base_var`) is propagated
            - the image scaling (:attr:`img_scale`) is averaged
            - the bit mask (:attr:`fullmask`) is joined (using an or operation)
            - if it's the same for both images, the spectrograph name
              (:attr:`PYP_SPEC`) is propagated
            - if it's the same for both images, the images units (:attr:`units`)
              is propagated
            - if both images provide source file names, the file lists are
              concatenated
            - the detector from the first image (``self``) is used for the
              returned image and the detector for the ``other`` image is
              *ignored*

        Args:
            other (:class:`~pypeit.images.pypeitimage.PypeItImage`):
                The image to subtract

        Returns:
            :class:`~pypeit.images.pypeitimage.PypeItImage`: The result of
            subtracting ``other`` from this image.
        """
        if not isinstance(other, PypeItImage):
            msgs.error('Image to subtract must be of type PypeItImage.')

        # Subtract the image
        newimg = self.image - other.image

        # Propagate the variance, if available
        if self.ivar is not None or other.ivar is not None:
            new_ivar = np.zeros(newimg.shape) if self.ivar is None else utils.inverse(self.ivar)
            if other.ivar is not None:
                new_ivar += utils.inverse(other.ivar)
            new_ivar = utils.inverse(new_ivar)
        else:
            new_ivar = None

        # Only create a new nimg if it's present in *both* images
        new_nimg = self.nimg + other.nimg if self.nimg is not None and other.nimg is not None \
                        else None

        # RN2
        if self.rn2img is not None or other.rn2img is not None:
            new_rn2 = np.zeros(newimg.shape) if self.rn2img is None else self.rn2img
            if other.rn2img is not None:
                new_rn2 += other.rn2img
        else:
            new_rn2 = None

        # Base variance (for the noise modeling, it's more important that this
        # is propagated compared to rn2)
        if self.base_var is not None or other.base_var is not None:
            new_base = np.zeros(newimg.shape) if self.base_var is None else self.base_var
            if other.base_var is not None:
                new_base += other.base_var
        else:
            new_base = None

        # Image scaling
        # TODO: This is bogus.  Maybe we should just set this to 1?  Either way,
        # trying to model the noise using our current approach won't be
        # mathematically correct for two subtracted images.  This will be worse
        # for more significant scaling.
        if self.img_scale is not None or other.img_scale is not None:
            new_img_scale = np.ones(newimg.shape) if self.img_scale is None else self.img_scale
            if other.img_scale is not None:
                new_img_scale = (new_img_scale + other.img_scale) / 2
        else:
            new_img_scale = None

        # Mask
        if self.fullmask is not None or other.fullmask is not None:
            new_fullmask = ImageBitMaskArray(newimg.shape) if self.fullmask is None \
                                else self.fullmask
            if other.fullmask is not None:
                new_fullmask |= other.fullmask
        else:
            new_fullmask = None

        # PYP_SPEC
        # TODO: Instead raise an error if they're not the same
        new_spec = self.PYP_SPEC if self.PYP_SPEC == other.PYP_SPEC else None

        # units
        # TODO: Instead raise an error if they're not the same
        new_units = self.units if self.units == other.units else None
        
        # spatial flexure
        # TODO: for now we assume the flexure is the same in both images
        # since doing anything else is not really feasible
        if self.spat_flexure is not None:
            new_spat_flexure = self.spat_flexure
        else:
            new_spat_flexure = None

        # Create the new image
        # TODO: We should instead *copy* the detector object; otherwise, it's
        # possible that it will be shared between multiple images.  Nominally,
        # this should be okay because the detector data is meant to be static,
        # but we should fix this.
        new_pypeitImage = PypeItImage(newimg, ivar=new_ivar, nimg=new_nimg, rn2img=new_rn2,
                                      base_var=new_base, img_scale=new_img_scale,
                                      fullmask=new_fullmask, detector=self.detector,
                                      PYP_SPEC=new_spec, units=new_units,
                                      spat_flexure=new_spat_flexure)

        # Files
        if self.files is not None and other.files is not None:
            new_pypeitImage.files = self.files + other.files

        return new_pypeitImage

    def show(self):
        """
        Show the image in a ginga viewer.
        """
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


