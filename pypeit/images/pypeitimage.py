""" Module for the PypeItImage include its Mask

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import numpy as np
import os
import inspect

from pypeit import msgs
from pypeit.images import detector_container, imagebitmask
from pypeit.core import procimg
from pypeit.display import display
from pypeit import datamodel
from pypeit import utils
from pypeit import masterframe

from IPython import embed


class PypeItImage(datamodel.DataContainer):
    """
    Class to hold a single image from a single detector in PypeIt
    and its related images (e.g. ivar, mask).

    Oriented in its spec,spat format

    The intent is to keep this object as light-weight as possible.

    Args:
        image (`numpy.ndarray`_ or None):
            See datamodel for description
        ivar (`numpy.ndarray`_, optional):
        rn2img (`numpy.ndarray`_, optional):
        bpm (`numpy.ndarray`_, optional):
        crmask (`numpy.ndarray`_, optional):
        fullmask (`numpy.ndarray`_, optional):
        detector (:class:`pypeit.images.data_container.DataContainer`):
        spat_flexure (:obj:`float`, optional):

    Attributes:
        head0 (astropy.io.fits.Header):
        detector (:class:`pypeit.images.detector_container.DetectorContainer`):
        files (list):
        rawheadlst (list):
            List containing headers of the raw image file
        master_key (str):
            Master key, only for Master frames
        master_dir (str):
            Master key, only for Master frames

    """
    version = '1.1.0'
    """Datamodel version number"""

    datamodel = {'image': dict(otype=np.ndarray, atype=np.floating, descr='Primary image data'),
                 'ivar': dict(otype=np.ndarray, atype=np.floating,
                              descr='Inverse variance image'),
                 'nimg': dict(otype=np.ndarray, atype=np.integer,
                              descr='If a combination of multiple images, this is the number of '
                                    'images that contributed to each pixel'),
                 'rn2img': dict(otype=np.ndarray, atype=np.floating,
                                descr='Read noise squared image'),
                 'base_var': dict(otype=np.ndarray, atype=np.floating,
                                  descr='Base-level image variance, excluding count shot-noise'),
                 'img_scale': dict(otype=np.ndarray, atype=np.floating,
                                   descr='Image count scaling applied (e.g., 1/flat-field)'),
                 'bpm': dict(otype=np.ndarray, atype=np.integer, descr='Bad pixel mask'),
                 'crmask': dict(otype=np.ndarray, atype=np.bool_, descr='CR mask image'),
                 'fullmask': dict(otype=np.ndarray, atype=np.integer, descr='Full image bitmask'),
                 'detector': dict(otype=detector_container.DetectorContainer,
                                  descr='Detector DataContainer'),
                 'PYP_SPEC': dict(otype=str, descr='PypeIt spectrograph name'),
                 # TODO: Use BUNIT instead?  This has a specific meaning in the
                 # FITS standard, so I'm not sure.
                 'units': dict(otype=str, descr='(Unscaled) Pixel units (e- or ADU)'),
                 # TODO: Consider forcing exptime to be a float.
                 'exptime': dict(otype=(int, float), descr='Effective exposure time (s)'),
                 'noise_floor': dict(otype=float, descr='Noise floor included in variance'),
                 'shot_noise': dict(otype=bool, descr='Shot-noise included in variance'),
                 'spat_flexure': dict(otype=float,
                                      descr='Shift, in spatial pixels, between this image '
                                            'and SlitTrace'),
                 'imgbitm': dict(otype=str, descr='List of BITMASK keys from ImageBitMask')}
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
            _d[key] = pypeitImage[key]
        # Instantiate
        slf = cls(**_d)
        # Internals
        slf.master_dir = pypeitImage.master_dir
        slf.master_key = pypeitImage.master_key
        # Return
        return slf

    # This needs to contain all datamodel items.
    # TODO: Not really. You don't have to pass everything to the
    # super().__init__ call...
    def __init__(self, image=None, ivar=None, nimg=None, rn2img=None, base_var=None,
                 img_scale=None, bpm=None, crmask=None, fullmask=None, detector=None,
                 spat_flexure=None, PYP_SPEC=None, units=None, exptime=None, noise_floor=None,
                 shot_noise=None, imgbitm=None):

        # Setup the DataContainer. Dictionary elements include
        # everything but self in the instantiation call.
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        _d = {k: values[k] for k in args[1:]}
        # Init
        super(PypeItImage, self).__init__(d=_d)

    def _init_internals(self):
        self.head0 = None
        self.process_steps = None
        self.files = None
        self.rawheadlist = None
        # Master stuff
        self.master_key = None
        self.master_dir = None

    def _validate(self):
        """
        Validate the slit traces.
        """
        if self.imgbitm is None:
            self.imgbitm = ','.join(list(self.bitmask.keys()))
        else:
            # Validate
            if self.imgbitm != ','.join(list(self.bitmask.keys())):
                msgs.error("Input BITMASK keys differ from current data model!")

    def _bundle(self):
        """
        Over-write default _bundle() method to write one
        HDU per image.  Any extras are in the HDU header of
        the primary image.

        Returns:
            :obj:`list`: A list of dictionaries, each list element is
            written to its own fits extension. See the description
            above.
        """
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
        return () if self.image is None else self.image.shape

    def build_crmask(self, par, subtract_img=None):
        """
        Generate the CR mask frame

        Mainly a wrapper to :func:`pypeit.core.procimg.lacosmic`

        Args:
            par (:class:`pypeit.par.pypeitpar.ProcessImagesPar`):
                Parameters that dictate the processing of the images.  See
                :class:`pypeit.par.pypeitpar.ProcessImagesPar` for the
                defaults.
            subtract_img (`numpy.ndarray`_, optional):
                If provided, subtract this from the image prior to CR detection

        Returns:
            `numpy.ndarray`_: Copy of self.crmask (boolean)

        """
        var = utils.inverse(self.ivar)
        use_img = self.image if subtract_img is None else self.image - subtract_img
        # Run LA Cosmic to get the cosmic ray mask
        self.crmask = procimg.lacosmic(use_img,
                                       self.detector['saturation'],
                                       self.detector['nonlinear'],
                                       varframe=var,
                                       maxiter=par['lamaxiter'],
                                       grow=par['grow'],
                                       remove_compact_obj=par['rmcompact'],
                                       sigclip=par['sigclip'],
                                       sigfrac=par['sigfrac'],
                                       objlim=par['objlim'])
        # Return
        return self.crmask.copy()

    def build_mask(self, saturation=None, mincounts=None, slitmask=None):
        """
        Construct the bit value mask used during extraction.

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
            saturation (:obj:`float`, optional):
                Saturation limit in counts or ADU (needs to match the input image)
                Defaults to self.detector['saturation']
            slitmask (`numpy.ndarray`_, optional):
                Slit mask image;  Pixels not in a slit are masked
            mincounts (:obj:`float`, optional):
                Defaults to self.detector['mincounts']
        """
        _mincounts = self.detector['mincounts'] if mincounts is None else mincounts
        _saturation = self.detector['saturation'] if saturation is None else saturation
        # Instatiate the mask
        self.reinit_mask()

        # Bad pixel mask
        if self.bpm is not None:
            self.update_mask('BPM', indx=self.bpm.astype(bool))

        # Cosmic rays
        if self.crmask is not None:
            self.update_mask('CR', indx=self.crmask.astype(bool))

        # Saturated pixels
        self.update_mask('SATURATION', indx=self.image>=_saturation)

        # Minimum counts
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
        # Pixels excluded from any slit.
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
        self.fullmask = np.zeros(self.image.shape, dtype=self.bitmask.minimum_dtype(asuint=True))

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


