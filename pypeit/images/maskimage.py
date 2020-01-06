""" Module for image mask related items """
from collections import OrderedDict
import numpy as np

from pypeit.bitmask import BitMask
from pypeit.core import procimg

from IPython import embed


class ImageBitMask(BitMask):
    """
    Define a bitmask used to set the reasons why each pixel in a science
    image was masked.
    """

    def __init__(self):
        # TODO:
        #   - Can IVAR0 and IVAR_NAN be consolidated into a single bit?
        #   - Is EXTRACT ever set?
        # TODO: This needs to be an OrderedDict for now to ensure that
        # the bits assigned to each key is always the same. As of python
        # 3.7, normal dict types are guaranteed to preserve insertion
        # order as part of its data model. When/if we require python
        # 3.7, we can remove this (and other) OrderedDict usage in favor
        # of just a normal dict.
        mask_dict = OrderedDict([
            ('BPM', 'Component of the instrument-specific bad pixel mask'),
            ('CR', 'Cosmic ray detected'),
            ('SATURATION', 'Saturated pixel'),
            ('MINCOUNTS', 'Pixel below the instrument-specific minimum counts'),
            ('OFFSLITS', 'Pixel does not belong to any slit'),
            ('IS_NAN', 'Pixel value is undefined'),
            ('IVAR0', 'Inverse variance is undefined'),
            ('IVAR_NAN', 'Inverse variance is NaN'),
            ('EXTRACT', 'Pixel masked during local skysub and extraction')
        ])
        super(ImageBitMask, self).__init__(list(mask_dict.keys()), descr=list(mask_dict.values()))

class ImageMask(object):
    """
    Class to handle masks associated with an Image

    Args:
        bpm (np.ndarray):
            Bad pixel mask
        crmask (np.ndarray, optional):
            Cosmic Ray mask (boolean)

    Attributes:
        mask (np.ndarray):
            The bitmask values for the full mask
        pars (dict):
            Used to hold parameters used when creating masks
    """

    bitmask = ImageBitMask()

    def __init__(self, bpm, crmask=None):

        self.bpm = bpm
        self.crmask = crmask

        # Internals
        self.mask = None

        # Data model
        self.mask_attributes = ('bpm', 'crmask', 'mask')

    def build_crmask(self, spectrograph, det, par, image, rawvarframe, subtract_img=None):
        """
        Generate the CR mask frame

        Mainly a wrapper to procimg.lacosmic

        Args:
            spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
                Spectrograph used to take the data.
            det (:obj:`int`, optional):
                The 1-indexed detector number to process.
            par (:class:`pypeit.par.pypeitpar.ProcessImagesPar`):
                Parameters that dictate the processing of the images.  See
                :class:`pypeit.par.pypeitpar.ProcessImagesPar` for the
                defaults.
            image (np.ndarray):
                Image to identify CR's in
            rawvarframe (np.ndarray):
                Variance image
            subtract_img (np.ndarray, optional):
                If provided, subtract this from the image prior to CR detection

        Returns:
            np.ndarray: Copy of self.crmask (boolean)

        """
        use_img = image if subtract_img is None else image - subtract_img
        # Run LA Cosmic to get the cosmic ray mask
        self.crmask = procimg.lacosmic(det, use_img,
                                  spectrograph.detector[det-1]['saturation'],
                                  spectrograph.detector[det-1]['nonlinear'],
                                  varframe=rawvarframe,
                                  maxiter=par['lamaxiter'],
                                  grow=par['grow'],
                                  remove_compact_obj=par['rmcompact'],
                                  sigclip=par['sigclip'],
                                  sigfrac=par['sigfrac'],
                                  objlim=par['objlim'])
        # Return
        return self.crmask.copy()

    def build_mask(self, image, ivar, saturation=1e10,
                   mincounts=-1e10, slitmask=None):
        """
        Return the bit value mask used during extraction.

        The mask keys are defined by :class:`ScienceImageBitMask`.  Any
        pixel with mask == 0 is valid, otherwise the pixel has been
        masked.  To determine why a given pixel has been masked::

            bitmask = ScienceImageBitMask()
            reasons = bm.flagged_bits(mask[i,j])

        To get all the pixel masked for a specific set of reasons::

            indx = bm.flagged(mask, flag=['CR', 'SATURATION'])

        Args:
            image (np.ndarray):
                Image
            ivar (np.ndarray or None):
                Inverse variance of the input image
            saturation (float, optional):
                Saturation limit in counts or ADU (needs to match the input image)
            slitmask (np.ndarray, optional):
                Slit mask image;  Pixels not in a slit are masked
            mincounts (float, optional):

        Returns:
            numpy.ndarray: Copy of the bit value mask for the science image.
        """
        # Instatiate the mask
        mask = np.zeros_like(image, dtype=self.bitmask.minimum_dtype(asuint=True))

        # Bad pixel mask
        indx = self.bpm.astype(bool)
        mask[indx] = self.bitmask.turn_on(mask[indx], 'BPM')

        # Cosmic rays
        if self.crmask is not None:
            indx = self.crmask.astype(bool)
            mask[indx] = self.bitmask.turn_on(mask[indx], 'CR')

        # Saturated pixels
        indx = image >= saturation
        mask[indx] = self.bitmask.turn_on(mask[indx], 'SATURATION')

        # Minimum counts
        indx = image <= mincounts
        mask[indx] = self.bitmask.turn_on(mask[indx], 'MINCOUNTS')

        # Undefined counts
        indx = np.invert(np.isfinite(image))
        mask[indx] = self.bitmask.turn_on(mask[indx], 'IS_NAN')

        if ivar is not None:
            # Bad inverse variance values
            indx = np.invert(ivar > 0.0)
            mask[indx] = self.bitmask.turn_on(mask[indx], 'IVAR0')

            # Undefined inverse variances
            indx = np.invert(np.isfinite(ivar))
            mask[indx] = self.bitmask.turn_on(mask[indx], 'IVAR_NAN')

        if slitmask is not None:
            indx = slitmask == -1
            mask[indx] = self.bitmask.turn_on(mask[indx], 'OFFSLITS')

        self.mask = mask
        return self.mask.copy()

    def update_mask_slitmask(self, slitmask):
        """
        Update a mask using the slitmask

        Args:
            slitmask (np.ndarray):
                Slitmask with -1 values pixels *not* in a slit

        """
        # Pixels excluded from any slit.
        indx = slitmask == -1
        # Finish
        self.mask[indx] = self.bitmask.turn_on(self.mask[indx], 'OFFSLITS')

    def update_mask_cr(self, crmask_new):
        """
        Update the mask bits for cosmic rays

        The original are turned off and the new
        ones are turned on.

        Args:
            crmask_new (np.ndarray):
                New CR mask
        """
        '''
        # Unset the CR bit from all places where it was set
        CR_old = (self.bitmask.unpack(self.mask, flag='CR'))[0]
        mask_new = np.copy(self.mask)
        mask_new[CR_old] = self.bitmask.turn_off(mask_new[CR_old], 'CR')
        # Now set the CR bit using the new crmask
        indx = crmask_new.astype(bool)
        mask_new[indx] = self.bitmask.turn_on(mask_new[indx], 'CR')
        # Save
        self.mask = mask_new.copy()
        '''
        self.mask = self.bitmask.turn_off(self.mask, 'CR')
        indx = crmask_new.astype(bool)
        self.mask[indx] = self.bitmask.turn_on(self.mask[indx], 'CR')



