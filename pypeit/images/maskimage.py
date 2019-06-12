""" Module for image mask related items """
from collections import OrderedDict
import numpy as np

from pypeit.bitmask import BitMask
from pypeit.core import procimg
from pypeit import utils


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

    bitmask = ImageBitMask()

    def __init__(self, bpm, crmask=None):

        self.bpm = bpm
        self.crmask = crmask

        # Internals
        self.mask = None

    def build_crmask(self, spectrograph, det, par, image, rawvarframe):
        """
        Generate the CR mask frame

        Wrapper to procimg.lacosmic

        Requires self.rawvarframe to exist

        Returns:
            np.ndarray: Copy of self.crmask

        """
        # Run LA Cosmic to get the cosmic ray mask
        self.crmask = procimg.lacosmic(det, image,
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

    def build_mask(self, image, sciivar, saturation=1e10,
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
            saturation (float, optional):
                Saturation limit in ADU
            mincounts (float, optional):
            slitmask (np.ndarray, optional):
                Slit mask image;  Pixels not in a slit are masked
            bpm (np.ndarray, optional):
                Bad pixel mask

        Returns:
            numpy.ndarray: Copy of the bit value mask for the science image.
        """
        # Instatiate the mask
        mask = np.zeros_like(image, dtype=self.bitmask.minimum_dtype(asuint=True))

        # Bad pixel mask
        indx = self.bpm.astype(bool)
        mask[indx] = self.bitmask.turn_on(mask[indx], 'BPM')

        # Cosmic rays
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

        # Bad inverse variance values
        indx = np.invert(sciivar > 0.0)
        mask[indx] = self.bitmask.turn_on(mask[indx], 'IVAR0')

        # Undefined inverse variances
        indx = np.invert(np.isfinite(sciivar))
        mask[indx] = self.bitmask.turn_on(mask[indx], 'IVAR_NAN')

        if slitmask is not None:
            indx = slitmask == -1
            mask[indx] = self.bitmask.turn_on(mask[indx], 'OFFSLITS')

        return mask



