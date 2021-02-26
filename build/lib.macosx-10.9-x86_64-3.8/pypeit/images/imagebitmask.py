""" Module for image mask related items """
import numpy as np

from pypeit.bitmask import BitMask
from collections import OrderedDict

from IPython import embed

# I am keeping this here to avoid circular imports as even core routines need occasional access

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


