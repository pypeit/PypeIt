""" Module for image mask related items """
from IPython import embed

import numpy as np

from pypeit.bitmask import BitMask


class ImageBitMask(BitMask):
    """
    Define a bitmask used to set the reasons why each pixel in a science
    image was masked.
    """

    def __init__(self):
        # TODO:
        #   - Can IVAR0 and IVAR_NAN be consolidated into a single bit?
        #   - Is EXTRACT ever set?
        mask_bits = dict(BPM='Component of the instrument-specific bad pixel mask',
                         CR='Cosmic ray detected',
                         SATURATION='Saturated pixel',
                         MINCOUNTS='Pixel below the instrument-specific minimum counts',
                         OFFSLITS='Pixel does not belong to any slit',
                         IS_NAN='Pixel value is undefined',
                         IVAR0='Inverse variance is undefined',
                         IVAR_NAN='Inverse variance is NaN',
                         EXTRACT='Pixel masked during local skysub and extraction',
                         BADSCALE='Bad image rescaling operation (e.g., flat value <= 0)',
                         STCKMASK='All pixels masked in image stack')
        super(ImageBitMask, self).__init__(list(mask_bits.keys()), descr=list(mask_bits.values()))


