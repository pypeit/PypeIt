""" Object to hold + process a single image"""

from pypeit import msgs
from pypeit import ginga

import numpy as np
from pypeit.images import maskimage

from IPython import embed


class PypeItImage(maskimage.ImageMask):
    """
    Class to hold a single image from a single detector in PypeIt
    Oriented in its spec,spat format

    Args:
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            Spectrograph used to take the data.
        det (:obj:`int`, optional):
            The 1-indexed detector number to process.

    Attributes:
        image (np.ndarray):
        datasec_img (np.ndarray):
            Used for the amplifiers
        head0 (astropy.io.fits.Header):
        orig_shape (tuple):
        binning_raw (tuple):  Binning in the raw image orientation (NAXIS1, NAXIS2)
        binning (tuple): Binning the PypeIt orientation (spec, spat)
        exptime (float): Exposure time of the image

    """

    def __init__(self, image, ivar=None, rn2img=None, bpm=None, state=None, binning=None):

        maskimage.ImageMask.__init__(self, bpm)

        # Required parameters
        self.image = image

        # Optional Attributes
        self.ivar = ivar
        self.rn2img = rn2img
        self.state = state
        self.binning = binning

    def show(self):
        """
        Simple show method
        """
        if self.image is None:
            msgs.warn("No image to show!")
            return
        ginga.show_image(self.image, chname='image')


