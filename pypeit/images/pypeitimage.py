""" Object to hold + process a single image"""

from pypeit import msgs
from pypeit import ginga

import numpy as np

from IPython import embed

class PypeItImage(object):
    """
    Class to hold a single image from a single detector in PypeIt

    Args:
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            Spectrograph used to take the data.
        det (:obj:`int`, optional):
            The 1-indexed detector number to process.

    Attributes:
        image (np.ndarray):
        head0 (astropy.io.fits.Header):
        orig_shape (tuple):
        binning_raw (tuple):  Binning in the raw image orientaion (NAXIS1, NAXIS2)
        binning (tuple): Binning the PypeIt orientation (spec, spat)
        exptime (float): Exposure time of the image

    """

    def __init__(self, spectrograph, det):

        # Required parameters
        self.spectrograph = spectrograph
        self.det = det

        # Attributes
        self.image = None
        self.head0 = None           # Image header
        self.orig_shape = None       # Shape of the image when loaded
        self.binning_raw = None     # Binning in the raw image orientation;  e.g. bin_1, bin_2 (for NAXIS1, NAXIS2)
        self.binning = None          # Binning in PypeIt orientation (spec, spat)
        self.exptime = None          # Required to generate variance image

    @property
    def shape(self):
        return () if self.image is None else self.image.shape

    def show(self):
        """
        Show the image in a ginga viewer.
        """
        if self.image is None:
            # TODO: This should fault.
            msgs.warn("No image to show!")
            return
        ginga.show_image(self.image, chname='image')

    def __repr__(self):
        txt = '<{:s}:'.format(self.__class__.__name__)
        if self.filename is not None:
            txt += ' file={}'.format(self.filename)
        txt += '>'

        return txt

