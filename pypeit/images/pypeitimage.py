""" Object to hold + process a single image"""

from pypeit import msgs

import numpy as np


class PypeItImage(object):

    def __init__(self, filename, spectrograph, det):

        # Required parameters
        self.spectrograph = spectrograph
        self.det = det
        self.filename = filename

        # Attributes
        self.image = None
        self.header = None           # Image header
        self.orig_shape = None       # Shape of the image when loaded

    @property
    def bpm(self):
        """
        Generate and return the bad pixel mask for this image
        Warning:  BPM masks are for processed (e.g. trimmed, rotated) images only!

        Returns:
            np.ndarray:  Bad pixel mask with a bad pixel = 1

        """
        bpm = self.spectrograph.bpm(shape=self.image.shape,
                                        filename=self.filename,
                                        det=self.det)
        return bpm

    @property
    def amps(self):
        """
        Return a list of the amplifier indices, 1-indexed

        Returns:
            list
        """
        amps = np.unique(self.datasec_img[self.datasec_img > 0]).tolist()
        # Return
        return amps

    def load(self):
        """
        Load the image from disk using the Spectrograph method load_raw_frame()

        Returns:
            np.ndarray, fits.Header:

        """
        self.image, self.header \
            = self.spectrograph.load_raw_frame(self.filename, det=self.det)
        # Shape
        self.orig_shape = self.image.shape
        #
        return self.image, self.header


