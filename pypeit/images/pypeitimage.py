""" Object to hold + process a single image"""

from pypeit import msgs
from pypeit import ginga

import numpy as np


class PypeItImage(object):

    def __init__(self, spectrograph, det):

        # Required parameters
        self.spectrograph = spectrograph
        self.det = det

        # Attributes
        self.image = None
        self.header = None           # Image header
        self.orig_shape = None       # Shape of the image when loaded
        self.filename = None       # Shape of the image when loaded

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

    def load_image(self, filename):
        """
        Load the image from disk using the Spectrograph method load_raw_frame()

        Returns:
            np.ndarray, fits.Header:

        """
        self.filename = filename
        self.image, self.header \
            = self.spectrograph.load_raw_frame(filename, det=self.det)
        # Shape
        self.orig_shape = self.image.shape
        #
        return self.image, self.header


    def show(self):
        viewer, ch = ginga.show_image(self.image, chname='image')

