""" Object to hold + process a single image"""

from pypeit import msgs
from pypeit import ginga

import numpy as np

from IPython import embed

class PypeItImage(object):

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
        self.filename = None         # Filename of the image
        self.exptime = None          # Required to generate variance image

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

    def load_rawimage(self, filename):
        """
        Load a raw image from disk using the Spectrograph method load_raw_frame()

        Returns:
            np.ndarray, fits.Header:

        """
        self.filename = filename
        self.image, self.head0 \
            = self.spectrograph.load_raw_frame(filename, det=self.det)
        # Shape
        self.orig_shape = self.image.shape
        # Exposure time
        self.exptime = self.spectrograph.get_meta_value(filename, 'exptime')
        # Binning
        self.binning = self.spectrograph.get_meta_value(filename, 'binning')
        embed(header='58 of pypeitimage')
        if self.spectrograph.detector[self.det-1]['specaxis'] == 1:
            self.binning_raw = (',').join(self.binning.split(',')[::-1])
        else:
            self.binning_raw = self.binning
        # Return
        return self.image, self.head0


    def show(self):
        """
        Simple show method
        """
        ginga.show_image(self.image, chname='image')

