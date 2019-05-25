""" Object to hold + process a single image"""

from pypeit import msgs

from pypeit.core import procimg
from astropy.io import fits


class PypeItImage(object):

    def __init__(self, filename, spectrograph, det):

        # Required parameters
        self.spectrograph = spectrograph
        self.det = det
        self.filename = filename

        # Attributes
        self.image = None
        self.header = None           # Image header

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
    def datasec_img(self):
        dimg = self.spectrograph.get_datasec_img(self.filename, self.det)
        return dimg

    @property
    def oscansec_img(self):
        oimg = self.spectrograph.get_oscansec_img(self.filename, self.det)
        return oimg

    def load(self):
        """
        Load the image from disk using the Spectrograph method load_raw_frame()

        Returns:
            np.ndarray, fits.Header:

        """
        self.image, self.header \
            = self.spectrograph.load_raw_frame(self.filename, det=self.det)
        #
        return self.image, self.header


