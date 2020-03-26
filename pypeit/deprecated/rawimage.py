""" Lightweight object to load and hold the RawImage and several other key
items for its processing."""

import numpy as np

from IPython import embed

# TODO -- Turn this into a DataContainer??
class OldRawImage(object):
    """
    Class to load and hold a raw image

    Args:
        filename (:obj:`str` or None):
            Filename
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            Spectrograph used to take the data.
        det (:obj:`int`):
            The 1-indexed detector number to process.

    Attributes:
        raw_image (np.ndarray):
            Raw image as a numpy array
        rawdatasec_img (np.ndarray):
            Holds the datasec_img which specifies the amp for each pixel in the
            raw image
        oscansec_img (np.ndarray):
            Holds the oscansec_img once loaded.  This is in the raw frame
        hdu (fits.HDUList):
            HDUList of the file
        exptime (float):
            Exposure time
        detector(:class:`pypeit.images.detector_container.Detector`):
        binning (:obj:`str`):
    """
    def __init__(self, filename, spectrograph, det):

        # Init me
        self.filename = filename
        self.spectrograph = spectrograph
        self.det = det

        # Load the raw image and the other items of interest
        self.detector, self.raw_image, self.hdu, self.exptime, self.rawdatasec_img, self.oscansec_img = self.spectrograph.get_rawimage(
            self.filename, self.det)

    def __repr__(self):
        return ('<{:s}: file={}> spectrograph={}'.format(
            self.__class__.__name__, self.filename, self.spectrograph.spectrograph))

