""" Lightweight object to load and hold the RawImage and a few additional bits and pieces """

import numpy as np

from IPython import embed


class RawImage(object):
    """
    Class to process a raw image

    Args:
        filename (:obj:`str` or None):
            Filename
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            Spectrograph used to take the data.
        det (:obj:`int`, optional):
            The 1-indexed detector number to process.

    Attributes:
        steps (dict):
            Dict describing the steps performed on the image
        _bpm (np.ndarray):
            Holds the bad pixel mask once loaded
        rawdatasec_img (np.ndarray):
            Holds the datasec_img which specifies the amp for each pixel in the
            raw image
        oscansec_img (np.ndarray):
            Holds the oscansec_img once loaded.  This is in the raw frame
        hdu (fits.HDUList):
            HDUList of the file
    """
    def __init__(self, filename, spectrograph, det):

        # Init me
        self.filename = filename
        self.spectrograph = spectrograph
        self.det = det

        # Load
        self.raw_image, self.hdu, self.exptime, self.rawdatasec_img, self.oscansec_img = self.spectrograph.get_rawimage(
            self.filename, self.det)

    @property
    def amps(self):
        """
        Return a list of the amplifier indices, 1-indexed

        Returns:
            list
        """
        return np.unique(self.rawdatasec_img[self.rawdatasec_img > 0]).tolist()

    def __repr__(self):
        return ('<{:s}: file={}> spectrograph={}'.format(
            self.__class__.__name__, self.filename, self.spectrograph.spectrograph))

