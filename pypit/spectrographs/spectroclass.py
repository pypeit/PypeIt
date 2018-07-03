""" Module to define the Spectrograph class
"""

import numpy as np

from astropy.io import fits

class Spectrograph(object):

    def __init__(self):

        self.spectrograph = 'generic'

    def load_raw_frame(self, raw_file, dataext=None, disp_dir=0, det=None):

        # Load the raw image
        raw_img, head0 = self.load_raw_img_head(raw_file, dataext=dataext, det=det)

        # Turn to float
        img = raw_img.astype(np.float)
        # Transpose?
        if disp_dir == 1:
            img = img.T
        # Return
        return img, head0

    def load_raw_img_head(self, raw_file, dataext=None, **null_kwargs):

        hdulist = fits.open(raw_file)
        raw_img = hdulist[dataext].data
        head0 = hdulist[0].header
        # Return
        return raw_img, head0
