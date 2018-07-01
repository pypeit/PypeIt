""" Module to define the Spectrograph class
"""

import numpy as np

from astropy.io import fits

class Spectrograph(object):

    def __init__(self):

        self.spectrograph = 'generic'

    def load_raw_frame(self, raw_file, dataext=None, disp_dir=0, det=None):
        hdulist = fits.open(raw_file)
        raw_img = hdulist[dataext].data
        head0 = hdulist[0].header

        # Return
        return self.fuss_with_raw(raw_img, disp_dir=disp_dir), head0

    def fuss_with_raw(self, raw_img, disp_dir=0):
        # Turn to float
        temp = raw_img.astype(np.float)
        # Transpose?
        if disp_dir == 1:
            temp = temp.T
        # Return
        return temp


