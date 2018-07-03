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

    def get_datasec(self, filename, det, settings_det):

        datasec, oscansec, naxis0, naxis1 = [], [], 0, 0
        for i in range(settings_det['numamplifiers']):
            sdatasec = "datasec{0:02d}".format(i+1)
            datasec.append(settings_det[sdatasec])
            soscansec = "oscansec{0:02d}".format(i+1)
            oscansec.append(settings_det[soscansec])

        # Read the image for the shape (just in case)
        temp, _ = self.load_raw_img_head(filename, det=det, dataext=settings_det['dataext01'],
                                    disp_dir=settings_det['dispaxis'])
        # Need naxis0, naxis1 too
        naxis0 = temp.shape[0]
        naxis1 = temp.shape[1]

        # Return
        return datasec, oscansec, naxis0, naxis1


