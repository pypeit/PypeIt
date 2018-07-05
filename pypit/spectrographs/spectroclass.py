""" Module to define the Spectrograph class
"""
from __future__ import absolute_import, division, print_function

import numpy as np

from astropy.io import fits

from abc import ABCMeta

class Spectrograph(object):
    """
    Generic class for spectrograph-specific codes
    """
    __metaclass__ = ABCMeta

    def __init__(self):

        self.spectrograph = 'generic'

    def load_raw_frame(self, raw_file, disp_dir, dataext=None, det=None):
        """
        Load the image (converted to np.float) and primary header of the input file

        The image is transposed, as needed, so that the spectral dimension
        runs along the columns

        Args:
            raw_file:  str, filename
            disp_dir: int
              This is an example where JFH would rant that this is known
              given the instrument so should not be a settting.  He is
              probably right.
            dataext: int, optional
              Extension in the FITS list for the data
            det: int, optional
              Desired detector

        Returns:
            img: ndarray
              Converted to np.float and transpsed if necessary
            head0: Header

        """

        # Load the raw image
        raw_img, head0 = self.load_raw_img_head(raw_file, dataext=dataext, det=det)

        # Turn to float
        img = raw_img.astype(np.float)
        # Transpose?
        if disp_dir == 1:
            img = img.T
        # Return
        return img, head0

    def load_raw_img_head(self, raw_file, dataext, **null_kwargs):
        """
        Generic raw image reader

        Args:
            raw_file: str
            dataext: int
            **null_kwargs:
              Captured and never used

        Returns:
            raw_img: ndarray
              Raw image;  likely unsigned int
            head0: Header

        """
        # Open and go
        hdulist = fits.open(raw_file)
        from IPython import embed
        raw_img = hdulist[dataext].data
        head0 = hdulist[0].header
        # Return
        return raw_img, head0

    def get_datasec(self, filename, det, settings_det):
        """
        Load up the datasec and oscansec and also naxis0 and naxis1

        Args:
            filename: str
              data filename
            det: int
              Detector specification
            settings_det: ParSet
              numamplifiers
              dispaxis
              dataext01

        Returns:
            datasec: list
            oscansec: list
            naxis0: int
            naxis1: int
        """

        datasec, oscansec, naxis0, naxis1 = [], [], 0, 0
        for i in range(settings_det['numamplifiers']):
            sdatasec = "datasec{0:02d}".format(i+1)
            datasec.append(settings_det[sdatasec])
            soscansec = "oscansec{0:02d}".format(i+1)
            oscansec.append(settings_det[soscansec])

        # Read the image for the shape (just in case)
        temp, _ = self.load_raw_frame(filename, settings_det['dispaxis'],
                                      det=det, dataext=settings_det['dataext'])
        # Need naxis0, naxis1 too
        naxis0 = temp.shape[0]
        naxis1 = temp.shape[1]

        # Return
        return datasec, oscansec, naxis0, naxis1

    def bpm(self, shape=None, **null_kwargs):
        """
        Generate a generic (empty) BPM

        Args:
            shape: tuple, REQUIRED
            **null_kwargs:

        Returns:
            bpm: ndarray, int
              0=not masked; 1=masked

        """
        bpm = np.zeros((shape[0], shape[1]), dtype=int)
        #
        return bpm

    def setup_arcparam(self, **null_kwargs):
        modify_dict = None
        return modify_dict
