""" Module to define the Spectrograph class
"""
from __future__ import absolute_import, division, print_function

import numpy as np

from astropy.io import fits

from abc import ABCMeta

from ..par.pypitpar import DetectorPar

class Spectrograph(object):
    """
    Generic class for spectrograph-specific codes
    """
    __metaclass__ = ABCMeta

    def __init__(self):
        self.spectrograph = 'generic'
        self.detector = None

    def _check_detector(self):
        # Check the detector
        if self.detector is None:
            raise ValueError('Must first define spectrograph detector parameters!')
        for d in self.detector:
            if not isinstance(d, DetectorPar):
                raise TypeError('Detector parameters must be specified using DetectorPar.')

    def load_raw_frame(self, raw_file, det=None):
        """
        Load the image (converted to np.float) and primary header of the input file

        The image is transposed, as needed, so that the spectral dimension
        runs along the columns

        Args:
            raw_file:  str, filename
            dataext: int, optional
              Extension in the FITS list for the data
            det: int, optional
              Desired detector

        Returns:
            img: ndarray
              Converted to np.float and transposed if necessary
            head0: Header

        """
        self._check_detector()
        _det = 1 if det is None else det

        # Load the raw image
        raw_img, head0 = self.load_raw_img_head(raw_file, dataext=self.detector[_det-1]['dataext'],
                                                det=_det)

        # Turn to float
        img = raw_img.astype(np.float)
        # Transpose?
        if self.detector[_det-1]['dispaxis'] == 1:
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
        raw_img = hdulist[dataext].data
        head0 = hdulist[0].header
        # Return
        return raw_img, head0

    def get_datasec(self, filename, det):
        """
        Load up the datasec and oscansec and also naxis0 and naxis1

        Args:
            filename (str):
                data filename
            det (int):
                Detector number

        Returns:
            datasec: list
            oscansec: list
            naxis0: int
            naxis1: int
        """
        self._check_detector()
        # TODO: This seems like a lot of effort to get the size of the
        # image.
        # Read the image for the shape (just in case)
        temp, _ = self.load_raw_frame(filename, self.detector['dispaxis'], det=det,
                                      dataext=self.detector['dataext'])
        return (self.detector['datasec'], self.detector['oscansec']) + temp.shape

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
        return np.zeros((shape[0], shape[1]), dtype=int)

    def setup_arcparam(self, **null_kwargs):
        return None

    @property
    def ndet(self):
        """Return the number of detectors."""
        return 0 if self.detector is None else len(self.detector)

