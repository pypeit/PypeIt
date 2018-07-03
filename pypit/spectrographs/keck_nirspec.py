""" Module for Keck/NIRSPEC specific codes
"""
from __future__ import absolute_import, division, print_function

try:
    basestring
except NameError:  # For Python 3
    basestring = str

import glob

import numpy as np
from astropy.io import fits

from pypit import msgs
from pypit import ardebug as debugger
from. import spectroclass

class KeckNIRSPECSpectrograph(spectroclass.Spectrograph):
    """
    Child to handle Keck/LRIS specific code
    """

    def __init__(self):

        # Get it started
        spectroclass.Spectrograph.__init__(self)
        self.spectrograph = 'keck_nirspec'

    def bpm(self, shape=None, **null_kwargs):
        """ Generate a BPM

        Parameters
        ----------
        shape : tuple, REQUIRED

        Returns
        -------
        badpix : ndarray

        """
        # Edges of the detector are junk
        msgs.info("Custom bad pixel mask for NIRSPEC")
        self.bpm = np.zeros((self.shape[0], self.shape[1]))
        self.bpm[:, :20] = 1.
        self.bpm[:, 1000:] = 1.

