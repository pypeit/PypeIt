""" Module for Keck/NIRES specific codes
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

class KeckNIRESpectrograph(spectroclass.Spectrograph):
    """
    Child to handle Keck/NIRES specific code
    """

    def __init__(self):

        # Get it started
        spectroclass.Spectrograph.__init__(self)
        self.spectrograph = 'keck_nires'

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
        msgs.info("Custom bad pixel mask for NIRES")
        self.bpm = np.zeros((self.shape[0], self.shape[1]))
        #self.bpm[:, :20] = 1.
        #self.bpm[:, 1000:] = 1.

    def setup_arcparam(self, arcparam, fitstbl=None, arc_idx=None,
                       msarc_shape=None, **null_kwargs):
        """

        Args:
            arcparam:
            disperser:
            fitstbl:
            arc_idx:
            msarc_shape:
            **null_kwargs:

        Returns:

        """
        arcparam = dict(llist='',
                        disp=2.,  # Ang/unbinned pixel
                        b1=0.,  # Pixel fit term (binning independent)
                        b2=0.,  # Pixel fit term
                        lamps=['OH_triplespec'],  # Line lamps on
                        wv_cen=0.,  # Estimate of central wavelength
                        wvmnx=[9000., 25000.],  # Guess at wavelength range
                        disp_toler=0.1,  # 10% tolerance
                        match_toler=3.,  # Matching tolerance (pixels)
                        min_ampl=1000.,  # Minimum amplitude
                        func='legendre',  # Function for fitting
                        n_first=1,  # Order of polynomial for first fit
                        n_final=3,  # Order of polynomial for final fit
                        nsig_rej=2.,  # Number of sigma for rejection
                        nsig_rej_final=2.0,  # Number of sigma for rejection (final fit)
                        Nstrong=13)  # Number of lines for auto-analysis

