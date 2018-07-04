""" Module for TNG/Dolores
"""
from __future__ import absolute_import, division, print_function

try:
    basestring
except NameError:  # For Python 3
    basestring = str

import glob

import numpy as np

from pypit import msgs
from pypit import ardebug as debugger
from pypit.spectrographs import spectroclass

class TngDoloresSpectrograph(spectroclass.Spectrograph):
    """
    Child to handle Shane/Kast specific code
    """

    def __init__(self):

        # Get it started
        spectroclass.Spectrograph.__init__(self)
        self.spectrograph = 'NULL'

    def setup_arcparam(self, arcparam, disperser=None, msarc_shape=None,
                       binspectral=None, **null_kwargs):
        """
        Setup the arc parameters

        Args:
            arcparam: dict
            disperser: str, REQUIRED
            **null_kwargs:
              Captured and never used

        Returns:
            arcparam is modified in place
            modify_dict: dict

        """
        arcparam['lamps'] = ['NeI', 'HgI']
        if disperser == 'LR-R':
            arcparam['n_first'] = 2  # Too much curvature for 1st order
            arcparam['disp'] = 2.61  # Ang per pixel (unbinned)
            arcparam['disp_toler'] = 0.1  # Ang per pixel (unbinned)
            arcparam['wvmnx'][0] = 4470.0
            arcparam['wvmnx'][1] = 10073.0
            arcparam['wv_cen'] = 7400.
            arcparam['b1'] = 1. / arcparam['disp'] / msarc_shape[0] / binspectral
        else:
            msgs.error('Not ready for this disperser {:s}!'.format(disperser))

