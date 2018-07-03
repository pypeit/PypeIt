""" Module for Shane/Kast specific codes
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

class WhtIsisSpectrograph(spectroclass.Spectrograph):
    """
    Child to handle Shane/Kast specific code
    """

    def __init__(self):

        # Get it started
        spectroclass.Spectrograph.__init__(self)
        self.spectrograph = 'NULL'

class WhtIsisBlueSpectrograph(WhtIsisSpectrograph):
    """
    Child to handle WHT/ISIS blue specific code
    """
    def __init__(self):

        # Get it started
        spectroclass.Spectrograph.__init__(self)
        self.spectrograph = 'wht_isis_blue'

    def setup_arcparam(self, arcparam, disperser=None, fitstbl=None,
                       arc_idx=None, msarc_shape=None, **null_kwargs):
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
        modify_dict = dict(NeI={'min_wave': 3000.,'min_intensity': 299,
                                'min_Aki': 0.},ArI={'min_intensity': 399.})
        arcparam['lamps']=['CuI','NeI','ArI']
        if fitstbl["dichroic"][arc_idx].strip() == '5300':
            arcparam['wvmnx'][1] = 6000.
        else:
            msgs.error('Not ready for this dichroic {:s}!'.format(disperser))
        if disperser == 'R300B':
            arcparam['n_first']=1  #
            arcparam['disp']=0.80  # Ang per pixel (unbinned)
            arcparam['b1']= 1./arcparam['disp']/msarc_shape[0]
        else:
            msgs.error('Not ready for this disperser {:s}!'.format(disperser))
        #
        return modify_dict

