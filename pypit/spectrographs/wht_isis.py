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
from ..par.pypitpar import DetectorPar
from . import spectrograph

from pypit import ardebug as debugger

class WhtIsisSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Shane/Kast specific code
    """

    def __init__(self):
        super(WhtIsisSpectrograph, self).__init__()
        self.spectrograph = 'NULL'

class WhtIsisBlueSpectrograph(WhtIsisSpectrograph):
    """
    Child to handle WHT/ISIS blue specific code
    """
    def __init__(self):
        # Get it started
        super(WhtIsisBlueSpectrograph, self).__init__()
        self.spectrograph = 'wht_isis_blue'
        self.detector = [
                # Detector 1
                DetectorPar(dataext=1, dispaxis=0, xgap=0., ygap=0., ysize=1., platescale=0.225,
                            darkcurr=0.0, saturation=65535., nonlinear=0.76, numamplifiers=1,
                            gain=1.2, ronoise=5.0, datasec='[:,2:4030]', suffix='_blue')
            ]
        

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

