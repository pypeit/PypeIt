""" Module for TNG/Dolores
"""
from __future__ import absolute_import, division, print_function


import numpy as np

from pypeit import msgs
from pypeit.par.pypeitpar import DetectorPar
from pypeit.spectrographs import spectrograph
from pypeit import telescopes

from pypeit import debugger

class TngDoloresSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Shane/Kast specific code
    """

    def __init__(self):
        super(TngDoloresSpectrograph, self).__init__()
        self.spectrograph = 'tng_dolores'
        self.telescope = telescopes.TNGTelescopePar()
        self.camera = 'DOLORES'
        self.detector = [
                # Detector 1
                DetectorPar(dataext         = 0,
                            dispaxis        = 1,
                            dispflip=False,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.252,
                            darkcurr        = 0.0,
                            saturation      = 65535.,
                            nonlinear       = 0.76,
                            numamplifiers   = 1,
                            gain            = 0.97,
                            ronoise         = 9.0,
                            datasec         = '[51:,1:2045]',
                            oscansec        = '[51:,2054:]',
                            suffix          = '_lrr'
                            )
            ]
        # Uses default primary_hdrext
        self.timeunit = 'isot'
        # self.sky_file = ?

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
        arcparam['nonlinear_counts'] = self.detector[0]['nonlinear']*self.detector[0]['saturation']

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

