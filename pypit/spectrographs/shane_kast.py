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
from .. import telescopes

from pypit import ardebug as debugger


class ShaneKastSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Shane/Kast specific code
    """
    def __init__(self):
        # Get it started
        super(ShaneKastSpectrograph, self).__init__()
        self.spectrograph = 'shane_kast'
        self.telescope = telescopes.ShaneTelescopePar()
        self.timeunit = 's'


class ShaneKastBlueSpectrograph(ShaneKastSpectrograph):
    """
    Child to handle Shane/Kast blue specific code
    """
    def __init__(self):
        # Get it started
        super(ShaneKastBlueSpectrograph, self).__init__()
        self.spectrograph = 'shane_kast_blue'
        self.camera = 'KASTb'
        self.detector = [
                # Detector 1
                DetectorPar(dataext         = 0,
                            dispaxis        = 1,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.43,
                            darkcurr        = 0.0,
                            saturation      = 65535.,
                            nonlinear       = 0.76,
                            numamplifiers   = 2,
                            gain            = [1.2, 1.2],
                            ronoise         = [3.7, 3.7],
                            datasec         = [ '[1:1024,:]', '[1025:2048,:]'],
                            oscansec        = [ '[2050:2080,:]', '[2081:2111,:]'],
                            suffix          = '_blue'
                            )
            ]
        # Uses timeunit from parent class
        # Uses default primary_hdrext
        self.sky_file = 'sky_kastb_600.fits'
        

    def setup_arcparam(self, arcparam, disperser=None, **null_kwargs):
        """
        Setup the arc parameters

        Args:
            arcparam: dict
            disperser: str, REQUIRED
            **null_kwargs:
              Captured and never used

        Returns:
            arcparam is modified in place

        """
        arcparam['lamps'] = ['CdI','HgI','HeI']
        if disperser == '600/4310':
            arcparam['disp']=1.02
            arcparam['b1']=6.88935788e-04
            arcparam['b2']=-2.38634231e-08
            arcparam['wvmnx'][1] = 6000.
            arcparam['wv_cen'] = 4250.
        else:
            msgs.error('Not ready for this disperser {:s}!'.format(disperser))


class ShaneKastRedSpectrograph(ShaneKastSpectrograph):
    """
    Child to handle Shane/Kast red specific code
    """
    def __init__(self):
        # Get it started
        super(ShaneKastRedSpectrograph, self).__init__()
        self.spectrograph = 'shane_kast_red'
        self.camera = 'KASTr'
        self.detector = [
                # Detector 1
                DetectorPar(dataext         = 0,
                            dispaxis        = 0,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.43,
                            darkcurr        = 0.0,
                            saturation      = 65535.,
                            nonlinear       = 0.76,
                            numamplifiers   = 2,
                            gain            = [1.9, 1.9],
                            ronoise         = [3.8, 3.8],
                            datasec         = ['[2:511,:]', '[513:525,:]'],
                            oscansec        = ['[527:625,:]', '[627:725,:]'],
                            suffix          = '_red'
                            )
            ]
        # Uses timeunit from parent class
        # Uses default primary_hdrext
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

        """
        arcparam['lamps'] = ['NeI','HgI','HeI','ArI']
        if disperser == '600/7500':
            arcparam['disp']=1.30
            arcparam['b1']= 1./arcparam['disp']/msarc_shape[0] / binspectral
            arcparam['wvmnx'][0] = 5000.
            arcparam['n_first']=2 # Should be able to lock on
        elif disperser == '1200/5000':
            arcparam['disp']=0.63
            arcparam['b1']= 1./arcparam['disp']/msarc_shape[0] / binspectral
            arcparam['wvmnx'][0] = 5000.
            arcparam['n_first']=2 # Should be able to lock on
            arcparam['wv_cen'] = 6600.
        else:
            msgs.error('Not ready for this disperser {:s}!'.format(disperser))

class ShaneKastRedRetSpectrograph(ShaneKastSpectrograph):
    """
    Child to handle Shane/Kast red specific code
    """
    def __init__(self):
        # Get it started
        super(ShaneKastRedRetSpectrograph, self).__init__()
        self.spectrograph = 'shane_kast_red_ret'
        # WARNING: This is not unique wrt ShaneKastRed...
        self.camera = 'KASTr'
        self.detector = [
                # Detector 1
                DetectorPar(dataext         = 0,
                            dispaxis        = 1,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.774,
                            darkcurr        = 0.0,
                            saturation      = 65535.,
                            nonlinear       = 0.76,
                            numamplifiers   = 1,
                            gain            = 3.0,
                            ronoise         = 12.5,
                            oscansec        = '[1203:1232,:]',
                            suffix          = '_red'
                            )
            ]
        # Uses timeunit from parent class
        # Uses default primary_hdrext
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

        """
        arcparam['lamps'] = ['NeI','HgI','HeI','ArI']
        if disperser == '600/7500':
            arcparam['disp']=2.35
            arcparam['b1']= 1./arcparam['disp']/msarc_shape[0] / binspectral
            arcparam['wvmnx'][0] = 5000.
            arcparam['n_first']=2 # Should be able to lock on
        elif disperser == '1200/5000':
            arcparam['disp']=1.17
            arcparam['b1']= 1./arcparam['disp']/msarc_shape[0] / binspectral
            arcparam['wvmnx'][0] = 5000.
            arcparam['n_first']=2 # Should be able to lock on
        else:
            msgs.error('Not ready for this disperser {:s}!'.format(disperser))

