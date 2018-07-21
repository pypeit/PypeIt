# Module for guiding construction of the Wavelength Image
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np
import os

#from importlib import reload

from pypeit import msgs
from pypeit import arutils
from pypeit import masterframe
from pypeit import ginga

from pypeit import ardebug as debugger

class WaveImage(masterframe.MasterFrame):
    """Class to generate the Wavelength Image

    Parameters
    ----------
    tilts : ndarray
      Tilt image
    wv_calib : dict
      1D wavelength solutions
    settings : dict
    setup : str
    maskslits : ndarray
      True = skip this slit
    slitpix : ndarray
      Specifies locations of pixels in the slits

    Attributes
    ----------
    frametype : str
      Hard-coded to 'wave'
    wave : ndarray
      Wavelength image

    steps : list
      List of the processing steps performed
    """
    # Frametype is a class attribute
    frametype = 'wave'

    def __init__(self, slitpix, tilts, wv_calib, setup=None, directory_path=None, mode=None, 
                 maskslits=None):

        # MasterFrame
        masterframe.MasterFrame.__init__(self, self.frametype, setup,
                                         directory_path=directory_path, mode=mode)

        # Required parameters (but can be None)
        self.slitpix = slitpix
        self.tilts = tilts
        self.wv_calib = wv_calib

        # Optional parameters
        self.maskslits = maskslits

        # Attributes
        self.steps = [] # List to hold ouptut from inspect about what module create the image?

        # Main output
        self.wave = None

    def _build_wave(self):
        """
        Main algorithm to build the wavelength image

        Returns
        -------
        self.wave : ndarray
          Wavelength image

        """
        # TODO: self.maskslits cannot be None
        # Loop on slits
        ok_slits = np.where(~self.maskslits)[0]
        self.wave = np.zeros_like(self.tilts)
        for slit in ok_slits:
            iwv_calib = self.wv_calib[str(slit)]
            tmpwv = arutils.func_val(iwv_calib['fitc'], self.tilts, iwv_calib['function'],
                                     minv=iwv_calib['fmin'], maxv=iwv_calib['fmax'])
            word = np.where(self.slitpix == slit+1)
            self.wave[word] = tmpwv[word]
        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.wave

    def show(self, item='wave'):
        """
        Show the image

        Parameters
        ----------
        item : str, optional

        Returns
        -------

        """
        if item == 'wave':
            if self.wave is not None:
                ginga.show_image(self.wave)
        else:
            msgs.warn("Not able to show this type of image")

    def __repr__(self):
        # Generate sets string
        txt = '<{:s}: >'.format(self.__class__.__name__)
        return txt


