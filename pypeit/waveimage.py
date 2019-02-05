# Module for guiding construction of the Wavelength Image
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np
import os

#from importlib import reload

from pypeit import msgs
from pypeit import utils
from pypeit import masterframe
from pypeit import ginga
from pypeit.core import pixels

from pypeit import debugger


class WaveImage(masterframe.MasterFrame):
    """Class to generate the Wavelength Image

    Parameters
    ----------
    tilts : ndarray
      Tilt image
    wv_calib : dict
      wavelength solution dictionary
    settings : dict
    master_key : str
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

    def __init__(self, tslits_dict, tilts, wv_calib, spectrograph, binning = None, master_key=None, master_dir=None,
                 reuse_masters=False, maskslits=None):

        # MasterFrame
        masterframe.MasterFrame.__init__(self, self.frametype, master_key,
                                         master_dir=master_dir, reuse_masters=reuse_masters)

        # Required parameters (but can be None)
        self.tslits_dict = tslits_dict
        self.tilts = tilts
        self.wv_calib = wv_calib
        self.spectrograph = spectrograph
        self.binning = binning
        self.slitmask = pixels.tslits2mask(self.tslits_dict)
        self.par = wv_calib['par']

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
        nspec =self.slitmask.shape[0]

        # Error checking on the wv_calib
        #if (nspec-1) != int(self.wv_calib[str(0)]['fmax']):
        #    msgs.error('Your wavelength fits used inconsistent normalization. Something is wrong!')

        # Ff this is echelle print out a status message and do some error checking
        if self.par['echelle']:
            msgs.info('Evaluating 2-d wavelength solution for echelle....')
            if len(self.wv_calib['fit2d']['orders']) != len(ok_slits):
                msgs.error('wv_calib and ok_slits do not line up. Something is very wrong!')

        # Unpack some 2-d fit parameters if this is echelle
        for slit in ok_slits:
            thismask = (self.slitmask == slit)
            if self.par['echelle']:
                order = self.spectrograph.slit2order(slit)
                # evaluate solution
                tmpwv = utils.func_val(self.wv_calib['fit2d']['coeffs'], self.tilts[thismask], self.wv_calib['fit2d']['func2d'],
                                       x2=np.ones_like(self.tilts[thismask])*order,
                                       minx=self.wv_calib['fit2d']['min_spec'], maxx=self.wv_calib['fit2d']['max_spec'],
                                       minx2=self.wv_calib['fit2d']['min_order'], maxx2=self.wv_calib['fit2d']['max_order'])/order
            else:
                iwv_calib = self.wv_calib[str(slit)]
                tmpwv = utils.func_val(iwv_calib['fitc'], self.tilts[thismask], iwv_calib['function'],
                                       minx=iwv_calib['fmin'], maxx=iwv_calib['fmax'])
            self.wave[thismask] = tmpwv
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


