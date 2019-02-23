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

    Args:
        tslits_dict (dict): dict from TraceSlits class (e.g. slitpix)
        tilts (np.ndarray): Tilt image
        wv_calib (dict): wavelength solution dictionary
            Parameters are read from wv_calib['par']
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            The `Spectrograph` instance that sets the
            instrument used to take the observations.  Used to set
            :attr:`spectrograph`.
        master_key (:obj:`str`, optional):
            The string identifier for the instrument configuration.  See
            :class:`pypeit.masterframe.MasterFrame`.
        master_dir (str, optional): Path to master frames
        maskslits (ndarray, optional): True = skip this slit
        reuse_masters (bool, optional):  Load from disk if possible

    Attributes:
        frametype : str
          Hard-coded to 'wave'
        wave (ndarray): Wavelength image
        steps (list): List of the processing steps performed

    """
    # Frametype is a class attribute
    frametype = 'wave'

    def __init__(self, tslits_dict, tilts, wv_calib, spectrograph, maskslits,
                 master_key=None, master_dir=None, reuse_masters=False):

        # MasterFrame
        masterframe.MasterFrame.__init__(self, self.frametype, master_key,
                                         master_dir, reuse_masters=reuse_masters)

        # Required parameters (but can be None)
        self.tslits_dict = tslits_dict
        self.tilts = tilts
        self.wv_calib = wv_calib
        self.spectrograph = spectrograph
        self.slitmask = pixels.tslits2mask(self.tslits_dict) if tslits_dict is not None else None
        self.par = wv_calib['par'] if wv_calib is not None else None

        # Optional parameters
        self.maskslits = maskslits

        # Attributes
        self.steps = [] # List to hold ouptut from inspect about what module create the image?

        # Main output
        self.wave = None

    def _build_wave(self):
        """
        Main algorithm to build the wavelength image

        Returns:
            ndarray: Wavelength image

        """
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

        Args:
            item (str, optional):

        Returns:

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




def load_waveimage(filename):
    """
    Utility function which enables one to load the waveimage from a master file in one line of code without
    instantiating the class.

    Args:
        filename (str): Master file name

    Returns:
        dict:  The trace slits dict

    """

    waveImage = WaveImage(None, None, None, None, None)
    waveimage, _ = waveImage.load_master(filename)

    return waveimage
