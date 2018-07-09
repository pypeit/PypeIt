# Module for generating the BPM image
from __future__ import absolute_import, division, print_function

import numpy as np
import os

from pypit import msgs
from pypit.core import arprocimg

from .spectrographs.spectrograph import Spectrograph
from .spectrographs.util import load_spectrograph

from pypit import ardebug as debugger

try:
    basestring
except NameError:
    basestring = str

# Does not need to be global, but I prefer it
frametype = 'bpm'

class BPMImage(object):
    """
    This class is primarily designed to generate an Bad Pixel Image
      The master() method will return the image

    Should provide both shape and filename to ensure that the
    spectograph can construct the bpm, if reduce_badpix is False.

    There are several ways to build a BPM:
       1. keck_lris_red
          Provide binning and detector number to generate this instrument specific BPM
       2. keck_deimos
          Provide detector number to generate this instrument specific BPM
       3. From a bias image (not well tested)
          Set reduce_badpix='bias'
          Provide msbias image
       4. Dummy image from shape

    Parameters
    ----------
    spectrograph : str (optional)
       Used to specify properties of the detector (for processing)
       Attempt to set with settings['run']['spectrograph'] if not input
    settings : dict (optional)
      Settings for trace slits
    det : int, optional
      Detector index
      Required for LRISr and DEIMOS
    binning : Table (optional)
        Could remove if we put binning in settings['detector']
    shape : tuple
      Image shape;  used to construct a dummy BPM if all else fails
    msbias : ndarray, optional
      Used to construct the BPM if reduce_badpix=='bias'
    reduce_badpix : str, optional
      'bias' -- Build from bias images

    Attributes
    ----------
    bpm : ndarray

    """
    # Keep order same as processimages (or else!)
    def __init__(self, spectrograph, shape=None, filename=None, det=None, msbias=None, trim=True):

        # Spectrograph is required
        if isinstance(spectrograph, basestring):
            self.spectrograph = load_spectrograph(spectrograph=spectrograph)
        elif isinstance(spectrograph, Spectrograph):
            self.spectrograph = spectrograph
        else:
            raise TypeError('Must provide a name or instance for the Spectrograph.')

        # Used to construct the BPM using the spectrograph class
        self.shape = shape
        self.filename = filename
        self.det = det

        # Used to construct the BPM from the bias
        self.msbias = msbias
        self.trim = trim

        # Attributes (set after init)
        self.frametype = frametype

        # Output
        self.bpm_img = None

    def build(self):
        """
        Generate the BPM Image

        Returns
        -------
        self.bpm

        """
        # Will raise an exception if could not construct the BPM
        if self.msbias is None:
            # TODO: Is this trimmed?
            self.bpm_img = self.spectrograph.bpm(shape=self.shape, filename=self.filename,
                                                 det=self.det)
            return self.bpm_img

        # Get the data sections
        datasec, one_indexed, include_end, transpose \
                = self.spectrograph.get_image_section(self.filename, self.det, section='datasec')
        datasec = [ arparse.sec2slice(sec, one_indexed=one_indexed, include_end=include_end,
                                      require_dim=2, transpose=transpose) for sec in datasec ]
       
        # Identify the bad pixels
        self.bpm_img = arprocimg.find_bad_pixels(self.msbias,
                                    self.spectrograph.detector[self.det-1]['numamplifiers'],
                                    datasec, trim=self.trim)
        return self.bpm_img


