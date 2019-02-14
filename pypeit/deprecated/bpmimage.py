# Module for generating the BPM image
from __future__ import absolute_import, division, print_function

import numpy as np
import os

from pypeit import msgs
from pypeit.core import procimg
from pypeit.core import parse

from pypeit.spectrographs.util import load_spectrograph

from pypeit import debugger


class BPMImage(object):
    """
    This class is primarily designed to generate an Bad Pixel Image
      The master() method will return the image

    Should provide both shape and filename to ensure that the
    spectrograph can construct the bpm, if reduce_badpix is False.

    There are several ways to build a BPM:
       1. keck_lris_red
          Provide binning and detector number to generate this instrument specific BPM
       2. keck_deimos
          Provide detector number to generate this instrument specific BPM
       3. From a bias image (not well tested)
          Set reduce_badpix='bias'
          Provide msbias image
       4. Dummy image from shape

    Args:
        spectrograph (str or :class:pypeit.spectrographs.spectrograph.Spectrograph):
           Used to specify properties of the detector (for processing)
           Attempt to set with settings['run']['spectrograph'] if not input
        shape (tuple):
          Image shape;  used to construct a dummy BPM if all else fails
          spec, spat This is the TRIMMED size of the raw images.
        det (int, optional):
          Detector index
          Required for LRISr and DEIMOS

    Attributes:
        bpm_img (np.ndarray): BPM image

    """

    # Frametype is a class attribute
    frametype = 'bpm'

    # Keep order same as processimages (or else!)
    def __init__(self, spectrograph, shape, det=None):

        # This function interprets both strings and spectrograph
        # instances now
        self.spectrograph = load_spectrograph(spectrograph)

        # Used to construct the BPM using the spectrograph class
        self.shape = shape
        #self.filename = filename
        self.det = det

        # Used to construct the BPM from the bias
        #self.msbias = msbias JFH This is not yet supported anyway
        #self.trim = trim -- Killing this option for now

        # spectrograph or msbias must be defined
        if self.spectrograph is None:
            msgs.error('BPMImage instantiation incomplete.  Must provide spectrograph.')


        # Output
        self.bpm_img = None

    def build(self, filename=None):
        """
        Generate the BPM Image

        Simmpl

        Args:
            datasec:
            filename:

        Returns:

        """
        self.bpm_img = self.spectrograph.bpm(shape=self.shape, det=self.det, filename=filename)

