# Module for generating the BPM image
from __future__ import absolute_import, division, print_function

import numpy as np
import os

from pypit import msgs
from pypit.core import arprocimg

from .spectrographs.util import load_spectrograph

from pypit import ardebug as debugger

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
    def __init__(self, spectrograph, shape=None, filename=None, det=None, msbias=None):

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
        self.bpm_img = self.spectrograph.bpm(shape=self.shape, filename=self.filename,
                                             det=self.det) if self.msbias is None else \
                            arprocimg.badpix(self.msbias,
                                    self.spectrograph.detector[self.det-1]['numamplifiers'],
                                    self.spectrograph.detector[self.det-1]['datasec'])
        return self.bpm_img
#            # Get all of the bias frames for this science frame
#            if self.msbias is None:
#                msgs.warn("No bias frame provided!")
#                msgs.info("Not preparing a bad pixel mask")
#                return False
#            # TODO -- Deal better with this datasec kludge
#            datasec = []
#            for i in range(self.settings['detector']['numamplifiers']):
#                sdatasec = "datasec{0:02d}".format(i+1)
#                datasec.append(self.settings['detector'][sdatasec])


