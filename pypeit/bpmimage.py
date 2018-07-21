# Module for generating the BPM image
from __future__ import absolute_import, division, print_function

import numpy as np
import os

from pypit import msgs
from pypit.core import arprocimg

from pypit.spectrographs.spectrograph import Spectrograph
from pypit.spectrographs.util import load_spectrograph

from pypit import ardebug as debugger

try:
    basestring
except NameError:
    basestring = str

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

        This is the UNTRIMMED size of the raw images.

    msbias : ndarray, optional
      Used to construct the BPM if reduce_badpix=='bias'
    reduce_badpix : str, optional
      'bias' -- Build from bias images

    Attributes
    ----------
    bpm : ndarray

    """

    # Frametype is a class attribute
    frametype = 'bpm'

    # Keep order same as processimages (or else!)
    def __init__(self, spectrograph=None, shape=None, filename=None, det=None, msbias=None,
                 trim=True):

        # Spectrograph is required
        if isinstance(spectrograph, basestring):
            self.spectrograph = load_spectrograph(spectrograph=spectrograph)
        elif isinstance(spectrograph, Spectrograph):
            self.spectrograph = spectrograph
        else:
            self.spectrograph = None

        # Used to construct the BPM using the spectrograph class
        self.shape = shape
        self.filename = filename
        self.det = det

        # Used to construct the BPM from the bias
        self.msbias = msbias
        self.trim = trim

        # spectrograph or msbias must be defined
        if self.spectrograph is None and msbias is None:
            msgs.error('BPMImage instantiation incomplete.  Must provide spectrograph or msbias.')

        # Output
        self.bpm_img = None

    def build(self, datasec=None):
        """
        Generate the BPM Image

        Returns
        -------
        self.bpm

        """
        # Will raise an exception if could not construct the BPM
        if self.msbias is None:
            # WARNING: This assumes shape is the untrimmed size of the
            # image!
            self.bpm_img = self.spectrograph.bpm(shape=self.shape, filename=self.filename,
                                                 det=self.det)
            if self.trim:
                mask = self.spectrograph.get_datasec_img(filename=self.filename, det=self.det) < 1

        else:
            _datasec = datasec
            if self.spectrograph is None and _datasec is None:
                _datasec = [[slice(None),slice(None)]]
                _numamplifiers = 1
            if _datasec is None:
                # Get the data sections from the spectrograph definition
                _datasec, one_indexed, include_end, transpose \
                        = self.spectrograph.get_image_section(self.filename, self.det,
                                                              section='datasec')
                _datasec = [ arparse.sec2slice(sec, one_indexed=one_indexed,
                                               include_end=include_end, require_dim=2,
                                               transpose=transpose)
                                    for sec in _datasec ]
                _numamplifiers = self.spectrograph.detector[self.det-1]['numamplifiers']

            # Identify the bad pixels.  WARNING: When self.trim is True,
            # this assumes that msbias is not trimmed.
            self.bpm_img = arprocimg.find_bad_pixels(self.msbias, _numamplifiers, _datasec)

            if self.trim:
                mask = np.ones_like(self.bpm_img, dtype=bool)
                for d in _datasec:
                    mask[d] = False

        # Trim
        if self.trim:
            self.bpm_img = arprocimg.trim_frame(self.bpm_img, mask)
        return self.bpm_img


