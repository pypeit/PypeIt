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

    Parameters
    ----------
    spectrograph : str or Spectrograph
       Used to specify properties of the detector (for processing)
       Attempt to set with settings['run']['spectrograph'] if not input
    settings : dict (optional)
      Settings for trace slits
    det : int, optional
      Detector index
      Required for LRISr and DEIMOS
    filename : str, optional
      Used primarily for binning...
    shape : tuple
      Image shape;  used to construct a dummy BPM if all else fails

        This is the TRIMMED size of the raw images.

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
    def __init__(self, spectrograph, shape=None, det=None):

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

    def build(self, datasec=None):
        """
        Generate the BPM Image

        Returns
        -------
        self.bpm

        """
        # Will raise an exception if could not construct the BPM
        if self.shape is not None:
            # WARNING: This assumes shape is the untrimmed size of the
            # image!
            self.bpm_img = self.spectrograph.bpm(shape=self.shape, det=self.det)
            #if self.trim:
            #    mask = self.spectrograph.get_datasec_img(filename=self.filename, det=self.det) < 1
        else:
            # JFH This is all deprecated here
            debugger.set_trace() # TOO EXPERIMENTAL
            _datasec = datasec
            if self.spectrograph is None and _datasec is None:
                _datasec = [(slice(None),slice(None))]
                _numamplifiers = 1
            if _datasec is None:
                # Get the data sections from the spectrograph definition
                _datasec, one_indexed, include_end, transpose \
                        = self.spectrograph.get_image_section(self.filename, self.det,
                                                              section='datasec')
                _datasec = [ parse.sec2slice(sec, one_indexed=one_indexed,
                                               include_end=include_end, require_dim=2,
                                               transpose=transpose)
                                    for sec in _datasec ]
                _numamplifiers = self.spectrograph.detector[self.det-1]['numamplifiers']

            # Identify the bad pixels.  WARNING: When self.trim is True,
            # this assumes that msbias is not trimmed.
            self.bpm_img = procimg.find_bad_pixels(self.msbias, _numamplifiers, _datasec)

            if self.trim:
                mask = np.ones_like(self.bpm_img, dtype=bool)
                for d in _datasec:
                    mask[d] = False

        # Trim
        #if self.trim:
        #    self.bpm_img = procimg.trim_frame(self.bpm_img, mask)
        return self.bpm_img


