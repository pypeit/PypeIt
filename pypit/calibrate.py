""" Class for guiding calibration object generation in PYPIT
"""
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np
import os


from pypit import msgs
from pypit import ardebug as debugger
from pypit import arcimage
from pypit import biasframe


# For out of PYPIT running
if msgs._debug is None:
    debug = debugger.init()
    debug['develop'] = True
    msgs.reset(debug=debug, verbosity=2)

# Does not need to be global, but I prefer it
frametype = 'bias'

# Place these here or elsewhere?
#  Wherever they be, they need to be defined, described, etc.
#  These are settings beyond those in the Parent class (ProcessImages)
additional_default_settings = {frametype: {'useframe': 'none'}}


class Calibrate(object):
    """
    This class is primarily designed to guide the generation of calibration images
    and objects in PYPIT

    Parameters
    ----------

    Attributes
    ----------

    Inherited Attributes
    --------------------
    """
    def __init__(self, fitstbl):

        # Parameters unique to this Object
        self.fitstbl = fitstbl

        # Attributes
        self.calib_dict = {}
        self.det = None
        self.sci_ID = None
        self.settings = None
        self.setup = None

        # Internals
        self.msarc = None
        self.msbias = None

    def set(self, setup, det, sci_ID, settings):
        self.setup = setup
        self.det = det
        self.sci_ID = sci_ID
        self.settings = settings.copy()
        #
        if self.setup not in self.calib_dict.keys():
            self.calib_dict[self.setup] = {}

    def get_arc(self, bias=None):
        # Checks
        if bias is not None:
            self.msbias = bias
        if self.msbias is None:
            msgs.error("msbias needs to be set prior to arc")
        #
        if 'arc' in self.calib_dict[self.setup].keys():
            self.msarc = self.calib_dict[self.setup]['arc']
        else:
            # Grab it -- msarc will be a 2D image
            self.msarc, self.arcImage = arcimage.get_msarc(self.det, self.setup, self.sci_ID,
                                          self.fitstbl, self.settings, self.msbias)
            # Save
            self.calib_dict[self.setup]['arc'] = self.msarc
        # Return
        return self.msarc

    def get_bias(self):
        for item in ['setup', 'det', 'sci_ID', 'settings']:
            if getattr(self, item) is None:
                msgs.error("Use self.set to specify '{:s}' prior to generating the bias".format(item))

        if 'bias' in self.calib_dict[self.setup].keys():
            self.msbias = self.calib_dict[self.setup]['bias']
        else:
            # Grab it
            #   Bias will either be an image (ndarray) or a command (str, e.g. 'overscan') or none
            self.msbias, self.biasFrame = biasframe.get_msbias(
                self.det, self.setup, self.sci_ID, self.fitstbl, self.settings)
            # Save
            self.calib_dict[self.setup]['bias'] = self.msbias
        # Return
        return self.msbias

    def full_calibrate(self):
        self.msbias = self.get_bias()
        self.msarc = self.get_arc()
