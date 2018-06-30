""" Class for guiding calibration object generation in PYPIT
"""
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np
import os


from pypit import msgs
from pypit import ardebug as debugger
from pypit import arpixels

from pypit import arcimage
from pypit import biasframe
from pypit.spectrographs import bpmimage

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
        self.mbpm = None
        self.pixlocn = None

    def set(self, setup, det, sci_ID, settings):
        self.setup = setup
        self.det = det
        self.sci_ID = sci_ID
        self.settings = settings.copy()
        # Setup the calib_dict
        if self.setup not in self.calib_dict.keys():
            self.calib_dict[self.setup] = {}

    def get_arc(self, bias=None):
        # Checks
        if bias is not None:
            self.msbias = bias
        if self.msbias is None:
            msgs.error("msbias needs to be set prior to arc")
        # Check
        self._chk_set(['setup', 'det', 'sci_ID', 'settings'])
        #
        if 'arc' in self.calib_dict[self.setup].keys():
            self.msarc = self.calib_dict[self.setup]['arc']
        else:
            # Instantiate with everything needed to generate the image (in case we do)
            self.arcImage = arcimage.ArcImage([], spectrograph=self.fitstbl['instrume'][0],
                                settings=self.settings, det=self.det, setup=self.setup,
                                sci_ID=self.sci_ID, msbias=self.msbias, fitstbl=self.fitstbl)
            # Load the MasterFrame (if it exists and is desired)?
            self.msarc = self.arcImage.master()
            if self.msarc is None:  # Otherwise build it
                msgs.info("Preparing a master {0:s} frame".format(self.arcImage.frametype))
                self.msarc = self.arcImage.build_image()
                # Save to Masters
                self.arcImage.save_master(self.msarc, raw_files=self.arcImage.file_list, steps=self.arcImage.steps)
            # Return
            return self.msarc
            '''
            # Grab it -- msarc will be a 2D image
            self.msarc, self.arcImage = arcimage.get_msarc(self.det, self.setup, self.sci_ID,
                                          self.fitstbl, self.settings, self.msbias)
            '''
            # Save
            self.calib_dict[self.setup]['arc'] = self.msarc
        # Return
        return self.msarc

    def get_bias(self):
        self._chk_set(['setup', 'det', 'sci_ID', 'settings'])
        if 'bias' in self.calib_dict[self.setup].keys():
            self.msbias = self.calib_dict[self.setup]['bias']
        else:
            # Init
            self.biasFrame = biasframe.BiasFrame(settings=self.settings, setup=self.setup,
                                  det=self.det, fitstbl=self.fitstbl, sci_ID=self.sci_ID)
            # Load the MasterFrame (if it exists and is desired) or the command (e.g. 'overscan')
            msbias = self.biasFrame.master()
            if msbias is None:  # Build it and save it
                msbias = self.biasFrame.build_image()
                self.biasFrame.save_master(msbias, raw_files=self.biasFrame.file_list,
                                      steps=self.biasFrame.steps)
            # Return
            return msbias#, biasFrame
            '''
            # Grab it
            #   Bias will either be an image (ndarray) or a command (str, e.g. 'overscan') or none
            self.msbias, self.biasFrame = biasframe.get_msbias(
                self.det, self.setup, self.sci_ID, self.fitstbl, self.settings)
            # Save
            self.calib_dict[self.setup]['bias'] = self.msbias
            '''
        # Return
        return self.msbias

    def get_bpm(self, arc=None, bias=None):
        if arc is not None:
            self.msarc = arc
        if bias is not None:
            self.msbias = bias
        # Check me
        self._chk_set(['det', 'settings'])
        ###############################################################################
        # Generate a bad pixel mask (should not repeat)
        if 'bpm' in self.calib_dict[self.setup].keys():
            self.msbpm = self.calib_dict[self.setup]['bpm']
        else:
            scidx = np.where((self.fitstbl['sci_ID'] == self.sci_ID) & self.fitstbl['science'])[0][0]
            bpmImage = bpmimage.BPMImage(spectrograph=self.fitstbl['instrume'][0],
                            settings=self.settings, det=self.det,
                            shape=self.msarc.shape,
                            binning=self.fitstbl['binning'][scidx],
                            reduce_badpix=self.settings['reduce']['badpix'],
                            msbias=self.msbias)
            self.msbpm = bpmImage.build()
            '''
            # Grab it -- msbpm is a 2D image
            scidx = np.where((self.fitstbl['sci_ID'] == self.sci_ID) & self.fitstbl['science'])[0][0]
            self.msbpm, _ = bpmimage.get_mspbm(self.det, self.fitstbl['instrume'][0],
                                          self.settings, self.msarc.shape,
                                          binning=self.fitstbl['binning'][scidx],
                                          reduce_badpix=self.settings['reduce']['badpix'],
                                          msbias=self.msbias)
            '''
            # Save
            self.calib_dict[self.setup]['bpm'] = self.msbpm
        # Return
        return self.msbpm

    def get_slits(self):
        if 'trace' in self.calib_dict[self.setup].keys():  # Internal
            tslits_dict = self.calib_dict[self.setup]['trace']
        else:
            # Setup up the settings (will be Refactored with settings)
            # Get it -- Key arrays are in the tslits_dict
            tslits_dict, _ = traceslits.get_tslits_dict(
                det, setup, spectrograph, sci_ID, ts_settings, tsettings, fitstbl, pixlocn,
                msbias, msbpm, datasec_img, trim=settings.argflag['reduce']['trim'])
            # Save in calib
            calib_dict[setup]['trace'] = tslits_dict

        ###############################################################################
        # Initialize maskslits
        nslits = tslits_dict['lcen'].shape[1]
        maskslits = np.zeros(nslits, dtype=bool)

    def make_pixlocn(self, arc=None):
        if arc is not None:
            self.msarc = arc
        # Check
        self._chk_set('settings')
        #
        xgap = self.settings['detector']['xgap']
        ygap = self.settings['detector']['ygap']
        ysize = self.settings['detector']['ysize']
        self.pixlocn = arpixels.core_gen_pixloc(self.msarc.shape, xgap=xgap, ygap=ygap, ysize=ysize)
        # Return
        return self.pixlocn

    def _chk_set(self, items):
        for item in items:
            if getattr(self, item) is None:
                msgs.error("Use self.set to specify '{:s}' prior to generating the bias".format(item))


    def full_calibrate(self):
        self.msbias = self.get_bias()
        self.msarc = self.get_arc(self.msbias)
        self.msbpm = self.get_bpm(self.msarc, self.msbias)
        self.pixlocn = self.make_pixlocn(self.msarc)
