# Module for guiding Slit/Order tracing
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np
import os

from astropy.io import fits

from pypit import msgs
from pypit import ardebug as debugger
from pypit import armasters
from pypit import processimages
from pypit import masterframe
from pypit import ginga


# For out of PYPIT running
if msgs._debug is None:
    debug = debugger.init()
    debug['develop'] = True
    msgs.reset(debug=debug, verbosity=2)

# Place these here or elsewhere?
#  Wherever they be, they need to be defined, described, etc.
#  These are settings beyond those in the Parent class (ProcessImages)
frametype = 'bias'
additional_default_settings = {frametype: {'useframe': 'none'}}


class BiasFrame(processimages.ProcessImages, masterframe.MasterFrame):
    """
    This class is primarily designed to generate a Bias frame for bias subtraction
      It also contains I/O methods for the Master frames of PYPIT
      And the run() method will return a simple command (str) if that is the specified setting
      in settings['bias']['useframe']

    Parameters
    ----------
    file_list : list (optional)
      List of filenames
    settings : dict (optional)
      Settings for trace slits
    setup : str (optional)
      Setup tag
    det : int, optional
      Detector index, starts at 1
    ind : list (optional)
      Indices for bias frames (if a Bias image may be generated)
    fitsdict : dict (optional)
      FITS info (mainly for filenames)

    Attributes
    ----------
    images : list
    stack : ndarray
    frametype : str
      Set to 'bias'
    """
    def __init__(self, settings=None, setup=None, file_list=[], det=1, ind=[], fitsdict=None):
        # Hard-coded
        self.frametype = frametype

        # Start us up
        processimages.ProcessImages.__init__(self, file_list, settings=settings)


        # Parameters
        self.det = det
        self.setup = setup
        self.ind = ind
        self.fitsdict = fitsdict

        # Settings
        # The copy allows up to update settings with user settings without changing the original
        if settings is None:
            # Defaults
            self.settings = processimages.default_settings.copy().update(additional_default_settings)
        else:
            self.settings = settings.copy()
            # The following is somewhat kludgy and the current way we do settings may
            #   not touch all the options (not sure .update() would help)
            if 'combine' not in settings.keys():
                self.settings['combine'] = settings['bias']['combine']

        # MasterFrames -- Is this required?
        debugger.set_trace()
        masterframe.MasterFrame.__init__(self, self.frametype, self.setup, self.settings)

        # Child-specific Internals
        #    See ProcessImages for the rest

    def combine(self, overwrite=False):
        # Over-write?
        if (inspect.stack()[0][3] in self.steps) & (not overwrite):
            msgs.warn("Images already combined.  Use overwrite=True to do it again.")
            return

        # Allow for one-stop-shopping
        if 'load_images' not in self.steps:
            self.load_images()

        # Create proc_images from raw_images if need be
        self.proc_images = np.zeros((self.raw_images[0].shape[0],
                                         self.raw_images[0].shape[1],
                                         self.nloaded))
        for kk,image in enumerate(self.raw_images):
                self.proc_images[:,:,kk] = image
        # Combine
        self.stack = self._combine(frametype=self.frametype)
        # Step
        self.steps.append(inspect.stack()[0][3])
        return self.stack.copy()

    def build_master(self):
        # Generate a bias or dark image (or load a pre-made Master by PYPIT)?
        if self.settings['bias']['useframe'] in ['bias', 'dark']:
            # Load the MasterFrame if it exists
            msbias, header, raw_files = armasters.core_load_master_frame("bias", self.setup, self.mdir)
            if msbias is None:
                msgs.info("Preparing a master {0:s} frame".format(self.settings[self.frametype]['useframe']))
                # Get all of the bias frames for this science frame
                if self.nfiles == 0:
                    for i in range(len(self.ind)):
                        self.file_list.append(self.fitsdict['directory'][self.ind[i]]+self.fitsdict['filename'][self.ind[i]])
                # Combine
                msbias = self.combine()
                # Save to Masters
                self.save_master(msbias, raw_files=self.file_list, steps=self.steps)
            else:
                # Prevent over-writing the master bias when it is time to save
                self.settings['reduce']['masters']['loaded'].append(self.frametype)
        # Simple command?
        elif self.settings[self.frametype]['useframe'] in ['overscan', 'none']:
            if self.settings[self.frametype]['useframe'] == 'none':
                msgs.info("Will not perform bias/dark subtraction")
            return self.settings[self.frametype]['useframe']
        # It must be a user-specified file the user wishes to load
        else:
            msbias_name = self.settings['run']['directory']['master']+u'/'+self.settings['bias']['useframe']
            msbias, head = armasters.load_master(msbias_name, frametype=self.frametype)
            self.settings['reduce']['masters']['loaded'].append(self.frametype+self.setup)

        self.stack = msbias
        return msbias.copy()

