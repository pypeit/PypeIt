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
    # Keep order same as processimages (or else!)
    def __init__(self, file_list=[], settings=None, det=1, setup=None, ind=[], fitsdict=None):

        # Parameters unique to this Object
        self.frametype = frametype
        self.ind = ind
        self.fitsdict = fitsdict

        # Start us up
        processimages.ProcessImages.__init__(self, file_list, settings=settings, det=det)

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
                self.settings['combine'] = settings[self.frametype]['combine']

        # MasterFrames -- Is this required?
        masterframe.MasterFrame.__init__(self, self.frametype, setup, self.settings)

        # Child-specific Internals
        #    See ProcessImages for the rest

    def process_bias(self, overwrite=False):
        # Wrapper
        self.stack = self.process(bias_subtract=None, trim=False, overwrite=overwrite)
        return self.stack.copy()

    def build_master(self):
        # Generate a bias or dark image (or load a pre-made Master by PYPIT)?
        if self.settings[self.frametype]['useframe'] in ['bias', 'dark']:
            # Load the MasterFrame if it exists
            msframe, header, raw_files = self.load_master_frame()
            if msframe is None:
                msgs.info("Preparing a master {0:s} frame".format(self.settings[self.frametype]['useframe']))
                # Get all of the bias frames for this science frame
                if self.nfiles == 0:
                    for i in range(len(self.ind)):
                        self.file_list.append(self.fitsdict['directory'][self.ind[i]]+self.fitsdict['filename'][self.ind[i]])
                # Combine
                msframe = self.process_bias()
                # Save to Masters
                self.save_master(msframe, raw_files=self.file_list, steps=self.steps)
            else:
                # Prevent over-writing the master frame when it is time to save
                self.settings['reduce']['masters']['loaded'].append(self.frametype)
        # Simple command?
        elif self.settings[self.frametype]['useframe'] in ['overscan', 'none']:
            if self.settings[self.frametype]['useframe'] == 'none':
                msgs.info("Will not perform bias/dark subtraction")
            return self.settings[self.frametype]['useframe']
        # It must be a user-specified file the user wishes to load
        else:
            msframe_name = self.settings['run']['directory']['master']+u'/'+self.settings[self.frametype]['useframe']
            msframe, head, _ = armasters._core_load(msframe_name, frametype=self.frametype)
            self.settings['reduce']['masters']['loaded'].append(self.frametype+self.setup)

        self.stack = msframe
        return msframe.copy()

