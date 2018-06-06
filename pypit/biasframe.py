""" Module for guiding Bias subtraction including generating a Bias image as desired
"""
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np
import os


from pypit import msgs
from pypit import ardebug as debugger
from pypit import armasters
from pypit.core import arsort
from pypit import processimages
from pypit import masterframe


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


class BiasFrame(processimages.ProcessImages, masterframe.MasterFrame):
    """
    This class is primarily designed to generate a Bias frame for bias subtraction
      It also contains I/O methods for the Master frames of PYPIT
      The build_master() method will return a simple command (str) if that is the specified setting
      in settings['bias']['useframe']

    Parameters
    ----------
    file_list : list (optional)
      List of filenames
    spectrograph : str (optional)
       Used to specify properties of the detector (for processing)
       Attempt to set with settings['run']['spectrograph'] if not input
    settings : dict (optional)
      Settings for trace slits
    setup : str (optional)
      Setup tag
    det : int, optional
      Detector index, starts at 1
    fitstbl : Table (optional)
      FITS info (mainly for filenames)
    sci_ID : int (optional)
      Science ID value
      used to match bias frames to the current science exposure

    Attributes
    ----------
    frametype : str
      Set to 'bias'

    Inherited Attributes
    --------------------
    stack : ndarray
    """
    # Keep order same as processimages (or else!)
    def __init__(self, file_list=[], spectrograph=None, settings=None, det=1, setup=None, fitstbl=None,
                 sci_ID=None):

        # Parameters unique to this Object
        self.fitstbl = fitstbl
        self.sci_ID = sci_ID

        # Start us up
        processimages.ProcessImages.__init__(self, file_list, spectrograph=spectrograph, settings=settings, det=det)

        # Attributes  (set after ProcessImages call)
        self.frametype = frametype

        # Settings
        # The copy allows up to update settings with user settings without changing the original
        if settings is None:
            # Defaults
            self.settings = processimages.default_settings().update(additional_default_settings)
        else:
            self.settings = settings.copy()
            # The following is somewhat kludgy and the current way we do settings may
            #   not touch all the options (not sure .update() would help)
            if 'combine' not in settings.keys():
                self.settings['combine'] = settings[self.frametype]['combine']

        # MasterFrames
        masterframe.MasterFrame.__init__(self, self.frametype, setup, self.settings)

        # Child-specific Internals
        #    See ProcessImages for the rest

    def process_bias(self, overwrite=False):
        """
        Process the input bias frames with ProcessImages.process()
          Avoid bias subtraction
          Avoid trim

        Parameters
        ----------
        overwrite : bool, optional

        Returns
        -------
        stack : ndarray

        """
        # Wrapper
        self.stack = self.process(bias_subtract=None, trim=False, overwrite=overwrite)
        return self.stack.copy()

    def build_image(self):
        """
        Generate the image

        Returns
        -------
        stack : ndarray

        """
        # Get all of the bias frames for this science frame
        if self.nfiles == 0:
            self.file_list = arsort.list_of_files(self.fitstbl, 'bias', self.sci_ID)
        # Combine
        self.stack = self.process_bias()
        #
        return self.stack

    def master(self):
        """
        Load the master frame from disk, as the settings allow
        or return the command
        or return None

        Returns
        -------
        msframe : ndarray or str or None

        """
        # Generate a bias or dark image (or load a pre-made Master by PYPIT)?
        if self.settings[self.frametype]['useframe'] in ['bias', 'dark']:
            # Load the MasterFrame if it exists and user requested one to load it
            msframe, header, raw_files = self.load_master_frame()
            if msframe is None:
                return None
            else:
                # Prevent over-writing the master frame when it is time to save
                self.settings['masters']['loaded'].append(self.frametype)
        # Simple command?
        elif self.settings[self.frametype]['useframe'] in ['overscan', 'none']:
            if self.settings[self.frametype]['useframe'] == 'none':
                msgs.info("Will not perform bias/dark subtraction")
            return self.settings[self.frametype]['useframe']
        # It must be a user-specified file the user wishes to load
        else:
            msframe_name = self.settings['run']['directory']['master']+u'/'+self.settings[self.frametype]['useframe']
            msframe, head, _ = armasters._core_load(msframe_name, frametype=self.frametype)
            self.settings['masters']['loaded'].append(self.frametype+self.setup)

        # Put in
        self.stack = msframe
        return msframe.copy()

