""" Module for generating a normalized pixel flat
NOT IMPLEMENTED YET
"""

from __future__ import absolute_import, division, print_function

import inspect
import numpy as np
import os

from pypit import msgs
from pypit import processimages
from pypit import masterframe
from pypit import arproc
from pypit.core import arsort

from pypit import ardebug as debugger

# For out of PYPIT running
if msgs._debug is None:
    debug = debugger.init()
    debug['develop'] = True
    msgs.reset(debug=debug, verbosity=2)

# Does not need to be global
frametype = 'pixelflat'


class NormPixelFlat(processimages.ProcessImages, masterframe.MasterFrame):
    """
    This class is primarily designed to generate a normalized pixel flat image
      The master() method will return the image (and build it if necessary)

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
      Setup tag;  required for MasterFrame functionality
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
      Set to 'trace_image'

    Inherited Attributes
    --------------------
    stack : ndarray
    """
    # Keep order same as processimages (or else!)
    def __init__(self, file_list=[], spectrograph=None, settings=None, det=1, setup=None, sci_ID=None,
                 msbias=None, datasec_img=None, fitstbl=None):

        # Parameters unique to this Object
        self.sci_ID = sci_ID
        self.msbias = msbias
        self.fitstbl = fitstbl
        self.datasec_img = datasec_img

        # Start us up
        processimages.ProcessImages.__init__(self, file_list, spectrograph=spectrograph, settings=settings, det=det)

        # Attributes (set after init)
        self.frametype = frametype
        self.setup = setup

        # Settings
        # The copy allows up to update settings with user settings without changing the original
        if settings is None:
            # Defaults
            self.settings = processimages.default_settings.copy()
        else:
            self.settings = settings.copy()
            # The following is somewhat kludgy and the current way we do settings may
            #   not touch all the options (not sure .update() would help)
            if 'combine' not in settings.keys():
                self.settings['combine'] = settings[self.frametype]['combine']

        # Child-specific Internals
        #    See ProcessImages for the rest

        # MasterFrames
        masterframe.MasterFrame.__init__(self, self.frametype, self.setup, self.settings)

    def build_image(self):
        # Get list of arc frames for this science frame
        if self.nfiles == 0:
            self.file_list = arsort.list_of_files(self.fitstbl, self.frametype, self.sci_ID)
        # Combine
        self.stack = self.process(bias_subtract=self.msbias)
        #
        return self.stack

    def master(self):
        """
        Build the master frame and save to disk
         OR
        Load the master frame from disk

        Returns
        -------
        msframe : ndarray

        """
        # Load the MasterFrame if it exists and user requested one to load it
        msframe, header, raw_files = self.load_master_frame()
        # Build?
        if msframe is None:
            msgs.info("Preparing a master {0:s} frame".format(self.frametype))
            # Stack the images
            #   Could save a little time by reusing mstrace, but no bother
            mspixelflat = self.build_image()
            # Apply gain (instead of ampsec scale)
            mspixelflat *= arproc.gain_frame(self.datasec_img, self.settings['detector']['numamplifiers'],
                                             self.settings['detector']['gain'])
            # Normalize the flat field
            msgs.info("Normalizing the pixel flat")
            slit_profiles, mstracenrm, msblaze, flat_ext1d, extrap_slit = \
                arproc.slit_profile(mspixelflat, self.det,
                                    ntcky=self.settings.argflag['reduce']['flatfield']['params'][0])

            # Save to Masters
            self.save_master(msframe, raw_files=self.file_list, steps=self.steps)
        else:
            # Prevent over-writing the master frame when it is time to save
            self.settings['reduce']['masters']['loaded'].append(self.frametype+self.setup)
            # Put in
            self.stack = msframe
        # Return
        return msframe.copy()
