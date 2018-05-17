# Module for generating the Trace image
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

# Does not need to be global, but I prefer it
frametype = 'trace_image'

# Place these here or elsewhere?
#  Wherever they be, they need to be defined, described, etc.
#  These are settings beyond those in the Parent class (ProcessImages)
additional_default_settings = {frametype: {'useframe': 'none'}}


class TraceImage(processimages.ProcessImages):
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
    ind : list (optional)
      Indices for bias frames (if a Bias image may be generated)
    fitstbl : Table (optional)
      FITS info (mainly for filenames)

    Attributes
    ----------
    frametype : str
      Set to 'trace_image'

    Inherited Attributes
    --------------------
    stack : ndarray
    """
    # Keep order same as processimages (or else!)
    def __init__(self, file_list, spectrograph=None, settings=None, det=1):

        # Parameters unique to this Object

        # Start us up
        processimages.ProcessImages.__init__(self, file_list, spectrograph=spectrograph, settings=settings, det=det)

        # Attributes (set after init)
        self.frametype = frametype

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
                self.settings['combine'] = settings['trace']['combine']

        # Child-specific Internals
        #    See ProcessImages for the rest


