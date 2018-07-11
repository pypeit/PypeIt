# Module for generating the Trace image
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np


from pypit import msgs
from pypit import ardebug as debugger
from pypit import processimages


# For out of PYPIT running
if msgs._debug is None:
    debug = debugger.init()
    debug['develop'] = True
    msgs.reset(debug=debug, verbosity=2)

# Does not need to be global, but I prefer it
# JFH Not sure why you prefer this, since the global scope is irrelvant here unless we add more code, and global is
# always less readable than local.

frametype = 'trace_image'

# ToDO This syntax for instantiaion differs slightlhy from ArcIMage, and I think they need to be homogenized.
# TODO I don't understand why TraceImage is its own class, since the traceimg is going to be the same as the flat image
# but the Flatimage does not have its own class but is instead embedded in the Flatfield class. Maybe the idea
# is to have the flexibilty to have the traceimage not be the flat image? But this could also apply to the tilts, since
# I might want to use arcs for the wavelength solution but use the science data itself for the traceimage. I think these
# image classes are just adding a bunch of extra code for no reason. I get that the motivation is to avoid making the msarc twice
# but perhaps we can pay that computational cost (which is small) for simpler code.
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
    #ToDO datasec_img is not documented!

    Attributes
    ----------
    frametype : str
      Set to 'trace_image'

    Inherited Attributes
    --------------------
    stack : ndarray
    """
    # Keep order same as processimages (or else!)

    # TODO user settings need to be added to these codes!
    def __init__(self, file_list, spectrograph=None, settings=None, det=1, datasec_img=None):

        # Parameters unique to this Object

        # Start us up
        processimages.ProcessImages.__init__(self, file_list, spectrograph=spectrograph,
                                             settings=settings, det=det, datasec_img=datasec_img)

        # Attributes (set after init)
        self.frametype = frametype

        # Settings
        # The copy allows up to update settings with user settings without changing the original
        if settings is None:
            # Defaults
            self.settings = processimages.default_settings()
        else:
            self.settings = settings.copy()
            # The following is somewhat kludgy and the current way we do settings may
            #   not touch all the options (not sure .update() would help)
            if 'combine' not in settings.keys():
                self.settings['combine'] = settings['trace']['combine']

        # Child-specific Internals
        #    See ProcessImages for the rest


