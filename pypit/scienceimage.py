# Module for the ScienceImage class
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np
import os

#from importlib import reload

from pypit import msgs
from pypit import processimages
from pypit import armasters
from pypit import masterframe
from pypit.core import arprocimg
from pypit.core import arflat
from pypit import ginga

from pypit import ardebug as debugger

# For out of PYPIT running
if msgs._debug is None:
    debug = debugger.init()
    debug['develop'] = True
    msgs.reset(debug=debug, verbosity=2)

# Does not need to be global, but I prefer it
frametype = 'science'


class ScienceImage(processimages.ProcessImages):
    """
    This class will generat the pixel-level FlatField
      The master() method returns the image

    Parameters
    ----------
    file_list : list
      List of raw files to produce the flat field
    spectrograph : str
    settings : dict-like
    msbias : ndarray or str or None
    tslits_dict : dict
      dict from TraceSlits class (e.g. slitpix)
    tilts : ndarray
      tilts from WaveTilts class
    det : int
    setup : str

    Attributes
    ----------
    frametype : str
      Set to 'pixelflat'
    mspixelflat : ndarray
      Stacked image
    mspixelflatnrm : ndarray
      Normalized flat
    extrap_slit
    msblaze : ndarray
      Blaze function fit to normalize
    blazeext :
    slit_profiles : ndarray
      Slit profile(s)
    self.ntckx : int
      Number of knots in the spatial dimension
    self.ntcky : int
      Number of knots in the spectral dimension

    """
    # Keep order same as processimages (or else!)
    def __init__(self, file_list=[], spectrograph=None, settings=None,
                 tslits_dict=None, tilts=None, det=None, setup=None, datasec_img=None,
                 bpm=None):

        # Parameters unique to this Object
        self.det = det
        self.setup = setup
        self.tslits_dict = tslits_dict
        self.tilts = tilts

        # Start us up
        processimages.ProcessImages.__init__(self, file_list, spectrograph=spectrograph,
                                             settings=settings, det=det, datasec_img=datasec_img,
                                             bpm=bpm)

        # Attributes (set after init)
        self.frametype = frametype

        # Key outputs
        self.mspixelflat = None
        self.mspixelflatnrm = None

        # Settings
        # The copy allows up to update settings with user settings without changing the original
        if settings is None:
            # Defaults
            self.settings = processimages.default_settings()
        else:
            self.settings = settings.copy()

        # Child-specific Internals
        #    See ProcessImages


    def show(self, attr, display='ginga'):
        """
        Show one of the internal images

        Parameters
        ----------
        attr : str
          mspixelflat -- Show the combined flat image, unnormalized
          norm -- Show the combined normalized flat image
        display : str, optional

        Returns
        -------

        """
        if attr == 'mspixelflat':
            if self.mspixelflat is not None:
                ginga.show_image(self.mspixelflat)
        elif attr == 'norm':
            if self.mspixelflatnrm is not None:
                ginga.show_image(self.mspixelflatnrm)
