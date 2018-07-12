# Module for generating the Trace image
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np

from pypit import msgs
from pypit import processimages
from pypit.par import pypitpar

from pypit import ardebug as debugger

# TODO: This syntax for instantiaion differs slightlhy from ArcIMage,
# and I think they need to be homogenized.

# TODO: I don't understand why TraceImage is its own class, since the
# traceimg is going to be the same as the flat image but the Flatimage
# does not have its own class but is instead embedded in the Flatfield
# class. Maybe the idea is to have the flexibilty to have the traceimage
# not be the flat image? But this could also apply to the tilts, since I
# might want to use arcs for the wavelength solution but use the science
# data itself for the traceimage. I think these image classes are just
# adding a bunch of extra code for no reason. I get that the motivation
# is to avoid making the msarc twice but perhaps we can pay that
# computational cost (which is small) for simpler code.

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
   
    # Frametype is a class attribute
    frametype = 'trace_image'

    def __init__(self, spectrograph, file_list=[], det=1, par=None):
        self.par = pypitpar.FrameGroupPar('trace') if par is None else par
        processimages.ProcessImages.__init__(self, spectrograph, file_list=file_list, det=det,
                                             overscan_par=self.par['overscan'],
                                             combine_par=self.par['combine'],
                                             lacosmic_par=self.par['lacosmic'])


