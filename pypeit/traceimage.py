# Module for generating the Trace image
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np

from pypeit import msgs
from pypeit import processimages
from pypeit.par import pypeitpar

from pypeit import debugger


class TraceImage(processimages.ProcessImages):
    """
    Generate the image for tracing the slits/orders, typically from a
    set of flat-field or twilight sky exposures

    Args:
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            The `Spectrograph` instance that sets the
            instrument used to take the observations.  Used to set
            :attr:`spectrograph`.
        files (list, optional): List of filenames
        det (:obj:`int`, optional):
            The 1-indexed detector number to process.
        par (:class:`pypeit.par.pypeitpar.FrameGroupPar`):
            The parameters used to type and process the arc frames.

    Attributes:
        frametype (str): Set to 'trace_image'

    """
   
    # Frametype is a class attribute
    frametype = 'trace_image'

    def __init__(self, spectrograph, files=None, det=1, par=None):
        self.par = pypeitpar.FrameGroupPar('trace') if par is None else par
        processimages.ProcessImages.__init__(self, spectrograph, self.par['process'],
                                             files=files, det=det)


