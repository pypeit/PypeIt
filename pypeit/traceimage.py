"""
Module for generating the Trace image
"""
import numpy as np

from pypeit import msgs
from pypeit.par import pypeitpar
from pypeit.images import combinedimage
from pypeit.core import procimg


class TraceImage(combinedimage.CombinedImage):
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
    """
    # Frametype is a class attribute
    frametype = 'trace_image'

    def __init__(self, spectrograph, files=None, det=1, par=None):
        self.par = pypeitpar.FrameGroupPar('trace') if par is None else par
        # Start us up
        combinedimage.CombinedImage.__init__(self, spectrograph, det, self.par['process'],
                                             files=files, frametype=self.frametype)

    def build_image(self, bias=None):
        # Load
        self.load_images()
        # Process + combine
        process_steps = procimg.init_process_steps(bias, self.proc_par)
        process_steps += ['trim']
        process_steps += ['apply_gain']
        self.process_images(process_steps, bias=bias)
        self.image = self.combine()
