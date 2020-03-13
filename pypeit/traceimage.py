"""
Module for generating the Trace image
"""
import numpy as np

from pypeit import msgs
from pypeit.par import pypeitpar
from pypeit.images import calibrationimage
from pypeit.images import pypeitimage
from pypeit.core import procimg
from IPython import embed

class TraceImage(calibrationimage.CalibrationImage):
    """
    Simple DataContainer for the Trace Image
    """

    # Peg the version of this class to that of PypeItImage
    version = pypeitimage.PypeItImage.version

    # I/O
    output_to_disk = ('TRACe_IMAGE', 'TRACE_FULLMASK', 'TRACE_DETECTOR')
    hdu_prefix = 'TRACE_'

    # Master fun
    master_type = 'Trace'
    frametype = 'trace'


class BuildTraceImage(calibrationimage.BuildCalibrationImage):
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

    def __init__(self, spectrograph, files=None, det=1, par=None, bias=None):

        self.par = pypeitpar.FrameGroupPar('trace') if par is None else par
        # Start us up
        calibrationimage.BuildCalibrationImage.__init__(self, spectrograph, det,
                                                        self.par['process'], files=files)
        # Processing steps
        self.process_steps = procimg.init_process_steps(bias, self.par['process'])
        self.process_steps += ['trim']
        self.process_steps += ['apply_gain']
        self.process_steps += ['orient']
        # TODO: CR masking for the trace images causes major issues with
        # the edge tracing because it identifies the edges of the slits
        # as cosmic rays.  CR parameters need to be optimized for this
        # to work.
#        self.process_steps += ['crmask']

