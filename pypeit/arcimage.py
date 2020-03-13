"""
Module for generating the Arc image.

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""
import os

from pypeit.par import pypeitpar
from pypeit.images import calibrationimage
from pypeit.images import pypeitimage
from pypeit.core import procimg

from IPython import embed


class ArcImage(calibrationimage.CalibrationImage):
    """
    Simple DataContainer for the Arc Image
    """
    # Peg the version of this class to that of PypeItImage
    version = pypeitimage.PypeItImage.version

    # I/O
    output_to_disk = ('ARC_IMAGE', 'ARC_FULLMASK', 'ARC_DETECTOR')
    hdu_prefix = 'ARC_'

    # Master fun
    master_type = 'Arc'
    frametype = 'arc'


class BuildArcImage(calibrationimage.BuildCalibrationImage):
    """
    Generate an ArcImage by processing and combining one or more arc frames.


    Args:
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            The `Spectrograph` instance that sets the
            instrument used to take the observations.  Used to set
            :attr:`spectrograph`.
        files (:obj:`list`, optional):
            The list of files to process.
            Can be an empty list or None
        det (:obj:`int`, optional):
            The 1-indexed detector number to process.
        par (:class:`pypeit.par.pypeitpar.FrameGroupPar`):
            The parameters used to type and process the arc frames.
        msbias (ndarray or str, optional):
            Guides bias subtraction

    Attributes:
        msbias (ndarray):
            Bias image or bias-subtraction method; see
            :func:`pypeit.processimages.ProcessImages.process`.
    """
    # Frametype is a class attribute
    frametype = 'arc'
    image_type = ArcImage

    def __init__(self, spectrograph, files=None, det=1, par=None, msbias=None):
    
        # Parameters unique to this Object
        self.msbias = msbias

        # Parameters
        self.par = pypeitpar.FrameGroupPar(self.frametype) if par is None else par

        # Start us up
        calibrationimage.BuildCalibrationImage.__init__(self, spectrograph, det, self.par['process'], files=files)

        # Process steps
        self.process_steps = procimg.init_process_steps(self.msbias, self.par['process'])
        self.process_steps += ['trim']
        self.process_steps += ['orient']
        self.process_steps += ['apply_gain']

