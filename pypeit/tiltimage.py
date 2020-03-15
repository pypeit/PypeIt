"""
Module for generating the Tilt image.

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""
import os
import inspect
import numpy as np

from pypeit import msgs
from pypeit import masterframe
from pypeit.par import pypeitpar
from pypeit.images import calibrationimage
from pypeit.images import pypeitimage
from pypeit.core import procimg

from IPython import embed


class TiltImage(calibrationimage.CalibrationImage):
    """
    Simple DataContainer for the Tilt Image
    """

    # Peg the version of this class to that of PypeItImage
    version = pypeitimage.PypeItImage.version

    # I/O
    output_to_disk = ('TILT_IMAGE', 'TILT_FULLMASK', 'TILT_DETECTOR')
    hdu_prefix = 'TILT_'

    # Master fun
    master_type = 'Tilt'
    frametype = 'tilt'


class BuildTiltImage(calibrationimage.BuildCalibrationImage):
    """
    Generate a Tilt Image by processing and combining one or more tilt frames.

    See :class:`pypeit.images.BuildCalibrationImage` for the __init__
    """

    # Frametype is a class attribute
    frametype = 'tilt'
    image_type = TiltImage

    postbias_process_steps = ['trim']
    postbias_process_steps += ['orient']
    postbias_process_steps += ['apply_gain']

