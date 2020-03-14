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
    output_to_disk = ('ARC_IMAGE', 'ARC_FULLMASK', 'ARC_DETECTOR')# 'ARC_DETECTOR_CONTAINER')
    hdu_prefix = 'ARC_'

    # Master fun
    master_type = 'Arc'
    frametype = 'arc'


class BuildArcImage(calibrationimage.BuildCalibrationImage):
    """
    Generate an ArcImage by processing and combining one or more arc frames.

    See :class:`pypeit.images.BuildCalibrationImage` for the __init__
    """
    # Define the processing steps *after* bias/overscan subtraction
    postbias_process_steps = ['trim']
    postbias_process_steps += ['orient']
    postbias_process_steps += ['apply_gain']

    image_type = ArcImage
