""" Uber object for calibration images, e.g. arc, flat """

import inspect

import os
import numpy as np

from abc import ABCMeta

from pypeit import msgs
from pypeit.par import pypeitpar
from pypeit.images import combineimage
from pypeit.images import pypeitimage
from pypeit.core import procimg

from IPython import embed


class ArcImage(pypeitimage.PypeItImage):
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
    file_format = 'fits'


class BiasImage(pypeitimage.PypeItImage):
    """
    Simple DataContainer for the Tilt Image
    """
    # Set the version of this class
    version = pypeitimage.PypeItImage.version

    # Output to disk
    output_to_disk = ('BIAS_IMAGE',)
    hdu_prefix = 'BIAS_'
    master_type = 'Bias'
    file_format = 'fits'


class TiltImage(pypeitimage.PypeItImage):
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
    file_format = 'fits'


class TraceImage(pypeitimage.PypeItImage):
    """
    Simple DataContainer for the Trace Image
    """

    # Peg the version of this class to that of PypeItImage
    version = pypeitimage.PypeItImage.version

    # I/O
    output_to_disk = ('TRACE_IMAGE', 'TRACE_FULLMASK', 'TRACE_DETECTOR')
    hdu_prefix = 'TRACE_'

    # Master fun
    master_type = 'Trace'
    frametype = 'trace'


'''
class BuildCalibrationImage(object):
    """
    Class to generate (and hold) a combined calibration image from a list of input images

    Args:
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            Spectrograph used to take the data.
        det (:obj:`int`):
            The 1-indexed detector number to process.
        proc_par (:class:`pypeit.par.pypeitpar.ProcessImagesPar`):
            Parameters that dictate the processing of the images.  See
            :class:`pypeit.par.pypeitpar.ProcessImagesPar` for the defaults.
        files (list):
            List of filenames to be processed and combined
        bias (`np.ndarray`, str, optional):
            Bias image or str or None describing how to bias subtract

    Attributes:
        pypeitImage (:class:`pypeit.images.pypeitimage.PypeItImage`):
        file_list (list):
            List of files to process
        process_steps (list):
            List of processing steps to be used

    """
    __metaclass__ = ABCMeta


    def __init__(self, spectrograph, det, frame_par, files, bias=None):

        # Required items
        self.spectrograph = spectrograph
        self.det = det

        # Assign the internal list of files
        self._set_files(files)

        # Required parameters
        if not isinstance(frame_par, pypeitpar.FrameGroupPar):
            msgs.error('Provided ParSet for must be type FrameGroupPar.')
        self.frame_par = frame_par  # This must be named this way as it is frequently a child

        # Init Process steps
        self.process_steps = procimg.set_process_steps(bias, self.frame_par)
        #self.process_steps += self.postbias_process_steps

        # Standard output
        self.pypeitImage = None

    @property
    def nfiles(self):
        """
        The number of calibration files
        """
        return len(self.file_list)

    def _set_files(self, files, check=False):
        """
        Assign the provided files to :attr:`files`.

        Args:
            files (None, :obj:`str`, :obj:`list`):
                The files to process.
            check (:obj:`bool`, optional):
                Check that the files exist on disk.

        Raises:
            PypeItError:
                Raised if the input objects have the wrong type.
        """
        if files is None:
            self.file_list = []
        elif isinstance(files, str):
            self.file_list = [files]
        elif isinstance(files, list):
            if not np.all([isinstance(f, str) for f in files]):
                msgs.error('File list elements must be strings.')
            self.file_list = files
        else:
            msgs.error('Provides files must be None, a string name, or a list of strings.')

        if check:
            self._check_files()

    def _check_files(self):
        """
        Check that files in :attr:`files` exist.

        Raises:
            PypeItError:
                Raised if any of the files don't exist.
        """
        for f in self.file_list:
            if not os.path.isfile(f):
                msgs.error('{0} does not exist!'.format(f))

    def build_image(self, bias=None, bpm=None, ignore_saturation=True):
        """
        Load, process and combine the images for a Calibration

        Args:
            bias (np.ndarray, optional):
                Bias image
            bpm (np.ndarray, optional):
                Bad pixel mask
            ignore_saturation (bool, optional):
                If True, turn off the saturation flag in the individual images before stacking
                This avoids having such values set to 0 which for certain images (e.g. flat calibrations)
                has unintended consequences.

        Returns:
            PypeItImage:
                Or a child (most likely)

        """
        combineImage = combineimage.CombineImage(self.spectrograph, self.det, self.frame_par['process'], self.file_list)
        pypeitImage = combineImage.run(self.process_steps, bias, bpm=bpm, ignore_saturation=ignore_saturation)

        # Decorate according to the type of calibration
        #   Primarily for handling MasterFrames
        #   WARNING, any internals in pypeitImage are lost here
        if self.frame_par['frametype'] == 'bias':
            finalImage = BiasImage.from_pypeitimage(pypeitImage)
        elif self.frame_par['frametype'] == 'arc':
            finalImage = ArcImage.from_pypeitimage(pypeitImage)
        elif self.frame_par['frametype'] == 'tilt':
            finalImage = TiltImage.from_pypeitimage(pypeitImage)
        else:
            embed(header='193 of calibrationimage')

        # Internals
        finalImage.process_steps = self.process_steps

        # Return
        return finalImage



    def __repr__(self):
        return ('<{:s}: nfiles={}, steps={}>'.format(
            self.__class__.__name__, self.nfiles, self.process_steps))
'''


def buildcalibrationimage(spectrograph, det, frame_par, files, bias=None, bpm=None, ignore_saturation=True):
    """

    Args:
        spectrograph:
        det:
        frame_par:
        files:
        bias:
        bpm:
        ignore_saturation:

    Returns:

    """

    if not isinstance(frame_par, pypeitpar.FrameGroupPar):
        msgs.error('Provided ParSet for must be type FrameGroupPar.')
    process_steps = procimg.set_process_steps(bias, frame_par)
    #
    combineImage = combineimage.CombineImage(spectrograph, det, frame_par['process'], files)
    pypeitImage = combineImage.run(process_steps, bias, bpm=bpm, ignore_saturation=ignore_saturation)
    #
    # Decorate according to the type of calibration
    #   Primarily for handling MasterFrames
    #   WARNING, any internals in pypeitImage are lost here
    if frame_par['frametype'] == 'bias':
        finalImage = BiasImage.from_pypeitimage(pypeitImage)
    elif frame_par['frametype'] == 'arc':
        finalImage = ArcImage.from_pypeitimage(pypeitImage)
    elif frame_par['frametype'] == 'tilt':
        finalImage = TiltImage.from_pypeitimage(pypeitImage)
    elif frame_par['frametype'] == 'trace':
        finalImage = TraceImage.from_pypeitimage(pypeitImage)
    elif frame_par['frametype'] in ['pixelflat']:
        finalImage = pypeitImage
    else:
        finalImage = None
        embed(header='193 of calibrationimage')

    # Internals
    finalImage.process_steps = process_steps

    # Return
    return finalImage

