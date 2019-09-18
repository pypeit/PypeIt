""" Uber object for calibration images, e.g. arc, flat """

import inspect

import os
import numpy as np


from pypeit import msgs
from pypeit.par import pypeitpar
from pypeit.images import combineimage

from IPython import embed


class CalibrationImage(object):
    """
    Class to generate (and hold) a combined calibration image from a list of input images or
    simply to hold a previously generated calibration image.

    Args:
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            Spectrograph used to take the data.
        det (:obj:`int`, optional):
            The 1-indexed detector number to process.
        proc_par (:class:`pypeit.par.pypeitpar.ProcessImagesPar` or None):
            Parameters that dictate the processing of the images.  See
            :class:`pypeit.par.pypeitpar.ProcessImagesPar` for the
            defaults.

        files (list, optional):
            List of filenames to be combined

    Attributes:
        pypeitImage (:class:`pypeit.images.pypeitimage.PypeItImage`):
        file_list (list):
            List of files to process
        process_steps (list):
            List of processing steps to be used

    """
    def __init__(self, spectrograph, det, proc_par, files=None):

        # Required items
        self.spectrograph = spectrograph
        self.det = det

        # Assign the internal list of files
        self._set_files(files)

        # Required parameters
        if proc_par is not None:
            if not isinstance(proc_par, pypeitpar.ProcessImagesPar):
                msgs.error('Provided ParSet for must be type ProcessImagesPar.')
        self.proc_par = proc_par  # This must be named this way as it is frequently a child

        # Process steps
        self.process_steps = []

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

        """
        combineImage = combineimage.CombineImage(self.spectrograph, self.det, self.proc_par, self.file_list)
        self.pypeitImage = combineImage.run(self.process_steps, bias, bpm=bpm, ignore_saturation=ignore_saturation)
        # Return
        return self.pypeitImage


    def __repr__(self):
        return ('<{:s}: nfiles={}, steps={}>'.format(
            self.__class__.__name__, self.nfiles, self.process_steps))





