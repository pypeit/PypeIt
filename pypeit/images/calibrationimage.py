""" Object to hold + process a single image"""

import inspect

import os
import numpy as np


from pypeit import msgs

from pypeit.core import combine
from pypeit.par import pypeitpar

from pypeit.images import pypeitimage
from pypeit.images import processrawimage


from IPython import embed


class CalibrationImage(pypeitimage.PypeItImage):
    """
    Class to generate a combined calibration image from a list of input images or
    simply to hold a previously generated calibration image.

    Args:
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            Spectrograph used to take the data.
        det (:obj:`int`, optional):
            The 1-indexed detector number to process.
        proc_par (:class:`pypeit.par.pypeitpar.ProcessImagesPar`):
            Parameters that dictate the processing of the images.  See
            :class:`pypeit.par.pypeitpar.ProcessImagesPar` for the
            defaults.

        files (list, optional):
            List of filenames to be combined
        frametype (str, optional): Frame type

    Attributes:
        image (np.ndarray):
        file_list (list): List of files to process
        process_steps (list): List of processing steps to be used

    """


    def __init__(self, spectrograph, det, proc_par, files=None):

        # Init me
        pypeitimage.PypeItImage.__init__(self, spectrograph, det)

        # Assign the internal list of files
        self._set_files(files)

        # Required parameters
        if not isinstance(proc_par, pypeitpar.ProcessImagesPar):
            msgs.error('Provided ParSet for must be type ProcessImagesPar.')
        self.proc_par = proc_par  # This musts be named this way as it is frequently a child

        # Internal images
        self.image = None

        # Process steps
        self.process_steps = []

    @property
    def shape(self):
        return () if self.image is None else self.image.shape

    @property
    def nfiles(self):
        """

        Returns:
            int: Number of files in the file_list

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

    def build_image(self, bias=None, bpm=None):
        """
        Load, process and combine the images for a Calibration

        Args:
            bias (np.ndarray, optional):
                Bias image
            bpm (np.ndarray, optional):
                Bad pixel mask

        Returns:
            np.ndarray: Copy of self.image

        """
        if self.nfiles == 0:
            msgs.warn("Need to provide a non-zero list of files")
            return
        # Load up image array
        image_arr = None
        for kk,file in enumerate(self.file_list):
            # Process raw file
            processrawImage = processrawimage.ProcessRawImage(file, self.spectrograph,
                                                           self.det, self.proc_par)
            image = processrawImage.process(self.process_steps, bias=bias, bpm=bpm)
            # Instantiate the image stack
            if image_arr is None:
                image_arr = np.zeros((image.shape[0], image.shape[1], self.nfiles))
            # Hold
            image_arr[:,:,kk] = image

        # Combine
        if self.nfiles == 1:
            self.image = image_arr[:,:,0]
        else:
            self.image = combine.comb_frames(image_arr,
                                         saturation=self.spectrograph.detector[self.det-1]['saturation'],
                                         method=self.proc_par['combine'],
                                         satpix=self.proc_par['satpix'],
                                         cosmics=self.proc_par['sigrej'],
                                         n_lohi=self.proc_par['n_lohi'],
                                         sig_lohi=self.proc_par['sig_lohi'],
                                         replace=self.proc_par['replace'])
        # Return
        return self.image.copy()

    def __repr__(self):
        return ('<{:s}: nfiles={}, steps={}>'.format(
            self.__class__.__name__, self.nfiles, self.process_steps))





