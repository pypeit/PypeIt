""" Object to hold + process a single image"""

import inspect

import os
import numpy as np


from pypeit import msgs

from pypeit import utils
from pypeit.core import procimg
from pypeit.core import combine
from pypeit.par import pypeitpar

from pypeit.images import pypeitimage
from pypeit.images import processimage


from IPython import embed

# REMOVE THIS
from importlib import reload
reload(procimg)



class CombinedImage(pypeitimage.PypeItImage):
    """
    Class to generate a combined image from a list of input images or
    simply to hold a previously generated combined image.

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
        pimages (list): List of ProcessImage objects
        image (np.ndarray):
        file_list (list): List of files to process

    """


    def __init__(self, spectrograph, det, proc_par, files=None, frametype=None):

        # Init me
        pypeitimage.PypeItImage.__init__(self, spectrograph, det)

        # Assign the internal list of files
        self._set_files(files)

        # Required parameters
        if not isinstance(proc_par, pypeitpar.ProcessImagesPar):
            msgs.error('Provided ParSet for must be type ProcessImagesPar.')
        self.proc_par = proc_par  # This musts be named this way as it is frequently a child

        # Optional parameters
        self.frametype = frametype

        # Internal images
        self.pimages = []
        self.image = None

    @property
    def nimages(self):
        """

        Returns:
            int: Number of ProcessImage objects loaded

        """
        return len(self.pimages)

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

    def build_stack(self, bpm=None, reject_cr=False):
        """
        Generate a set of stack arrays useful for image processing

        These are generated from the internal ProcessImage objects
        held in self.pimages which need to have been previously
        generated/loaded

        Note:  To save on memory, the image attributes of each
        ProcessImage in self.pimages is reset after it contributes
        to the stack arrays.

        Args:
            bpm (np.ndarray, optional):
                Bad pixel mask image
            reject_cr (bool, optional):
                If true, generate a CR image

        Returns:
            tuple: np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray
                sciimg_stack, sciivar_stack, rn2img_stack, crmask_stack, mask_stack

        """

        # For now, this winds up getting 2 copies for everything.
        #  One in the self.pimages list and one in the stacks

        # Get it ready
        shape = (self.nimages, self.pimages[0].image.shape[0], self.pimages[0].image.shape[1])
        sciimg_stack = np.zeros(shape)
        sciivar_stack= np.zeros(shape)
        rn2img_stack = np.zeros(shape)
        crmask_stack = np.zeros(shape, dtype=bool)
        mask_stack = np.zeros(shape, self.pimages[0].bitmask.minimum_dtype(asuint=True))

        # Loop on the ProcessImage objects
        for kk, pimage in enumerate(self.pimages):
            # Construct raw variance image and turn into inverse variance
            rawvarframe = pimage.build_rawvarframe()
            sciivar_stack[kk, :, :] = utils.calc_ivar(rawvarframe)
            # Mask cosmic rays
            if reject_cr:
                crmask_stack[kk, :, :] = pimage.build_crmask()
            sciimg_stack[kk,:,:] = pimage.image
            # Build read noise squared image
            rn2img_stack[kk, :, :] = pimage.build_rn2img()
            # Final mask for this image
            mask_stack[kk, :, :] = pimage.build_mask(
                bpm=bpm,
                saturation=self.spectrograph.detector[self.det - 1]['saturation'],
                mincounts = self.spectrograph.detector[self.det - 1]['mincounts'])

            # Delete/reset the image internals now
            pimage._reset_internals()

        # Return
        return sciimg_stack, sciivar_stack, rn2img_stack, crmask_stack, mask_stack

    def combine(self):
        """
        Combine the images held in self.pimages

        Wrapper to combine.comb_frames

        Returns:
            np.ndarray: Copy of self.image

        """
        if self.nimages == 1:
            self.image = self.pimages[0].image
        else:
            #
            # Build the image stack
            image_arr = np.zeros((self.pimages[0].image.shape[0],
                                         self.pimages[0].image.shape[1],
                                         self.nimages))
            # This is a bit memory expensive...
            for kk,iimage in enumerate(self.pimages):
                image_arr[:,:,kk] = iimage.image

            # Do it
            self.image = combine.comb_frames(image_arr,
                                             frametype=self.frametype,
                                             saturation=self.spectrograph.detector[self.det-1]['saturation'],
                                             method=self.proc_par['combine'],
                                             satpix=self.proc_par['satpix'],
                                             cosmics=self.proc_par['sigrej'],
                                             n_lohi=self.proc_par['n_lohi'],
                                             sig_lohi=self.proc_par['sig_lohi'],
                                             replace=self.proc_par['replace'])
        # Return
        return self.image.copy()

    def load_images(self, reload=False):
        """
        Load up the images in self.file_list into a
        list of ProcessImages in self.pimages

        Args:
            reload (bool, optional):

        """
        if (not reload) and (self.nimages > 0):
            msgs.warn("Images already loaded.  Use reload if you wish")
            return
        for file in self.file_list:
            # Instantiate
            processImage = processimage.ProcessImage(file, self.spectrograph, self.det, self.proc_par)
            # Load
            processImage.load_rawimage(file)
            # Append
            self.pimages.append(processImage)

    def process_images(self, process_steps, pixel_flat=None, illum_flat=None,
                       bias=None, bpm=None):
        """
        Process the images held in self.pimages

        Note:  The processing steps are currently 'frozen' as is.
          We may choose to allow optional ordering of the steps

        Args:
            process_steps (list):
                List of processing steps
            pixel_flat (np.ndarray, optional):
                Pixel flat image
            illum_flat (np.ndarray, optional):
                Illumination flat
            bias (np.ndarray, optional):
                Bias image
            bpm (np.ndarray, optional):
                Bad pixel mask image

        """

        for pimage in self.pimages:
            # Standard order
            #   -- May need to allow for other order some day..
            if 'subtract_bias' in process_steps:
                pimage.subtract_bias(bias)
            if 'subtract_overscan' in process_steps:
                pimage.subtract_overscan()
            if 'trim' in process_steps:
                pimage.trim()
            if 'apply_gain' in process_steps:
                pimage.apply_gain()
            # Always orient
            pimage.orient()
            # Flat field
            if 'flatten' in process_steps:
                #embed(header='190 of combinedimage')
                pimage.flatten(pixel_flat, illum_flat=illum_flat, bpm=bpm)






