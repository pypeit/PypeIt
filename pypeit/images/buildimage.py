""" Object to hold + process a single image"""

import inspect

import os
import numpy as np


from pypeit import msgs

from pypeit.core import coadd2d
from pypeit.par import pypeitpar
from pypeit import utils

from pypeit.images import pypeitimage
from pypeit.images import processrawimage
from pypeit.images import rawimage
from pypeit.images import maskimage

from IPython import embed


class BuildImage(object):
    """
    Class to generate an image from one or more files (and other pieces).

    The core processing steps are handled by ProcessRawImage
    This object is mainly for combining multiple images

    Args:
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            Spectrograph used to take the data.
        det (:obj:`int`, optional):
            The 1-indexed detector number to process.
        par (:class:`pypeit.par.pypeitpar.ProcessImagesPar`):
            Parameters that dictate the processing of the images.  See
            :class:`pypeit.par.pypeitpar.ProcessImagesPar` for the
            defaults.
        bpm (np.ndarray):
            Bad pixel mask.  Held in ImageMask

        frametype (str, optional): Frame type
        files (list, optional):
            List of filenames that went into the loaded image

    Attributes:
        ivar (np.narray):
            Inverse variance image
        rn2img (np.narray):
            Read noise**2 image
        filename (str):
            Required to build from a Raw image
    """
    def __init__(self, spectrograph, det, par, files):

        # Required parameters
        self.spectrograph = spectrograph
        self.det = det
        if not isinstance(par, pypeitpar.ProcessImagesPar):
            msgs.error('Provided ParSet for must be type ProcessImagesPar.')
        self.par = par  # This musts be named this way as it is frequently a child
        self.files = files

    def process_one(self, filename, process_steps, bias, pixel_flat, illum_flat=None, bpm=None):
        """
        Instantiate from a single file

        This will also generate the ivar, crmask, rn2img and mask

        Args:
            spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
                Spectrograph used to take the data.
            det (:obj:`int`, optional):
                The 1-indexed detector number to process.
            par (:class:`pypeit.par.pypeitpar.ProcessImagesPar`):
                Parameters that dictate the processing of the images.  See
                :class:`pypeit.par.pypeitpar.ProcessImagesPar` for the
                defaults.
            bpm (np.ndarray):
                Bad pixel mask.  Held in ImageMask
            filename (str):
                Filename
            bias (np.ndarray or None):
                Bias image
            pixel_flat (np.ndarray):
                Flat image
            illum_flat (np.ndarray, optional):
                Illumination image

        Returns:
            ScienceImage:

        """
        # Load raw image
        rawImage = rawimage.RawImage(filename, self.spectrograph, self.det)

        # Process
        processrawImage = processrawimage.ProcessRawImage(rawImage, self.par, bpm=bpm)
        processedImage = processrawImage.process(process_steps, bias=bias, pixel_flat=pixel_flat,
                                                 illum_flat=illum_flat)

        # Return
        return processedImage

    def run(self, process_steps, bias, pixel_flat=None, illum_flat=None,
                       sigma_clip=False, bpm=None, sigrej=None, maxiters=5):
        """
        Instantiate from file list

        This may also generate the ivar, crmask, rn2img and mask

        Args:
            process_steps (list):
            bpm (np.ndarray, optional):
                Bad pixel mask.  Held in ImageMask
            bias (np.ndarray or None):
                Bias image or instruction
            pixel_flat (np.ndarray, optional):
                Flat image
            illum_flat (np.ndarray, optional):
                Illumination image
            sigrej (int or float, optional): Rejection threshold for sigma clipping.
                 Code defaults to determining this automatically based on the number of images provided.
            maxiters (int, optional):

        Returns:
            ScienceImage:

        """
        # Loop on the files
        for kk, ifile in enumerate(self.files):
            # Process a single image
            pypeitImage = self.process_one(ifile, process_steps, bias, pixel_flat, illum_flat=illum_flat, bpm=bpm)
            # Are we all done?
            if len(self.files) == 1:
                return pypeitImage
            elif kk == 0:
                # Get ready
                nimages = len(self.files)
                shape = (nimages, pypeitImage.bpm.shape[0], pypeitImage.bpm.shape[1])
                img_stack = np.zeros(shape)
                ivar_stack= np.zeros(shape)
                rn2img_stack = np.zeros(shape)
                crmask_stack = np.zeros(shape, dtype=bool)
                # Mask
                bitmask = maskimage.ImageBitMask()
                mask_stack = np.zeros(shape, bitmask.minimum_dtype(asuint=True))
            # Process
            img_stack[kk,:,:] = pypeitImage.image
            # Construct raw variance image and turn into inverse variance
            if pypeitImage.ivar is not None:
                ivar_stack[kk, :, :] = pypeitImage.ivar
            else:
                ivar_stack[kk, :, :] = 1.
            # Mask cosmic rays
            if pypeitImage.crmask is not None:
                crmask_stack[kk, :, :] = pypeitImage.crmask
            # Read noise squared image
            if pypeitImage.rn2img is not None:
                rn2img_stack[kk, :, :] = pypeitImage.rn2img
            # Final mask for this image
            mask_stack[kk, :, :] = pypeitImage.mask

        # Coadd them
        weights = np.ones(nimages)/float(nimages)
        img_list = [img_stack]
        var_stack = utils.inverse(ivar_stack)
        var_list = [var_stack, rn2img_stack]
        img_list_out, var_list_out, outmask, nused = coadd2d.weighted_combine(
            weights, img_list, var_list, (mask_stack == 0),
            sigma_clip=sigma_clip, sigma_clip_stack=img_stack, sigrej=sigrej, maxiters=maxiters)

        # Build the last one
        #slf = ScienceImage.__init__(spectrograph, det, par, img_list_out[0],
        #                               utils.inverse(var_list_out[0]), bpm,
        #                               rn2img=var_list_out[1], crmask=np.invert(outmask), files=file_list)
        pypeitImage = pypeitimage.PypeItImage(img_list_out[0],
                                              ivar=utils.inverse(var_list_out[0]),
                                              bpm=pypeitImage.bpm,
                                              rn2img=var_list_out[1],
                                              crmask=np.invert(outmask))
        pypeitImage.build_mask(pypeitImage.image, pypeitImage.ivar,
                               saturation=self.spectrograph.detector[self.det-1]['saturation'],
                               mincounts=self.spectrograph.detector[self.det-1]['mincounts'])
        # Return
        return pypeitImage

    @property
    def nfiles(self):
        """
        Number of files in the files attribute

        Returns:
            int

        """
        if isinstance(self.files, list):
            return len(self.files)
        else:
            return 0


