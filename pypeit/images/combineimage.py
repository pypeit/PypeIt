""" Class to generate an image from one or more files (and other pieces).
"""

import inspect

import os
import numpy as np


from pypeit import msgs

from pypeit.core import combine
from pypeit.par import pypeitpar
from pypeit import utils

from pypeit.images import pypeitimage
from pypeit.images import processrawimage
from pypeit.images import rawimage
from pypeit.images import maskimage

from IPython import embed


class CombineImage(object):
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

    """
    def __init__(self, spectrograph, det, par, files):

        # Required parameters
        self.spectrograph = spectrograph
        self.det = det
        if not isinstance(par, pypeitpar.ProcessImagesPar):
            msgs.error('Provided ParSet for must be type ProcessImagesPar.')
        self.par = par  # This musts be named this way as it is frequently a child
        self.files = files
        if self.nfiles == 0:
            msgs.error('Combineimage requires a list of files to instantiate')

    def process_one(self, filename, process_steps, bias, pixel_flat=None, illum_flat=None, bpm=None):
        """
        Process a single image

        Args:
            filename (str):
                File to process
            process_steps (list):
                List of processing steps
            bias (np.ndarray or None):
                Bias image
            pixel_flat (np.ndarray, optional):
                Flat image
            illum_flat (np.ndarray, optional):
                Illumination image
            bpm (np.ndarray, optional):
                Bad pixel mask

        Returns:
            :class:`pypeit.images.pypeitimage.PypeItImage`:

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
            ignore_saturation=False, sigma_clip=True, bpm=None, sigrej=None, maxiters=5):
        """
        Generate a PypeItImage from a list of images

        Mainly a wrapper to coadd2d.weighted_combine()

        This may also generate the ivar, crmask, rn2img and mask

        Args:
            process_steps (list):
            bias (np.ndarray or None):
                Bias image or instruction
            pixel_flat (np.ndarray, optional):
                Flat image
            illum_flat (np.ndarray, optional):
                Illumination image
            sigma_clip (bool, optional):
                Perform sigma clipping
            sigrej (int or float, optional): Rejection threshold for sigma clipping.
                 Code defaults to determining this automatically based on the number of images provided.
            maxiters (int, optional):
                Number of iterations for the clipping
            bpm (np.ndarray, optional):
                Bad pixel mask.  Held in ImageMask
            ignore_saturation (bool, optional):
                If True, turn off the saturation flag in the individual images before stacking
                This avoids having such values set to 0 which for certain images (e.g. flat calibrations)
                can have unintended consequences.

        Returns:
            :class:`pypeit.images.pypeitimage.PypeItImage`:

        """
        # Loop on the files
        nimages = len(self.files)
        for kk, ifile in enumerate(self.files):
            # Process a single image
            pypeitImage = self.process_one(ifile, process_steps, bias, pixel_flat=pixel_flat,
                                           illum_flat=illum_flat, bpm=bpm)
            # Are we all done?
            if len(self.files) == 1:
                return pypeitImage
            elif kk == 0:
                # Get ready
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
            # TODO This seems kludgy to me. Why not just pass ignore_saturation to process_one and ignore the saturation
            # when the mask is actually built, rather than untoggling the bit here
            if ignore_saturation:  # Important for calibrations as we don't want replacement by 0
                indx = pypeitImage.bitmask.flagged(pypeitImage.mask, flag=['SATURATION'])
                pypeitImage.mask[indx] = pypeitImage.bitmask.turn_off(pypeitImage.mask[indx], 'SATURATION')
            mask_stack[kk, :, :] = pypeitImage.mask

        # Coadd them
        weights = np.ones(nimages)/float(nimages)
        img_list = [img_stack]
        var_stack = utils.inverse(ivar_stack)
        var_list = [var_stack, rn2img_stack]
        img_list_out, var_list_out, outmask, nused = combine.weighted_combine(
            weights, img_list, var_list, (mask_stack == 0),
            sigma_clip=sigma_clip, sigma_clip_stack=img_stack, sigrej=sigrej, maxiters=maxiters)

        # Build the last one
        final_pypeitImage = pypeitimage.PypeItImage(img_list_out[0],
                                                    ivar=utils.inverse(var_list_out[0]),
                                                    bpm=pypeitImage.bpm,
                                                    rn2img=var_list_out[1],
                                                    crmask=np.invert(outmask),
                                                    binning=pypeitImage.binning)
        nonlinear_counts = self.spectrograph.nonlinear_counts(self.det,
                                                              apply_gain='apply_gain' in process_steps)
        final_pypeitImage.build_mask(final_pypeitImage.image, final_pypeitImage.ivar,
                               saturation=nonlinear_counts, #self.spectrograph.detector[self.det-1]['saturation'],
                               mincounts=self.spectrograph.detector[self.det-1]['mincounts'])
        # Return
        return final_pypeitImage

    @property
    def nfiles(self):
        """
        Number of files in the files attribute

        Returns:
            int

        """
        return len(self.files) if isinstance(self.files, (np.ndarray, list)) else 0


