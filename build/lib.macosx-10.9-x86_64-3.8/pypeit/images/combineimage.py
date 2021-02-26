""" Class to generate an image from one or more files (and other pieces).

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import inspect

import os
import numpy as np


from pypeit import msgs

from pypeit.core import combine
from pypeit.par import pypeitpar
from pypeit import utils

from pypeit.images import pypeitimage
from pypeit.images import rawimage
from pypeit.images import imagebitmask

from IPython import embed


class CombineImage:
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

    def run(self, bias=None, flatimages=None, ignore_saturation=False, sigma_clip=True,
            bpm=None, sigrej=None, maxiters=5, slits=None, dark=None, combine_method='weightmean'):
        """
        Generate a PypeItImage from a list of images

        This may also generate the ivar, crmask, rn2img and mask

        Args:
            bias (:class:`pypeit.images.buildimage.BiasImage`, optional): Bias image
            flatimages (:class:`pypeit.flatfield.FlatImages`, optional):  For flat fielding
            dark (:class:`pypeit.images.buildimage.DarkImage`, optional): Dark image
            slits (:class:`pypeit.slittrace.SlitTraceSet`, optional): Slit object
            sigma_clip (bool, optional):
                Perform sigma clipping
            sigrej (tuple, optional): Rejection thresholds (lo/hi) for sigma clipping.
                 Code defaults to determining this automatically based on the number of images provided.
            maxiters (int, optional):
                Number of iterations for the clipping
            bpm (`numpy.ndarray`_, optional):
                Bad pixel mask.  Held in ImageMask
            ignore_saturation (:obj:`bool`, optional):
                If True, turn off the saturation flag in the individual images before stacking
                This avoids having such values set to 0 which for certain images (e.g. flat calibrations)
                can have unintended consequences.
            combine_method (str):
                Method to combine images
                Allowed options are 'weightmean', 'median'

        Returns:
            :class:`pypeit.images.pypeitimage.PypeItImage`:

        """
        # Loop on the files
        nimages = len(self.files)
        lampstat = []
        for kk, ifile in enumerate(self.files):
            # Load raw image
            rawImage = rawimage.RawImage(ifile, self.spectrograph, self.det)
            # Process
            pypeitImage = rawImage.process(self.par, bias=bias, bpm=bpm, dark=dark,
                                           flatimages=flatimages, slits=slits)
            #embed(header='96 of combineimage')
            # Are we all done?
            if nimages == 1:
                return pypeitImage
            elif kk == 0:
                # Get ready
                shape = (nimages, pypeitImage.image.shape[0], pypeitImage.image.shape[1])
                img_stack = np.zeros(shape)
                ivar_stack= np.zeros(shape)
                rn2img_stack = np.zeros(shape)
                crmask_stack = np.zeros(shape, dtype=bool)
                # Mask
                bitmask = imagebitmask.ImageBitMask()
                mask_stack = np.zeros(shape, bitmask.minimum_dtype(asuint=True))
            # Grab the lamp status
            lampstat += [self.spectrograph.get_lamps_status(pypeitImage.rawheadlist)]
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
                indx = pypeitImage.bitmask.flagged(pypeitImage.fullmask, flag=['SATURATION'])
                pypeitImage.fullmask[indx] = pypeitImage.bitmask.turn_off(
                    pypeitImage.fullmask[indx], 'SATURATION')
            mask_stack[kk, :, :] = pypeitImage.fullmask

        # Check that the lamps being combined are all the same:
        if not lampstat[1:] == lampstat[:-1]:
            msgs.warn("The following files contain different lamp status")
            # Get the longest strings
            maxlen = max([len("Filename")]+[len(os.path.split(x)[1]) for x in self.files])
            maxlmp = max([len("Lamp status")]+[len(x) for x in lampstat])
            strout = "{0:" + str(maxlen) + "}  {1:s}"
            # Print the messages
            print(msgs.indent() + '-'*maxlen + "  " + '-'*maxlmp)
            print(msgs.indent() + strout.format("Filename", "Lamp status"))
            print(msgs.indent() + '-'*maxlen + "  " + '-'*maxlmp)
            for ff, file in enumerate(self.files):
                print(msgs.indent() + strout.format(os.path.split(file)[1], " ".join(lampstat[ff].split("_"))))
            print(msgs.indent() + '-'*maxlen + "  " + '-'*maxlmp)

        # Coadd them
        weights = np.ones(nimages)/float(nimages)
        img_list = [img_stack]
        var_stack = utils.inverse(ivar_stack)
        var_list = [var_stack, rn2img_stack]
        if combine_method == 'weightmean':
            img_list_out, var_list_out, gpm, nused = combine.weighted_combine(
                weights, img_list, var_list, (mask_stack == 0),
                sigma_clip=sigma_clip, sigma_clip_stack=img_stack, sigrej=sigrej, maxiters=maxiters)
        elif combine_method == 'median':
            img_list_out = [np.median(img_stack, axis=0)]
            var_list_out = [np.median(var_stack, axis=0)]
            var_list_out += [np.median(rn2img_stack, axis=0)]
            gpm = np.ones_like(img_list_out[0], dtype='bool')
        else:
            msgs.error("Bad choice for combine.  Allowed options are 'median', 'weightmean'.")

        # Build the last one
        final_pypeitImage = pypeitimage.PypeItImage(img_list_out[0],
                                                    ivar=utils.inverse(var_list_out[0]),
                                                    bpm=pypeitImage.bpm,
                                                    rn2img=var_list_out[1],
                                                    crmask=np.logical_not(gpm),
                                                    detector=pypeitImage.detector,
                                                    PYP_SPEC=pypeitImage.PYP_SPEC)
        # Internals
        final_pypeitImage.rawheadlist = pypeitImage.rawheadlist
        final_pypeitImage.process_steps = pypeitImage.process_steps

        nonlinear_counts = self.spectrograph.nonlinear_counts(pypeitImage.detector,
                                                              apply_gain=self.par['apply_gain'])
        final_pypeitImage.build_mask(saturation=nonlinear_counts)
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

