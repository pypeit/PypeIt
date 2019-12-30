""" Object to hold + process a single image.
This module also includes the build_from_list() method
which is how the ScienceImage is most frequently generated. """

import inspect

import os
import numpy as np


from pypeit import msgs

from pypeit.core import procimg
from pypeit.par import pypeitpar
from pypeit import utils

from pypeit.images import pypeitimage
from pypeit.images import combineimage

from IPython import embed


class ScienceImage(pypeitimage.PypeItImage):
    """
    Class to generate and hold a science image

    Child of PypeItImage

    Args:
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            Spectrograph used to take the data.
        det (:obj:`int`):
            The 1-indexed detector number to process.
        par (:class:`pypeit.par.pypeitpar.ProcessImagesPar`):
            Parameters that dictate the processing of the images.  See
            :class:`pypeit.par.pypeitpar.ProcessImagesPar` for the
            defaults.
        image (np.ndarray):
        ivar (np.ndarray):
        bpm (np.ndarray):
            Bad pixel mask.  Held in ImageMask
        rn2img (np.ndarray, optional):
        crmask (np.ndarray, optional):
        mask (np.ndarray, optional):
        files (list, optional):
            List of filenames that went into the loaded image

    """
    frametype = 'science'

    def __init__(self, spectrograph, det, par, image, ivar, bpm, rn2img=None,
                 crmask=None, mask=None, files=[]):

        # Init me
        pypeitimage.PypeItImage.__init__(self, image, ivar=ivar, rn2img=rn2img,
                                         bpm=bpm, crmask=crmask, mask=mask)

        # Required attribs
        self.spectrograph = spectrograph
        if not isinstance(par, pypeitpar.ProcessImagesPar):
            msgs.error('Provided ParSet for must be type ProcessImagesPar.')
        self.par = par
        self.det = det

        # Not required
        self.files = files

    def build_crmask(self, subtract_img=None):
        """
        Call to ImageMask.build_crmask which will
        generate the cosmic ray mask

        Args:
            subtract_img (np.ndarray, optional):
                Image to be subtracted off of self.image prior to CR evaluation

        Returns:
            np.ndarray: Boolean array of self.crmask

        """
        return super(ScienceImage, self).build_crmask(self.spectrograph, self.det,
                                                      self.par, self.image,
                                                      utils.inverse(self.ivar),
                                                      subtract_img=subtract_img).copy()

    def build_mask(self, saturation=1e10, mincounts=-1e10, slitmask=None):
        """
        Call to ImageMask.build_mask()
        This generates the full Image mask

        Args:
            saturation (float, optional):
            mincounts (float, optional):
            slitmask (np.ndarray, optional):

        Returns:
            np.ndarray:  The full mask, held in self.mask

        """
        super(ScienceImage, self).build_mask(self.image, self.ivar,
                                             saturation=saturation,
                                             mincounts=mincounts,
                                             slitmask=slitmask)
        return self.mask.copy()

    def update_mask_cr(self, subtract_img=None):
        """
        Updates the CR mask values in self.mask
        through a call to ImageMask.build_crmask which
        generates a new CR mask and then a call to
        ImageMask.update_mask_cr() which updates self.mask

        Args:
            subtract_img (np.ndarray, optional):
                If provided, this is subtracted from self.image prior to
                CR masking
        """
        # Generate the CR mask (and save in self.crmask)
        super(ScienceImage, self).build_crmask(self.spectrograph, self.det,
                                               self.par, self.image,
                                               utils.inverse(self.ivar),
                                               subtract_img=subtract_img).copy()
        # Now update the mask
        super(ScienceImage, self).update_mask_cr(self.crmask)

    def __sub__(self, other):
        """
        Subtract a ScienceImage object from another
        Extras (e.g. ivar, masks) are included if they are present

        Args:
            other (ScienceImage):

        Returns:
            ScienceImage:

        """
        if not isinstance(other, ScienceImage):
            msgs.error("Misuse of the subtract method")
        # Images
        newimg = self.image - other.image

        # Mask time
        outmask_comb = (self.mask == 0) & (other.mask == 0)

        # Variance
        if self.ivar is not None:
            new_ivar = utils.inverse(utils.inverse(self.ivar) + utils.inverse(other.ivar))
            new_ivar[np.invert(outmask_comb)] = 0
        else:
            new_ivar = None

        # RN2
        if self.rn2img is not None:
            new_rn2 = self.rn2img + other.rn2img
        else:
            new_rn2 = None

        # Files
        new_files = self.files + other.files

        # Instantiate
        new_sciImg = ScienceImage(self.spectrograph, self.det, self.par,
            newimg, new_ivar, self.bpm, rn2img=new_rn2, files=new_files)
        #TODO: KW properly handle adding the bits
        crmask_diff = new_sciImg.build_crmask()
        # crmask_eff assumes evertything masked in the outmask_comb is a CR in the individual images
        new_sciImg.crmask = crmask_diff | np.invert(outmask_comb)
        # Note that the following uses the saturation and mincounts held in
        # self.spectrograph.detector[self.det-1]
        new_sciImg.build_mask()

        return new_sciImg

    def __repr__(self):
        repr = '<{:s}: files={}'.format(self.__class__.__name__, self.files)
        # Image
        rdict = {}
        for attr in ['image', 'ivar', 'rn2img', 'crmask', 'mask']:
            if getattr(self, attr) is not None:
                rdict[attr] = True
            else:
                rdict[attr] = False
        repr += ' images={}'.format(rdict)
        repr = repr + '>'
        return repr


def build_from_file_list(spectrograph, det, par, bpm,
                   file_list, bias, pixel_flat, illum_flat=None,
                   sigma_clip=False, sigrej=None, maxiters=5):
    """
    Build a ScienceImage from a file list
    using a default set of process steps

    This will also generate the ivar, crmask, rn2img and mask

    Args:
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            Spectrograph used to take the data.
        det (:obj:`int`):
            The 1-indexed detector number to process.
        par (:class:`pypeit.par.pypeitpar.ProcessImagesPar`):
            Parameters that dictate the processing of the images.  See
            :class:`pypeit.par.pypeitpar.ProcessImagesPar` for the
            defaults.
        bpm (np.ndarray):
            Bad pixel mask.  Held in ImageMask
        file_list (list):
            List of files
        bias (np.ndarray or None):
            Bias image
        pixel_flat (np.ndarray):
            Flat image
        illum_flat (np.ndarray, optional):
            Illumination image
        sigrej (int or float, optional): Rejection threshold for sigma clipping.
             Code defaults to determining this automatically based on the numberr of images provided.
        maxiters (int, optional):

    Returns:
        ScienceImage:

    """
    # Process steps
    process_steps = procimg.init_process_steps(bias, par)
    process_steps += ['trim', 'apply_gain', 'orient']
    if (pixel_flat is not None) or (illum_flat is not None):
        process_steps += ['flatten']
    process_steps += ['extras']
    if par['cr_reject']:
        process_steps += ['crmask']

    combineImage = combineimage.CombineImage(spectrograph, det, par, file_list)
    pypeitImage = combineImage.run(process_steps, bias, bpm=bpm, pixel_flat=pixel_flat,
                                 illum_flat=illum_flat, sigma_clip=sigma_clip,
                                 sigrej=sigrej, maxiters=maxiters)

    # Instantiate
    slf = ScienceImage(spectrograph, det, par, pypeitImage.image, pypeitImage.ivar,
                                pypeitImage.bpm, rn2img=pypeitImage.rn2img,
                                crmask=pypeitImage.crmask, mask=pypeitImage.mask,
                                files=file_list)
    # Return
    return slf


