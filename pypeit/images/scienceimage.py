""" Object to hold + process a single image"""

import inspect

import os
import numpy as np


from pypeit import msgs

from pypeit.core import procimg
from pypeit.core import coadd2d
from pypeit.par import pypeitpar
from pypeit import utils

from pypeit.images import pypeitimage
from pypeit.images import processrawimage
from pypeit.images import maskimage

from IPython import embed


class ScienceImage(pypeitimage.PypeItImage, maskimage.ImageMask):
    """
    Class to hold and process a science image

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
    frametype = 'science'

    def __init__(self, spectrograph, det, par, bpm):

        # Init me
        maskimage.ImageMask.__init__(self, bpm)
        pypeitimage.PypeItImage.__init__(self, spectrograph, det)

        # Required parameters
        if not isinstance(par, pypeitpar.ProcessImagesPar):
            msgs.error('Provided ParSet for must be type ProcessImagesPar.')
        self.par = par  # This musts be named this way as it is frequently a child

        # Internal images
        self.ivar = None
        self.rn2img = None

        # Other internals
        self.filename = None
        self.files = []

    @classmethod
    def from_images(cls, spectrograph, det, par, bpm,
                    image, ivar, rn2img, crmask=None, mask=None,
                    files=None):
        """
        Instantiate from a set of images

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
            image (np.ndarray):
            ivar (np.ndarray):
                Inverse variance image
            rn2img (np.ndarray):
                Read noise**2 image
            crmask (np.ndarray, optional):
                CR mask
            mask (np.ndarray, optional):
                Full mask
            files (list, optional):
                Image list

        Returns:
            ScienceImage:

        """
        # Init
        slf = cls(spectrograph, det, par, bpm)
        # Other images
        slf.image = image
        slf.ivar = ivar
        slf.rn2img = rn2img
        # Masks
        slf.crmask = crmask
        slf.mask = mask
        # Files
        if files is not None:
            slf.files = files
        # Return
        return slf

    @classmethod
    def from_single_file(cls, spectrograph, det, par, bpm,
                         filename, bias, pixel_flat, illum_flat=None):
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
        slf = cls(spectrograph, det, par, bpm)
        # Build up
        prawImage = processrawimage.ProcessRawImage(filename, slf.spectrograph,
                                                    slf.det, slf.par)
        # Save a bit
        slf.filename = filename
        slf.exptime = prawImage.exptime
        # Process steps
        process_steps = procimg.init_process_steps(bias, slf.par)
        process_steps += ['trim', 'apply_gain', 'orient']
        if (pixel_flat is not None) or (illum_flat is not None):
            process_steps += ['flatten']
        # Do it
        slf.image = prawImage.process(process_steps, pixel_flat=pixel_flat,
                                       bias=bias, illum_flat=illum_flat,
                                       bpm=slf.bpm, debug=True)
        # Build the rest
        slf.build_ivar()
        slf.build_crmask()
        slf.build_rn2img()
        slf.build_mask(saturation=slf.spectrograph.detector[slf.det-1]['saturation'],
                       mincounts=slf.spectrograph.detector[slf.det-1]['mincounts'])
        # Return
        return slf

    @classmethod
    def from_file_list(cls, spectrograph, det, par, bpm,
                       file_list, bias, pixel_flat, illum_flat=None,
                       sigma_clip=False, sigrej=None, maxiters=5):
        """
        Instantiate from file list

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
        # Single file?
        if len(file_list) == 1:
            return cls.from_single_file(spectrograph, det, par, bpm,
                         file_list[0], bias, pixel_flat, illum_flat=illum_flat)

        # Continue with an actual list
        # Get it ready
        nimages = len(file_list)
        shape = (nimages, bpm.shape[0], bpm.shape[1])
        sciimg_stack = np.zeros(shape)
        sciivar_stack= np.zeros(shape)
        rn2img_stack = np.zeros(shape)
        crmask_stack = np.zeros(shape, dtype=bool)

        # Mask
        bitmask = maskimage.ImageBitMask()
        mask_stack = np.zeros(shape, bitmask.minimum_dtype(asuint=True))

        # Loop on the files
        for kk, ifile in enumerate(file_list):
            # Instantiate
            sciImage = ScienceImage.from_single_file(spectrograph, det, par, bpm,
                ifile, bias, pixel_flat, illum_flat=illum_flat)
            # Process
            sciimg_stack[kk,:,:] = sciImage.image
            # Construct raw variance image and turn into inverse variance
            sciivar_stack[kk, :, :] = sciImage.ivar
            # Mask cosmic rays
            crmask_stack[kk, :, :] = sciImage.crmask
            # Build read noise squared image
            rn2img_stack[kk, :, :] = sciImage.build_rn2img()
            # Final mask for this image
            mask_stack[kk, :, :] = sciImage.mask

        # Coadd them
        weights = np.ones(nimages)/float(nimages)
        img_list = [sciimg_stack]
        var_stack = utils.inverse(sciivar_stack, positive=True)
        var_list = [var_stack, rn2img_stack]
        img_list_out, var_list_out, outmask, nused = coadd2d.weighted_combine(
            weights, img_list, var_list, (mask_stack == 0),
            sigma_clip=sigma_clip, sigma_clip_stack=sciimg_stack, sigrej=sigrej, maxiters=maxiters)

        # Build the last one
        slf = ScienceImage.from_images(spectrograph, det, par, bpm, img_list_out[0],
                                       utils.inverse(var_list_out[0], positive=True),
                                       var_list_out[1], np.invert(outmask), files=file_list)
        slf.build_mask(saturation=slf.spectrograph.detector[slf.det-1]['saturation'],
                       mincounts=slf.spectrograph.detector[slf.det-1]['mincounts'])
        # Return
        return slf


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
                                                      utils.inverse(self.ivar, positive=True),
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

    def build_ivar(self):
        """
        Generate the Inverse Variance frame

        Uses procimg.variance_frame

        Returns:
            np.ndarray: Copy of self.ivar

        """
        msgs.info("Generating raw variance frame (from detected counts [flat fielded])")
        # Convenience
        detector = self.spectrograph.detector[self.det-1]
        # For RN variations by Amp
        datasec_img = self.spectrograph.get_datasec_img(filename=self.filename,
                                                     det=self.det)
        # Generate
        rawvarframe = procimg.variance_frame(datasec_img, self.image,
                                                  detector['gain'], detector['ronoise'],
                                                  numamplifiers=detector['numamplifiers'],
                                                  darkcurr=detector['darkcurr'],
                                                  exptime=self.exptime)
        # Ivar
        self.ivar = utils.inverse(rawvarframe, positive=True)
        # Return
        return self.ivar.copy()

    # TODO sort out dark current here. Need to pass exposure time for that.
    def build_rn2img(self):
        """
        Generate the model read noise squared image

        Currently only used by ScienceImage.

        Wrapper to procimg.rn_frame

        Returns:
            np.ndarray: Copy of the read noise squared image

        """
        msgs.info("Generating read noise image from detector properties and amplifier layout)")
        # Convenience
        detector = self.spectrograph.detector[self.det-1]
        datasec_img = self.spectrograph.get_datasec_img(filename=self.filename,
                                                        det=self.det)
        # Build it
        self.rn2img = procimg.rn_frame(datasec_img,
                                       detector['gain'],
                                       detector['ronoise'],
                                       numamplifiers=detector['numamplifiers'])
        # Return
        return self.rn2img.copy()


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
                                                      utils.inverse(self.ivar, positive=True),
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
        # Images
        newimg = self.image - other.image

        # Mask time
        outmask_comb = (self.mask == 0) & (other.mask == 0)

        # Variance
        if self.ivar is not None:
            new_ivar = utils.inverse(utils.inverse(self.ivar, positive=True) + utils.inverse(other.ivar, positive=True))
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
        new_sciImg = ScienceImage.from_images(self.spectrograph, self.det, self.par, self.bpm,
                                              newimg, new_ivar, new_rn2, files=new_files)
        #TODO: Check this out JFH
        #embed(header='279 in sciImg')
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




