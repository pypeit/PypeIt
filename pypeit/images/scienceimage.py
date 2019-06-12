""" Object to hold + process a single image"""

import inspect

import os
import numpy as np


from pypeit import msgs

from pypeit.core import procimg
from pypeit.par import pypeitpar
from pypeit import utils

from pypeit.images import pypeitimage
from pypeit.images import processrawimage
from pypeit.images import maskimage


from IPython import embed

# REMOVE THIS
from importlib import reload
reload(procimg)


class ScienceImage(pypeitimage.PypeItImage, maskimage.ImageMask):
    """
    Class to hold and process a science image

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
        steps (list): List of steps used

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

        #
        self.spectrograph = spectrograph
        self.det = det

        # Internal images
        self.image = None
        self.rawvarframe = None
        self.rn2img = None

        # Other internals
        self.filename = None

    @classmethod
    def from_images(cls, spectrograph, det, par, bpm,
                    image, ivar, rn2img, crmask=None, mask=None):
        # Init
        slf = cls(spectrograph, det, par, bpm)
        # Other images
        slf.image = image
        slf.ivar = ivar
        slf.rn2img = rn2img
        # Masks
        slf.crmask = crmask
        slf.mask = mask
        # Return
        return slf


    def build_crmask(self, subtract_img=None):
        return super(ScienceImage, self).build_crmask(self.spectrograph, self.det,
                                                      self.par, self.image,
                                                      self.rawvarframe,
                                                      subtract_img=subtract_img).copy()

    def build_mask(self, saturation=1e10, mincounts=-1e10, slitmask=None):
        return super(ScienceImage, self).build_mask(self.image, self.ivar,
                                                    saturation=saturation,
                                                    mincounts=mincounts,
                                                    slitmask=slitmask)

    def build_ivar(self):
        """
        Generate the Raw Variance frame
        Currently only used by ScienceImage.

        Wrapper to procimg.variance_frame

        Returns:
            np.ndarray: Copy of self.rawvarframe

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
        self.ivar = utils.calc_ivar(rawvarframe)
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

    def process_raw(self, filename, bias, pixel_flat, illum_flat=None):
        # Build up
        prawImage = processrawimage.ProcessRawImage(filename, self.spectrograph,
                                                    self.det, self.par,
                                                    frametype=self.frametype)
        # Save a bit
        self.filename = filename
        self.exptime = prawImage.exptime
        #
        process_steps = procimg.init_process_steps(bias, self.par)
        process_steps += ['trim', 'apply_gain', 'orient']
        if (pixel_flat is not None) or (illum_flat is not None):
            process_steps += ['flatten']

        self.image = prawImage.process(process_steps, pixel_flat=pixel_flat,
                                       bias=bias, illum_flat=illum_flat, debug=True)
        return self.image.copy()

    def update_mask_cr(self, subtract_img=None):
        # Generate the CR mask (and save in self.crmask)
        _ = super(ScienceImage, self).build_crmask(self.spectrograph, self.det,
                                                      self.par, self.image,
                                                      self.rawvarframe,
                                                      subtract_img=subtract_img).copy()
        # Now update the mask
        _ = super(ScienceImage, self).update_mask_cr(self.crmask)


