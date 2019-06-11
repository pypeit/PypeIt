""" Object to hold + process a single image"""

import inspect

import os
import numpy as np
from collections import OrderedDict


from pypeit import msgs

from pypeit import utils
from pypeit.core import procimg
from pypeit.core import combine
from pypeit.par import pypeitpar

from pypeit.images import pypeitimage
from pypeit.images import processrawimage
from pypeit.bitmask import BitMask


from IPython import embed

# REMOVE THIS
from importlib import reload
reload(procimg)

class ProcessImagesBitMask(BitMask):
    """
    Define a bitmask used to set the reasons why each pixel in a science
    image was masked.
    """

    def __init__(self):
        # TODO:
        #   - Can IVAR0 and IVAR_NAN be consolidated into a single bit?
        #   - Is EXTRACT ever set?
        # TODO: This needs to be an OrderedDict for now to ensure that
        # the bits assigned to each key is always the same. As of python
        # 3.7, normal dict types are guaranteed to preserve insertion
        # order as part of its data model. When/if we require python
        # 3.7, we can remove this (and other) OrderedDict usage in favor
        # of just a normal dict.
        mask = OrderedDict([
            ('BPM', 'Component of the instrument-specific bad pixel mask'),
            ('CR', 'Cosmic ray detected'),
            ('SATURATION', 'Saturated pixel'),
            ('MINCOUNTS', 'Pixel below the instrument-specific minimum counts'),
            ('OFFSLITS', 'Pixel does not belong to any slit'),
            ('IS_NAN', 'Pixel value is undefined'),
            ('IVAR0', 'Inverse variance is undefined'),
            ('IVAR_NAN', 'Inverse variance is NaN'),
            ('EXTRACT', 'Pixel masked during local skysub and extraction')
        ])
        super(ProcessImagesBitMask, self).__init__(list(mask.keys()), descr=list(mask.values()))


class ScienceImage(processrawimage.ProcessRawImage):
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
    bitmask = ProcessImagesBitMask()
    frametype = 'science'

    def __init__(self, spectrograph, det, par, filename=None):

        # Init me
        processrawimage.ProcessRawImage.__init__(self, filename, spectrograph, det, par, frametype=self.frametype)

        # Required parameters
        if not isinstance(par, pypeitpar.ProcessImagesPar):
            msgs.error('Provided ParSet for must be type ProcessImagesPar.')
        self.par = par  # This musts be named this way as it is frequently a child

        # Internal images
        self.image = None

        self.rawvarframe = None
        self.ivar = None
        self.crmask = None
        self.mask = None
        self.rn2img = None

    def build_crmask(self):
        """
        Generate the CR mask frame

        Wrapper to procimg.lacosmic

        Requires self.rawvarframe to exist

        Returns:
            np.ndarray: Copy of self.crmask

        """
        if self.rawvarframe is None:
            msgs.error("Need to generate the rawvariance frame first!")
        # Run LA Cosmic to get the cosmic ray mask
        self.crmask = procimg.lacosmic(self.det, self.image,
                                  self.spectrograph.detector[self.det-1]['saturation'],
                                  self.spectrograph.detector[self.det-1]['nonlinear'],
                                  varframe=self.rawvarframe,
                                  maxiter=self.par['lamaxiter'],
                                  grow=self.par['grow'],
                                  remove_compact_obj=self.par['rmcompact'],
                                  sigclip=self.par['sigclip'],
                                  sigfrac=self.par['sigfrac'],
                                  objlim=self.par['objlim'])
        # Return
        return self.crmask.copy()

    def build_mask(self, slitmask=None, saturation=1e10, mincounts=-1e10):
        """
        Return the bit value mask used during extraction.

        Wrapper to procimg.build_mask

        Args:
            saturation (float, optional):
                Saturation limit in ADU
            mincounts (float, optional):
            slitmask (np.ndarray, optional):
                Slit mask image;  Pixels not in a slit are masked

        Returns:
            numpy.ndarray: Copy of the bit value mask for the science image.

        """
        sciivar = utils.calc_ivar(self.rawvarframe)
        self.mask = procimg.build_mask(self.bitmask, self.image, sciivar,
                                       self.bpm, self.crmask,
                                       saturation=saturation, slitmask=slitmask,
                                       mincounts=mincounts)
        return self.mask.copy()

    def build_rawvarframe(self):
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
        # Generate
        self.rawvarframe = procimg.variance_frame(self.datasec_img, self.image,
                                                  detector['gain'], detector['ronoise'],
                                                  numamplifiers=detector['numamplifiers'],
                                                  darkcurr=detector['darkcurr'],
                                                  exptime=self.exptime)
        # Return
        return self.rawvarframe.copy()

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
        # Build it
        self.rn2img = procimg.rn_frame(self.datasec_img,
                                       detector['gain'],
                                       detector['ronoise'],
                                       numamplifiers=detector['numamplifiers'])
        # Return
        return self.rn2img.copy()

    def process_raw(self, filename, bias, pixel_flat, bpm, illum_flat=None):
        self.filename = filename
        self.load()
        process_steps = procimg.init_process_steps(bias, self.par)
        process_steps += ['trim', 'apply_gain', 'orient']
        if (pixel_flat is not None) or (illum_flat is not None):
            process_steps += ['flatten']

        self.process(process_steps, pixel_flat=pixel_flat, bias=bias,
                     bpm=bpm, illum_flat=illum_flat)
        return self.image.copy()


