""" Object to hold + process a single image"""

import inspect

import numpy as np
from collections import OrderedDict

from pypeit import msgs
from pypeit import utils

from pypeit.core import procimg
from pypeit.core import flat

from pypeit.images import pypeitimage
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


class ProcessImage(pypeitimage.PypeItImage):

    bitmask = ProcessImagesBitMask()  # The bit mask interpreter

    def __init__(self, filename, spectrograph, det, par, frametype=None):

        # Init me
        pypeitimage.PypeItImage.__init__(self, spectrograph, det)
        # Required parameters
        self.par = par  # ProcessImagesPar
        self.filename = filename

        # Optional parameters
        self.frametype = frametype

        # Attributes
        self._reset_internals()

        # All possible processing steps
        #  Note these have to match the method names below
        self.steps = dict(subtract_bias=False,
                          subtract_overscan=False,
                          subtract_dark=False,
                          trim=False,
                          apply_gain=False,
                          orient=False,
                          flatten=False,
                          )
    @property
    def amps(self):
        """
        Return a list of the amplifier indices, 1-indexed

        Returns:
            list
        """
        amps = np.unique(self.rawdatasec_img[self.rawdatasec_img > 0]).tolist()
        # Return
        return amps

    @property
    def rawdatasec_img(self):
        rdimg = self.spectrograph.get_rawdatasec_img(self.filename, self.det)
        return rdimg

    @property
    def datasec_img(self):
        dimg = self.rawdatasec_img
        # Fuss
        if self.steps['trim']:
            dimg = procimg.trim_frame(dimg, dimg < 1)
        if self.steps['orient']:
            dimg = self.spectrograph.orient_image(dimg, self.det)
        # Return
        return dimg

    @property
    def oscansec_img(self):
        oimg = self.spectrograph.get_oscansec_img(self.filename, self.det)
        return oimg

    def _reset_steps(self):
        for key in self.steps.keys():
            self.steps[key] = False

    def _reset_internals(self):
        """
        Init or free up memory
        """
        self.rawvarframe = None
        self.crmask = None
        self.mask = None
        self.rn2img = None

    def apply_gain(self, force=False):
        step = inspect.stack()[0][3]
        # Check if already trimmed
        if self.steps[step] and (not force):
            msgs.warn("Gain was already applied. Returning")
            return self.image.copy()

        gain = self.spectrograph.detector[self.det-1]['gain']
        # Enforce gain being a list
        if not isinstance(gain, list):
            gain = [gain]
        # Apply
        self.image *= procimg.gain_frame(self.rawdatasec_img, gain, trim=self.steps['trim'])
        self.steps[step] = True
        # Return
        return self.image.copy()

    def build_crmask(self):
        """
        Generate the CR mask frame

        Wrapper to procimg.lacosmic

        Requires self.rawvarframe to exist

        Returns:

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

    def build_mask(self, bpm=None, saturation=1e10, mincounts=-1e10, slitmask=None):
        """
        Return the bit value mask used during extraction.

        The mask keys are defined by :class:`ScienceImageBitMask`.  Any
        pixel with mask == 0 is valid, otherwise the pixel has been
        masked.  To determine why a given pixel has been masked::

            bitmask = ScienceImageBitMask()
            reasons = bm.flagged_bits(mask[i,j])

        To get all the pixel masked for a specific set of reasons::

            indx = bm.flagged(mask, flag=['CR', 'SATURATION'])

        Returns:
            numpy.ndarray: The bit value mask for the science image.
        """
        # Init
        sciivar = utils.calc_ivar(self.rawvarframe)
        if bpm is None:
            bpm = self.bpm

        # Instatiate the mask
        mask = np.zeros_like(self.image, dtype=self.bitmask.minimum_dtype(asuint=True))

        # Bad pixel mask
        indx = bpm.astype(bool)
        mask[indx] = self.bitmask.turn_on(mask[indx], 'BPM')

        # Cosmic rays
        indx = self.crmask.astype(bool)
        mask[indx] = self.bitmask.turn_on(mask[indx], 'CR')

        # Saturated pixels
        indx = self.image >= saturation
        mask[indx] = self.bitmask.turn_on(mask[indx], 'SATURATION')

        # Minimum counts
        indx = self.image <= mincounts
        mask[indx] = self.bitmask.turn_on(mask[indx], 'MINCOUNTS')

        # Undefined counts
        indx = np.invert(np.isfinite(self.image))
        mask[indx] = self.bitmask.turn_on(mask[indx], 'IS_NAN')

        # Bad inverse variance values
        indx = np.invert(sciivar > 0.0)
        mask[indx] = self.bitmask.turn_on(mask[indx], 'IVAR0')

        # Undefined inverse variances
        indx = np.invert(np.isfinite(sciivar))
        mask[indx] = self.bitmask.turn_on(mask[indx], 'IVAR_NAN')

        if slitmask is not None:
            indx = slitmask == -1
            mask[indx] = self.bitmask.turn_on(mask[indx], 'OFFSLITS')

        # Save
        self.mask = mask

        return mask.copy()

    def build_rawvarframe(self):
        """
        Generate the Raw Variance frame
        Currently only used by ScienceImage.

        Wrapper to procimg.variance_frame

        Args:
            trim:

        Returns:
            np.ndarray: self.rawvarframe

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

        Returns
        -------
        self.rn2img : ndarray

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

    def flatten(self, pixel_flat, illum_flat=None, bpm=None, force=False):
        step = inspect.stack()[0][3]
        # Check if already trimmed
        if self.steps[step] and (not force):
            msgs.warn("Image was already flat fielded.  Returning the current image")
            return self.image.copy()
        # BPM
        if bpm is None:
            bpm = self.bpm
        # Do it
        self.image = flat.flatfield(self.image, pixel_flat, bpm, illum_flat=illum_flat)
        self.steps[step] = True
        # Return
        return self.image.copy()

    def orient(self, force=False):
        step = inspect.stack()[0][3]
        # Orient the image to have blue/red run bottom to top
        # Check if already oriented
        if self.steps[step] and (not force):
            msgs.warn("Image was already oriented.  Returning current image")
            return self.image.copy()
        # Orient me
        self.image = self.spectrograph.orient_image(self.image, self.det)
        self.steps[step] = True
        # Return
        return self.image.copy()

    def subtract_bias(self, bias_image, force=False):
        step = inspect.stack()[0][3]
        # Check if already trimmed
        if self.steps[step] and (not force):
            msgs.warn("Image was already bias subtracted.  Returning the current image")
            return self.image.copy()
        # Do it
        self.image -= bias_image
        self.steps[step] = True
        # Return
        return self.image.copy()

    def subtract_overscan(self, force=False):
        step = inspect.stack()[0][3]
        # Check if already trimmed
        if self.steps[step] and (not force):
            msgs.warn("Image was already trimmed")

        temp = procimg.new_subtract_overscan(self.image, self.rawdatasec_img, self.oscansec_img,
                                         method=self.par['overscan'],
                                         params=self.par['overscan_par'])
        # Fill
        self.steps[step] = True
        self.image = temp
        # Return
        return self.image.copy()

    def trim(self, force=False):
        step = inspect.stack()[0][3]
        # Check input image matches the original
        if self.orig_shape is not None:
            if self.image.shape != self.orig_shape:
                msgs.warn("Image shape does not match original.  Returning current image")
                return self.image.copy()
        # Check if already trimmed
        if self.steps[step] and (not force):
            msgs.warn("Image was already trimmed.  Returning current image")
            return self.image
        # Do it
        trim_image = procimg.trim_frame(self.image, self.rawdatasec_img < 1)
        # Overwrite
        self.image = trim_image
        self.steps[step] = True
        # Return
        return self.image.copy()

    def __repr__(self):
        return ('<{:s}: file={}, steps={}>'.format(
            self.__class__.__name__, self.filename, self.steps))

