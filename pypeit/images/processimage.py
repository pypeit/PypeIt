""" Object to hold + process a single image"""

import inspect

import numpy as np

from pypeit import msgs

from pypeit.core import procimg
from pypeit.core import flat

from pypeit.images import pypeitimage

from IPython import embed

# REMOVE THIS
from importlib import reload
reload(procimg)


class ProcessImage(pypeitimage.PypeItImage):

    def __init__(self, filename, spectrograph, det, par, frametype=None):

        # Init me
        pypeitimage.PypeItImage.__init__(self, filename, spectrograph, det)
        # Required parameters
        self.par = par  # ProcessImagesPar

        # Optional parameters
        self.frametype = frametype

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
    def datasec_img(self):
        dimg = self.spectrograph.get_datasec_img(self.filename, self.det)
        return dimg

    @property
    def oscansec_img(self):
        oimg = self.spectrograph.get_oscansec_img(self.filename, self.det)
        return oimg

    def reset_steps(self):
        for key in self.steps.keys():
            self.steps[key] = False

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
        self.image *= procimg.gain_frame(self.datasec_img, gain, trim=self.steps['trim'])
        self.steps[step] = True
        # Return
        return self.image.copy()

    def flatten(self, pixel_flat, illum_flat=None, bpm=None, force=False):
        step = inspect.stack()[0][3]
        # Check if already trimmed
        if self.steps[step] and (not force):
            msgs.warn("Image was already flat fielded.  Returning the current image")
            return self.image
        # BPM
        if bpm is None:
            bpm = self.bpm
        # Do it
        self.image = flat.flatfield(self.image, pixel_flat, self.bpm,
                                    illum_flat=illum_flat)
        self.steps[step] = True
        # Return
        return self.image.copy()

    def subtract_bias(self, bias_image, force=False):
        step = inspect.stack()[0][3]
        # Check if already trimmed
        if self.steps[step] and (not force):
            msgs.warn("Image was already bias subtracted.  Returning the current image")
            return self.image
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

        temp = procimg.new_subtract_overscan(self.image, self.datasec_img, self.oscansec_img,
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
                return self.image
        # Check if already trimmed
        if self.steps[step] and (not force):
            msgs.warn("Image was already trimmed.  Returning current image")
            return self.image
        # Do it
        trim_image = procimg.trim_frame(self.image, self.datasec_img < 1)
        # Overwrite
        self.image = trim_image
        self.steps[step] = True
        # Return
        return self.image.copy()

    def orient(self, force=False):
        step = inspect.stack()[0][3]
        # Orient the image to have blue/red run bottom to top
        # Check if already oriented
        if self.steps[step] and (not force):
            msgs.warn("Image was already oriented.  Returning current image")
            return self.image
        # Transpose?
        if self.spectrograph.detector[self.det-1]['specaxis'] == 1:
            self.image = self.image.T
        # Flip spectgral axis?
        if self.spectrograph.detector[self.det-1]['specflip'] is True:
            self.image = np.flip(self.image, axis=0)
        # Flip spatial axis?
        if self.spectrograph.detector[self.det-1]['spatflip'] is True:
            self.image = np.flip(self.image, axis=1)

        self.steps[step] = True
        # Return
        return self.image.copy()


