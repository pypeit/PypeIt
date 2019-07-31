""" Object process a single raw image"""

import inspect

import numpy as np

from pypeit import msgs
from pypeit.core import procimg
from pypeit.core import flat
from pypeit.par import pypeitpar
from pypeit.images import pypeitimage
from pypeit import utils

from IPython import embed


class ProcessRawImage(object):
    """
    Class to process a raw image

    Args:
        filename (:obj:`str` or None):
            Filename
        par (:class:`pypeit.par.pypeitpar.ProcessImagesPar`):
            Parameters that dictate the processing of the images.  See
            :class:`pypeit.par.pypeitpar.ProcessImagesPar` for the
            defaults.
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            Spectrograph used to take the data.
        det (:obj:`int`, optional):
            The 1-indexed detector number to process.

    Attributes:
        steps (dict):
            Dict describing the steps performed on the image
        datasec_img (np.ndarray):
            Holds the datasec_img which specifies the amp for each pixel in the
            current self.image image
    """
    def __init__(self, rawImage, par, bpm=None):

        # Required parameters
        if not isinstance(par, pypeitpar.ProcessImagesPar):
            msgs.error("Bad par input")
        self.par = par  # ProcessImagesPar
        self._bpm = bpm

        # Grab a few things from rawImage (for convenience and for processing)
        #   Could just keep rawImage in the object, if preferred
        self.spectrograph = rawImage.spectrograph
        self.det = rawImage.det
        self.filename = rawImage.filename
        self.datasec_img = rawImage.rawdatasec_img.copy()
        self.oscansec_img = rawImage.oscansec_img
        self.image = rawImage.raw_image.copy()
        self.orig_shape = rawImage.raw_image.shape
        self.exptime = rawImage.exptime

        # Binning of the processed image
        self.binning = self.spectrograph.get_meta_value(self.filename, 'binning')

        # Attributes
        self.ivar = None
        self.rn2img = None

        # All possible processing steps
        #  Note these have to match the method names below
        self.steps = dict(subtract_bias=False,
                          subtract_overscan=False,
                          subtract_dark=False,
                          trim=False,
                          orient=False,
                          apply_gain=False,
                          flatten=False,
                          )

    @property
    def bpm(self):
        """
        Generate and return the bad pixel mask for this image
        Warning:  BPM masks are for processed (e.g. trimmed, rotated) images only!

        Returns:
            np.ndarray:  Bad pixel mask with a bad pixel = 1

        """
        if self._bpm is None:
            self._bpm = self.spectrograph.bpm(shape=self.image.shape,
                                    filename=self.filename,
                                    det=self.det)
        return self._bpm

    def _reset_steps(self):
        """
        Reset all the processing steps to False

        Should consider setting the Image to None too..
        """
        for key in self.steps.keys():
            self.steps[key] = False

    def apply_gain(self, force=False):
        """
        Apply the Gain values to self.image

        Args:
            force (bool, optional):

        Returns:
            np.ndarray:  copy of self.image

        """
        step = inspect.stack()[0][3]
        # Check if already trimmed
        if self.steps[step] and (not force):
            msgs.warn("Gain was already applied. Returning")
            return self.image.copy()

        gain = np.atleast_1d(self.spectrograph.detector[self.det - 1]['gain']).tolist()
        # Apply
        self.image *= procimg.gain_frame(self.datasec_img, gain) #, trim=self.steps['trim'])
        self.steps[step] = True
        # Return
        return self.image.copy()

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
        #embed(header='333 of sciencimage')
        #datasec_img = self.spectrograph.get_datasec_img(filename=self.filename,
        #                                             det=self.det)
        # Generate
        rawvarframe = procimg.variance_frame(self.datasec_img, self.image,
                                             detector['gain'], detector['ronoise'],
                                             numamplifiers=detector['numamplifiers'],
                                             darkcurr=detector['darkcurr'],
                                             exptime=self.exptime)
        # Ivar
        self.ivar = utils.inverse(rawvarframe)
        # Return
        return self.ivar.copy()

    # TODO sort out dark current here. Need to use exposure time for that.
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

    def process(self, process_steps, pixel_flat=None, illum_flat=None,
                       bias=None, debug=False):
        """
        Process the image

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

        Returns:
            pypeitimage.PypeItImage:

        """
        # Get started
        # Standard order
        #   -- May need to allow for other order some day..
        if 'subtract_overscan' in process_steps:
            self.subtract_overscan()
        if 'trim' in process_steps:
            self.trim()
        if 'orient' in process_steps:
            self.orient()
        if 'subtract_bias' in process_steps: # Bias frame, if it exists, is trimmed and oriented
            self.subtract_bias(bias)
        if 'apply_gain' in process_steps:
            self.apply_gain()
        # Flat field
        if 'flatten' in process_steps:
            self.flatten(pixel_flat, illum_flat=illum_flat, bpm=self.bpm)

        # Fresh BPM
        bpm = self.spectrograph.bpm(self.filename, self.det, shape=self.image.shape)

        # Extras
        if 'extras' in process_steps:
            self.build_ivar()
            self.build_rn2img()

        # Return a PypeItImage
        pypeitImage = pypeitimage.PypeItImage(self.image, state=self.steps, binning=self.binning,
                                                       ivar=self.ivar, rn2img=self.rn2img, bpm=bpm)
        return pypeitImage

    def flatten(self, pixel_flat, illum_flat=None, bpm=None, force=False):
        """
        Flat field the proc image

        Wrapper to flat.flatfield

        Args:
            pixel_flat (np.ndarray):
                Pixel flat image
            illum_flat (np.ndarray, optional):
                Illumination flat image
            bpm (np.ndarray, optional):
                Bad pixel mask image;  if provided, over-rides internal one
            force (bool, optional):
                Force the processing even if the image was already processed

        """
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

    '''
    def load_rawframe(self):
        """
        Load a raw image from disk using the Spectrograph method load_raw_frame()

        Also loads up the binning, exposure time, and header of the Primary image
        And the HDUList in self.hdu

        """
        # Load
        #self.image, self.hdu, \
        #    = self.spectrograph.load_raw_frame(self.filename, det=self.det)
        #self.exptime, self.datasec_img, self.oscansec_img, self.binning_raw \
        #    = self.spectrograph.load_raw_extras(self.hdu, self.det)
        # Load itup
        self.image, self.hdu, self.exptime, self.datasec_img, self.oscansec_img = self.spectrograph.get_rawimage(
            self.filename, self.det)

        self.head0 = self.hdu[0].header
        # Shape
        self.orig_shape = self.image.shape
    '''

    def orient(self, force=False):
        """
        Orient the image in the PypeIt format with spectra running blue (down)
        to red (up).

        Args:
            force (bool, optional):
                Force the processing even if the image was already processed

        """
        step = inspect.stack()[0][3]
        # Orient the image to have blue/red run bottom to top
        # Check if already oriented
        if self.steps[step] and not force:
            msgs.warn("Image was already oriented.  Returning current image")
            return self.image.copy()
        # Orient me
        self.image = self.spectrograph.orient_image(self.image, self.det)
        self.datasec_img = self.spectrograph.orient_image(self.datasec_img, self.det)
        #
        self.steps[step] = True

    def subtract_bias(self, bias_image, force=False):
        """
        Perform bias subtraction

        Args:
            bias_image (PypeItImage):
                Bias image
            force (bool, optional):
                Force the processing even if the image was already processed
        """
        step = inspect.stack()[0][3]
        # Check if already bias subtracted
        if self.steps[step] and (not force):
            msgs.warn("Image was already bias subtracted.  Returning the current image")
            return self.image.copy()
        # Do it
        self.image -= bias_image.image
        self.steps[step] = True

    def subtract_overscan(self, force=False):
        """
        Analyze and subtract the overscan from the image

        Args:
            force (bool, optional):
                Force the processing even if the image was already processed

        """
        step = inspect.stack()[0][3]
        # Check if already overscan subtracted
        if self.steps[step] and (not force):
            msgs.warn("Image was already overscan subtracted!")

        temp = procimg.subtract_overscan(self.image, self.datasec_img, self.oscansec_img,
                                         method=self.par['overscan'],
                                         params=self.par['overscan_par'])
        # Fill
        self.steps[step] = True
        self.image = temp

    def trim(self, force=False):
        """
        Trim the image to include only the science data

        Args:
            force (bool, optional):
                Force the processing even if the image was already processed

        """
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
        self.image = procimg.trim_frame(self.image, self.datasec_img < 1)
        self.datasec_img = procimg.trim_frame(self.datasec_img, self.datasec_img < 1)
        #
        self.steps[step] = True

    def __repr__(self):
        return ('<{:s}: file={}, steps={}>'.format(
            self.__class__.__name__, self.filename, self.steps))

