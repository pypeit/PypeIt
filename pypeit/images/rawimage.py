""" Object to load and process a single raw image

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import inspect
from copy import deepcopy
import numpy as np
from astropy import stats

from pypeit import msgs
from pypeit.core import procimg
from pypeit.core import flat
from pypeit.core import flexure
from pypeit.images import pypeitimage
from pypeit import utils
from pypeit.display import display

from IPython import embed

class RawImage(object):
    """
    Class to load and process a raw image

    Args:
        ifile (:obj:`str`):
        bpm (`numpy.ndarray`_, optional):
            Bad pixel mask

    Attributes:
        rawimage (`numpy.ndarray`_):
        steps (dict):
            Dict describing the steps performed on the image
        datasec_img (`numpy.ndarray`_):
            Holds the datasec_img which specifies the amp for each pixel in the
            current self.image image.  This is modified as the image is, i.e.
            orientation and trimming.
        spat_flexure_shift (float):
            Holds the spatial flexure shift, if calculated
        image (`numpy.ndarray`_):
    """
    def __init__(self, ifile, spectrograph, det):

        # Required parameters
        self.filename = ifile
        self.spectrograph = spectrograph  # One for convenience
        self.det = det

        # Load
        # Load the raw image and the other items of interest
        self.detector, self.rawimage, self.hdu, self.exptime, self.rawdatasec_img, \
            self.oscansec_img = self.spectrograph.get_rawimage(self.filename, self.det)

        # Grab items from rawImage (for convenience and for processing)
        #   Could just keep rawImage in the object, if preferred
        self.headarr = deepcopy(self.spectrograph.get_headarr(self.hdu))

        # Key attributes
        self.rawimage = self.rawimage.copy()
        self.image = self.rawimage.copy()
        self.datasec_img = self.rawdatasec_img.copy()
        self.ronoise = self.detector['ronoise']

        # Attributes
        self.par = None
        self.ivar = None
        self.rn2img = None
        self.spat_flexure_shift = None

        # All possible processing steps
        #  Note these have to match the method names below
        self.steps = dict(subtract_bias=False,
                          subtract_overscan=False,
                          subtract_dark=False,
                          subtract_pattern=False,
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
            `numpy.ndarray`_:  Bad pixel mask with a bad pixel = 1

        """
        if self._bpm is None:
            self._bpm = self.spectrograph.bpm(shape=self.image.shape,
                                    filename=self.filename,
                                    det=self.det)
        return self._bpm

    def apply_gain(self, force=False):
        """
        Apply the Gain values to self.image

        Args:
            force (bool, optional):

        Returns:
            `numpy.ndarray`:  copy of self.image

        """
        step = inspect.stack()[0][3]
        # Check if already trimmed
        if self.steps[step] and (not force):
            msgs.warn("Gain was already applied. Returning")
            return self.image.copy()

        gain = np.atleast_1d(self.detector['gain']).tolist()
        # Apply
        self.image *= procimg.gain_frame(self.datasec_img, gain)
        self.steps[step] = True
        # Return
        return self.image.copy()

    def build_ivar(self):
        """
        Generate the Inverse Variance frame

        Uses procimg.variance_frame

        Returns:
            `numpy.ndarray`_: Copy of self.ivar

        """
        # Generate
        rawvarframe = procimg.variance_frame(self.datasec_img, self.image,
                                             self.detector['gain'], self.detector['ronoise'],
                                             darkcurr=self.detector['darkcurr'],
                                             exptime=self.exptime,
                                             rnoise=self.rn2img)
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
            `numpy.ndarray`_: Copy of the read noise squared image

        """
        # Check if we need to determine the read noise directly from the overscan region from any amplifier
        numamps = len(self.ronoise)
        for amp in range(numamps):
            if self.ronoise[amp] <= 0.0:
                biaspix = self.rawimage[self.oscansec_img==amp+1] * self.detector['gain'][amp]
                _, _, stddev = stats.sigma_clipped_stats(biaspix, sigma=5)
                self.ronoise[amp] = stddev
                msgs.info("Read noise of amplifier {0:d} = {1:.3f} e-".format(amp+1, self.ronoise[amp]))

        # Build it
        self.rn2img = procimg.rn_frame(self.datasec_img,
                                       self.detector['gain'],
                                       self.ronoise)
        # Return
        return self.rn2img.copy()

    def process(self, par, bpm=bpm, flatimages=None, bias=None,
                slits=None, debug=False, dark=None):
        """
        Process the image

        Note:  The processing step order is currently 'frozen' as is.
          We may choose to allow optional ordering

        Here are the allowed steps, in the order they will be applied:
            subtract_overscan -- Analyze the overscan region and subtract from the image
            trim -- Trim the image down to the data (i.e. remove the overscan)
            orient -- Orient the image in the PypeIt orientation (spec, spat) with blue to red going down to up
            subtract_bias -- Subtract a bias image
            apply_gain -- Convert to counts, amp by amp
            flatten -- Divide by the pixel flat and (if provided) the illumination flat
            extras -- Generate the RN2 and IVAR images
            crmask -- Generate a CR mask

        Args:
            par (:class:`pypeit.par.pypeitpar.ProcessImagesPar`):
                Parameters that dictate the processing of the images.  See
                :class:`pypeit.par.pypeitpar.ProcessImagesPar` for the
                defaults.
            bpm (`numpy.ndarray`_, optional):
            flatimages (:class:`pypeit.flatfield.FlatImages`):
            bias (`numpy.ndarray`_, optional):
                Bias image
            slits (:class:`pypeit.slittrace.SlitTraceSet`, optional):
                Used to calculate spatial flexure between the image and the slits

        Returns:
            :class:`pypeit.images.pypeitimage.PypeItImage`:

        """
        self.par = par
        self._bpm = bpm

        # Get started
        # Standard order
        #   -- May need to allow for other order some day..
        if par['use_pattern']:  # Note, this step *must* be done before use_overscan
            self.subtract_pattern()
        if par['use_overscan']:
            self.subtract_overscan()
        if par['trim']:
            self.trim()
        if par['orient']:
            self.orient()
        if par['use_biasimage']:  # Bias frame, if it exists, is *not* trimmed nor oriented
            self.subtract_bias(bias)
        if par['use_darkimage']:  # Dark frame, if it exists, is TODO:: check: trimmed, oriented (and oscan/bias subtracted?)
            self.subtract_dark(dark)
        if par['apply_gain']:
            self.apply_gain()

        # This needs to come after trim, orient
        # Calculate flexure -- May not be used, but always calculated when slits are provided
        if slits is not None and self.par['spat_flexure_correct']:
            self.spat_flexure_shift = flexure.spat_flexure_shift(self.image, slits)

        # Generate the illumination flat, as needed
        illum_flat = None
        if self.par['use_illumflat']:
            if flatimages is None:
                msgs.error("Cannot illumflatten, no such image generated. Add one or more illumflat images to your PypeIt file!!")
            if slits is None:
                msgs.error("Need to provide slits to create illumination flat")
            illum_flat = flatimages.fit2illumflat(slits, flexure_shift=self.spat_flexure_shift)
            if debug:
                left, right = slits.select_edges(flexure=self.spat_flexure_shift)
                viewer, ch = display.show_image(illum_flat, chname='illum_flat')
                display.show_slits(viewer, ch, left, right)  # , slits.id)
                #
                orig_image = self.image.copy()
                viewer, ch = display.show_image(orig_image, chname='orig_image')
                display.show_slits(viewer, ch, left, right)  # , slits.id)

        # Apply the relative spectral illumination
        spec_illum = 1.0
        if self.par['use_specillum']:
            if flatimages is None or flatimages.get_spec_illum() is None:
                msgs.error("Spectral illumination correction desired but not generated/provided.")
            else:
                spec_illum = flatimages.get_spec_illum().copy()

        # Flat field -- We cannot do illumination flat without a pixel flat (yet)
        if self.par['use_pixelflat'] or self.par['use_illumflat']:
            if flatimages is None or flatimages.get_pixelflat() is None:
                msgs.error("Flat fielding desired but not generated/provided.")
            else:
                self.flatten(flatimages.get_pixelflat()/spec_illum, illum_flat=illum_flat, bpm=self.bpm)

        # Fresh BPM
        bpm = self.spectrograph.bpm(self.filename, self.det, shape=self.image.shape)

        # Extras
        self.build_rn2img()
        self.build_ivar()

        # Generate a PypeItImage
        pypeitImage = pypeitimage.PypeItImage(self.image, ivar=self.ivar, rn2img=self.rn2img,
                                              bpm=bpm, detector=self.detector,
                                              spat_flexure=self.spat_flexure_shift,
                                              PYP_SPEC=self.spectrograph.name)
        pypeitImage.rawheadlist = self.headarr
        pypeitImage.process_steps = [key for key in self.steps.keys() if self.steps[key]]

        # Mask(s)
        if par['mask_cr']:
            pypeitImage.build_crmask(self.par)
        #
        nonlinear_counts = self.spectrograph.nonlinear_counts(self.detector,
                                                              apply_gain=self.par['apply_gain'])
        # Build
        pypeitImage.build_mask(saturation=nonlinear_counts)

        # Return
        return pypeitImage

    def flatten(self, pixel_flat, illum_flat=None, bpm=None, force=False):
        """
        Flat field the proc image

        Wrapper to flat.flatfield

        Args:
            pixel_flat (`numpy.ndarray`_):
                Pixel flat image
            illum_flat (`numpy.ndarray`_, optional):
                Illumination flat image
            bpm (`numpy.ndarray`_, optional):
                Bad pixel mask image;  if provided, over-rides internal one
            force (bool, optional):
                Force the processing even if the image was already processed

        """
        step = inspect.stack()[0][3]
        # Check if already trimmed
        if self.steps[step] and (not force):
            msgs.warn("Image was already flat fielded.  Returning the current image")
            return self.image.copy()
        # Check
        if pixel_flat is None:
            msgs.error("You requested flattening but provided no pixel flat!")
        # BPM
        if bpm is None:
            bpm = self.bpm
        # Do it
        self.image = flat.flatfield(self.image, pixel_flat, bpm, illum_flat=illum_flat)
        self.steps[step] = True

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
        self.image = self.spectrograph.orient_image(self.detector, self.image)#, self.det)
        self.datasec_img = self.spectrograph.orient_image(self.detector, self.datasec_img)#, self.det)
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
        if self.steps[step] and not force:
            msgs.warn("Image was already bias subtracted.  Returning the current image")
            return self.image.copy()
        # Do it
        self.image -= bias_image.image
        self.steps[step] = True

    def subtract_dark(self, dark_image, force=False):
        """
        Perform dark image subtraction

        Args:
            dark_image (PypeItImage):
                Dark image
        """
        step = inspect.stack()[0][3]
        # Check if already bias subtracted
        if self.steps[step] and not force:
            msgs.warn("Image was already dark subtracted.  Returning the current image")
            return self.image.copy()
        # Do it
        self.image -= dark_image.image
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
                                         method=self.par['overscan_method'],
                                         params=self.par['overscan_par'])
        # Fill
        self.steps[step] = True
        self.image = temp

    def subtract_pattern(self):
        """
        Analyze and subtract the pattern noise from the image
        """
        step = inspect.stack()[0][3]
        # Check if already overscan subtracted
        if self.steps[step]:
            msgs.warn("Image was already pattern subtracted!")

        # Grab the frequency, if it exists in the header
        # For some instruments, PYPEITFRQ is added to the header in get_rawimage() in the spectrograph file
        # See keck_kcwi.py for an example
        frequency = []
        try:
            # Grab a list of all the amplifiers
            amps = np.sort(np.unique(self.oscansec_img[np.where(self.oscansec_img > 0)]))
            for amp in amps:
                frequency.append(self.hdu[0].header['PYPFRQ{0:02d}'.format(amp)])
            # Final check to make sure the list isn't empty (which it shouldn't be, anyway)
            if len(frequency) == 0:
                frequency = None
        except KeyError:
            frequency = None
        # Generate a new image with the pattern removed
        temp = procimg.subtract_pattern(self.image, self.datasec_img, self.oscansec_img, frequency=frequency)

        # Fill
        self.steps[step] = True
        self.image = temp.copy()
        self.rawimage = temp.copy()

    def trim(self, force=False):
        """
        Trim the image to include only the science data

        Args:
            force (bool, optional):
                Force the processing even if the image was already processed

        """
        step = inspect.stack()[0][3]
        # Check input image matches the original
        if self.rawimage.shape is not None and self.image.shape != self.rawimage.shape:
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
        return ('<{:s}: file={}, steps={}>'.format(self.__class__.__name__, self.filename,
                                                   self.steps))

