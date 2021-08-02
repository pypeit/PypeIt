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

# TODO: I don't understand why we have some of these attributes.  E.g., why do
# we need both hdu and headarr?
class RawImage:
    """
    Class to load and process a raw image

    Args:
        ifile (:obj:`str`):
            File with the data.
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            The spectrograph from which the data was collected.
        det (:obj:`int`):
            The detector to load/process.

    Attributes:
        filename (:obj:`str`):
            Original file name with the data.
        spectrograph (:class:`~pypeit.spectrograph.spectrographs.Spectrograph`):
            Spectrograph instance with the instrument-specific properties and
            methods.
        det (:obj:`int`):
            1-indexed detector number
        detector (:class:`~pypeit.images.detector_container.DetectorContainer`):
            Detector characteristics
        rawimage (`numpy.ndarray`_):
            The raw, not trimmed or reoriented, image data for the detector.
        hdu (`astropy.io.fits.HDUList`_):
            The full list of HDUs provided by :attr:`filename`.
        exptime (:obj:`float`):
            Frame exposure time in seconds.
        rawdatasec_img (`numpy.ndarray`_):
            The original, not trimmed or reoriented, image identifying which
            amplifier was used to read each section of the raw image.
        oscansec_img (`numpy.ndarray`_):
            The original, not trimmed or reoriented, image identifying the
            overscan regions in the raw image read by each amplifier.
        headarr (:obj:`list`):
            A list of `astropy.io.fits.Header`_ objects with the headers for all
            extensions in :attr:`hdu`.
        image (`numpy.ndarray`_):
            The processed image.  This starts as identical to :attr:`rawimage`
            and then altered by the processing steps; see :func:`process`.
        ronoise (:obj:`list`):
            The readnoise (in e-/ADU) for each of the detector amplifiers.
        par (:class:`~pypeit.par.pypeitpar.ProcessImagesPar`):
            Parameters that dictate the processing of the images.
        ivar (`numpy.ndarray`_):
            The inverse variance of :attr:`image`, the processed image.
        rn2img (`numpy.ndarray`_):
            The readnoise variance image.
        steps (:obj:`dict`):
            Dictionary containing a set of booleans that track the processing
            steps that have been performed.
        datasec_img (`numpy.ndarray`_):
            Image identifying which amplifier was used to read each section of
            the *processed* image.
        spat_flexure_shift (:obj:`float`):
            The spatial flexure shift in pixels, if calculated
    """
    def __init__(self, ifile, spectrograph, det):

        # Required parameters
        self.filename = ifile
        self.spectrograph = spectrograph
        self.det = det

        # Load the raw image and the other items of interest
        self.detector, self.rawimage, self.hdu, self.exptime, self.rawdatasec_img, \
                self.oscansec_img = self.spectrograph.get_rawimage(self.filename, self.det)

        # Grab items from rawImage (for convenience and for processing)
        #   Could just keep rawImage in the object, if preferred
        self.headarr = deepcopy(self.spectrograph.get_headarr(self.hdu))

        # Key attributes
        self.image = self.rawimage.copy()
        self.datasec_img = self.rawdatasec_img.copy()
        # NOTE: Prevent estimate_readnoise() from altering self.detector using
        # deepcopy
        self.ronoise = deepcopy(self.detector['ronoise'])

        # Attributes
        self.par = None
        self.ivar = None
        self.var = None
        self.rn2img = None
        self.spat_flexure_shift = None
        self._bpm = None

        # All possible processing steps.  NOTE: These have to match the
        # method/function names.  Their order here matches there execution order
        # in self.process(), but that's not necessary.
        self.steps = dict(subtract_pattern=False,
                          subtract_overscan=False,
                          trim=False,
                          orient=False,
                          subtract_bias=False,
                          apply_gain=False,
                          subtract_dark=False,
                          spatial_flexure_shift=False,
                          flatten=False)

    @property
    def bpm(self):
        """
        Generate and return the bad pixel mask for this image.

        .. warning::

            BPMs are for processed (e.g. trimmed, rotated) images only!

        Returns:
            `numpy.ndarray`_:  Bad pixel mask with a bad pixel = 1

        """
        if self._bpm is None:
            # TODO: Was this ever being called?
#            self._bpm = self.spectrograph.bpm(shape=self.image.shape, filename=self.filename,
#                                              det=self.det)
            # TODO: Pass bias if there is one?
            self._bpm = self.spectrograph.bpm(self.filename, self.det, shape=self.image.shape)
        return self._bpm

    @property
    def use_flat(self):
        """
        Return a flag setting if the flat data should be used in the image
        processing.
        """
        if self.par is None:
            return False
        # TODO: At the moment, we can only perform any of the flat-field
        # corrections if we are applying the pixel-flat correction.
#        return self.par['use_pixelflat'] or self.par['use_specillum'] or self.par['use_illumflat']
        return self.par['use_pixelflat']

    def apply_gain(self, force=False):
        """
        Use the gain to convert images from ADUs to electrons/counts.

        Conversion applied to :attr:`image, :attr:`var`, and :attr:`rn2img`.

        Args:
            force (:obj:`bool`, optional):
                Force the gain to be applied to the image, even if the step log
                (:attr:`steps`) indicates that it already has been.
        """
        step = inspect.stack()[0][3]
        if self.steps[step] and not force:
            # Already applied
            msgs.warn('Gain was already applied.')
            return
        gain = procimg.gain_frame(self.datasec_img, np.atleast_1d(self.detector['gain']))
        self.image *= gain
        if self.var is not None:
            # Convert the variance
            self.var *= gain**2
        if self.rn2img is not None:
            # Convert the detector variance
            self.rn2img *= gain**2
        self.steps[step] = True

#    def build_ivar(self):
#        """
#        Generate the inverse variance in the image.
#
#        This is a simple wrapper for :func:`~pypeit.core.procimg.variance_frame`.
#
#        Returns:
#            `numpy.ndarray`_: The inverse variance in the image.
#        """
#        # Generate
#        rawvarframe = procimg.variance_frame(self.datasec_img, self.image,
#                                             self.detector['gain'], self.ronoise,
#                                             # TODO: dark current is expected to
#                                             # be in e-/s.
#                                             darkcurr=self.detector['darkcurr'],
#                                             # TODO: exptime is expected to be
#                                             # in seconds
#                                             exptime=self.exptime,
#                                             rnoise=self.rn2img)
#        return utils.inverse(rawvarframe)

    def estimate_readnoise(self):
        """
        Estimate the readnoise (in electrons) based on the overscan regions of
        the image.

        If the readnoise is not known for any of the amplifiers (i.e., if
        :attr:`ronoise` is :math:`\leq 0`), the function estimates it using the
        standard deviation in the overscan region.

        .. warning::

            This function edits :attr:`ronoise` in place.
        """
        # TODO: For what instruments do we use this?  For those instruments,
        # does this provide a robust measure of the readnoise?  If we calculate
        # the readnoise in this way, we shouldn't be including additional
        # digitization noise in detector variance image...
        for amp in range(len(self.ronoise)):
            if self.ronoise[amp] > 0:
                continue
            # TODO: What happens where there is no overscan region?  Do all the
            # supported instruments include overscan sections?
            biaspix = self.rawimage[self.oscansec_img==amp+1] * self.detector['gain'][amp]
            self.ronoise[amp] = stats.sigma_clipped_stats(biaspix, sigma=5)[-1]
            msgs.info(f'Estimated readnoise of amplifier {amp+1} = {self.ronoise[amp]:.3f} e-')

    def build_rn2img(self, units='e-', digitization=True):
        """
        Generate the model detector variance image (:attr:`rn2img`).

        This is primarily a wrapper for :func:`~pypeit.core.procimg.rn2_frame`.

        Args:
            units (:obj:`str`, optional):
                Units for the output variance.  Options are ``'e-'`` for
                variance in square electrons (counts) or ``'ADU'`` for square
                ADU.

        Returns:
            `numpy.ndarray`_: Readnoise variance image.
        """
        if not np.all(self.ronoise > 0):
            # TODO: Consider just calling estimate_readnoise here...
            msgs.error('Some readnoise values <=0; first call estimate_readnoise.')
        # Compute and return the readnoise variance image 
        return procimg.rn2_frame(self.datasec_img, self.detector['gain'], self.ronoise,
                                 units=units, digitization=digitization)

    def process(self, par, bpm=None, flatimages=None, bias=None, slits=None, debug=False,
                dark=None):
        """
        Process the image.

        See further discussion of :ref:`image_proc` in ``PypeIt``.

        .. note::

            The processing step order is currently 'frozen' as is.  However, in
            the future, we may choose to allow optional ordering.

        The processing steps used (depending on the parameter toggling in
        :attr:`par`), in the order they will be applied are:

            #. :func:`subtract_pattern`: Analyze and subtract the pattern noise
               from the image.

            #. :func:`build_rn2img`: Construct the detector variance image,
               which includes readnoise and digitization error.  If any of the
               amplifiers on the detector do not have a measured readnoise, the
               readnoise is estimated using :func:`estimate_readnoise`.
            
            #. :func:`subtract_overscan`: Analyze the overscan region and
               subtract from the image

            #. :func:`trim`: Trim the image down to the data (i.e. remove the
               overscan)

            #. :func:`orient`: Orient the image in the PypeIt orientation ---
               spectral coordinates ordered along the first axis and spatial
               coordinates ordered along the second, ``(spec, spat)`` --- with
               blue to red going from small pixel numbers to large.

            #. :func:`subtract_bias`: Subtract a bias image.  The shape and
               orientation of the bias image must match the *processed* image.
               I.e., if you trim and orient this image, you must also have
               trimmed and oriented the bias frames.  The units of the bias
               image *must* be ADU.

            #. :func:`apply_gain`: Convert from ADU to electrons, amp by amp.

            #. :func:`subtract_dark`: Subtract a dark image.  The shape and
               orientation of the dark image must match the *processed* image.
               I.e., if you trim and orient this image, you must also have
               trimmed and oriented the dark frames.  The units of the dark
               frame *must* be electrons (counts).

            #. :func:`add_shot_noise`: Include shot-noise from the image
               variance.

            #. :func:`impose_noise_floor`: Add uncertainty to the variance image
               such that the S/N per pixel is never reaches higher than
               ``1/noise_floor``.

            #. :func:`spatial_flexure_shift`: Measure any spatial shift due to flexure.

            #. :func:`flatten`: Divide by the pixel-to-pixel, spatial and
               spectral response functions.

            #. :func:`~pypeit.images.pypeitimage.PypeItImage.build_crmask`:
               Generate a cosmic-ray mask

        Args:
            par (:class:`~pypeit.par.pypeitpar.ProcessImagesPar`):
                Parameters that dictate the processing of the images.  See
                :class:`pypeit.par.pypeitpar.ProcessImagesPar` for the
                defaults.
            bpm (`numpy.ndarray`_, optional):
                The bad-pixel mask.  This is used to *overwrite* the default
                bad-pixel mask for this spectrograph.  The shape must match a
                trimmed and oriented processed image.
            flatimages (:class:`~pypeit.flatfield.FlatImages`, optional):
                Flat-field images used to apply flat-field corrections.
            bias (:class:`~pypeit.images.pypeitimage.PypeItImage`, optional):
                Bias image for bias subtraction.
            slits (:class:`~pypeit.slittrace.SlitTraceSet`, optional):
                Used to calculate spatial flexure between the image and the
                slits, if requested via the ``spat_flexure_correct`` parameter
                in :attr:`par`; see
                :func:`~pypeit.core.flexure.spat_flexure_shift`.  Also used to
                construct the slit illumination profile, if requested via the
                ``use_illumflat`` parameter in :attr:`par`; see
                :func:`~pypeit.flatfield.FlatImages.fit2illumflat`.
            debug (:obj:`bool`, optional):
                Run in debug mode.
            dark (:class:`~pypeit.images.pypeitimage.PypeItImage`):
                Dark image

        Returns:
            :class:`~pypeit.images.pypeitimage.PypeItImage`:

        """
        # TODO: Reset steps to all be false at the beginning of the function?
        # If we're careful with book-keeping of the attributes (i.e., resetting
        # image to rawimage, etc), calling process multiple times could be a way
        # for developers to *re*process an image multiple times from scratch
        # with different parameters to test changes.

        # Set input to attributes
        self.par = par
        if bpm is not None:
            self._bpm = bpm

        # Check the input
        if self.par['use_biasimage'] and bias is None:
            msgs.error('No bias available for bias subtraction!')
        if self.par['use_darkimage'] and dark is None:
            msgs.error('No dark available for dark subtraction!')
        if self.par['spat_flexure_correct'] and slits is None:
            msgs.error('Spatial flexure correction requested but no slits provided.')
        if self.use_flat and flatimages is None:
            msgs.error('Flat-field corrections requested but no flat-field images generated '
                       'or provided.  Make sure you have flat-field images in your PypeIt file!')

        # Apply additive corrections.  The order matters and is fixed.
        #   - Subtract any fixed pattern defined for the instrument.  NOTE: This
        #     step *must* be done before use_overscan
        if self.par['use_pattern']:
            self.subtract_pattern()

        #   - Estimate the readnoise using the overscan regions.  NOTE: I'm
        #     assuming this needs to be done *after* the sine pattern has been
        #     subtracted.
        if not np.all(self.ronoise > 0):
            self.estimate_readnoise()

        #   - Get detector variance image.  Starts with the shape and
        #     orientation of the raw image, and the units are in ADU to match
        #     the current units of the image being processed.  Note that the
        #     readnoise variance image is kept separate because it is used again
        #     when adding the shot noise.
        self.rn2img = self.build_rn2img(units='ADU')
        self.var = self.rn2img.copy()

        #   - Subtract the overscan.  Uncertainty from the overscan subtraction
        #     is added to the variance.
        if self.par['use_overscan']:
            self.subtract_overscan()

        # TODO: Do we need to keep trim and orient as optional?
        #   - Trim to the data region.  This trims the image, rn2img, var, and
        #     datasec_img.
        if self.par['trim']:
            self.trim()

        #   - Re-orient to PypeIt convention. This re-orients the image, rn2img,
        #     var, and datasec_img.
        if self.par['orient']:
            self.orient()
            
        #   - Check the shape of the bpm
        if self.bpm.shape != self.image.shape:
            msgs.error(f'Bad-pixel mask has incorrect shape: found {self.bpm.shape}, expected '
                       f'{self.image.shape}.  The shape must match a trimmed and oriented '
                       'PypeIt-processed image.')
            
        #   - Subtract master bias
        if self.par['use_biasimage']:
            # Bias frame.  Shape and orientation must match *processed* image,
            # and the units *must* be ADU.  Uncertainty from the bias
            # subtraction is added to the variance.
            self.subtract_bias(bias)
            
        #   - Convert from ADU to electrons/counts.  This converts the image,
        #     the variance, and the detector variance (because the latter is
        #     used by `add_shot_noise`).  TODO: I don't think there's any
        #     particular reason why this could not be done early, but it's
        #     important that it happen before we subtract the dark so that the
        #     dark can have a different exposure time.
        if self.par['apply_gain']:
            self.apply_gain()
            
        #   - Subtract master dark.
        if self.par['use_darkimage']:
            # Dark frame.  Shape and orientation must match *processed* image,
            # and the units *must* be in electrons/counts.  Uncertainty from the
            # dark subtraction is added to the variance.
            # TODO: With exptime now included in PypeItImage, we can explicitly
            # check and/or account for the difference in exposure time between
            # the dark image and the image being processed.
            self.subtract_dark(dark)
        
        # Use the bias-subtracted and dark-subtracted counts to include shot
        # noise in the error budget.
        if self.par['shot_noise']:
            self.add_shot_noise()
        
        # Impose a noise floor.  This does nothing if self.par['noise_floor'] is
        # not greater than 0.
        self.impose_noise_floor()

        # TODO: Shouldn't the cosmic-ray detection/masking happen here, before
        # flat-fielding?

        # Calculate flexure, if slits are provided and the correction is
        # requested.  NOTE: This step must come after trim, orient (just like
        # bias and dark subtraction) and before field flattening.
        self.spat_flexure_shift = None if slits is None or not self.par['spat_flexure_correct'] \
                                    else self.spatial_flexure_shift(slits)

        if self.use_flat:
            # Flat-field the data.  This propagates the flat-fielding
            # corrections to the variance.
            self.flatten(flatimages, slits=slits, debug=debug)

        # Calculate the inverse variance
        self.ivar = utils.inverse(self.var)

        # Generate a PypeItImage
        pypeitImage = pypeitimage.PypeItImage(self.image, ivar=self.ivar, rn2img=self.rn2img,
                                              bpm=self.bpm, detector=self.detector,
                                              spat_flexure=self.spat_flexure_shift,
                                              PYP_SPEC=self.spectrograph.name,
                                              units='e-' if self.par['apply_gain'] else 'ADU',
                                              exptime=self.exptime)
        pypeitImage.rawheadlist = self.headarr
        pypeitImage.process_steps = [key for key in self.steps.keys() if self.steps[key]]

        # Mask(s)
        if self.par['mask_cr']:
            pypeitImage.build_crmask(self.par)
        nonlinear_counts = self.spectrograph.nonlinear_counts(self.detector,
                                                              apply_gain=self.par['apply_gain'])
        pypeitImage.build_mask(saturation=nonlinear_counts)

        # Return
        return pypeitImage

    def spatial_flexure_shift(self, slits, force=False):
        """
        Calculate a spatial shift in the edge traces due to flexure.

        This is a simple wrapper for
        :func:`~pypeit.core.flexure.spat_flexure_shift`.

        Args:
            slits (:class:`~pypeit.slittrace.SlitTraceSet`, optional):
                Slit edge traces
            force (:obj:`bool`, optional):
                Force the image to be field flattened, even if the step log
                (:attr:`steps`) indicates that it already has been.

        """
        step = inspect.stack()[0][3]
        if self.steps[step] and not force:
            # Already field flattened
            msgs.warn('Spatial flexure shift already calculated.')
            return
        self.spat_flexure_shift = flexure.spat_flexure_shift(self.image, slits)
        self.steps[step] = True

    def flatten(self, flatimages, slits=None, force=False, debug=False):
        """
        Field flatten the processed image.

        This method calculates a slit-illumination correction and a spectral
        illumination correction, as well as the pixel-to-pixel correction, and
        multiplicatively removes them from the current image.  If available, the
        calculation is propagated to the variance image; however, no uncertainty
        in the flat-field corrections are included.

        .. warning::

            If you want the spatial flexure to be accounted for, you must first
            calculate the shift using
            :func:`~pypeit.images.rawimage.RawImage.spatial_flexure_shift`.

        Args:
            flatimages (:class:`~pypeit.flatfield.FlatImages`):
                Flat-field images used to apply flat-field corrections.
            slits (:class:`~pypeit.slittrace.SlitTraceSet`, optional):
                Used to construct the slit illumination profile, and only 
                required if this is to be calculated and normalized out.  See
                :func:`~pypeit.flatfield.FlatImages.fit2illumflat`.
            force (:obj:`bool`, optional):
                Force the image to be field flattened, even if the step log
                (:attr:`steps`) indicates that it already has been.
            debug (:obj:`bool`, optional):
                Run in debug mode.
        """
        step = inspect.stack()[0][3]
        if self.steps[step] and not force:
            # Already field flattened
            msgs.warn('Image was already flat fielded.')
            return

        # Check input
        if flatimages.pixelflat_norm is None:
            # We cannot do any flat-field correction without a pixel flat (yet)
            msgs.error("Flat fielding desired but not generated/provided.")
        if self.par['use_illumflat'] and slits is None:
            msgs.error('Need to provide slits to create illumination flat.')
        if self.par['use_specillum'] and flatimages.pixelflat_spec_illum is None:
            msgs.error("Spectral illumination correction desired but not generated/provided.")

        # Generate the illumination flat, as needed
        illum_flat = 1.0
        if self.par['use_illumflat']:
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
        spec_illum = 1.
        if self.par['use_specillum']:
            # TODO: Why is the convention for spec_illum different from the
            # other two flat-field corrections?  I.e., is it worth rethinking
            # the definition of spec_illum to be the illumination profile itself
            # instead of the scale factor required to remove the illumination
            # profile? cf. pypeit.flatfield.FlatField.spectral_illumination
            spec_illum = utils.inverse(flatimages.pixelflat_spec_illum)

        # Apply flat-field correction
        ret = flat.flatfield(self.image, flatimages.pixelflat_norm * illum_flat * spec_illum,
                             varframe=self.var)
        if self.var is None:
            self.image = ret
        else:
            self.image, self.var = ret

        # TODO: Include image pixels that are 0 in the bpm?

        # TODO: Recreating the bpm can be an expensive operation for some
        # spectrographs because it has to open the file and get the
        # trimmed/oriented shape.  We should avoid having to recreate it.  So I
        # changed the call to self.flatten so that the bpm is no longer passed
        # to it (it wasn't used, only altered), but we need to keep track of
        # this if/when edits are made to this class.
#        # Fresh BPM
#        self._bpm = self.spectrograph.bpm(self.filename, self.det, shape=self.image.shape)
        self.steps[step] = True

    def orient(self, force=False):
        """
        Orient :attr:`image`, :attr:`var` (if it exists), :attr:`rn2img` (if it
        exists), and :attr:`datasec_img` such that they follow the ``PypeIt``
        convention with spectra running blue (down) to red (up) and with orders
        decreasing from high (left) to low (right).

        Args:
            force (:obj:`bool`, optional):
                Force the image to be re-oriented, even if the step log
                (:attr:`steps`) indicates that it already has been.
        """
        step = inspect.stack()[0][3]
        # Check if already oriented
        if self.steps[step] and not force:
            msgs.warn('Image was already oriented.')
            return
        # Orient the image to have blue/red run bottom to top
        self.image = self.spectrograph.orient_image(self.detector, self.image)
        if self.var is not None:
            self.var = self.spectrograph.orient_image(self.detector, self.var)
        if self.rn2img is not None:
            self.rn2img = self.spectrograph.orient_image(self.detector, self.rn2img)
        self.datasec_img = self.spectrograph.orient_image(self.detector, self.datasec_img)
        self.steps[step] = True

    def subtract_bias(self, bias_image, force=False):
        """
        Subtract a bias image.
        
        If the ``bias_image`` object includes an inverse variance image and if
        :attr:`var` is available, the error in the bias is propagated to the
        bias-subtracted image.

        Args:
            bias_image (:class:`~pypeit.images.pypeitimage.PypeItImage`):
                Bias image
            force (:obj:`bool`, optional):
                Force the image to be subtracted, even if the step log
                (:attr:`steps`) indicates that it already has been.
        """
        step = inspect.stack()[0][3]
        if self.steps[step] and not force:
            # Already bias subtracted
            msgs.warn('Image was already bias subtracted.')
            return
        self.image -= bias_image.image
        # TODO: Also incorporate the mask?
        if self.var is not None and bias_image.ivar is not None:
            self.var += utils.inverse(bias_image.ivar)
        self.steps[step] = True

    # TODO: expscale is currently not a parameter that the user can control.
    # Should it be?
    def subtract_dark(self, dark_image, expscale=True, force=False):
        """
        Subtract a dark image.

        The dark image is expected to have been bias subtracted and be in
        electon counts.  The dark image is automatically rescaled to the
        exposure time of this image; see ``expscale``.

        If the ``dark_image`` object includes an inverse variance image and if
        :attr:`var` is available, the error in the dark is propagated to the
        dark-subtracted image.

        Args:
            dark_image (:class:`~pypeit.images.pypeitimage.PypeItImage`):
                Dark image
            expscale (:obj:`bool`, optional):
                Scale the dark image by the ratio of the exposure times so that
                the counts per second represented by the dark image are
                appropriately subtracted from image being processed.
            force (:obj:`bool`, optional):
                Force the image to be subtracted, even if the step log
                (:attr:`steps`) indicates that it already has been.
        """
        step = inspect.stack()[0][3]
        if self.steps[step] and not force:
            # Already bias subtracted
            msgs.warn('Image was already dark subtracted.')
            return
        scale = self.exptime / dark_image.exptime if expscale else 1.
        self.image -= dark_image.image*scale
        # TODO: Also incorporate the mask?
        if self.var is not None and dark_image.ivar is not None:
            self.var += scale**2 * utils.inverse(dark_image.ivar)
        self.steps[step] = True

    def subtract_overscan(self, force=False):
        """
        Analyze and subtract the overscan from the image

        Args:
            force (:obj:`bool`, optional):
                Force the image to be overscan subtracted, even if the step log
                (:attr:`steps`) indicates that it already has been.
        """
        step = inspect.stack()[0][3]
        if self.steps[step] and not force:
            # Already overscan subtracted
            msgs.warn("Image was already overscan subtracted!")
            return
        ret = procimg.subtract_overscan(self.image, self.datasec_img, self.oscansec_img,
                                        method=self.par['overscan_method'],
                                        params=self.par['overscan_par'], var=self.var)
        if self.var is None:
            self.image = ret
        else:
            self.image, self.var = ret
        self.steps[step] = True

    def subtract_pattern(self):
        """
        Analyze and subtract the pattern noise from the image.
        """
        step = inspect.stack()[0][3]
        if self.steps[step]:
            # Already pattern subtracted
            msgs.warn("Image was already pattern subtracted!")
            return

        # Grab the frequency, if it exists in the header.  For some instruments,
        # PYPEITFRQ is added to the header in get_rawimage() in the spectrograph
        # file.  See keck_kcwi.py for an example.
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
        # TODO: Why does this overwrite the raw image?  Is it so that the
        # readnoise estimation doesn't include the sine pattern?  Does the new
        # order of forcing the calculation of the readnoise directly after
        # subtracting the pattern get around this?
        # TODO: Since subtract_pattern effectively must happen first, does it
        # make sense to pass self.rawimage instead of self.image?
        # TODO: Is there an error estimate in the pattern subtraction that we
        # can propagate to the error budget?
        self.rawimage = procimg.subtract_pattern(self.image, self.datasec_img, self.oscansec_img,
                                                 frequency=frequency)
        self.image = self.rawimage.copy()
        self.steps[step] = True

    def trim(self, force=False):
        """
        Trim :attr:`image`, :attr:`var` (if it exists), :attr:`rn2img` (if it
        exists), and :attr:`datasec_img` to include only the science data.

        Args:
            force (:obj:`bool`, optional):
                Force the image to be trimmed, even if the step log
                (:attr:`steps`) indicates that it already has been.
        """
        step = inspect.stack()[0][3]
        if self.rawimage.shape is not None and self.image.shape != self.rawimage.shape:
            # Image *must* have been trimmed already because shape does not
            # match raw image
            self.steps[step] = True
            msgs.warn('Image shape does not match original.')
            return
        if self.steps[step] and not force:
            # Already trimmed
            msgs.warn('Image was already trimmed.')
            return
        self.image = procimg.trim_frame(self.image, self.datasec_img < 1)
        if self.var is not None:
            self.var = procimg.trim_frame(self.var, self.datasec_img < 1)
        if self.rn2img is not None:
            self.rn2img = procimg.trim_frame(self.rn2img, self.datasec_img < 1)
        # NOTE: This must be done last because the untrimmed image is used for
        # deciding what to trim!
        self.datasec_img = procimg.trim_frame(self.datasec_img, self.datasec_img < 1)
        self.steps[step] = True

    def add_shot_noise(self):
        r"""
        Use the image to calculate and add shot noise to the error budget.

        The method assumes the image is bias- and dark-subtracted, and that the
        units are electrons/counts.  The calculation of the shot noise includes
        a term that adjusts the variance at low count levels to account for the
        Gaussian approximation of a Poisson distribution throughout the rest of
        the code base (*need a reference for this*).  If availabile, the
        calculation also includes shot noise from the dark current, using the
        exposure time and the dark current (e-/s) provided by the spectrograph
        detector (:class:`~pypeit.images.detector_container.DetectorContainer`)
        object.  In detail, the added shot noise variance is:

        .. math::

            V = | C - \sqrt(2 V_{\rm rn}) | + C_{\rm dark} t_{\rm exp}

        where :math:`V_{\rm rn}` is the readnoise variance (see
        :func:`~pypeit.core.procimg.rn2_frame`), :math:`C` is the number of
        bias- and dark-subtracted counts, :math:`C_{\rm dark}` is the number of
        dark counts per second, and :math:`t_{\rm exp}` is the exposure time
        (:attr:`exptime`) in seconds.

        .. note::

            The variance attribute (:attr:`var`) is edited in-place.

        """
        self.var += np.absolute(self.image - np.sqrt(2*self.rn2img))
        if self.detector['darkcurr'] > 0:
            # TODO: This approach assumes that the dark-current provided/known
            # for each detector is more accurate and precise than using the
            # counts in the dark image itself (assuming one is available) to
            # calculate the dark-current shot noise.
            self.var += self.detector['darkcurr'] * self.exptime

    def impose_noise_floor(self):
        """
        Impose a noise floor by inflating the error budget.

        If the :attr:`var` is not None and ``noise_floor`` value in :attr:`par`
        is greater than 0, the quantity ``(noise_floor*self.image)**2`` is added
        to :attr:`var`.  This effectively sets a maximum S/N per pixel.
        """
        if self.var is not None and self.par['noise_floor'] > 0:
            self.var += (self.par['noise_floor'] * self.image)**2

    def __repr__(self):
        return ('<{:s}: file={}, steps={}>'.format(self.__class__.__name__, self.filename,
                                                   self.steps))


