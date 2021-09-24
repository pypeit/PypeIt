""" Object to load and process a single raw image

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import os
import inspect
from copy import deepcopy

from IPython import embed

import numpy as np

from astropy import stats

from pypeit import msgs
from pypeit.core import parse
from pypeit.core import procimg
from pypeit.core import flat
from pypeit.core import flexure
from pypeit.images import pypeitimage
from pypeit import utils
from pypeit.display import display

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
        proc_var (`numpy.ndarray`_):
            The sum of the variance components added by the image processing;
            i.e., the error in the overscan subtraction, bias subtraction, etc.
        base_var (`numpy.ndarray`_):
            The base-level variance in the processed image.  See
            :func:`~pypeit.core.procimg.base_variance`.
        var (`numpy.ndarray`_):
            The aggregate variance in the processed image.  This is the primary
            array used during :func:`process` to track uncertainties;
            :attr:`ivar` is created by inverting this at the end of the
            processing method.
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
        self.rn2img = None
        self.proc_var = None
        self.base_var = None
        self.spat_flexure_shift = None
        self.img_scale = None
        self._bpm = None

        # All possible processing steps.  NOTE: These have to match the
        # method/function names.  Their order here matches there execution order
        # in self.process(), but that's not necessary.
        self.steps = dict(apply_gain=False,
                          subtract_pattern=False,
                          subtract_overscan=False,
                          trim=False,
                          orient=False,
                          subtract_bias=False,
                          subtract_dark=False,
                          spatial_flexure_shift=False,
                          flatfield=False)

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

        # Behavior needs to be different if the image has been trimmed or not
        gain = procimg.gain_frame(self.datasec_img, np.atleast_1d(self.detector['gain']))
        # TODO: Could also check if steps['trim'] is true.  Is either better or worse?
        if self.rawimage.shape is not None and self.image.shape == self.rawimage.shape:
            # Image is raw, so need to include overscan sections
            # TODO: This only works assuming that there is *no* overlap between
            # oscansec_img and datasec_img.  There shouldn't be, but the code
            # doesn't check this...
            gain += procimg.gain_frame(self.oscansec_img, np.atleast_1d(self.detector['gain']))
        self.image *= gain
        # NOTE: In ``process``, ``apply_gain`` is called first, meaning that all
        # the variance arrays should be None.
        self.steps[step] = True

    def build_ivar(self):
        """
        Generate the inverse variance in the image.

        This is a simple wrapper for :func:`~pypeit.core.procimg.base_variance`
        and :func:`~pypeit.core.procimg.variance_model`.

        Returns:
            `numpy.ndarray`_: The inverse variance in the image.
        """
        # TODO: This approach assumes that the dark-current provided/known for
        # each detector is more accurate and precise than using the counts in
        # the dark image itself (assuming one is available) to calculate the
        # dark-current shot noise.  Instead, we could pass an image of the dark
        # current, if available.
        npix = np.prod(parse.parse_binning(self.detector.binning))
        darkcurr = npix*self.detector['darkcurr'] \
                    if self.detector['darkcurr'] > 0 and self.par['shot_noise'] else None
        self.base_var = procimg.base_variance(self.rn2img, darkcurr=darkcurr, exptime=self.exptime,
                                              proc_var=self.proc_var, count_scale=self.img_scale)
        var = procimg.variance_model(self.base_var,
                                     counts=self.image if self.par['shot_noise'] else None,
                                     count_scale=self.img_scale,
                                     noise_floor=self.par['noise_floor'])
        return utils.inverse(var)

    def estimate_readnoise(self):
        """
        Estimate the readnoise (in electrons) based on the overscan regions of
        the image.

        If the readnoise is not known for any of the amplifiers (i.e., if
        :attr:`ronoise` is :math:`\leq 0`) or if explicitly requested using the
        ``empirical_rn`` parameter, the function estimates it using the standard
        deviation in the overscan region.

        .. warning::

            This function edits :attr:`ronoise` in place.
        """
        if self.oscansec_img.shape != self.image.shape:
            msgs.error('Must estimate readnoise before trimming the image.')
        for amp in range(len(self.ronoise)):
            if self.ronoise[amp] > 0 and not self.par['empirical_rn']:
                # Skip if the readnoise is defined and the empirical readnoise
                # estimate was not explicitly requested.
                continue
            if not np.any(self.oscansec_img==amp+1):
                msgs.error(f'Cannot estimate readnoise for amplifier {amp+1}; no overscan region!')
            gain = 1. if self.steps['apply_gain'] else self.detector['gain'][amp]
            biaspix = self.image[self.oscansec_img==amp+1] * gain
            self.ronoise[amp] = stats.sigma_clipped_stats(biaspix, sigma=5)[-1]
            msgs.info(f'Estimated readnoise of amplifier {amp+1} = {self.ronoise[amp]:.3f} e-')

    def build_rn2img(self, units='e-', digitization=False):
        """
        Generate the model readnoise variance image (:attr:`rn2img`).

        This is primarily a wrapper for :func:`~pypeit.core.procimg.rn2_frame`.

        Args:
            units (:obj:`str`, optional):
                Units for the output variance.  Options are ``'e-'`` for
                variance in square electrons (counts) or ``'ADU'`` for square
                ADU.
            digitization (:obj:`bool`, optional):
                Include digitization error in the calculation.

        Returns:
            `numpy.ndarray`_: Readnoise variance image.
        """
        if not np.all(self.ronoise > 0):
            # TODO: Consider just calling estimate_readnoise here...
            msgs.error('Some readnoise values <=0; first call estimate_readnoise.')

        # Compute and return the readnoise variance image 
        rn2 = procimg.rn2_frame(self.datasec_img, self.ronoise, units=units,
                                gain=self.detector['gain'], digitization=digitization)
        # TODO: Could also check if steps['trim'] is true.  Is either better or worse?
        if self.rawimage.shape is not None and self.image.shape == self.rawimage.shape \
                and np.any(self.oscansec_img > 0):
            # Image is raw, so need to include overscan sections
            # TODO: This only works assuming that there is *no* overlap between
            # oscansec_img and datasec_img.  There shouldn't be, but the code
            # doesn't check this...
            rn2 += procimg.rn2_frame(self.oscansec_img, self.ronoise, units=units,
                                     gain=self.detector['gain'], digitization=digitization)
        return rn2

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

            #. :func:`apply_gain`: The first step is to convert the image units
               from ADU to electrons, amp by amp, using the gain provided by the
               :class:`~pypeit.images.detector_container.DetectorContainer`
               instance(s) for each
               :class:`~pypeit.spectrographs.spectrograph.Spectrograph`
               subclass.

            #. :func:`subtract_pattern`: Analyze and subtract sinusoidal pattern
               noise from the image; see
               :func:`~pypeit.core.procimg.subtract_pattern`.

            #. :func:`build_rn2img`: Construct the readnoise variance image,
               which includes readnoise and digitization error.  If any of the
               amplifiers on the detector do not have a measured readnoise or if
               explicitly requested using the ``empirical_rn`` parameter, the
               readnoise is estimated using :func:`estimate_readnoise`.
            
            #. :func:`subtract_overscan`: Use the detector overscan region to
               measure and subtract the frame-dependent bias level along the
               readout direction.

            #. :func:`trim`: Trim the image to include the data regions only
               (i.e. remove the overscan).

            #. :func:`orient`: Orient the image in the PypeIt orientation ---
               spectral coordinates ordered along the first axis and spatial
               coordinates ordered along the second, ``(spec, spat)`` --- with
               blue to red going from small pixel numbers to large.

            #. :func:`subtract_bias`: Subtract the master bias image.  The shape
               and orientation of the bias image must match the *processed*
               image.  I.e., if you trim and orient this image, you must also
               have trimmed and oriented the bias frames.

            #. :func:`subtract_dark`: Subtract the master dark image.  The shape
               and orientation of the dark image must match the *processed*
               image.  I.e., if you trim and orient this image, you must also
               have trimmed and oriented the dark frames.  The dark image is
               *automatically* scaled by the ratio of the exposure times to
               ensure the counts/s in the dark are removed from the image being
               processed.

            #. :func:`spatial_flexure_shift`: Measure any spatial shift due to
               flexure.

            #. :func:`flatfield`: Divide by the pixel-to-pixel, spatial and
               spectral response functions.

            #. :func:`build_ivar`: Construct a model estimate of the variance in
               the image based in the readnoise, errors from the additive
               processing steps, shot-noise from the observed counts (see the
               ``shot_noise`` parameter), a rescaling due to the flat-field
               correction, and a noise floor that sets a maximum S/N per pixel
               (see the ``noise_floor`` parameter); see
               :func:`~pypeit.core.procimg.variance_model`.
            
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
            :class:`~pypeit.images.pypeitimage.PypeItImage`: The processed
            image.
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

        #   - Convert from ADU to electron counts.
        if self.par['apply_gain']:
            self.apply_gain()

        # Apply additive corrections.  The order matters and is fixed.
        #
        #   - Subtract any fixed pattern defined for the instrument.  NOTE: This
        #     step *must* be done before subtract_overscan
        if self.par['use_pattern']:
            self.subtract_pattern()

        #   - Estimate the readnoise using the overscan regions.  NOTE: This
        #     needs to be done *after* any pattern subtraction.
        if self.par['empirical_rn'] or not np.all(self.ronoise > 0):
            self.estimate_readnoise()

        #   - Initialize the variance images.  Starts with the shape and
        #     orientation of the raw image.  Note that the readnoise variance
        #     and "process" variance images are kept separate because they are
        #     used again when modeling the noise during object extraction.
        self.rn2img = self.build_rn2img()
        self.proc_var = np.zeros(self.rn2img.shape, dtype=float)

        #   - Subtract the overscan.  Uncertainty from the overscan subtraction
        #     is added to the variance.
        if self.par['use_overscan']:
            self.subtract_overscan()

        # TODO: Do we need to keep trim and orient as optional?

        #   - Trim to the data region.  This trims the image, rn2img, proc_var,
        #     var, and datasec_img.
        if self.par['trim']:
            self.trim()

        #   - Re-orient to PypeIt convention. This re-orients the image, rn2img,
        #     proc_var, var, and datasec_img.
        if self.par['orient']:
            self.orient()

        #   - Check the shape of the bpm
        if self.bpm.shape != self.image.shape:
            # The BPM is the wrong shape.  Assume this is because the
            # calibrations were taken with a different binning than the science
            # data, and assume this is because a BPM was provided as an
            # argument.  First erase the "protected" attribute, then access it
            # again using the @property method, which will recreated it based on
            # the expected shape for this frame.
            bpm_shape = self.bpm.shape                  # Save the current shape for the warning
            self._bpm = None                            # This erases the current bpm attribute
            if self.bpm.shape != self.image.shape:      # This recreates it
                # This should only happen because of a coding error, not a user error
                msgs.error(f'CODING ERROR: From-scratch BPM has incorrect shape!')
            # If the above was successful, the code can continue, but first warn
            # the user that the code ignored the provided bpm.
            msgs.warn(f'Bad-pixel mask has incorrect shape: found {bpm_shape}, expected '
                      f'{self.image.shape}.  Assuming this is because different binning used for '
                      'various frames.  Recreating BPM specifically for this frame '
                      f'({os.path.basename(self.filename)}) and assuming the difference in the '
                      'binning will be handled later in the code.')
            
        #   - Subtract master bias
        if self.par['use_biasimage']:
            # Bias frame.  Shape and orientation must match *processed* image,.
            # Uncertainty from the bias subtraction is added to the variance.
            self.subtract_bias(bias)

        #   - Subtract dark current.  The tabulated dark current is *always*
        #     subtracted regardless of whether or not use_darkimage is true.  If
        #     a dark frame is provided and subtracted, its shape and orientation
        #     must match *processed* image, and the units *must* be in
        #     electrons/counts.  If available, uncertainty from the dark
        #     subtraction is added to the variance.
        self.subtract_dark(dark_image=dark if self.par['use_darkimage'] else None)

        # TODO: We should be performing the cosmic-ray detection/masking here,
        # *before* flat-fielding, right?.

        # Calculate flexure, if slits are provided and the correction is
        # requested.  NOTE: This step must come after trim, orient (just like
        # bias and dark subtraction) and before field flattening.
        self.spat_flexure_shift = None if slits is None or not self.par['spat_flexure_correct'] \
                                    else self.spatial_flexure_shift(slits)

        # Flat-field the data.  This propagates the flat-fielding corrections to
        # the variance.  The returned bpm is propagated to the PypeItImage
        # bitmask below.
        flat_bpm = self.flatfield(flatimages, slits=slits, debug=debug) if self.use_flat else None

        # Calculate the inverse variance
        # TODO: This approach assumes that the dark-current provided/known for
        # each detector is more accurate and precise than using the counts in
        # the dark image itself (assuming one is available) to calculate the
        # dark-current shot noise.  We could pass ``darkcurr=dark`` if ``dark``
        # is available instead...
        self.ivar = self.build_ivar()

        # Generate a PypeItImage.
        # NOTE: To reconstruct the variance model, you need base_var, image,
        # img_scale, noise_floor, and shot_noise.
        pypeitImage = pypeitimage.PypeItImage(self.image, ivar=self.ivar, rn2img=self.rn2img,
                                              base_var=self.base_var, img_scale=self.img_scale,
                                              bpm=self.bpm, detector=self.detector,
                                              spat_flexure=self.spat_flexure_shift,
                                              PYP_SPEC=self.spectrograph.name,
                                              units='e-' if self.par['apply_gain'] else 'ADU',
                                              exptime=self.exptime,
                                              noise_floor=self.par['noise_floor'],
                                              shot_noise=self.par['shot_noise'])
        pypeitImage.rawheadlist = self.headarr
        pypeitImage.process_steps = [key for key in self.steps.keys() if self.steps[key]]

        # Mask(s)
        if self.par['mask_cr']:
            pypeitImage.build_crmask(self.par)
        nonlinear_counts = self.spectrograph.nonlinear_counts(self.detector,
                                                              apply_gain=self.par['apply_gain'])
        pypeitImage.build_mask(saturation=nonlinear_counts)
        if flat_bpm is not None:
            pypeitImage.update_mask('BADSCALE', indx=flat_bpm)

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

    def flatfield(self, flatimages, slits=None, force=False, debug=False):
        """
        Field flatten the processed image.

        This method uses the results of the flat-field modeling code (see
        :class:`~pypeit.flatfield.FlatField`) and any measured spatial shift due
        to flexure to construct slit-illumination, spectral response, and
        pixel-to-pixel response corrections, and multiplicatively removes them
        from the current image.  If available, the calculation is propagated to
        the variance image; however, no uncertainty in the flat-field
        corrections are included.

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

        Returns:
            `numpy.ndarray`_: Returns a boolean array flagging pixels were the
            total applied flat-field value (i.e., the combination if the
            pixelflat and illumination corrections) was <=0.
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
        spec_illum = flatimages.pixelflat_spec_illum if self.par['use_specillum'] else 1.

        # Apply flat-field correction
        # NOTE: Using flat.flatfield to effectively multiply image*img_scale is
        # a bit overkill...
        total_flat = flatimages.pixelflat_norm * illum_flat * spec_illum
        self.img_scale = utils.inverse(total_flat)
        self.image, flat_bpm = flat.flatfield(self.image, total_flat)
        self.steps[step] = True
        return flat_bpm

    def orient(self, force=False):
        """
        Orient image attributes such that they follow the ``PypeIt`` convention
        with spectra running blue (down) to red (up) and with orders decreasing
        from high (left) to low (right).
        
        This edits :attr:`image`, :attr:`rn2img` (if it exists),
        :attr:`proc_var` (if it exists), and :attr:`datasec_img` in place.

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
        if self.rn2img is not None:
            self.rn2img = self.spectrograph.orient_image(self.detector, self.rn2img)
        if self.proc_var is not None:
            self.proc_var = self.spectrograph.orient_image(self.detector, self.proc_var)
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
        if bias_image.ivar is not None and self.proc_var is not None:
            self.proc_var += utils.inverse(bias_image.ivar)
        self.steps[step] = True

    # TODO: expscale is currently not a parameter that the user can control.
    # Should it be?
    def subtract_dark(self, dark_image=None, expscale=True, force=False):
        """
        Subtract detector dark current.

        This method subtracts both the tabulated dark-current for the detector
        and a dark image.  For this to be appropriate, the dark image (if
        provided) *must* also have had the tabulated dark-current value
        subtracted from it.

        Also, the processing of the dark image (if provided) should match the
        processing of the image being processed.  For example, if this image has
        been bias subtracted, so should be the dark image.

        If the ``dark_image`` object includes an inverse variance estimate and
        if :attr:`proc_var` is available, the uncertainty in the dark will be
        propagated to the dark-subtracted image.

        .. warning::

            The current default automatically scales the dark frame to match the
            exposure time of the image being processed.  Typically dark frames
            should have the same exposure time as the image being processed, so
            this will have no effect.  However, beware if that's not the case,
            and make sure this scaling is appropriate.  Use ``expscale`` to
            turn it off.

        Args:
            dark_image (:class:`~pypeit.images.pypeitimage.PypeItImage`, optional):
                The dark image to subtract.  If None, only the tabulated
                dark-current is subtracted from *all* pixels.
            expscale (:obj:`bool`, optional):
                Scale the dark image (if provided) by the ratio of the exposure
                times so that the counts per second represented by the dark
                image are appropriately subtracted from image being processed.
            force (:obj:`bool`, optional):
                Force the dark to be subtracted, even if the step log
                (:attr:`steps`) indicates that it already has been.
        """
        step = inspect.stack()[0][3]
        if self.steps[step] and not force:
            # Already bias subtracted
            msgs.warn('Image was already dark subtracted.')
            return
        # TODO: Is the dark-current amplifier dependent?  Also, this usage means
        # that darkcurr cannot be None.
        # Tabulated dark current is in e-/pixel/hour and exptime is in s, the
        # 3600 factor converts the dark current to e-/pixel/s.
        npix = np.prod(parse.parse_binning(self.detector.binning))
        dark_count = npix*self.detector['darkcurr'] * self.exptime / 3600.
        if dark_image is not None:
            # TODO: Include a warning when the scaling is "signficant"?
            scale = self.exptime / dark_image.exptime if expscale else 1.
            # Warn the user if the counts in the dark image are significantly
            # different from the tabulated value.  The 50% difference is not
            # well justified but, for now, hard-coded.
            med_dark = np.median(scale * dark_image.image)
            if np.absolute(med_dark) > 0.5*dark_count:
                msgs.warn(f'Dark-subtracted dark frame has significant signal remaining.  Median '
                          f'count is {med_dark:.2f}; warning threshold = +/-{0.5*dark_count:.2f}.')
            # NOTE: This converts dark_count from a scalar to an array
            dark_count += scale * dark_image.image
            if dark_image.ivar is not None and self.proc_var is not None:
                # Include the variance in the dark image
                # TODO: Also incorporate the mask?
                self.proc_var += scale**2 * utils.inverse(dark_image.ivar)
        self.image -= dark_count
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
        # Subtract the overscan.  var is the variance in the overscan
        # subtraction
        self.image, var = procimg.subtract_overscan(self.image, self.datasec_img,
                                                    self.oscansec_img,
                                                    method=self.par['overscan_method'],
                                                    params=self.par['overscan_par'],
                                                    var=self.rn2img)
        # Parse the returned value
        if self.rn2img is not None:
            # Include the variance in the overscan subtraction in the
            # "processing variance"
            self.proc_var += var
        self.steps[step] = True

    def subtract_pattern(self):
        """
        Analyze and subtract the pattern noise from the image.

        This is primarily a wrapper for
        :func:`~pypeit.core.procimg.subtract_pattern`.

        """
        step = inspect.stack()[0][3]
        if self.steps[step]:
            # Already pattern subtracted
            msgs.warn("Image was already pattern subtracted!")
            return

        # The image must have an overscan region for this to work.
        if not np.any(self.oscansec_img > 0):
            msgs.error('Image has no overscan region.  Pattern noise cannot be subtracted.')

        # The image cannot have already been trimmed
        if self.oscansec_img.shape != self.image.shape:
            msgs.error('Must estimate readnoise before trimming the image.')

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
        # Subtract the pattern and overwrite the current image
        self.image = procimg.subtract_pattern(self.image, self.datasec_img, self.oscansec_img,
                                              frequency=frequency)
        self.steps[step] = True

    def trim(self, force=False):
        """
        Trim image attributes to include only the science data.
        
        This edits :attr:`image`, :attr:`rn2img` (if it exists),
        :attr:`proc_var` (if it exists), and :attr:`datasec_img` in place.

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
        if self.rn2img is not None:
            self.rn2img = procimg.trim_frame(self.rn2img, self.datasec_img < 1)
        if self.proc_var is not None:
            self.proc_var = procimg.trim_frame(self.proc_var, self.datasec_img < 1)
        # NOTE: This must be done last because the untrimmed image is used for
        # deciding what to trim!
        self.datasec_img = procimg.trim_frame(self.datasec_img, self.datasec_img < 1)
        self.steps[step] = True

    def __repr__(self):
        return ('<{:s}: file={}, steps={}>'.format(self.__class__.__name__, self.filename,
                                                   self.steps))


