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
from pypeit.core import arc
from pypeit.core import parse
from pypeit.core import procimg
from pypeit.core import flat
from pypeit.core import flexure
from pypeit.core.mosaic import build_image_mosaic
from pypeit.images import pypeitimage
from pypeit import utils
from pypeit.display import display


# TODO: I don't understand why we have some of these attributes.  E.g., why do
# we need both hdu and headarr?
class RawImage:
    """
    Class to load and process raw images.

    Generally speaking, this class should only be used as follows:

    .. code-block:: python

        # Load the raw data and prepare the object
        rawImage = RawImage(file, spectrograph, det)
        pypeitImage = rawImage.process(par)

    modulo details of the keyword arguments in
    :func:`~pypeit.images.rawimage.RawImage.process`.  The class has many
    methods that handle each step of the processing but the order of these steps
    matters, meaning they're not guaranteed to succeed if done out of order.
    This is most relevant when processing multiple detector images into an image
    mosaic; see :func:`process`.

    Args:
        ifile (:obj:`str`):
            File with the data.
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            The spectrograph from which the data was collected.
        det (:obj:`int`, :obj:`tuple`):
            1-indexed detector(s) to read.  An image mosaic is selected using a
            :obj:`tuple` with the detectors in the mosaic, which must be one of
            the allowed mosaics returned by
            :func:`~pypeit.spectrographs.spectrograph.Spectrograph.allowed_mosaics`.

    Attributes:
        filename (:obj:`str`):
            Original file name with the data.
        spectrograph (:class:`~pypeit.spectrograph.spectrographs.Spectrograph`):
            Spectrograph instance with the instrument-specific properties and
            methods.
        det (:obj:`int`, :obj:`tuple`):
            1-indexed detector number(s); see class argument.
        detector (:class:`~pypeit.images.detector_container.DetectorContainer`, :class:`~pypeit.images.mosaic.Mosaic`):
            Mosaic/Detector characteristics
        rawimage (`numpy.ndarray`_):
            The raw, not trimmed or reoriented, image data for the detector(s).
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

        # NOTE: The binning is expected to be the same for all images in a
        # mosaic, but it's left to the raw image reader for each spectrograph.
        # See, e.g., gemini_gmos.

        # Number of loaded images (needed in case the raw image is used to
        # create a mosaic)
        self.nimg = 1 if self.rawimage.ndim == 2 else self.rawimage.shape[0]

        # Re-package so that, independent of whether or not the image is a
        # detector mosaic, the attributes are all organized with the number of
        # detectors along their first axis.
        if self.nimg > 1:
            self.mosaic = self.detector
            self.detector = self.mosaic.detectors
        else:
            self.mosaic = None
            self.detector = np.array([self.detector])
            self.rawimage = np.expand_dims(self.rawimage, 0)
            self.rawdatasec_img = np.expand_dims(self.rawdatasec_img, 0)
            self.oscansec_img = np.expand_dims(self.oscansec_img, 0)

        # Grab items from rawImage (for convenience and for processing)
        #   Could just keep rawImage in the object, if preferred
        # TODO: Why do we need to deepcopy this?
        self.headarr = deepcopy(self.spectrograph.get_headarr(self.hdu))

        # Key attributes
        self.image = self.rawimage.copy()
        self.datasec_img = self.rawdatasec_img.copy()
        # NOTE: Prevent estimate_readnoise() from altering self.detector using
        # deepcopy
        self.ronoise = np.array([deepcopy(d['ronoise']) for d in self.detector])

        # Attributes
        self.par = None
        self.ivar = None
        self.rn2img = None
        self.dark = None
        self.dark_var = None
        self.proc_var = None
        self.base_var = None
        self.spat_flexure_shift = None
        self.img_scale = None
        self.det_img = None
        self._bpm = None

        # All possible processing steps.  NOTE: These have to match the
        # method/function names.  Their order here matches there execution order
        # in self.process(), but that's not necessary.
        self.steps = dict(apply_gain=False,
                          subtract_pattern=False,
                          subtract_overscan=False,
                          subtract_continuum=False,
                          trim=False,
                          orient=False,
                          subtract_bias=False,
                          subtract_dark=False,
                          build_mosaic=False,
                          spatial_flexure_shift=False,
                          flatfield=False)

    @property
    def shape(self):
        return () if self.image is None else self.image.shape

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
            # TODO: Pass msbias if there is one?  Only if `bpm_usebias` is
            # true in the calibrations parameter set, but we don't have access
            # to that here.  Add it as a parameter of ProcessImagesPar?
            self._bpm = self.spectrograph.bpm(self.filename, self.det, shape=self.image.shape[1:])
            if self.nimg == 1:
                self._bpm = np.expand_dims(self._bpm, 0)
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

    @property
    def use_slits(self):
        """
        Return a flag setting if the slit-edge traces should be used in the
        image processing.  The slits are required if a spatial flexure
        correction is requested and/or when the slit-illumination profile is
        removed.
        """
        if self.par is None:
            return False
        return self.par['spat_flexure_correct'] or (self.use_flat and self.par['use_illumflat'])

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

        # Have the images been trimmed?
        not_trimmed = self.rawimage.shape is not None and self.image.shape == self.rawimage.shape

        # Construct an image with the gain in each pixel
        gain = [None]*self.nimg
        for i in range(self.nimg):
            # Behavior needs to be different if the image has been trimmed or not
            gain[i] = procimg.gain_frame(self.datasec_img[i],
                                         np.atleast_1d(self.detector[i]['gain']))
            if not_trimmed:
                # Image is raw, so need to include overscan sections
                # TODO: This only works assuming that there is *no* overlap between
                # oscansec_img and datasec_img.  There shouldn't be, but the code
                # doesn't check this...
                gain[i] += procimg.gain_frame(self.oscansec_img[i],
                                              np.atleast_1d(self.detector[i]['gain']))
            # Set gain to 1 outside of the datasec and oscansec sections.
            gain[i][gain[i]==0] = 1
        # Convert from DN to counts
        self.image *= np.array(gain)

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
        if self.dark is None and self.par['shot_noise']:
            msgs.error('Dark image has not been created!  Run build_dark.')
        _dark = self.dark if self.par['shot_noise'] else None
        _counts = self.image if self.par['shot_noise'] else None
        # NOTE: self.dark is expected to be in *counts*.  This means that
        # procimg.base_variance should be called with exptime=None.  If the
        # exposure time is provided, the units of the dark current are expected
        # to be in e-/hr!
        self.base_var = procimg.base_variance(self.rn2img, darkcurr=_dark, #exptime=self.exptime,
                                              proc_var=self.proc_var, count_scale=self.img_scale)
        var = procimg.variance_model(self.base_var, counts=_counts, count_scale=self.img_scale,
                                     noise_floor=self.par['noise_floor'])
        return utils.inverse(var)

    def estimate_readnoise(self):
        """ Estimate the readnoise (in electrons) based on the overscan regions of
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
        for i in range(self.nimg):
            for amp in range(len(self.ronoise[i])):
                if self.ronoise[i,amp] > 0 and not self.par['empirical_rn']:
                    # Skip if the readnoise is defined and the empirical readnoise
                    # estimate was not explicitly requested.
                    continue
                if not np.any(self.oscansec_img[i]==amp+1):
                    msgs.error(f'Cannot estimate readnoise for amplifier {amp+1}.  Raw image '
                               'has no overscan region!')
                gain = 1. if self.steps['apply_gain'] else self.detector[i]['gain'][amp]
                biaspix = self.image[i,self.oscansec_img[i]==amp+1] * gain
                self.ronoise[i,amp] = stats.sigma_clipped_stats(biaspix, sigma=5)[-1]
                msgs.info(f'Estimated readnoise of amplifier {amp+1} = '
                          f'{self.ronoise[i,amp]:.3f} e-')

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

        # Have the images been trimmed?
        not_trimmed = self.rawimage.shape is not None and self.image.shape == self.rawimage.shape

        rn2 = [None]*self.nimg
        for i in range(self.nimg):
            # Compute and return the readnoise variance image 
            rn2[i] = procimg.rn2_frame(self.datasec_img[i], self.ronoise[i], units=units,
                                       gain=self.detector[i]['gain'], digitization=digitization)
            # TODO: Could also check if steps['trim'] is true.  Is either better or worse?
            if not_trimmed and np.any(self.oscansec_img[i] > 0):
                # Image is raw, so need to include overscan sections
                # TODO: This only works assuming that there is *no* overlap between
                # oscansec_img and datasec_img.  There shouldn't be, but the code
                # doesn't check this...
                rn2[i] += procimg.rn2_frame(self.oscansec_img[i], self.ronoise[i], units=units,
                                            gain=self.detector[i]['gain'],
                                            digitization=digitization)
        return np.array(rn2)

    def process(self, par, bpm=None, flatimages=None, bias=None, slits=None, dark=None,
                mosaic=False, debug=False):
        """
        Process the data.

        See further discussion of :ref:`image_proc` in ``PypeIt``.

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

            #. :func:`subtract_bias`: Subtract the processed bias image.  The
               shape and orientation of the bias image must match the
               *processed* image.  I.e., if you trim and orient this image, you
               must also have trimmed and oriented the bias frames.

            #. :func:`build_dark`: Create dark-current images using both the
               tabulated dark-current value for each detector and any directly
               observed dark images.  The shape and orientation of the observed
               dark image must match the *processed* image.  I.e., if you trim
               and orient this image, you must also have trimmed and oriented
               the dark frames.  To scale the dark image by the ratio of the
               exposure times to ensure the counts/s in the dark are removed
               from the image being processed, set the ``dark_expscale``
               parameter to true.

            #. :func:`subtract_dark`: Subtract the processed dark image and
               propagate any error.

            #. :func:`build_mosaic`: If data from multiple detectors are being
               processed as components of a detector mosaic, this resamples the
               individual images into a single image mosaic.  The current
               "resampling" scheme is restricted to nearest grid-point
               interpolation; see .  The placement of this step is important in that
               all of the previous corrections (overscan, trim, orientation,
               bias- and dark-subtraction) are done on the individual detector
               images.  However, after this point, we potentially need the slits
               and flat-field images which are *only defined in the mosaic
               frame*.  Because of this, bias and dark frames should *never* be
               reformatted into a mosaic.

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
            dark (:class:`~pypeit.images.pypeitimage.PypeItImage`):
                Dark image
            mosaic (:obj:`bool`, optional):
                When processing multiple detectors, resample the images into a
                mosaic.  If flats or slits are provided (and used), this *must*
                be true because these objects are always defined in the mosaic
                frame.
            debug (:obj:`bool`, optional):
                Run in debug mode.

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
            if self._bpm.ndim == 2:
                self._bpm = np.expand_dims(self._bpm, 0)

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
        if self.use_slits and slits is None:
            # TODO: I think this should only happen as a developer error, not a
            # user error, but I'm not sure.
            msgs.error('Processing steps requested that require slit-edge traces, but they were '
                       'not provided!')
        if self.nimg > 1 and not mosaic and (self.use_flat or self.use_slits):
            msgs.error('Mosaicing must be performed if multiple detectors are processed and '
                       'either flat-fielding or spatial flexure corrections are applied.')
        if self.nimg == 1 and mosaic:
            msgs.warn('Only processing a single detector; mosaicing is ignored.')

        msgs.info(f'Performing basic image processing on {os.path.basename(self.filename)}.')
        # TODO: Checking for bit saturation should be done here.

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

            # TODO: The logic of whether or not the BPM uses msbias to identify
            # bad pixels is difficult to follow.  Where and how the bpm is
            # created, and whether or not it uses msbias should be more clear.

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
            
        #   - Subtract processed bias
        if self.par['use_biasimage']:
            # Bias frame.  Shape and orientation must match *processed* image,.
            # Uncertainty from the bias subtraction is added to the variance.
            self.subtract_bias(bias)

        # TODO: Checking for count (well-depth) saturation should be done here.

        #   - Create the dark current image(s).  The dark-current image *always*
        #     includes the tabulated dark current and the call below ensures
        #     that this is *always* subtracted from the image being processed,
        #     regardless of whether or not use_darkimage is true.  If a dark
        #     frame is provided and subtracted, its shape and orientation must
        #     match the *processed* image, and the units *must* be in
        #     electrons/counts.
        self.build_dark(dark_image=dark if self.par['use_darkimage'] else None,
                        expscale=self.par['dark_expscale'])

        #   - Subtract dark current.  This simply subtracts the dark current
        #     from the image being processed.  If available, uncertainty from
        #     the dark subtraction is added to the processing variance.
        self.subtract_dark()

        # Create the mosaic images.  This *must* come after trimming and
        # orienting and before determining the spatial flexure shift and
        # applying any flat-fielding.
        if self.nimg > 1 and mosaic:
            self.build_mosaic()

        # Calculate flexure, if slits are provided and the correction is
        # requested.  NOTE: This step must come after trim, orient (just like
        # bias and dark subtraction) and before field flattening.  Also the
        # function checks that the slits exist if running the spatial flexure
        # correction, so no need to do it again here.
        self.spat_flexure_shift = self.spatial_flexure_shift(slits) \
                                    if self.par['spat_flexure_correct'] else None

        # Flat-field the data.  This propagates the flat-fielding corrections to
        # the variance.  The returned bpm is propagated to the PypeItImage
        # bitmask below.
        flat_bpm = self.flatfield(flatimages, slits=slits, debug=debug) if self.use_flat else None

        # Calculate the inverse variance
        self.ivar = self.build_ivar()

        #   - Subtract continuum level
        if self.par['subtract_continuum']:
            # Calculate a simple smooth continuum image, and subtract this from the frame
            self.subtract_continuum()

        # Generate a PypeItImage.
        # NOTE: To reconstruct the variance model, you need base_var, image,
        # img_scale, noise_floor, and shot_noise.
        _det, _image, _ivar, _datasec_img, _det_img, _rn2img, _base_var, _img_scale, _bpm \
                = self._squeeze()
        # NOTE: BPM MUST BE A BOOLEAN!
        pypeitImage = pypeitimage.PypeItImage(_image, ivar=_ivar, amp_img=_datasec_img,
                                              det_img=_det_img, rn2img=_rn2img, base_var=_base_var,
                                              img_scale=_img_scale, detector=_det,
                                              spat_flexure=self.spat_flexure_shift,
                                              PYP_SPEC=self.spectrograph.name,
                                              units='e-' if self.par['apply_gain'] else 'ADU',
                                              exptime=self.exptime,
                                              noise_floor=self.par['noise_floor'],
                                              shot_noise=self.par['shot_noise'],
                                              bpm=_bpm.astype(bool))

        pypeitImage.rawheadlist = self.headarr
        pypeitImage.process_steps = [key for key in self.steps.keys() if self.steps[key]]

        # Mask(s)
        if self.par['mask_cr']:
            # TODO: CR rejection of the darks was failing for HIRES for some reason...
            pypeitImage.build_crmask(self.par)

        pypeitImage.build_mask(saturation='default', mincounts='default')
        if flat_bpm is not None:
            pypeitImage.update_mask('BADSCALE', indx=flat_bpm)

        # Return
        return pypeitImage

    def _squeeze(self):
        """
        Convenience method for preparing attributes for construction of a
        :class:`~pypeit.images.pypeitimage.PypeItImage`.
        
        The issue is that :class:`~pypeit.images.rawimage.RawImage` image arrays
        are *always* 3D, even if there's only one image.  This is acceptable
        because use of :class:`~pypeit.images.rawimage.RawImage` is relatively
        self-contained.  It's really a namespace used for the image processing
        that disappears as soon as the image processing is done.

        :class:`~pypeit.images.pypeitimage.PypeItImage`, on the other hand, is a
        core class that is shared by many subclasses and used throughout the
        code base, meaning that it doesn't make sense to keep single images in
        3D arrays.

        This method "squeezes" (see `numpy.squeeze`_) the arrays used to
        construct a :class:`~pypeit.images.pypeitimage.PypeItImage` so that they
        are 3D only if they have to be.

        Returns:
            :obj:`tuple`: Returns the
            :class:`pypeit.images.detector_container.DetectorContainer` or
            :class:`pypeit.images.mosaic.Mosaic` instance, and the reshaped
            arrays with the image flux, inverse variance, amplifier number,
            detector number, readnoise-squared image, base-level variance, image
            scaling factor, and bad-pixel mask.
        """
        _det = self.detector[0] if self.mosaic is None else self.mosaic
        if self.nimg == 1:
            return _det, self.image[0], \
                   None if self.ivar is None else self.ivar[0], \
                   None if self.datasec_img is None else self.datasec_img[0], \
                   None if self.det_img is None else self.det_img[0], \
                   None if self.rn2img is None else self.rn2img[0], \
                   None if self.base_var is None else self.base_var[0], \
                   None if self.img_scale is None else self.img_scale[0], \
                   None if self.bpm is None else self.bpm[0]
        return _det, self.image, self.ivar, self.datasec_img, self.det_img, self.rn2img, \
                self.base_var, self.img_scale, self.bpm

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

        Return:
            float: The calculated flexure correction

        """
        step = inspect.stack()[0][3]
        if self.steps[step] and not force:
            # Already field flattened
            msgs.warn('Spatial flexure shift already calculated.')
            return
        if self.nimg > 1:
            msgs.error('CODING ERROR: Must use a single image (single detector or detector '
                       'mosaic) to determine spatial flexure.')
        self.spat_flexure_shift = flexure.spat_flexure_shift(self.image[0], slits)
        self.steps[step] = True
        # Return
        return self.spat_flexure_shift

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
        if self.nimg > 1:
            msgs.error('CODING ERROR: Can only apply flat field to a single image (single '
                       'detector or detector mosaic).')

        # Generate the illumination flat, as needed
        illum_flat = 1.0
        if self.par['use_illumflat']:
            # TODO :: We don't have tilts here yet... might be ever so slightly better, especially on very tilted slits
            illum_flat = flatimages.fit2illumflat(slits, spat_flexure=self.spat_flexure_shift, finecorr=False)
            illum_flat *= flatimages.fit2illumflat(slits, spat_flexure=self.spat_flexure_shift, finecorr=True)
            if debug:
                left, right = slits.select_edges(flexure=self.spat_flexure_shift)
                viewer, ch = display.show_image(illum_flat, chname='illum_flat')
                display.show_slits(viewer, ch, left, right)  # , slits.id)
                #
                viewer, ch = display.show_image(self.image[0], chname='orig_image')
                display.show_slits(viewer, ch, left, right)  # , slits.id)

        # Apply the relative spectral illumination
        spec_illum = flatimages.pixelflat_spec_illum if self.par['use_specillum'] else 1.

        # Apply flat-field correction
        # NOTE: Using flat.flatfield to effectively multiply image*img_scale is
        # a bit overkill...
        total_flat = flatimages.pixelflat_norm * illum_flat * spec_illum
        self.img_scale = np.expand_dims(utils.inverse(total_flat), 0)
        self.image[0], flat_bpm = flat.flatfield(self.image[0], total_flat)
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
        self.image = np.array([self.spectrograph.orient_image(d, i) 
                               for d,i in zip(self.detector, self.image)])
        if self.rn2img is not None:
            self.rn2img = np.array([self.spectrograph.orient_image(d, i)
                                    for d,i in zip(self.detector, self.rn2img)])
        if self.proc_var is not None:
            self.proc_var = np.array([self.spectrograph.orient_image(d, i)
                                      for d,i in zip(self.detector, self.proc_var)])
        self.datasec_img = np.array([self.spectrograph.orient_image(d, i)
                                     for d,i in zip(self.detector, self.datasec_img)])
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
        _bias = bias_image.image if self.nimg > 1 else np.expand_dims(bias_image.image, 0)
        if self.image.shape != _bias.shape:
            msgs.error('Shape mismatch with bias image!')
        self.image -= _bias
        # TODO: Also incorporate the mask?
        if bias_image.ivar is not None and self.proc_var is not None:
            self.proc_var += utils.inverse(bias_image.ivar if self.nimg > 1 
                                           else np.expand_dims(bias_image.ivar, 0))
        self.steps[step] = True

    # TODO: expscale is currently not a parameter that the user can control.
    # Should it be?
    def build_dark(self, dark_image=None, expscale=False):
        """
        Build the dark image data used for dark subtraction and error
        propagation.

        If a dark image is not provided, the dark image is simply the tabulated
        value and the error is set to None.  Otherwise, the dark is the
        combination of the tabulated dark-current for the detector and a dark
        image.  For this to be appropriate, the dark image (if provided) *must*
        also have had the tabulated dark-current value subtracted from it.

        Also, the processing of the dark image (if provided) should match the
        processing of the image being processed.  For example, if this image has
        been bias subtracted, so should be the dark image.

        If the ``dark_image`` object includes an inverse variance estimate, this
        is used to set the dark-current error.

        .. warning::

            Typically dark frames should have the same exposure time as the
            image being processed.  However, beware if that's not the case, and 
            make sure any use of exposure time scaling of the counts (see
            ``expscale``) is appropriate!

        Args:
            dark_image (:class:`~pypeit.images.pypeitimage.PypeItImage`, optional):
                The *observed* dark image in counts (not counts/s).  If None,
                only the tabulated dark-current are used to construct the dark
                image(s).
            expscale (:obj:`bool`, optional):
                Scale the dark image (if provided) by the ratio of the exposure
                times so that the counts per second represented by the dark
                image are correct.
        """
        # TODO: Is the dark-current amplifier dependent?  Also, this usage means
        # that darkcurr cannot be None.

        # Tabulated dark current is in e-/pixel/hour and exptime is in s, the
        # 3600 factor converts the dark current to e-/pixel, and the
        # multiplication by the number of binned pixels converts to e- (per
        # binned pixel).  The dark counts per binned pixel is determined
        # separately for each detector.
        self.dark = np.array([np.prod(parse.parse_binning(d.binning))*d.darkcurr
                                * self.exptime / 3600. for d in self.detector])

        if dark_image is None:
            # Reformat the tabulated dark counts into an image for each
            # detector.
            self.dark = np.repeat(self.dark, np.prod(self.image.shape[1:])
                                  ).reshape(self.image.shape)
            return

        # Include the deviations from the tabulated value as determined by an
        # observed dark image.
        _dark = dark_image.image if self.nimg > 1 else np.expand_dims(dark_image.image, 0)
        if self.image.shape != _dark.shape:
            # Shapes must match
            msgs.error(f'Dark image shape mismatch; expected {self.image.shape}, '
                       f'found {_dark.shape}.')

        # Scale the observed dark counts by the ratio of the exposure times.
        # TODO: Include a warning when the scaling is "significant"?
        scale = self.exptime / dark_image.exptime if expscale else 1.

        # Warn the user if the counts in the dark image are significantly
        # different from the tabulated value.  The 50% difference is not
        # well justified but, for now, hard-coded.
        med_dark = np.median(scale * _dark, axis=(1,2))
        large_signal = np.absolute(med_dark) > 0.5*self.dark
        if any(large_signal):
            med_str = np.array2string(med_dark, formatter={'float_kind':lambda x: "%.2f" % x},
                                      separator=',')
            drk_str = np.array2string(0.5*self.dark, formatter={'float_kind':lambda x: "%.2f" % x},
                                      separator=',')
            msgs.warn(f'Dark-subtracted dark frame has significant signal remaining.  Median '
                        f'counts are {med_str}; warning threshold = +/- {drk_str}.')

        # Combine the tabulated and observed dark values
        # NOTE: This converts self.dark from a vector to an array
        self.dark = scale * _dark + self.dark[:,None,None]
        if dark_image.ivar is not None:
            _dark_ivar = dark_image.ivar if self.nimg > 1 else np.expand_dims(dark_image.ivar, 0)
            # Include the scaling in the error, if available
            self.dark_var = scale**2 * utils.inverse(_dark_ivar)

    def subtract_dark(self, force=False):
        """
        Subtract detector dark current.

        The :attr:`dark` and :attr:`dark_ivar` arrays must have already been
        constructed using :func:`build_dark`.  If they aren't, a warning is
        thrown and nothing is done.

        Args:
            force (:obj:`bool`, optional):
                Force the dark to be subtracted, even if the step log
                (:attr:`steps`) indicates that it already has been.
        """
        step = inspect.stack()[0][3]
        if self.steps[step] and not force:
            # Already bias subtracted
            msgs.warn('Image was already dark subtracted.')
            return

        if self.dark is None:
            msgs.error('Dark image has not been created!  Run build_dark.')

        self.image -= self.dark
        if self.dark_var is not None:
            # Include the variance in the dark image
            self.proc_var += self.dark_var
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

        # NOTE: procimg.subtract_overscan checks that the provided images all
        # have the same shape.  This will fault if the images have already been
        # trimmed.
        _os_img = [None]*self.nimg
        var = [None]*self.nimg
        for i in range(self.nimg):
            # Subtract the overscan.  var is the variance in the overscan
            # subtraction.
            _os_img[i], var[i] = procimg.subtract_overscan(self.image[i], self.datasec_img[i],
                                                           self.oscansec_img[i],
                                                           method=self.par['overscan_method'],
                                                           params=self.par['overscan_par'],
                                                           var=None if self.rn2img is None 
                                                                else self.rn2img[i])
        self.image = np.array(_os_img)
        # Parse the returned value
        if self.rn2img is not None:
            # Include the variance in the overscan subtraction in the
            # "processing variance"
            self.proc_var += np.array(var)
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

        # The image cannot have already been trimmed
        if self.oscansec_img.shape != self.image.shape:
            msgs.error('Must estimate readnoise before trimming the image.')

        # Calculate the slit image
        _ps_img = [None]*self.nimg
        for i in range(self.nimg):
            # The image must have an overscan region for this to work.
            if not np.any(self.oscansec_img[i] > 0):
                msgs.error('Image has no overscan region.  Pattern noise cannot be subtracted.')

            patt_freqs = self.spectrograph.calc_pattern_freq(self.image[i], self.datasec_img[i],
                                                             self.oscansec_img[i], self.hdu)
            # Final check to make sure the list isn't empty (which it shouldn't be, anyway)
            if len(patt_freqs) == 0:
                patt_freqs = None
            # Subtract the pattern and overwrite the current image
            _ps_img[i] = procimg.subtract_pattern(self.image[i], self.datasec_img[i],
                                                  self.oscansec_img[i], frequency=patt_freqs)
        self.image = np.array(_ps_img)
        self.steps[step] = True

    def subtract_continuum(self, force=False):
        """
        Subtract the continuum level from the image.

        Args:
            force (:obj:`bool`, optional):
                Force the continuum to be subtracted, even if the step log
                (:attr:`steps`) indicates that it already has been.
        """
        step = inspect.stack()[0][3]
        if self.steps[step] and not force:
            # Already bias subtracted
            msgs.warn('Image was already continuum subtracted.')
            return

        # Generate the continuum image
        for ii in range(self.nimg):
            cont = np.zeros((self.image.shape[1], self.image.shape[2]))
            for rr in range(self.image.shape[2]):
                cont_now, cont_mask = arc.iter_continuum(self.image[ii, :, rr])
                cont[:,rr] = cont_now
            self.image[ii,:,:] -= cont
        #cont = ndimage.median_filter(self.image, size=(1,101,3), mode='reflect')
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
            msgs.warn('Image shape does not match raw image.  Assuming it was already trimmed.')
            return
        if self.steps[step] and not force:
            # Already trimmed
            msgs.warn('Image was already trimmed.')
            return
        self.image = np.array([procimg.trim_frame(i, d < 1)
                               for i, d in zip(self.image, self.datasec_img)])
        if self.rn2img is not None:
            self.rn2img = np.array([procimg.trim_frame(r, d < 1)
                                    for r, d in zip(self.rn2img, self.datasec_img)])
        if self.proc_var is not None:
            self.proc_var = np.array([procimg.trim_frame(p, d < 1)
                                      for p, d in zip(self.proc_var, self.datasec_img)])
        # NOTE: This must be done last because the untrimmed image is used for
        # deciding what to trim!
        self.datasec_img = np.array([procimg.trim_frame(d, d < 1) for d in self.datasec_img])
        self.steps[step] = True

    def build_mosaic(self):
        """
        When processing multiple detectors, this remaps the detector data to a mosaic.

        This is largely a wrapper for multiple calls to
        :func:`~pypeit.core.mosaic.build_image_mosaic`.  Resampling is currently
        restricted to nearest grid-point interpolation (``order=0``).

        Construction of the mosaic(s) must be done *after* the images have been
        trimmed and oriented to follow the ``PypeIt`` convention.

        This function remaps ``image``, ``datasec_img``, ``rn2img``, ``dark``,
        ``dark_var``, ``proc_var``, and ``base_var``.  These are all originally
        calculated in the native detector frame.  Because ``img_scale`` is only
        related to the flat-field images, it is not remapped because these
        images are always processed in the mosaic frame.
        """
        if self.nimg == 1:
            # NOTE: This also catches cases where the mosaicing has already been
            # performed.
            msgs.warn('There is only one image, so there is nothing to mosaic!')
            return

        # Check that the mosaicing is allowed
        if not self.steps['trim'] or not self.steps['orient']:
            msgs.error('Images must be trimmed and PypeIt-oriented before mosaicing.')

        # Create images that will track which detector contributes to each pixel
        # in the mosaic.  These images are created here first *before*
        # `self.image` is mosaiced below. NOTE: This assumes there is no overlap
        # in the detector mosaic (which should be true).
        self.det_img = np.array([np.full(img.shape, d.det, dtype=int)
                                    for img,d in zip(self.image, self.detector)])

        # Transform the image data to the mosaic frame.  This call determines
        # the shape of the mosaic image and adjusts the relative transforms to
        # the absolute mosaic frame.
        self.image, _, _img_npix, _tforms = build_image_mosaic(self.image, self.mosaic.tform, order=self.mosaic.msc_order)
        shape = self.image.shape
        # Maintain dimensionality
        self.image = np.expand_dims(self.image, 0)

        # For the remaining transforms, we can use the first call to skip over
        # the need to determine the output mosaic shape and the adjusted
        # transforms.

        # Transform the BPM and maintain its type
        bpm_type = self.bpm.dtype
        self._bpm = build_image_mosaic(self.bpm.astype(float), _tforms, mosaic_shape=shape, order=self.mosaic.msc_order)[0]
        # Include pixels that have no contribution from the original image in
        # the bad pixel mask of the mosaic.
        self._bpm[_img_npix < 1] = 1
        # np.round helps to deal with cases where the interpolation is performed
        # and values of adjacent pixels are combined
        self._bpm = np.expand_dims(np.round(self._bpm).astype(bpm_type), 0)

        # NOTE: The bitmask is set by a combination of pixels without any
        # contributions when creating the image mosaic and when mosaicing the
        # detector bpms.  The `_img_npix` equivalent for the other mosaics
        # (e.g., the rn2img mosaic) should be identical.

        # Get the pixels associated with each amplifier
        self.datasec_img = build_image_mosaic(self.datasec_img.astype(float), _tforms,
                                              mosaic_shape=shape, order=self.mosaic.msc_order)[0]
        self.datasec_img = np.expand_dims(np.round(self.datasec_img).astype(int), 0)

        # Get the pixels associated with each detector
        self.det_img = build_image_mosaic(self.det_img.astype(float), _tforms,
                                          mosaic_shape=shape, order=self.mosaic.msc_order)[0]
        self.det_img = np.expand_dims(np.round(self.det_img).astype(int), 0)

        # Transform all the variance arrays, as necessary
        if self.rn2img is not None:
            self.rn2img = build_image_mosaic(self.rn2img, _tforms, mosaic_shape=shape, order=self.mosaic.msc_order)[0]
            self.rn2img = np.expand_dims(self.rn2img, 0)
        if self.dark is not None:
            self.dark = build_image_mosaic(self.dark, _tforms, mosaic_shape=shape, order=self.mosaic.msc_order)[0]
            self.dark = np.expand_dims(self.dark, 0)
        if self.dark_var is not None:
            self.dark_var = build_image_mosaic(self.dark_var, _tforms, mosaic_shape=shape, order=self.mosaic.msc_order)[0]
            self.dark_var = np.expand_dims(self.dark_var, 0)
        if self.proc_var is not None:
            self.proc_var = build_image_mosaic(self.proc_var, _tforms, mosaic_shape=shape, order=self.mosaic.msc_order)[0]
            self.proc_var = np.expand_dims(self.proc_var, 0)
        if self.base_var is not None:
            self.base_var = build_image_mosaic(self.base_var, _tforms, mosaic_shape=shape, order=self.mosaic.msc_order)[0]
            self.base_var = np.expand_dims(self.base_var, 0)

        # TODO: Mosaicing means that many of the internals are no longer
        # valid/useful.  Specifically, I set ronoise, rawimage, rawdatsec_img,
        # and oscansec_img to None, which will force the code to barf if some of
        # the methods that use these are called after the mosaicing.
        self.ronoise = None
        self.rawimage = None
        self.rawdatasec_img = None
        self.oscansec_img = None
        self.nimg = 1
        self.steps[inspect.stack()[0][3]] = True

    def __repr__(self):
        return f'<{self.__class__.__name__}: file={self.filename}, nimg={self.nimg}, ' \
               f'steps={self.steps}>'


