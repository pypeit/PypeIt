""" Class to generate an image from one or more files (and other pieces).

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import os

from IPython import embed

import numpy as np

from pypeit import msgs
from pypeit.core import combine
from pypeit.core import procimg
from pypeit.par import pypeitpar
from pypeit import utils
from pypeit.images import pypeitimage
from pypeit.images import rawimage
from pypeit.images import imagebitmask

class CombineImage:
    """
    Process and combine detector images.

    All core processing steps for each image are handled by
    :class:`~pypeit.images.rawimage.RawImage`.  This class can be used to
    process both single images, lists of images, and detector mosaics.

    Args:
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            Spectrograph used to take the data.
        det (:obj:`int`, :obj:`tuple`):
            The 1-indexed detector number(s) to process.  If a tuple, it must
            include detectors viable as a mosaic for the provided spectrograph;
            see :func:`~pypeit.spectrographs.spectrograph.Spectrograph.allowed_mosaics`.
        par (:class:`~pypeit.par.pypeitpar.ProcessImagesPar`):
            Parameters that dictate the processing of the images.
        files (:obj:`str`, array-like):
            A set of one or more images to process and combine.

    Attributes:
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            Spectrograph used to take the data.
        det (:obj:`int`, :obj:`tuple`):
            The 1-indexed detector number(s) to process.
        par (:class:`~pypeit.par.pypeitpar.ProcessImagesPar`):
            Parameters that dictate the processing of the images.
        files (:obj:`list`):
            A set of one or more images to process and combine.

    """
    def __init__(self, spectrograph, det, par, files):
        self.spectrograph = spectrograph
        self.det = det
        if not isinstance(par, pypeitpar.ProcessImagesPar):
            msgs.error('Provided ParSet for must be type ProcessImagesPar.')
        self.par = par  # This musts be named this way as it is frequently a child
        self.files = list(files) if hasattr(files, '__len__') else [files]
        # NOTE: nfiles is a property method.  Defining files above must come
        # before this check!
        if self.nfiles == 0:
            msgs.error('CombineImage requires a list of files to instantiate')

    def run(self, bias=None, flatimages=None, ignore_saturation=False, sigma_clip=True,
            bpm=None, sigrej=None, maxiters=5, slits=None, dark=None, combine_method='mean',
            mosaic=False):
        r"""
        Process and combine all images.

        All processing is performed by the
        :class:`~pypeit.images.rawimage.RawImage` class; see 
        :func:`~pypeit.images.rawimage.RawImage.process`.

        If there is only one file (see :attr:`files`), this simply processes the
        file and returns the result.
        
        If there are multiple files, all the files are processed and the
        processed images are combined based on the ``combine_method``, where the
        options are:

            - 'mean': If ``sigma_clip`` is True, this is a sigma-clipped mean;
              otherwise, this is a simple average.  The combination is done
              using :func:`~pypeit.core.combine.weighted_combine`.

            - 'median': This is a simple masked median (using
              `numpy.ma.median`_).

        The errors in the image are also propagated through the stacking
        procedure; however, this isn't a simple propagation of the inverse
        variance arrays.  The image processing produces arrays with individual
        components used to construct the variance model for an individual frame.
        See :ref:`image_proc` and :func:`~pypeit.procimg.variance_model` for a
        description of these arrays.  Briefly, the relevant arrays are the
        readnoise variance (:math:`V_{\rm rn}`), the "processing" variance
        (:math:`V_{\rm proc}`), and the image scaling (i.e., the flat-field
        correction) (:math:`s`).  The variance calculation for the stacked image
        directly propagates the error in these.  For example, the propagated
        processing variance (modulo the masking) is:

        .. math::

            V_{\rm proc,stack} = \frac{\sum_i s_i^2 V_{{\rm
            proc},i}}\frac{s_{\rm stack}^2}

        where :math:`s_{\rm stack}` is the combined image scaling array,
        combined in the same way as the image data are combined.  This ensures
        that the reconstruction of the uncertainty in the combined image
        calculated using :func:`~pypeit.procimg.variance_model` accurately
        includes, e.g., the processing uncertainty.

        The uncertainty in the combined image, however, recalculates the
        variance model, using the combined image (which should have less noise)
        to set the Poisson statistics.  The same parameters used when processing
        the individual frames are applied to the combined frame; see
        :func:`~pypeit.images.rawimage.RawImage.build_ivar`.  This calculation
        is then the equivalent of when the observed counts are replaced by the
        model object and sky counts during sky subtraction and spectral
        extraction.

        Bitmasks from individual frames in the stack are *not* propagated to the
        combined image, except to indicate when a pixel was masked for all
        images in the stack (cf., ``ignore_saturation``).  Additionally, the
        instrument-specific bad-pixel mask, see the
        :func:`~pypeit.spectrographs.spectrograph.Spectrograph.bpm` method for
        each instrument subclass, saturated-pixel mask, and other default mask
        bits (e.g., NaN and non-positive inverse variance values) are all
        propagated to the combined-image mask; see
        :func:`~pypeit.images.pypeitimage.PypeItImage.build_mask`.
        
        .. warning::

            All image processing of the data in :attr:`files` *must* result
            in images of the same shape.

        Args:
            bias (:class:`~pypeit.images.buildimage.BiasImage`, optional):
                Bias image for bias subtraction; passed directly to
                :func:`~pypeit.images.rawimage.RawImage.process` for all images.
            flatimages (:class:`~pypeit.flatfield.FlatImages`, optional):
                Flat-field images for flat fielding; passed directly to
                :func:`~pypeit.images.rawimage.RawImage.process` for all images.
            ignore_saturation (:obj:`bool`, optional):
                If True, turn off the saturation flag in the individual images
                before stacking.  This avoids having such values set to 0, which
                for certain images (e.g. flat calibrations) can have unintended
                consequences.
            sigma_clip (:obj:`bool`, optional):
                When ``combine_method='mean'``, perform a sigma-clip the data;
                see :func:`~pypeit.core.combine.weighted_combine`.
            bpm (`numpy.ndarray`_, optional):
                Bad pixel mask; passed directly to
                :func:`~pypeit.images.rawimage.RawImage.process` for all images.
            sigrej (:obj:`float`, optional):
                When ``combine_method='mean'``, this sets the sigma-rejection
                thresholds used when sigma-clipping the image combination.
                Ignored if ``sigma_clip`` is False.  If None and ``sigma_clip``
                is True, the thresholds are determined automatically based on
                the number of images provided; see
                :func:`~pypeit.core.combine.weighted_combine``.
            maxiters (:obj:`int`, optional):
                When ``combine_method='mean'``) and sigma-clipping
                (``sigma_clip`` is True), this sets the maximum number of
                rejection iterations.  If None, rejection iterations continue
                until no more data are rejected; see
                :func:`~pypeit.core.combine.weighted_combine``.
            slits (:class:`~pypeit.slittrace.SlitTraceSet`, optional):
                Slit edge trace locations; passed directly to
                :func:`~pypeit.images.rawimage.RawImage.process` for all images.
            dark (:class:`~pypeit.images.buildimage.DarkImage`, optional):
                Dark-current image; passed directly to
                :func:`~pypeit.images.rawimage.RawImage.process` for all images.
            combine_method (:obj:`str`, optional):
                Method used to combine images.  Must be ``'mean'`` or
                ``'median'``; see above.
            mosaic (:obj:`bool`, optional):
                If multiple detectors are being processes, mosaic them into a
                single image.  See
                :func:`~pypeit.images.rawimage.RawImage.process`.

        Returns:
            :class:`~pypeit.images.pypeitimage.PypeItImage`: The combination of
            all the processed images.
        """
        # Check the input (i.e., bomb out *before* it does any processing)
        if self.nfiles == 0:
            msgs.error('Object contains no files to process!')
        if self.nfiles > 1 and combine_method not in ['mean', 'median']:
            msgs.error(f'Unknown image combination method, {combine_method}.  Must be '
                       '"mean" or "median".')

        # Loop on the files
        for kk, ifile in enumerate(self.files):
            # Load raw image
            rawImage = rawimage.RawImage(ifile, self.spectrograph, self.det)
            # Process
            pypeitImage = rawImage.process(self.par, bias=bias, bpm=bpm, dark=dark,
                                           flatimages=flatimages, slits=slits, mosaic=mosaic)

            if self.nfiles == 1:
                # Only 1 file, so we're done
                pypeitImage.files = self.files
                return pypeitImage
            elif kk == 0:
                # Allocate arrays to collect data for each frame
                shape = (self.nfiles,) + pypeitImage.shape
                img_stack = np.zeros(shape, dtype=float)
                scl_stack = np.ones(shape, dtype=float)
                rn2img_stack = np.zeros(shape, dtype=float)
                basev_stack = np.zeros(shape, dtype=float)
                gpm_stack = np.zeros(shape, dtype=bool)
                lampstat = [None]*self.nfiles
                exptime = np.zeros(self.nfiles, dtype=float)

            # Save the lamp status
            # TODO: As far as I can tell, this is the *only* place rawheadlist
            # is used.  Is there a way we can get this from fitstbl instead?
            lampstat[kk] = self.spectrograph.get_lamps_status(pypeitImage.rawheadlist)
            # Save the exposure time to check if it's consistent for all images.
            exptime[kk] = pypeitImage.exptime
            # Processed image
            img_stack[kk] = pypeitImage.image
            # Get the count scaling
            if pypeitImage.img_scale is not None:
                scl_stack[kk] = pypeitImage.img_scale
            # Read noise squared image
            if pypeitImage.rn2img is not None:
                rn2img_stack[kk] = pypeitImage.rn2img * scl_stack[kk]**2
            # Processing variance image
            if pypeitImage.base_var is not None:
                basev_stack[kk] = pypeitImage.base_var * scl_stack[kk]**2
            # Final mask for this image
            # TODO: This seems kludgy to me. Why not just pass ignore_saturation
            # to process_one and ignore the saturation when the mask is actually
            # built, rather than untoggling the bit here?
            if ignore_saturation:  # Important for calibrations as we don't want replacement by 0
                pypeitImage.update_mask('SATURATION', action='turn_off')
            # Get a simple boolean good-pixel mask for all the unmasked pixels
            gpm_stack[kk] = pypeitImage.select_flag(invert=True)

        # Check that the lamps being combined are all the same:
        if not lampstat[1:] == lampstat[:-1]:
            msgs.warn("The following files contain different lamp status")
            # Get the longest strings
            maxlen = max([len("Filename")]+[len(os.path.split(x)[1]) for x in self.files])
            maxlmp = max([len("Lamp status")]+[len(x) for x in lampstat])
            strout = "{0:" + str(maxlen) + "}  {1:s}"
            # Print the messages
            print(msgs.indent() + '-'*maxlen + "  " + '-'*maxlmp)
            print(msgs.indent() + strout.format("Filename", "Lamp status"))
            print(msgs.indent() + '-'*maxlen + "  " + '-'*maxlmp)
            for ff, file in enumerate(self.files):
                print(msgs.indent()
                      + strout.format(os.path.split(file)[1], " ".join(lampstat[ff].split("_"))))
            print(msgs.indent() + '-'*maxlen + "  " + '-'*maxlmp)

        # Do a similar check for exptime
        if np.any(np.absolute(np.diff(exptime)) > 0):
            # TODO: This should likely throw an error instead!
            msgs.warn('Exposure time is not consistent for all images being combined!  '
                      'Using the average.')
            comb_texp = np.mean(exptime)
        else:
            comb_texp = exptime[0]

        # Coadd them
        if combine_method == 'mean':
            weights = np.ones(self.nfiles, dtype=float)/self.nfiles
            img_list_out, var_list_out, gpm, nframes \
                    = combine.weighted_combine(weights,
                                               [img_stack, scl_stack],  # images to stack
                                               [rn2img_stack, basev_stack], # variances to stack
                                               gpm_stack, sigma_clip=sigma_clip,
                                               sigma_clip_stack=img_stack,  # clipping based on img
                                               sigrej=sigrej, maxiters=maxiters)
            comb_img, comb_scl = img_list_out
            comb_rn2, comb_basev = var_list_out
            # Divide by the number of images that contributed to each pixel
            comb_scl[gpm] /= nframes[gpm]

        elif combine_method == 'median':
            bpm_stack = np.logical_not(gpm_stack)
            nframes = np.sum(gpm_stack, axis=0)
            gpm = nframes > 0
            comb_img = np.ma.median(np.ma.MaskedArray(img_stack, mask=bpm_stack),axis=0).filled(0.)
            # TODO: I'm not sure if this is right.  Maybe we should just take
            # the masked average scale instead?
            comb_scl = np.ma.median(np.ma.MaskedArray(scl_stack, mask=bpm_stack),axis=0).filled(0.)
            # First calculate the error in the sum.  The variance is set to 0
            # for pixels masked in all images.
            comb_rn2 = np.ma.sum(np.ma.MaskedArray(rn2img_stack, mask=bpm_stack),axis=0).filled(0.)
            comb_basev = np.ma.sum(np.ma.MaskedArray(basev_stack, mask=bpm_stack),axis=0).filled(0.)
            # Convert to standard error in the median (pi/2 factor relates standard variance
            # in mean (sum(variance_i)/n^2) to standard variance in median)
            comb_rn2[gpm] *= np.pi/2/nframes[gpm]**2
            comb_basev[gpm] *= np.pi/2/nframes[gpm]**2
            # Divide by the number of images that contributed to each pixel
            comb_scl[gpm] *= np.pi/2/nframes[gpm]
        else:
            # NOTE: Given the check at the beginning of the function, the code
            # should *never* make it here.
            msgs.error("Bad choice for combine.  Allowed options are 'median', 'mean'.")

        # Recompute the inverse variance using the combined image
        comb_var = procimg.variance_model(comb_basev,
                                          counts=comb_img if self.par['shot_noise'] else None,
                                          count_scale=comb_scl,
                                          noise_floor=self.par['noise_floor'])

        # Build the combined image
        comb = pypeitimage.PypeItImage(image=comb_img, ivar=utils.inverse(comb_var), nimg=nframes,
                                       amp_img=pypeitImage.amp_img, det_img=pypeitImage.det_img,
                                       rn2img=comb_rn2, base_var=comb_basev, img_scale=comb_scl,
                                       # NOTE: This *must* be a boolean.
                                       bpm=np.logical_not(gpm), 
                                       # NOTE: The detector is needed here so
                                       # that we can get the dark current later.
                                       detector=pypeitImage.detector,
                                       PYP_SPEC=self.spectrograph.name,
                                       units='e-' if self.par['apply_gain'] else 'ADU',
                                       exptime=comb_texp, noise_floor=self.par['noise_floor'],
                                       shot_noise=self.par['shot_noise'])

        # Internals
        # TODO: Do we need these?
        comb.files = self.files
        comb.rawheadlist = pypeitImage.rawheadlist
        comb.process_steps = pypeitImage.process_steps

        # Build the base level mask
        comb.build_mask(saturation='default', mincounts='default')

        # Flag all pixels with no contributions from any of the stacked images.
        comb.update_mask('STCKMASK', indx=np.logical_not(gpm))

        # Return
        return comb

    @property
    def nfiles(self):
        """
        The number of files in :attr:`files`.
        """
        return len(self.files) if isinstance(self.files, (np.ndarray, list)) else 0


