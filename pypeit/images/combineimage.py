""" Class to generate an image from one or more files (and other pieces).

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from IPython import embed

import numpy as np

from pypeit import msgs
from pypeit.core import combine
from pypeit.core import procimg
from pypeit.par import pypeitpar
from pypeit import utils
from pypeit.images import pypeitimage
from pypeit.images import imagebitmask


class CombineImage:
    """
    Process and combine detector images. 

    Args:
        rawImages (:obj:`list`, :class:`~pypeit.images.pypeitimage.PypeItImage`):
            Either a single :class:`~pypeit.images.pypeitimage.PypeItImage`
            object or a list of one or more of these objects to be combined into
            an image.
        par (:class:`~pypeit.par.pypeitpar.ProcessImagesPar`):
            Parameters that dictate the processing of the images.

    Attributes:
        det (:obj:`int`, :obj:`tuple`):
            The 1-indexed detector number(s) to process.
        par (:class:`~pypeit.par.pypeitpar.ProcessImagesPar`):
            Parameters that dictate the processing of the images.
        rawImages (:obj:`list`):
            A list of one or more :class:`~pypeit.images.rawimage.RawImage` objects 
            to be combined.             
    """
    def __init__(self, rawImages, par):
        if not isinstance(par, pypeitpar.ProcessImagesPar):
            msgs.error('Provided ParSet for must be type ProcessImagesPar.')
        self.rawImages = list(rawImages) if hasattr(rawImages, '__len__') else [rawImages]
        self.par = par  # This musts be named this way as it is frequently a child

        # NOTE: nimgs is a property method.  Defining rawImages above must come
        # before this check!
        if self.nimgs == 0:
            msgs.error('CombineImage requires a list of files to instantiate')


    def run(self, ignore_saturation=False, maxiters=5):
        r"""
        Process and combine all images.

        All processing is performed by the
        :class:`~pypeit.images.rawimage.RawImage` class; see 
        :func:`~pypeit.images.rawimage.RawImage.process`.

        If there is only one file (see :attr:`files`), this simply processes the
        file and returns the result.
        
        If there are multiple files, all the files are processed and the
        processed images are combined based on the ``par['combine']``, where the
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
            ignore_saturation (:obj:`bool`, optional):
                If True, turn off the saturation flag in the individual images
                before stacking.  This avoids having such values set to 0, which
                for certain images (e.g. flat calibrations) can have unintended
                consequences.
            maxiters (:obj:`int`, optional):
                When ``par['combine']='mean'``) and sigma-clipping
                (``sigma_clip`` is True), this sets the maximum number of
                rejection iterations.  If None, rejection iterations continue
                until no more data are rejected; see
                :func:`~pypeit.core.combine.weighted_combine``.

        Returns:
            :class:`~pypeit.images.pypeitimage.PypeItImage`: The combination of
            all the processed images.
        """
        
        
        # Check the input (i.e., bomb out *before* it does any processing)
        if self.nimgs == 0:
            msgs.error('Object contains no files to process!')
        if self.nimgs > 1 and self.par['combine'] not in ['mean', 'median']:
            msgs.error(f'Unknown image combination method, {self.par["combine"]}.  Must be '
                       '"mean" or "median".')
        file_list = []
        # Loop on the files
        for kk, rawImage in enumerate(self.rawImages):
            if self.nimgs == 1:
                # Only 1 file, so we're done
                rawImage.files = [rawImage.filename]
                return rawImage
            elif kk == 0:
                # Allocate arrays to collect data for each frame
                shape = (self.nimgs,) + rawImage.shape
                img_stack = np.zeros(shape, dtype=float)
                scl_stack = np.ones(shape, dtype=float)
                rn2img_stack = np.zeros(shape, dtype=float)
                basev_stack = np.zeros(shape, dtype=float)
                gpm_stack = np.zeros(shape, dtype=bool)
                exptime = np.zeros(self.nimgs, dtype=float)

            # Save the exposure time to check if it's consistent for all images.
            exptime[kk] = rawImage.exptime
            # Processed image
            img_stack[kk] = rawImage.image
            # Get the count scaling
            if rawImage.img_scale is not None:
                scl_stack[kk] = rawImage.img_scale
            # Read noise squared image
            if rawImage.rn2img is not None:
                rn2img_stack[kk] = rawImage.rn2img * scl_stack[kk]**2
            # Processing variance image
            if rawImage.base_var is not None:
                basev_stack[kk] = rawImage.base_var * scl_stack[kk]**2
            # Final mask for this image
            # TODO: This seems kludgy to me. Why not just pass ignore_saturation
            # to process_one and ignore the saturation when the mask is actually
            # built, rather than untoggling the bit here?
            if ignore_saturation:  # Important for calibrations as we don't want replacement by 0
                rawImage.update_mask('SATURATION', action='turn_off')
            # Get a simple boolean good-pixel mask for all the unmasked pixels
            gpm_stack[kk] = rawImage.select_flag(invert=True)
            file_list.append(rawImage.filename)

        # Check that all exposure times are consistent
        # TODO: JFH suggests that we move this to calibrations.check_calibrations
        if np.any(np.absolute(np.diff(exptime)) > 0):
            # TODO: This should likely throw an error instead!
            msgs.warn('Exposure time is not consistent for all images being combined!  '
                      'Using the average.')
            comb_texp = np.mean(exptime)
        else:
            comb_texp = exptime[0]

        # scale the images to their mean, if requested, before combining
        if self.par['scale_to_mean']:
            msgs.info("Scaling images to have the same mean before combining")
            # calculate the mean of the images
            [mean_img], _, mean_gpm, _ = combine.weighted_combine(np.ones(self.nimgs, dtype=float)/self.nimgs,
                                                                  [img_stack],
                                                                  [rn2img_stack],
                                                                  # var_list is added because it is
                                                                  # required by the function but not used
                                                                  gpm_stack, sigma_clip=self.par['clip'],
                                                                  sigma_clip_stack=img_stack,
                                                                  sigrej=self.par['comb_sigrej'], maxiters=maxiters)

            # scale factor
            # TODO: Chose the median over the whole frame to avoid outliers.  Is this the right choice?
            _mscale = np.nanmedian(mean_img[None, mean_gpm]/img_stack[:, mean_gpm], axis=1)
            # reshape the scale factor
            mscale = _mscale[:, None, None]
            # scale the images
            img_stack *= mscale
            # scale the scales
            scl_stack *= mscale

            # scale the variances
            rn2img_stack *= mscale**2
            basev_stack *= mscale**2

        # Coadd them
        if self.par['combine'] == 'mean':
            weights = np.ones(self.nimgs, dtype=float)/self.nimgs
            img_list_out, var_list_out, gpm, nframes \
                    = combine.weighted_combine(weights,
                                               [img_stack, scl_stack],  # images to stack
                                               [rn2img_stack, basev_stack], # variances to stack
                                               gpm_stack, sigma_clip=self.par['clip'],
                                               sigma_clip_stack=img_stack,  # clipping based on img
                                               sigrej=self.par['comb_sigrej'], maxiters=maxiters)
            comb_img, comb_scl = img_list_out
            comb_rn2, comb_basev = var_list_out
            # Divide by the number of images that contributed to each pixel
            comb_scl[gpm] /= nframes[gpm]

        elif self.par['combine'] == 'median':
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
                                       amp_img=rawImage.amp_img, det_img=rawImage.det_img,
                                       rn2img=comb_rn2, base_var=comb_basev, img_scale=comb_scl,
                                       # NOTE: This *must* be a boolean.
                                       bpm=np.logical_not(gpm), 
                                       # NOTE: The detector is needed here so
                                       # that we can get the dark current later.
                                       detector=rawImage.detector,
                                       PYP_SPEC=rawImage.PYP_SPEC,
                                       units='e-' if self.par['apply_gain'] else 'ADU',
                                       exptime=comb_texp, noise_floor=self.par['noise_floor'],
                                       shot_noise=self.par['shot_noise'])

        # Internals
        # TODO: Do we need these?
        comb.files = file_list
        comb.rawheadlist = rawImage.rawheadlist
        comb.process_steps = rawImage.process_steps

        # Build the base level mask
        comb.build_mask(saturation='default' if not ignore_saturation else None, mincounts='default')

        # Flag all pixels with no contributions from any of the stacked images.
        comb.update_mask('STCKMASK', indx=np.logical_not(gpm))

        # Return
        return comb

    @property
    def nimgs(self):
        """
        The number of files in :attr:`files`.
        """
        return len(self.rawImages) if isinstance(self.rawImages, (np.ndarray, list)) else 0


