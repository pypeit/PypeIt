""" Class to generate an image from one or more files (and other pieces).

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import os

from IPython import embed

import numpy as np

from pypeit import msgs
from pypeit.core import combine
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
    process both single images and lists of images.

    Args:
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            Spectrograph used to take the data.
        det (:obj:`int`):
            The 1-indexed detector number to process.
        par (:class:`~pypeit.par.pypeitpar.ProcessImagesPar`):
            Parameters that dictate the processing of the images.
        files (array-like):
            A set of one or more images to process and combine.

    Attributes:
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            Spectrograph used to take the data.
        det (:obj:`int`):
            The 1-indexed detector number to process.
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

    # TODO: calling the combination method "weightmean" is confusing because it
    # doesn't do any weighting.  It's just a sigma-clipped mean.
    def run(self, bias=None, flatimages=None, ignore_saturation=False, sigma_clip=True,
            bpm=None, sigrej=None, maxiters=5, slits=None, dark=None, combine_method='weightmean'):
        """
        Process and combine all images.

        All processing is performed by the
        :class:`~pypeit.images.rawimage.RawImage` class, and images are combined
        based on the ``combine_method``, where the options are:

            - 'weightmean': If ``sigma_clip`` is True, this is a sigma-clipped
              mean; otherwise, this is a simple average.  The combination is
              done using :func:`~pypeit.core.combine.weighted_combine`.

            - 'median': This is a simple masked median (using
              `numpy.ma.median`_).

        Depending on the processing performed, this may also generate the
        inverse variance, cosmic-ray mask, read-noise variance image, and
        processing mask.

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
                When ``combine_method='weightmean'``, perform a sigma-clip the
                data; see :func:`~pypeit.core.combine.weighted_combine`.
            bpm (`numpy.ndarray`_, optional):
                Bad pixel mask; passed directly to
                :func:`~pypeit.images.rawimage.RawImage.process` for all images.
            sigrej (:obj:`float`, optional):
                When ``combine_method='weightmean'``, this sets the
                sigma-rejection thresholds used when sigma-clipping the image
                combination.  Ignored if ``sigma_clip`` is False.  If None and
                ``sigma_clip`` is True, the thresholds are determined
                automatically based on the number of images provided; see
                :func:`~pypeit.core.combine.weighted_combine``.
            maxiters (:obj:`int`, optional):
                When ``combine_method='weightmean'``) and sigma-clipping
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
            combine_method (str):
                Method used to combine images.  Must be ``'weightmean'`` or
                ``'median'``; see above.

        Returns:
            :class:`~pypeit.images.pypeitimage.PypeItImage`: The combination of
            all the processed images.
        """
        # Check the input (i.e., bomb out *before* it does any processing)
        if self.nfiles == 0:
            msgs.error('Object contains no files to process!')
        if self.nfiles > 1 and combine_method not in ['weightmean', 'median']:
            msgs.error(f'Unknown image combination method, {combine_method}.  Must be '
                       '"weightmean" or "median".')

        # Instantiate the bitmask
        bitmask = imagebitmask.ImageBitMask()

        # Loop on the files
        for kk, ifile in enumerate(self.files):
            # Load raw image
            rawImage = rawimage.RawImage(ifile, self.spectrograph, self.det)
            # Process
            pypeitImage = rawImage.process(self.par, bias=bias, bpm=bpm, dark=dark,
                                           flatimages=flatimages, slits=slits)
            if self.nfiles == 1:
                # Only 1 file, so we're done
                return pypeitImage
            elif kk == 0:
                # Allocate arrays to collect data for each frame
                shape = (self.nfiles,) + pypeitImage.shape
                img_stack = np.zeros(shape, dtype=float)
                ivar_stack= np.ones(shape, dtype=float)
                rn2img_stack = np.zeros(shape, dtype=float)
                crmask_stack = np.zeros(shape, dtype=bool)
                mask_stack = np.zeros(shape, bitmask.minimum_dtype(asuint=True))
                lampstat = [None]*self.nfiles

            # Save the lamp status
            lampstat[kk] = self.spectrograph.get_lamps_status(pypeitImage.rawheadlist)
            # Process
            img_stack[kk] = pypeitImage.image
            # Construct raw variance image and turn into inverse variance
            if pypeitImage.ivar is not None:
                ivar_stack[kk] = pypeitImage.ivar
            # Mask cosmic rays
            if pypeitImage.crmask is not None:
                crmask_stack[kk] = pypeitImage.crmask
            # Read noise squared image
            if pypeitImage.rn2img is not None:
                rn2img_stack[kk] = pypeitImage.rn2img
            # Final mask for this image
            # TODO: This seems kludgy to me. Why not just pass ignore_saturation
            # to process_one and ignore the saturation when the mask is actually
            # built, rather than untoggling the bit here?
            if ignore_saturation:  # Important for calibrations as we don't want replacement by 0
                indx = pypeitImage.boolean_mask(flag='SATURATION')
                pypeitImage.fullmask[indx] \
                        = pypeitImage.bitmask.turn_off(pypeitImage.fullmask[indx], 'SATURATION')
            # TODO: Come back to this.
            mask_stack[kk] = pypeitImage.fullmask

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

        # Coadd them
        var_stack = utils.inverse(ivar_stack)
        if combine_method == 'weightmean':
            weights = np.ones(self.nfiles, dtype=float)/self.nfiles
            img_list = [img_stack]
            var_list = [var_stack, rn2img_stack]
            img_list_out, var_list_out, gpm, nused \
                    = combine.weighted_combine(weights, img_list, var_list, mask_stack==0,
                                               sigma_clip=sigma_clip, sigma_clip_stack=img_stack,
                                               sigrej=sigrej, maxiters=maxiters)
        elif combine_method == 'median':
            bpm_stack = mask_stack > 0
            nstack = np.sum(np.logical_not(bpm_stack).astype(int), axis=0)
            gpm = nstack > 0
            img_list_out = [np.ma.median(np.ma.MaskedArray(img_stack, mask=bpm_stack),
                                         axis=0).filled(0.)]
            # First calculate the error in the sum
            var_list_out = [np.ma.sum(np.ma.MaskedArray(var_stack, mask=bpm_stack),
                                      axis=0).filled(0.)]
            var_list_out += [np.ma.sum(np.ma.MaskedArray(rn2img_stack, mask=bpm_stack),
                                       axis=0).filled(0.)]
            # Convert to standard error in the median (pi/2 factor relates standard variance
            # in mean (sum(variance_i)/n^2) to standard variance in median)
            var_list_out[0][gpm] *= np.pi/2/nstack[gpm]**2
            var_list_out[1][gpm] *= np.pi/2/nstack[gpm]**2
        else:
            # NOTE: Given the check at the beginning of the function, the code
            # should *never* make it here.
            msgs.error("Bad choice for combine.  Allowed options are 'median', 'weightmean'.")

        # Build the last one
        final_pypeitImage = pypeitimage.PypeItImage(img_list_out[0],
                                                    ivar=utils.inverse(var_list_out[0]),
                            # TODO: Why is this `pypeitImage.bpm` and not
                            # `np.logical_not(gpm)`?
                                                    bpm=pypeitImage.bpm,
                                                    rn2img=var_list_out[1],
                            # TODO: The mask can be set for a bunch of reasons.
                            # Why identify all of the masked pixels as cosmic rays?
                                                    crmask=np.logical_not(gpm),
                                                    detector=pypeitImage.detector,
                                                    PYP_SPEC=pypeitImage.PYP_SPEC)
        # Internals
        final_pypeitImage.rawheadlist = pypeitImage.rawheadlist
        final_pypeitImage.process_steps = pypeitImage.process_steps

        nonlinear_counts = self.spectrograph.nonlinear_counts(pypeitImage.detector,
                                                              apply_gain=self.par['apply_gain'])
        final_pypeitImage.build_mask(saturation=nonlinear_counts)
        # Return
        return final_pypeitImage

    @property
    def nfiles(self):
        """
        The number of files in :attr:`files`.
        """
        return len(self.files) if isinstance(self.files, (np.ndarray, list)) else 0


