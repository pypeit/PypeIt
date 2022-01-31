""" Uber object for calibration images, e.g. arc, flat

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import os
import numpy as np

from pypeit import msgs
from pypeit.par import pypeitpar
from pypeit.images import combineimage
from pypeit.images import pypeitimage
from pypeit.core import procimg
from pypeit.core.framematch import valid_frametype
from pypeit import utils

from IPython import embed


class ArcImage(pypeitimage.PypeItImage):
    """
    Simple DataContainer for the Arc Image
    """
    # Peg the version of this class to that of PypeItImage
    version = pypeitimage.PypeItImage.version

    # I/O
    output_to_disk = ('ARC_IMAGE', 'ARC_FULLMASK', 'ARC_DETECTOR')
    hdu_prefix = 'ARC_'

    # Master fun
    master_type = 'Arc'
    master_file_format = 'fits'


class AlignImage(pypeitimage.PypeItImage):
    """
    Simple DataContainer for the Alignment Image
    """
    # Peg the version of this class to that of PypeItImage
    version = pypeitimage.PypeItImage.version

    # I/O
    output_to_disk = ('ALIGN_IMAGE', 'ALIGN_FULLMASK', 'ALIGN_DETECTOR')
    hdu_prefix = 'ALIGN_'

    # Master fun
    master_type = 'Align'
    master_file_format = 'fits'


class BiasImage(pypeitimage.PypeItImage):
    """
    Simple DataContainer for the Bias Image
    """
    # Set the version of this class
    version = pypeitimage.PypeItImage.version

    # Output to disk
    output_to_disk = ('BIAS_IMAGE', 'BIAS_IVAR', 'BIAS_DETECTOR')
    hdu_prefix = 'BIAS_'
    master_type = 'Bias'
    master_file_format = 'fits'


class DarkImage(pypeitimage.PypeItImage):
    """
    Simple DataContainer for the Dark Image
    """
    # Set the version of this class
    version = pypeitimage.PypeItImage.version

    # Output to disk
    output_to_disk = ('DARK_IMAGE', 'DARK_IVAR', 'DARK_DETECTOR')
    hdu_prefix = 'DARK_'
    master_type = 'Dark'
    master_file_format = 'fits'


class TiltImage(pypeitimage.PypeItImage):
    """
    Simple DataContainer for the Tilt Image
    """

    # Peg the version of this class to that of PypeItImage
    version = pypeitimage.PypeItImage.version

    # I/O
    output_to_disk = ('TILT_IMAGE', 'TILT_FULLMASK', 'TILT_DETECTOR')
    hdu_prefix = 'TILT_'

    # Master fun
    master_type = 'Tiltimg'
    master_file_format = 'fits'


class TraceImage(pypeitimage.PypeItImage):
    """
    Simple DataContainer for the Trace Image
    """

    # Peg the version of this class to that of PypeItImage
    version = pypeitimage.PypeItImage.version

    # I/O
    output_to_disk = ('TRACE_IMAGE', 'TRACE_FULLMASK', 'TRACE_DETECTOR')
    hdu_prefix = 'TRACE_'


class SkyRegions(pypeitimage.PypeItImage):
    """
    Simple DataContainer for the SkyRegions Image
    """
    # Peg the version of this class to that of PypeItImage
    version = pypeitimage.PypeItImage.version

    # I/O
    output_to_disk = ('SKYREG_IMAGE')
    hdu_prefix = 'SKYREG_'

    # Master fun
    master_type = 'SkyRegions'
    master_file_format = 'fits.gz'


# Convert frame type into an Image
frame_image_classes = dict(
    bias=BiasImage,
    dark=DarkImage,
    arc=ArcImage,
    tilt=TiltImage,
    trace=TraceImage,
    align=AlignImage)


def buildimage_fromlist(spectrograph, det, frame_par, file_list, bias=None, bpm=None, dark=None,
                        flatimages=None, maxiters=5, ignore_saturation=True, slits=None,
                        mosaic=None):
    """
    Perform basic image processing on a list of images and combine the results.

    .. warning::

        For image mosaics (when ``det`` is a tuple) the processing behavior is
        hard-coded such that bias and dark frames are *not* reformatted into a
        mosaic image.  They are saved in their native multi-image format.
        Bad-pixel masks are also expected to be in multi-image format.  See
        :class:`~pypeit.images.rawimage.RawImage`.

    Args:
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            Spectrograph used to take the data.
        det (:obj:`int`, :obj:`tuple`):
            The 1-indexed detector number(s) to process.  If a tuple, it must
            include detectors viable as a mosaic for the provided spectrograph;
            see :func:`~pypeit.spectrographs.spectrograph.Spectrograph.allowed_mosaics`.
        frame_par (:class:`~pypeit.par.pypeitpar.FramePar`):
            Parameters that dictate the processing of the images.  See
            :class:`~pypeit.par.pypeitpar.ProcessImagesPar` for the
            defaults.
        file_list (:obj:`list`):
            List of files
        bias (:class:`~pypeit.images.buildimage.BiasImage`, optional):
            Bias image for bias subtraction; passed directly to
            :func:`~pypeit.images.rawimage.RawImage.process` for all images.
        bpm (`numpy.ndarray`_, optional):
            Bad pixel mask; passed directly to
            :func:`~pypeit.images.rawimage.RawImage.process` for all images.
        dark (:class:`~pypeit.images.buildimage.DarkImage`, optional):
            Dark-current image; passed directly to
            :func:`~pypeit.images.rawimage.RawImage.process` for all images.
        flatimages (:class:`~pypeit.flatfield.FlatImages`, optional):
            Flat-field images for flat fielding; passed directly to
            :func:`~pypeit.images.rawimage.RawImage.process` for all images.
        maxiters (:obj:`int`, optional):
            When ``combine_method='mean'``) and sigma-clipping
            (``sigma_clip`` is True), this sets the maximum number of
            rejection iterations.  If None, rejection iterations continue
            until no more data are rejected; see
            :func:`~pypeit.core.combine.weighted_combine``.
        ignore_saturation (:obj:`bool`, optional):
            If True, turn off the saturation flag in the individual images
            before stacking.  This avoids having such values set to 0, which
            for certain images (e.g. flat calibrations) can have unintended
            consequences.
        slits (:class:`~pypeit.slittrace.SlitTraceSet`, optional):
            Edge traces for all slits.  These are used to calculate spatial
            flexure between the image and the slits, and for constructing the
            slit-illumination correction.  See
            :class:`pypeit.images.rawimage.RawImage.process`.

    Returns:
        :class:`~pypeit.images.pypeitimage.PypeItImage`:  The processed and
        combined image.
    """
    # Check
    if not isinstance(frame_par, pypeitpar.FrameGroupPar):
        msgs.error('Provided ParSet for must be type FrameGroupPar.')
    if not valid_frametype(frame_par['frametype'], quiet=True):
        # NOTE: This should not be necessary because FrameGroupPar explicitly
        # requires frametype to be valid
        msgs.error(f'{frame_par["frametype"]} is not a valid PypeIt frame type.')

    # Should the detectors be reformatted into a single image mosaic?
    if mosaic is None:
        mosaic = isinstance(det, tuple) and frame_par['frametype'] not in ['bias', 'dark']

    # Do it
    combineImage = combineimage.CombineImage(spectrograph, det, frame_par['process'], file_list)
    pypeitImage = combineImage.run(bias=bias, bpm=bpm, dark=dark, flatimages=flatimages,
                                   sigma_clip=frame_par['process']['clip'],
                                   sigrej=frame_par['process']['comb_sigrej'],
                                   maxiters=maxiters, ignore_saturation=ignore_saturation,
                                   slits=slits, combine_method=frame_par['process']['combine'],
                                   mosaic=mosaic)
    # Decorate according to the type of calibration, primarily as needed for
    # handling MasterFrames.  WARNING: Any internals (i.e., the ones defined by
    # the _init_internals method) in pypeitImage are lost here.
    return frame_image_classes[frame_par['frametype']].from_pypeitimage(pypeitImage) \
            if frame_par['frametype'] in frame_image_classes.keys() else pypeitImage



