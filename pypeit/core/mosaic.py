"""
Provide basic mosaicing functions.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from IPython import embed

import numpy as np
from scipy import ndimage

from pypeit import msgs
from pypeit.core import transform
from pypeit.utils import inverse


def build_image_mosaic_transform(shape, shift, rot, binning):
    """
    Build the affine transform of a binned image.

    The order of operaions is as follows.  Steps 1, 3, and 5 are only
    performed if the image is rotated.  Steps 2 and 4 are only done if
    the binning is not square.

        #. Shift the coordinate system to the center of the image,
           assuming the (0,0) coordinate is at the center of the first
           pixel.

        #. For binning that is different in each dimension, scale the
           size of each pixel back to the correct aspect ratio.

        #. Rotate the image.

        #. Undo the dimension scaling to account for the binning aspect
           ratio.

        #. Undo the shift to the center of the image.

        #. Apply the requested shift, accounting for the binning.

    These steps are compiled into a single transformation matrix using
    :func:`pypeit.core.transform.affine_transform_series`.

    Args:
        shape (:obj:`tuple`):
            A two-tuple with the shape of the **binned** image.
        shift (:obj:`tuple`):
            A two-tuple with the nominal shift of the *unbinned* image
            in the mosaic in each dimension.
        rot (:obj:`float`):
            The counter-clockwise rotation in degrees of the
            **unbinned** image in the mosaic.
        binning (:obj:`tuple`):
            The number of pixels binned in each dimension.  This only
            has an effect on the results when the binning is not the
            same in both dimensions.

    Returns:
        `numpy.ndarray`_: The single coordinate transformation matrix
        that applies all transformations.  See
        :func:`pypeit.core.transform.affine_transform_series`.
    """
    tform = []
    if np.absolute(rot) > 0:
        # Offset to the center of the image
        tform += [dict(translation=(-(shape[0]-1)/2, -(shape[1]-1)/2))]
        if binning[0] != binning[1]:
            # Rescale back to square pixels
            tform += [dict(scale=(1.,binning[1]/binning[0]))]
        # Apply the rotation
        tform += [dict(rotation=np.radians(rot))]
        if binning[0] != binning[1]:
            # Undo the bin scaling
            tform += [dict(scale=(1.,binning[0]/binning[1]))]
        # Undo the offset to the image center
        tform += [dict(translation=((shape[0]-1)/2, (shape[1]-1)/2))]
    # Apply the shift
    tform += [dict(translation=(shift[0]/binning[0], shift[1]/binning[1]))]
    # Compile into a single transformation and return
    return transform.affine_transform_series(tform)


def prepare_mosaic(shape, tforms, buffer=0, inplace=False):
    r"""
    Prepare to mosaic images by determining the shape of the mosaic
    image and adjusting the transformation coordinates to the mosaic
    pixel coordinates.

    Args:
        shape (:obj:`tuple`):
            A two-tuple with the shape of the images to mosaic.
        tforms (:obj:`list`):
            A list of :math:`3\times3` `numpy.ndarray`_ objects with the
            transformations to apply to each image.  These are adjusted
            as necessary to perform the transformations within the
            coordinate system of the mosaic image.
        buffer (:obj:`int`, optional):
            An added buffer in each dimension that frames the mosaic
            image and should not have any image data.  Buffer pixels are
            set to 0.  Buffer must be non-negative.
        inplace (:obj:`bool`, optional):
            If True, alter the provided ``tforms`` in-place.  Otherwise,
            the returned transforms are new arrays.
    
    Returns:
        :obj:`tuple`: Returns the shape for the mosaic image as a
        two-tuple and the new transformation matrices as a list of
        `numpy.ndarray`_ objects.
    """
    # Use the number of transforms to set the number of images
    nimg = len(tforms)

    # Get a box that bounds the transformed coordinates of all the mosaic
    # images.
    coo = np.array([[-0.5,-0.5], [shape[0]+0.5,-0.5], [shape[0]+0.5, shape[1]+0.5],
                    [-0.5, shape[1]+0.5]]).astype(float)
    box = None
    for i in range(nimg):
        tc = transform.coordinate_transform_2d(coo, tforms[i], inverse=False)
        if box is None:
            box = np.vstack((np.floor(np.amin(tc, axis=0)), np.ceil(np.amax(tc, axis=0))))
            continue
        box[0] = np.amin(np.vstack((tc,box[0])), axis=0)
        box[1] = np.amax(np.vstack((tc,box[1])), axis=0)

    # Set the mosaic image shape
    if buffer < 0:
        msgs.error('Mosaic image buffer must be >= 0.')
    mosaic_shape = tuple(*np.ceil(np.diff(box, axis=0) + 2*buffer - 1).astype(int))

    # Adjust the image transformations to be within the limits of the mosaic
    # image.

    # NOTE: There's a subtlety in this as follows.  In order to reproduce the
    # mosaicing results for Gemini GMOS from DRAGONS (see
    # pypeit/tests/test_mosaic.py), the offsets below need to be whole pixels.
    # But this leads to an asymmetry in the transforms such that one can't
    # recover the input image frame by inversing the transform.  For now, I've
    # chosen the approach that gets closest to the DRAGONS result.
    _tforms = tforms if inplace else [None]*nimg
    for i in range(nimg):
        _tforms[i] = transform.affine_transform_series([dict(translation=(-(box[0,0]-buffer), #-0.5,
                                                                          -(box[0,1]-buffer))) #-0.5))
                                                       ]) @ tforms[i]
    return mosaic_shape, _tforms


def build_image_mosaic(imgs, tforms, ivar=None, bpm=None, mosaic_shape=None, cval=0., order=0,
                       overlap='combine'):
    r"""
    Use the provided images and transformation matrices to construct an image
    mosaic.

    .. warning::

        Beware when using ``order > 0``!

        Bad-pixel masks are *always* mapped to the mosaic image using
        ``order=0`` (i.e., without interpolation).  However, masked pixels are
        not excluded from the input images during the transformation.  For
        higher order interpolations (``order > 0``), this means that the masked
        pixels can contribute to the interpolation for any given output pixel.
        Users should appropriately consider how these pixels will affect the
        mosaic pixels *before* calling this function.

        Similarly, error propagation from the input image to the mosaic image is
        only approximate when ``order > 0``.  Error propagaion is performed
        simply by applying the coordinate transform to each variance image with
        the same order as used for the input image, and then combining those
        variances as necessary in overlap regions.

        Tests show that this approach is also not invertable.  I.e., iteratively transforming the image back and forth between the native and mosaic frames lead to image drifts.

    Args:
        imgs (:obj:`list`, `numpy.ndarray`_):
            List of `numpy.ndarray`_ images to include in the mosaic.  The shape
            of all the input images must be identical if ``mosaic_shape`` is
            None.
        tforms (:obj:`list`, `numpy.ndarray`_):
            List of `numpy.ndarray`_ objects with the transformation matrices
            necessary to convert between image and mosaic coordinates.  See
            :func:`pypeit.core.mosaic.build_image_mosaic_transform`.  The number
            of transforms must match the number of images.  If ``mosaic_shape``
            is None, the transforms are considered in a relative sense.  That
            is, the shape of the output mosaic is determined by applying these
            transforms to the bounding boxes of each image and then determining
            the shape needed to retain all pixels in the input images.  The
            transforms are then adjusted appropriately to map to this shape; see
            :func:`~pypeit.core.mosaic.prepare_mosaic`.  If ``mosaic_shape`` is
            *not* None, these transforms are expected to map directly to the
            output mosaic coordinates.
        ivar (:obj:`list`, `numpy.ndarray`_, optional):
            List of `numpy.ndarray`_ images with the inverse variance of the
            image data.  The number of inverse-variance images must match the
            number of images in the mosaic.  If None, inverse variance is
            returned as None.
        bpm (:obj:`list`, `numpy.ndarray`_, optional):
            List of boolean `numpy.ndarray`_ objects with the bad-pixel mask for
            each image in the mosaic.  The number of bad-pixel masks must match
            the number of images in the mosaic.  If None, all input pixels are
            considered valid.
        mosaic_shape (:obj:`tuple`, optional):
            Shape for the output image.  If None, the shape is determined by
            :func:`pypeit.core.mosaic.prepare_mosaic` and the shape of all the
            input images *must* be identical.
        cval (:obj:`float`, optional):
            The value used to fill empty pixels in the mosaic.
        order (:obj:`int`, optional):
            Order of the interpolation of each input image onto the mosaic grid.
            This is passed directly to `scipy.ndimage.affine_transform`_.
        overlap (:obj:`str`, optional):
            Keyword that indicates how to handle pixels in the regions where
            multiple images overlap in the mosaic.  Options are:

                - ``'combine'``: Average the values of the pixels and, if the
                  inverse variance is provided, propagate the error.
                - ``'error'``: Raise an exception.  Largely provided for testing
                  under the expectation that *no* pixels should overlap in the
                  mosaic.

    Returns:
        :obj:`tuple`: Four objects are returned. The first three are
        `numpy.ndarray`_ objects with the mosaic image, its inverse variance
        (None if no inverse variance is provided), and an integer array with the
        number of input pixels in each output pixel.  The last contains the
        detailed transformation matrices applied to each image.  If
        ``mosaic_shape`` is provided, these are identical to the input
        ``tforms``; otherwise, these are the transforms adjusted from the
        relative frame to the absolute mosaic frame given its determined shape;
        see :func:`~pypeit.core.mosaic.prepare_mosaic`.
    """
    # Check the input
    nimg = len(imgs)
    if len(tforms) != nimg:
        msgs.error('Number of image transformations does not match number of images to mosaic.')
    if ivar is not None and len(ivar) != nimg:
        msgs.error('If providing any, must provide inverse-variance for each image in the mosaic.')
    if bpm is not None and len(bpm) != nimg:
        msgs.error('If providing any, must provide bad-pixel masks for each image in the mosaic.')
    if overlap not in ['combine', 'error']:
        msgs.error(f'Unknown value for overlap ({overlap}), must be "combine" or "error".')

    # Get the output shape, if necessary
    if mosaic_shape is None:
        shape = imgs[0].shape
        if not np.all([img.shape == shape for img in imgs]):
            msgs.error('If output mosaic shape is not provided, all input images must have the '
                       'same shape!')
        mosaic_shape, _tforms = prepare_mosaic(shape, tforms)
    else:
        _tforms = tforms

    msgs.info(f'Constructing image mosaic with {nimg} images and output shape {mosaic_shape}.')

    if ivar is not None:
        var = [inverse(_ivar) for _ivar in ivar]

    mosaic_npix = np.zeros(mosaic_shape, dtype=int)
    mosaic_data = np.zeros(mosaic_shape, dtype=float)
    # NOTE: "mosaic_ivar" is actually the variance until it's inverted just
    # before output
    mosaic_ivar = None if ivar is None else np.zeros(mosaic_shape, dtype=float)

    # TODO: These loops can end up creating and destroying lots of big arrays.
    # Is there a way to make this faster?
    for i in range(nimg):
        _inv_tform = np.linalg.inv(_tforms[i])
        # TODO: I'm not a fan of using cval=np.nan, but this might be the most
        # robust way to catch mosaic pixels that are not filled by the input
        # image.  I'm not sure if we lose efficiency by using NaN instead of a
        # real scalar.
        _tform_img = ndimage.affine_transform(imgs[i], _inv_tform, output_shape=mosaic_shape,
                                              cval=np.nan, order=order)
        filled = np.logical_not(np.isnan(_tform_img))
        if bpm is not None:
            _tform_gpm = ndimage.affine_transform(np.logical_not(bpm[i]).astype(float), _inv_tform,
                                                  output_shape=mosaic_shape, cval=0., order=0)
            filled &= _tform_gpm > 0.
        mosaic_data[filled] += _tform_img[filled]
        mosaic_npix[filled] += 1
        if ivar is not None:
            # NOTE: "mosaic_ivar" is actually the variance until it's inverted
            # just before output
            mosaic_ivar[filled] += ndimage.affine_transform(var[i], _inv_tform,
                                                            output_shape=mosaic_shape,
                                                            cval=0., order=order)[filled]

    # TODO: This test is crude.  Input and output pixel sizes should be
    # identical, but I don't know if the order=0 approach will ever lead to a
    # single mosaic pixel being filled by more than one input pixel.  If so,
    # `overlap='error'` will catch both those cases and when multiple input
    # images overlap.
    has_overlap = np.any(mosaic_npix > 1)
    if has_overlap and overlap == 'error':
        # Input images should not be allowed to overlap
        msgs.error('Mosaic has pixels with contributions by more than one input image!')

    filled = mosaic_npix > 0
    # Average the overlapping pixels
    if has_overlap:
        mosaic_data[filled] /= mosaic_npix[filled]
    if mosaic_ivar is not None:
        # Propagate the error by averaging the variances
        if has_overlap:
            mosaic_ivar[filled] /= mosaic_npix[filled]
        # Revert to inverse variance
        mosaic_ivar = inverse(mosaic_ivar)

    return mosaic_data, mosaic_ivar, mosaic_npix, _tforms


