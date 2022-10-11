""" Module for image combining

.. include:: ../include/links.rst
"""
import numpy as np

from astropy import stats

from pypeit import msgs
from pypeit import utils

from IPython import embed


# TODO make weights optional and do uniform weighting without.
def weighted_combine(weights, sci_list, var_list, inmask_stack,
                     sigma_clip=False, sigma_clip_stack=None, sigrej=None, maxiters=5):
    r"""
    Combine multiple sets of images, all using the same weights and mask.

    The multiple sets of images and variances to combine must have the same
    shape --- ``(nimgs, nspec, nspat)`` --- and this shape must match the
    provided *single* mask set (``inmask_stack``).  The provided weights are
    broadcast to the necessary shape (see below), where one can provide one
    weight per image, one weight per image spatial coordinate (i.e.,
    wavelength-dependent weights), or independent weights for each pixel.

    Optionally, the image stack can be sigma-clipped by setting
    ``sigma_clip=True``.  If sigma-clipping is requested and no sigma-rejection
    thresholds are provided (``sigrej`` is None), the sigma-rejection thresholds
    are set *automatically* depending on the number of images to combine.  The
    default rejection thresholds are 1.1, 1.3, 1.6, 1.9, or 2.0 for,
    respectively, 3, 4, 5, 6, or :math:`\geq 7` images.  Sigma-clipping cannot
    be performed if there are fewer than 3 images.  The pixel rejection is based
    on a *single* image stack provided by ``sigma_clip_stack``, which does not
    necessarily need to be any of the image stacks provided by ``sci_list``.
    Pixels rejected by sigma-clipping the ``sigma_clip_stack`` array are applied
    to all image stacks in ``sci_list``.

    The combined images are collected into the returned image list, where the
    order of the list is identical to the input ``sci_list``.  The returned mask
    and pixel accounting array is identical for all stacked images.
    
    Parameters
    ----------
    weights : `numpy.ndarray`_
        Weights to use. Options for the shape of weights are:

            - ``(nimgs,)``: a single weight per image in the stack

            - ``(nimgs, nspec)``: wavelength dependent weights per image in the
              stack

            - ``(nimgs, nspec, nspat)``: weights input with the shape of the
              image stack

        Note that the weights are distinct from the mask, which is dealt with
        via the ``inmask_stack`` argument, meaning there should not be any
        weights that are set to zero (although in principle this would still
        work).
    sci_list : :obj:`list`
        List of floating-point `numpy.ndarray`_ image groups to stack.  Each
        image group *must* have the same shape: ``(nimgs, nspec, nspat)``.
    var_list : :obj:`list`
        List of floating-point `numpy.ndarray`_ images providing the variance
        for each image group.  The number of image groups and the shape of each
        group must match ``sci_list``.  These are used to propagate the error in
        the combined images.
    inmask_stack : `numpy.ndarray`_, boolean, shape (nimgs, nspec, nspat)
        Good-pixel mask (True=Good, False=Bad) for the input image stacks.  This
        single group of good-pixel masks is applied to *all* input image groups.
    sigma_clip : :obj:`bool`, optional, default = False
        Combine with a mask by sigma clipping the image stack.  Stacks can only
        be sigma-clipped if there are 3 or more images.
    sigma_clip_stack : `numpy.ndarray`_, float, shape (nimgs, nspec, nspat), optional, default = None
        The image stack to be used for the sigma clipping. For example, if the
        list of images to be combined with the weights is ``[sciimg_stack,
        waveimg_stack, tilts_stack]`` and you want to clip based on
        ``sciimg_stack``, you would set ``sigma_clip_stack=sciimg_stack``.
    sigrej : :obj:`int`, :obj:`float`, optional, default = None
        Rejection threshold for sigma clipping.  If None and ``sigma_clip`` is
        True, the rejection threshold is set based on the number of images to
        combine; see above.  This value is passed directly to
        `astropy.stats.SigmaClip`_ as its ``sigma`` parameter.
    maxiters : :obj:`int`, optional, default=5
        Maximum number of rejection iterations; see `astropy.stats.SigmaClip`_.

    Returns
    -------
    sci_list_out : :obj:`list`
        The list of ndarray float combined images with shape ``(nspec, nspat)``.
    var_list_out : :obj:`list`
        The list of ndarray propagated variance images with shape ``(nspec,
        nspat)``.
    gpm : boolean `numpy.ndarray`_, shape (nspec, nspat)
        Good pixel mask for combined image (True=Good, False=Bad).
    nused : integer `numpy.ndarray`_, shape (nspec, nspat)
        Number of pixels combined at each location in the stacked images.
    """

    shape = img_list_error_check(sci_list, var_list)

    nimgs = shape[0]
    img_shape = shape[1:]
    #nspec = shape[1]
    #nspat = shape[2]

    if nimgs == 1:
        # If only one image is passed in, simply return the input lists of images, but reshaped
        # to be (nspec, nspat)
        msgs.warn('Cannot combine a single image. Returning input images')
        sci_list_out = []
        for sci_stack in sci_list:
            sci_list_out.append(sci_stack.reshape(img_shape))
        var_list_out = []
        for var_stack in var_list:
            var_list_out.append(var_stack.reshape(img_shape))
        gpm = inmask_stack.reshape(img_shape)
        nused = gpm.astype(int)
        return sci_list_out, var_list_out, gpm, nused

    if sigma_clip and nimgs >= 3:
        if sigma_clip_stack is None:
            msgs.error('You must specify sigma_clip_stack; sigma-clipping is based on this array '
                       'and propagated to the arrays to be stacked.')
        if sigrej is None:
            # NOTE: If these are changed, make sure to update the doc-string!
            if nimgs == 3:
                sigrej = 1.1
            elif nimgs == 4:
                sigrej = 1.3
            elif nimgs == 5:
                sigrej = 1.6
            elif nimgs == 6:
                sigrej = 1.9
            else:
                sigrej = 2.0
        # sigma clip if we have enough images
        # mask_stack > 0 is a masked value. numpy masked arrays are True for masked (bad) values
        data = np.ma.MaskedArray(sigma_clip_stack, mask=np.logical_not(inmask_stack))
        sigclip = stats.SigmaClip(sigma=sigrej, maxiters=maxiters, cenfunc='median',
                                  stdfunc=utils.nan_mad_std)
        data_clipped, lower, upper = sigclip(data, axis=0, masked=True, return_bounds=True)
        mask_stack = np.logical_not(data_clipped.mask)  # mask_stack = True are good values
    else:
        if sigma_clip and nimgs < 3:
            msgs.warn('Sigma clipping requested, but you cannot sigma clip with less than 3 '
                      'images.  Proceeding without sigma clipping')
        mask_stack = inmask_stack  # mask_stack = True are good values

    nused = np.sum(mask_stack, axis=0)
    weights_stack = broadcast_weights(weights, shape)
    weights_mask_stack = weights_stack*mask_stack.astype(float)

    weights_sum = np.sum(weights_mask_stack, axis=0)
    inv_w_sum = 1./(weights_sum + (weights_sum == 0.0))
    sci_list_out = []
    for sci_stack in sci_list:
        sci_list_out.append(np.sum(sci_stack * weights_mask_stack, axis=0) * inv_w_sum)
    var_list_out = []
    for var_stack in var_list:
        var_list_out.append(np.sum(var_stack * weights_mask_stack**2, axis=0) * inv_w_sum**2)
    # Was it masked everywhere?
    gpm = np.any(mask_stack, axis=0)

    return sci_list_out, var_list_out, gpm, nused



def img_list_error_check(sci_list, var_list):
    """
    Utility routine for dealing with lists of image stacks for
    :func:`weighted_combine`. This routine checks that the images sizes are
    correct and routines the shape of the image stacks.

    Parameters
    ----------
    sci_list : :obj:`list`
        List of  float `numpy.ndarray`_ images (each being an image stack with
        shape ``(nimgs, nspec, nspat)``) which are to be combined with the
        weights, mask, and possibly sigma clipping.
    var_list : :obj:`list`
        List of  float `numpy.ndarray`_ variance images (each being an image
        stack with shape ``(nimgs, nspec, nspat)``) which are to be combined
        with proper erorr propagation, i.e.  using the weights squared, mask,
        and possibly sigma clipping.

    Returns
    -------
    shape : :obj:`tuple`
        The shapes of the image stacks -- ``(nimgs, nspec, nspat)``

    """
    shape_sci_list = []
    for img in sci_list:
        shape_sci_list.append(img.shape)
        if img.ndim < 2:
            msgs.error('Dimensionality of an image in sci_list is < 2')

    shape_var_list = []
    for img in var_list:
        shape_var_list.append(img.shape)
        if img.ndim < 2:
            msgs.error('Dimensionality of an image in var_list is < 2')

    for isci in shape_sci_list:
        if isci != shape_sci_list[0]:
            msgs.error('An image in sci_list have different dimensions')
        for ivar in shape_var_list:
            if ivar != shape_var_list[0]:
                msgs.error('An image in var_list have different dimensions')
            if isci != ivar:
                msgs.error('An image in sci_list had different dimensions than an image in var_list')

    shape = shape_sci_list[0]

    return shape



def broadcast_weights(weights, shape):
    """
    Utility routine to broadcast weights to be the size of image stacks specified by shape

    Args:
        weights (`numpy.ndarray`_):
            Weights to use. Options for the shape of weights are:
                - (nimgs,) -- a single weight per image in the stack
                - (nimgs, nspec) -- wavelength dependent weights per
                  image
                - (nimgs, nspec, nspat) -- weights already have the
                  shape of the image stack and are simply returned
        shape (tuple):
            Shape of the image stacks for weighted coadding. This is either (nimgs, nspec) for 1d extracted spectra or
            (nimgs, nspec, nspat) for 2d spectrum images

    Returns:
        `numpy.ndarray`_:
            Weights for the stack images with output shape
            described in the Args above.

    """
    # Create the weights stack images from the wavelength dependent weights, i.e. propagate these
    # weights to the spatial direction
    if weights.ndim == 1:
        # One float per image
        if len(shape) == 2:
            weights_stack = np.einsum('i,ij->ij', weights, np.ones(shape))
        elif len(shape) == 3:
            weights_stack = np.einsum('i,ijk->ijk', weights, np.ones(shape))
        else:
            msgs.error('Image shape is not supported')
    elif weights.ndim == 2:
        # Wavelength dependent weights per image
        if len(shape) == 2:
            if weights.shape != shape:
                msgs.error('The shape of weights does not match the shape of the image stack')
            weights_stack = weights
        elif len(shape) == 3:
            weights_stack = np.einsum('ij,k->ijk', weights, np.ones(shape[2]))
    elif weights.ndim == 3:
        # Full image stack of weights
        if weights.shape != shape:
            msgs.error('The shape of weights does not match the shape of the image stack')
        weights_stack = weights
    else:
        msgs.error('Unrecognized dimensionality for weights')

    return weights_stack


def broadcast_lists_of_weights(weights, shapes):
    """
    Utility routine to broadcast weights to be the size of image stacks specified by shape

    Parameters
    ----------
    weights : :obj:`list`
        List containing the weights to use. The length of weights must be nimgs,
        the number of images that are being combined. The options for the date
        type/shape for the individual elements of weights are:

            - :obj:`float`: a single weight per image in the stack

            - `numpy.ndarray`_ with shape ``(nspec,)``: wavelength dependent
              weights per image in the stack

            - `numpy.ndarray`_ with shape ``(nspec, nspat)``: weights input with
              the shape of each image stack and will be simply be returned

    shapes : :obj:`list`
        List with length of nimgs containing the tuples which are the shapes
        ``(nspec, nspat)`` of each image in the stack that should have their
        weights broadcast.


    Returns
    -------

    weights_list : :obj:`list` of `numpy.ndarray`_ objects
        Weight images where each image in the list has shape that was input via
        the shapes input parameter.

    """
    # Create the weights stack images from the wavelength dependent weights, i.e. propagate these
    # weights to the spatial direction

    weights_list = []
    for weight, shape in zip(weights, shapes):
        if isinstance(weight, float):
            weights_list.append(np.ones(shape, dtype=float) * weight)
        elif isinstance(weight, np.ndarray):
            if weight.ndim == 1:
                weights_list.append(np.broadcast_to(weight[:, np.newaxis], shape))
            elif weight.ndim == 2:
                weights_list.append(weight)
            else:
                msgs.error('Weights must be a float or a 1D or 2D ndarray')

    return weights_list

