""" Module for image combining
"""
import numpy as np

from astropy import stats

from pypeit import msgs
from pypeit import utils

from IPython import embed


def masked_weightmean(a, maskvalue):
    """
    .. todo::
        Document this!
    """
    num = np.ma.MaskedArray(a.copy(), mask=(a==maskvalue))
    num[np.invert(num.mask) & (num <= 1.0)] = 0.0
    num = np.ma.sum(np.ma.sqrt(num)*num, axis=2)
    den = np.ma.MaskedArray(a.copy(), mask=(a==maskvalue))
    den[np.invert(den.mask) & (den <= 1.0)] = 1.0
    den = np.ma.sum(np.sqrt(den), axis=2)
    return np.ma.divide(num, den).filled(maskvalue)


def maxnonsat(array, saturated):
    """
    .. todo::
        Document this!
    """
    minimum = np.amin(np.clip(array, None, saturated), axis=2)
    _array = np.ma.MaskedArray(array, mask=np.invert((array > 0.0) & (array<saturated)))
    maximum = np.ma.amax(_array, axis=2)
    maximum[maximum.mask] = minimum[maximum.mask]
    return maximum.data



# TODO make weights optional and do uniform weighting without.
def weighted_combine(weights, sci_list, var_list, inmask_stack,
                     sigma_clip=False, sigma_clip_stack=None, sigrej=None, maxiters=5):
    """

    Args:
        weights (ndarray):
            Weights to use. Options for the shape of weights are:

                - (nimgs,) -- a single weight per image in the stack
                - (nimgs, nspec) -- wavelength dependent weights per
                  image in the stack
                - (nimgs, nspec, nspat) -- weights input with the shape
                  of the image stack

             Note that the weights are distinct from the mask which is
             dealt with via inmask_stack argument so there should not be
             any weights that are set to zero (although in principle
             this would still work).

        sci_list: list
            List of  float ndarray images (each being an image stack with shape (nimgs, nspec, nspat))
            which are to be combined with the  weights, inmask_stack, and possibly sigma clipping
        var_list: list
            List of  float ndarray variance images (each being an image stack with shape (nimgs, nspec, nspat))
            which are to be combined with proper erorr propagation, i.e.
            using the  weights**2, inmask_stack, and possibly sigma clipping
        inmask_stack: ndarray, boolean, shape (nimgs, nspec, nspat)
            Array of input masks for the images. True = Good, False=Bad
        sigma_clip: bool, default = False
            Combine with a mask by sigma clipping the image stack. Only valid if nimgs > 2
        sigma_clip_stack: ndarray, float, shape (nimgs, nspec, nspat), default = None
            The image stack to be used for the sigma clipping. For example if
            if the list of images to be combined with the weights is [sciimg_stack, waveimg_stack, tilts_stack] you
            would be sigma clipping with sciimg_stack, and would set sigma_clip_stack = sciimg_stack
        sigrej: int or float, default = None
            Rejection threshold for sigma clipping. Code defaults to determining this automatically based
            on the numberr of images provided.
        maxiters:
            Maximum number of iterations for sigma clipping using astropy.stats.SigmaClip

    Returns:
        tuple: Returns the following:
            - sci_list_out: list: The list of ndarray float combined
              images with shape (nspec, nspat)
            - var_list_out: list: The list of ndarray propagated
              variance images with shape (nspec, nspat)
            - gpm: bool ndarray, shape (nspec, nspat): Good pixel mask for
              combined image. True=Good, False=Bad
            - nused: int ndarray, shape (nspec, nspat): Image of
              integers indicating the number of images that contributed
              to each pixel
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
            msgs.error('You must specify sigma_clip_stack, i.e. which quantity to use for sigma clipping')
        if sigrej is None:
            if nimgs <= 2:
                sigrej = 100.0  # Irrelevant for only 1 or 2 files, we don't sigma clip below
            elif nimgs == 3:
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
        data = np.ma.MaskedArray(sigma_clip_stack, np.invert(inmask_stack))
        sigclip = stats.SigmaClip(sigma=sigrej, maxiters=maxiters, cenfunc='median', stdfunc=utils.nan_mad_std)
        data_clipped, lower, upper = sigclip(data, axis=0, masked=True, return_bounds=True)
        mask_stack = np.invert(data_clipped.mask)  # mask_stack = True are good values
    else:
        if sigma_clip and nimgs < 3:
            msgs.warn('Sigma clipping requested, but you cannot sigma clip with less than 3 images. '
                      'Proceeding without sigma clipping')
        mask_stack = inmask_stack  # mask_stack = True are good values

    nused = np.sum(mask_stack, axis=0)
    weights_stack = broadcast_weights(weights, shape)
    weights_mask_stack = weights_stack*mask_stack

    weights_sum = np.sum(weights_mask_stack, axis=0)
    sci_list_out = []
    for sci_stack in sci_list:
        sci_list_out.append(np.sum(sci_stack*weights_mask_stack, axis=0)/(weights_sum + (weights_sum == 0.0)))
    var_list_out = []
    for var_stack in var_list:
        var_list_out.append(np.sum(var_stack * weights_mask_stack**2, axis=0) / (weights_sum + (weights_sum == 0.0))**2)
    # Was it masked everywhere?
    gpm = np.any(mask_stack, axis=0)

    return sci_list_out, var_list_out, gpm, nused


def img_list_error_check(sci_list, var_list):
    """
    Utility routine for dealing dealing with lists of image stacks for rebin2d and weigthed_combine routines below. This
    routine checks that the images sizes are correct and routines the shape of the image stacks.

    Args:
        sci_list: list
            List of  float ndarray images (each being an image stack with shape (nimgs, nspec, nspat))
            which are to be combined with the  weights, inmask_stack, and possibly sigma clipping
        var_list: list
            List of  float ndarray variance images (each being an image stack with shape (nimgs, nspec, nspat))
            which are to be combined with proper erorr propagation, i.e.
            using the  weights**2, inmask_stack, and possibly sigma clipping

    Returns:
        tuple: The shapes of the image stacks, (nimgs, nspec, nspat)

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
        weights (ndarray):
            Weights to use. Options for the shape of weights are:

                - (nimgs,) -- a single weight per image in the stack
                - (nimgs, nspec) -- wavelength dependent weights per
                  image
                - (nimgs, nspec, nspat) -- weights already have the
                  shape of the image stack and are simply returned
        shape: tuple of integers
            Shape of the image stacks for weighted coadding. This is either (nimgs, nspec) for 1d extracted spectra or
            (nimgs, nspec, nspat) for 2d spectrum images

    Returns:
        np.ndarray:

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

