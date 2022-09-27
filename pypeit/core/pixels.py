""" Routines related to mapping pixels to physical positions

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
from IPython import embed

import numpy as np

from pypeit import msgs


def phys_to_pix(array, pixlocn, axis):
    """
    Generate an array of physical pixel coordinates

    Parameters
    ----------
    array : `numpy.ndarray`_
      An array of physical pixel locations
    pixlocn : `numpy.ndarray`_
      A 3D array containing the x center, y center, x width and y width of each pixel.
    axis : int
      The axis that the input array probes

    Returns
    -------
    pixarr : ndarray
      The pixel locations of the input array (as seen on a computer screen)
    """
    if len(array.shape) > 2:
        msgs.error('Input array must have two dimensions or less!')
    _array = np.atleast_2d(array)
    doravel = len(array.shape) != 2

    diff = pixlocn[:,0,0] if axis == 0 else pixlocn[0,:,1]

    pix = np.argmin(np.absolute(_array[:,:,None] - diff[None,None,:]), axis=2).astype(int)
    return pix.ravel() if doravel else pix



# ToDO rewrite this function to use images rather than loops as in flat_fit.py
# TODO: This is used by core/extract.py and core/skysub.py.
# Flat-fielding now uses SlitTraceSet methods instead. Merge the usage?
def ximg_and_edgemask(lord_in, rord_in, slitpix, trim_edg=(3,3), xshift=0.):
    """
    Generate the ximg and edgemask frames

    Parameters
    ----------
    lord_in : `numpy.ndarray`_
        Array containing the left trace. This can either be a 2-d array with shape (nspec, nTrace)
        for multiple traces, or simply a 1-d array with shape  (nspec) for a single trace.

    rord_in : `numpy.ndarray`_
        Array containing the right trace. This can either be a 2-d array with shape (nspec, nTrace)
        for multiple traces, or simply a 1-d array with shape  (nspec) for a single trace.

    slitpix : `numpy.ndarray`_
        Image with shape (nspec, nspat) specifying pixel locations. This is created by core_slit_pixels above.

    trim_edg : tuple 
        How much to trim off each edge of each slit in pixels.
        integer or floats

    xshift : float, optional
        Future implementation may need to shift the edges

    Returns
    -------
    ximg : `numpy.ndarray`_
      Specifies spatial location of pixel in its own slit
      Scaled from 0 to 1
    edgemask : ndarray, bool
      True = Masked because it is too close to the edge
    """
    #; Generate the output image
    ximg = np.zeros_like(slitpix, dtype=float)
    # Intermediary arrays
    pixleft = np.zeros_like(slitpix, dtype=float)
    pixright = np.zeros_like(slitpix, dtype=float)
    #

    # This little bit of code allows the input lord and rord to either be (nspec, nslit) arrays or a single
    # vectors of size (nspec)
    if lord_in.ndim == 2:
        nslit = lord_in.shape[1]
        lord = lord_in
        rord = rord_in
    else:
        nslit = 1
        lord = lord_in.reshape(lord_in.size,1)
        rord = rord_in.reshape(rord_in.size,1)


    #; Loop over each slit
    for islit in range(nslit):
        #; How many pixels wide is the slit at each Y?
        xsize = rord[:, islit] - lord[:, islit]
        badp = xsize <= 0.
        if np.any(badp):
            meds = np.median(xsize)
            msgs.warn('Something goofy in slit # {:d}'.format(islit))
            msgs.warn('Probably a bad slit (e.g. a star box)')
            msgs.warn('It is best to expunge this slit')
            msgs.warn('Proceed at your own risk, with a slit width of {}'.format(meds))
            msgs.warn('Or set meds to your liking')
            #rord[:, islit] = lord[:, islit] + meds

        # Loop down the slit
        for iy in range(lord.shape[0]):
            # Set the pixels
            ix1 = max(int(np.ceil(lord[iy, islit])),0)
            ix2 = min(int(rord[iy, islit]),ximg.shape[1]-1)
            ix1 = min(ix1,ximg.shape[1] - 1)
            ix2 = max(ix2, 0)
            #
            ximg[iy, ix1:ix2+1] = (np.arange(ix2 - ix1 + 1) + ix1 - lord[iy, islit]) / xsize[iy]
            pixleft[iy, ix1:ix2+1] = (np.arange(ix2 - ix1 + 1) + ix1 - lord[iy, islit])
            pixright[iy, ix1:ix2+1] = (rord[iy, islit] - ix2 + np.flip(np.arange(ix2-ix1+1),0))

    # Generate the edge mask
    edgemask = (slitpix > 0) & np.any([pixleft < trim_edg[0], pixright < trim_edg[1]], axis=0)
    # Return
    return ximg, edgemask


