""" Routines related to mapping pixels to physical positions
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np

from pypeit import msgs
from pypeit import debugger

try:
    from pypeit import ginga
except ImportError:
    pass


def gen_pixloc(frame_shape, xgap=0, ygap=0, ysize=1., gen=True):
    """
    Generate an array of physical pixel coordinates

    Parameters
    ----------
    frame : ndarray
      uniformly illuminated and normalized flat field frame
    xgap : int (optional)
    ygap : int (optional)
    ysize : float (optional)
    gen : bool, optional
      Only allows True right now

    Returns
    -------
    locations : ndarray
      A 3D array containing the x center, y center, x width and y width of each pixel.
      The returned array has a shape:   frame.shape + (4,)
    """
    #dnum = settings.get_dnum(det)
    msgs.info("Deriving physical pixel locations on the detector")
    locations = np.zeros((frame_shape[0],frame_shape[1],4))
    if gen:
        msgs.info("Pixel gap in the dispersion direction = {0:4.3f}".format(xgap))
        msgs.info("Pixel size in the dispersion direction = {0:4.3f}".format(1.0))
        xs = np.arange(frame_shape[0]*1.0)*xgap
        xt = 0.5 + np.arange(frame_shape[0]*1.0) + xs
        msgs.info("Pixel gap in the spatial direction = {0:4.3f}".format(ygap))
        msgs.info("Pixel size in the spatial direction = {0:4.3f}".format(ysize))
        ys = np.arange(frame_shape[1])*ygap*ysize
        yt = ysize*(0.5 + np.arange(frame_shape[1]*1.0)) + ys
        xloc, yloc = np.meshgrid(xt, yt)
#		xwid, ywid = np.meshgrid(xs,ys)
        msgs.info("Saving pixel locations")
        locations[:,:,0] = xloc.T
        locations[:,:,1] = yloc.T
        locations[:,:,2] = 1.0
        locations[:,:,3] = ysize
    else:
        msgs.error("Have not yet included an algorithm to automatically generate pixel locations")
    return locations


def phys_to_pix(array, pixlocn, axis):
    """
    Generate an array of physical pixel coordinates

    Parameters
    ----------
    array : ndarray
      An array of physical pixel locations
    pixlocn : ndarray
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


'''
def slit_pixels(slf, frameshape, det):
    """ Wrapper to the core_slit_pixels method
    May be Deprecated in a future Refactor

    Parameters
    ----------
    slf : class
      Science Exposure Class
    frameshape : tuple
      A two element tuple providing the shape of a trace frame.
    det : int
      Detector index

    Returns
    -------

    """
    return core_slit_pixels(slf._lordloc[det-1], slf._rordloc[det-1], frameshape,
                settings.argflag['trace']['slits']['pad'])
'''

# ToDo add this function to the instrument class to enable ordermasks and trimming or orders
# ToDO rewrite this function to use images rather than loops as in flat_fit.py
def slit_pixels(all_lordloc_in, all_rordloc_in, frameshape, pad):
    """ Generate an image indicating the slit/order associated with each pixel.

    Parameters
    ----------
    all_lordloc : ndarray
        Array containing the left trace. This can either be a 2-d array with shape (nspec, nTrace)
        for multiple traces, or simply a 1-d array with shape  (nspec) for a single trace.

    all_rordloc : ndarray
        Array containing the right trace. This can either be a 2-d array with shape (nspec, nTrace)
        for multiple traces, or simply a 1-d array with shape  (nspec) for a single trace.

    frameshape : tuple
      A two element tuple providing the shape of a trace frame.

    pad : int

    Returns
    -------
    msordloc : ndarray
      An image assigning each pixel to a slit number. A zero value indicates
      that this pixel does not belong to any slit.
    """

    # This little bit of code allows the input lord and rord to either be (nspec, nslit) arrays or a single
    # vectors of size (nspec)
    if all_lordloc_in.ndim == 2:
        nslits = all_lordloc_in.shape[1]
        all_lordloc = all_lordloc_in
        all_rordloc = all_rordloc_in
    else:
        nslits = 1
        all_lordloc = all_lordloc_in.reshape(all_lordloc_in.size,1)
        all_rordloc = all_rordloc_in.reshape(all_rordloc_in.size,1)

#    nslits = all_lordloc.shape[1]
    msordloc = np.zeros(frameshape)
    for o in range(nslits):
        lordloc = all_lordloc[:, o]
        rordloc = all_rordloc[:, o]
        ordloc = locate_order(lordloc, rordloc, frameshape[0], frameshape[1], int(pad))
        word = np.where(ordloc != 0)
        if word[0].size == 0:
            msgs.warn("There are no pixels in slit {0:d}".format(o + 1))
            continue
        msordloc[word] = o + 1
    return msordloc


def locate_order(lordloc, rordloc, sz_x, sz_y, pad):
    """
    Generate a boolean image that identifies which pixels belong to the
    slit associated with the supplied left and right slit edges.

    Parameters
    ----------
    lordloc : ndarray
      Location of the left slit edges of 1 slit
    rordloc : ndarray
      Location of the right slit edges of 1 slit
    sz_x : int
      The size of an image in the spectral (0th) dimension
    sz_y : int
      The size of an image in the spatial (1st) dimension
    pad : int
      Additional pixels to pad the left and right slit edges

    Returns
    -------
    orderloc : ndarray
      An image the same size as the input frame, containing values from 0-1.
      0 = pixel is not in the specified slit
      1 = pixel is in the specified slit
    """
    ow = (rordloc-lordloc)/2.0
    oc = (rordloc+lordloc)/2.0
    ymin = (oc-ow).astype(int)-pad
    ymax = (oc+ow).astype(int)+1+pad
    indx = np.invert((ymax < 0) | (ymin >= sz_y))
    ymin[ymin < 0] = 0
    ymax[ymax > sz_y-1] = sz_y-1
    indx &= (ymax > ymin)

    orderloc = np.zeros((sz_x,sz_y), dtype=int)
    for x in np.arange(sz_x)[indx]:
        orderloc[x,ymin[x]:ymax[x]] = 1
    return orderloc


def pix_to_amp(naxis0, naxis1, datasec, numamplifiers):
    """ Generate a frame that identifies each pixel to an amplifier,
    and then trim it to the data sections.
    This frame can be used to later identify which trimmed pixels
    correspond to which amplifier

    Parameters
    ----------
    naxis0 : int
    naxis1 : int
    datasec : list
    numamplifiers : int

    Returns
    -------
    retarr : ndarray
      Frame assigning pixels to amplifiers

    """
    # For convenience
    # Initialize the returned array
    retarr = np.zeros((naxis0, naxis1))
    for i in range(numamplifiers):
        #datasec = "datasec{0:02d}".format(i+1)
        #x0, x1 = settings.spect[dnum][datasec][0][0], settings.spect[dnum][datasec][0][1]
        #y0, y1 = settings.spect[dnum][datasec][1][0], settings.spect[dnum][datasec][1][1]
        x0, x1 = datasec[i][0][0], datasec[i][0][1]
        y0, y1 = datasec[i][1][0], datasec[i][1][1]
        if x0 < 0: x0 += naxis0
        if x1 <= 0: x1 += naxis0
        if y0 < 0: y0 += naxis1
        if y1 <= 0: y1 += naxis1
        # Fill in the pixels for this amplifier
        xv = np.arange(x0, x1)
        yv = np.arange(y0, y1)
        w = np.ix_(xv, yv)
        try:
            retarr[w] = i+1
        except IndexError:
            debugger.set_trace()
        # Save these locations for trimming
        if i == 0:
            xfin = xv.copy()
            yfin = yv.copy()
        else:
            xfin = np.unique(np.append(xfin, xv.copy()))
            yfin = np.unique(np.append(yfin, yv.copy()))
    # Construct and array with the rows and columns to be extracted
    w = np.ix_(xfin, yfin)
    return retarr[w]

# ToDO rewrite this function to use images rather than loops as in flat_fit.py
def ximg_and_edgemask(lord_in, rord_in, slitpix, trim_edg=(3,3), xshift=0.):
    """
    Generate the ximg and edgemask frames

    Parameters
    ----------
    lord_in : ndarray
        Array containing the left trace. This can either be a 2-d array with shape (nspec, nTrace)
        for multiple traces, or simply a 1-d array with shape  (nspec) for a single trace.

    rord_in : ndarray
        Array containing the right trace. This can either be a 2-d array with shape (nspec, nTrace)
        for multiple traces, or simply a 1-d array with shape  (nspec) for a single trace.

    slitpix : ndarray
      Image with shape (nspec, nspat) specifying pixel locations. This is created by core_slit_pixels above.

    trim_edg : tuple
      How much to trim off each edge of each slit

    xshift : float, optional
      Future implementation may need to shift the edges

    Returns
    -------
    ximg : ndarray
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
            debugger.set_trace()
            rord[:, islit] = lord[:, islit] + meds

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
