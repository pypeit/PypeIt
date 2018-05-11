""" Routines related to mapping pixels to physical positions
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np

from pypit import msgs
from pypit import arparse as settings

try:
    from pypit import ginga
except ImportError:
    pass


def gen_pixloc(frame, det, **kwargs):
    """ Now a simple wrapper to core_gen_pixloc

    Parameters
    ----------
    frame
    det
    kwargs

    Returns
    -------

    """
    dnum = settings.get_dnum(det)
    xgap = settings.spect[dnum]['xgap']
    ygap = settings.spect[dnum]['ygap']
    ysize = settings.spect[dnum]['ysize']
    # Do it
    return core_gen_pixloc(frame, xgap=xgap, ygap=ygap, ysize=ysize, **kwargs)


def core_gen_pixloc(frame, xgap=0, ygap=0, ysize=1., gen=True):
    """
    Generate an array of physical pixel coordinates

    Parameters
    ----------
    frame : ndarray
      uniformly illuminated and normalized flat field frame
    xgap : int (optional)
    ygap : int (optional)
    ysize : float (optional)

    Returns
    -------
    locations : ndarray
      A 3D array containing the x center, y center, x width and y width of each pixel.
      The returned array has a shape:   frame.shape + (4,)
    """
    #dnum = settings.get_dnum(det)
    msgs.info("Deriving physical pixel locations on the detector")
    locations = np.zeros((frame.shape[0],frame.shape[1],4))
    if gen:
        msgs.info("Pixel gap in the dispersion direction = {0:4.3f}".format(xgap))
        msgs.info("Pixel size in the dispersion direction = {0:4.3f}".format(1.0))
        xs = np.arange(frame.shape[0]*1.0)*xgap
        xt = 0.5 + np.arange(frame.shape[0]*1.0) + xs
        msgs.info("Pixel gap in the spatial direction = {0:4.3f}".format(ygap))
        msgs.info("Pixel size in the spatial direction = {0:4.3f}".format(ysize))
        ys = np.arange(frame.shape[1])*ygap*ysize
        yt = ysize*(0.5 + np.arange(frame.shape[1]*1.0)) + ys
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
    diff = pixlocn[:,0,0] if axis == 0 else pixlocn[0,:,1]

#    print('calling phys_to_pix')
#    t = time.clock()
#    _pixarr = arcytrace.phys_to_pix(np.array([array]).T, diff).flatten() \
#                if len(np.shape(array)) == 1 else arcytrace.phys_to_pix(array, diff)
#    print('Old phys_to_pix: {0} seconds'.format(time.clock() - t))
#    t = time.clock()
    pixarr = new_phys_to_pix(array, diff)
#    print('New phys_to_pix: {0} seconds'.format(time.clock() - t))
#    assert np.sum(_pixarr != pixarr) == 0, 'Difference between old and new phys_to_pix, pixarr'

    return pixarr


def new_phys_to_pix(array, diff):
    if len(array.shape) > 2:
        msgs.error('Input array must have two dimensions or less!')
    if len(diff.shape) != 1:
        msgs.error('Input difference array must be 1D!')
    _array = np.atleast_2d(array)
    doravel = len(array.shape) != 2
    pix = np.argmin(np.absolute(_array[:,:,None] - diff[None,None,:]), axis=2).astype(int)
    return pix.ravel() if doravel else pix


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


def core_slit_pixels(all_lordloc, all_rordloc, frameshape, pad):
    """ Generate an image indicating the slit/order associated with each pixel.

    Parameters
    ----------
    all_lordloc : ndarray
    all_rordloc : ndarray
    frameshape : tuple
      A two element tuple providing the shape of a trace frame.
    pad : int

    Returns
    -------
    msordloc : ndarray
      An image assigning each pixel to a slit number. A zero value indicates
      that this pixel does not belong to any slit.
    """

    nslits = all_lordloc.shape[1]
    msordloc = np.zeros(frameshape)
    for o in range(nslits):
        lordloc = all_lordloc[:, o]
        rordloc = all_rordloc[:, o]
#        print('calling locate_order')
#        t = time.clock()
#        _ordloc = arcytrace.locate_order(lordloc, rordloc, frameshape[0], frameshape[1],
#                                         settings.argflag['trace']['slits']['pad'])
#        print('Old locate_order: {0} seconds'.format(time.clock() - t))
#        t = time.clock()
        ordloc = new_locate_order(lordloc, rordloc, frameshape[0], frameshape[1], pad)
#        print('New locate_order: {0} seconds'.format(time.clock() - t))
#        assert np.sum(_ordloc != ordloc) == 0, \
#                    'Difference between old and new locate_order'
        word = np.where(ordloc != 0)
        if word[0].size == 0:
            msgs.warn("There are no pixels in slit {0:d}".format(o + 1))
            continue
        msordloc[word] = o + 1
    return msordloc


def new_locate_order(lordloc, rordloc, sz_x, sz_y, pad):
    """ Generate a boolean image that identifies which pixels
    belong to the slit associated with the supplied left and
    right slit edges.

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

