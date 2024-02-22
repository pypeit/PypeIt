""" Routines related to mapping pixels to physical positions
"""
import numpy as np

from pypeit import msgs
from pypeit import debugger

from IPython import embed



def tslits2mask(tslits_dict, pad=None):
    """ Generate an image indicating the slit/order associated with each pixel.

    Parameters
    ----------
    slit_left : ndarray
        Array containing the left trace. This can either be a 2-d array with shape (nspec, nTrace)
        for multiple traces, or simply a 1-d array with shape  (nspec) for a single trace.

    slit_righ : ndarray
        Array containing the right trace. This can either be a 2-d array with shape (nspec, nTrace)
        for multiple traces, or simply a 1-d array with shape  (nspec) for a single trace.

    nspat : tuple
      Spatial dimension of the frames

    pad : int or float
      Pad the mask in both dimensions by this amount.

    Returns
    -------
    slitmask : ndarray int
      An image assigning each pixel to a slit number. A value of -1 indicates
      that this pixel does not belong to any slit.
    """

    # This little bit of code allows the input lord and rord to either be (nspec, nslit) arrays or a single
    # vectors of size (nspec)

    slit_left = tslits_dict['slit_left']
    slit_righ = tslits_dict['slit_righ']
    nslits = tslits_dict['nslits']
    nspec = tslits_dict['nspec']
    nspat = tslits_dict['nspat']
    spec_min = tslits_dict['spec_min']
    if spec_min is None:
        spec_min = np.zeros(nslits, dtype=int)
    spec_max = tslits_dict['spec_max']
    if spec_max is None:
        spec_max = np.full(nslits, nspec-1, dtype=int)
    if pad is None:
        pad = tslits_dict['pad']

    slitmask = np.full((nspec, nspat),-1,dtype=int)
    spat_img, spec_img = np.meshgrid(np.arange(nspat), np.arange(nspec))

    for islit in range(nslits):
        left_trace_img = np.outer(slit_left[:,islit], np.ones(nspat))  # left slit boundary replicated spatially
        righ_trace_img = np.outer(slit_righ[:,islit], np.ones(nspat))  # left slit boundary replicated spatially
        thismask = (spat_img > (left_trace_img - pad)) & (spat_img < (righ_trace_img + pad)) & \
                   (spec_img >= spec_min[islit]) & (spec_img <= spec_max[islit])
        if not np.any(thismask):
            msgs.warn("There are no pixels in slit {:d}".format(islit))
            continue
        slitmask[thismask] = islit

    return slitmask


