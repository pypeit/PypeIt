# Licensed under a 3-clause BSD style license - see PYDL_LICENSE.rst
# -*- coding: utf-8 -*-
# Also cite https://doi.org/10.5281/zenodo.1095150 when referencing PYDL
import warnings

from IPython import embed

import numpy as np


def bspline_model(x, action, lower, upper, coeff, n, nord, npoly):
    """
    Calculate the bspline model.

    Args:
        x (`numpy.ndarray`_):
            The independent variable in the fit.
        action (`numpy.ndarray`_):
            Action matrix. See
            :func:`pypeit.bspline.bspline.bspline.action.` The shape
            of the array is expected to be ``nd`` by ``npoly*nord``.
        lower (`numpy.ndarray`_):
            Vector with the starting indices along the second axis of
            action used to construct the model.
        upper (`numpy.ndarray`_):
            Vector with the (inclusive) ending indices along the
            second axis of action used to construct the model.
        coeff (`numpy.ndarray`_):
            The model coefficients used for each action.
        n (:obj:`int`):
            Number of unmasked measurements included in the fit.
        nord (:obj:`int`):
            Fit order.
        npoly (:obj:`int`):
            Polynomial per fit order.

    Returns:
        `numpy.ndarray`: The best fitting bspline model at all
        provided :math:`x`.
    """
    # TODO: Can we save some of these objects to self so that we
    # don't have to recreate them?
    # TODO: x is always 1D right?
    # TODO: Used for testing bspline
#    np.savez_compressed('bspline_model.npz', x=x, action=action, lower=lower, upper=upper,
#                        coeff=coeff, n=n, nord=nord, npoly=npoly)
#    raise ValueError('Entered bspline_model')
    yfit = np.zeros(x.shape, dtype=x.dtype)
    spot = np.arange(npoly * nord, dtype=int)
    nowidth = np.invert(upper+1 > lower)
    for i in range(n-nord+1):
        if nowidth[i]:
            continue
        yfit[lower[i]:upper[i]+1] = np.dot(action[lower[i]:upper[i]+1,:],
                                           coeff.flatten('F')[i*npoly+spot])
    return yfit


def intrv(nord, breakpoints, x):
    """
    Find the segment between breakpoints which contain each value in
    the array x.

    The minimum breakpoint is nbkptord -1, and the maximum
    is nbkpt - nbkptord - 1.

    Parameters
    ----------
    nord : :obj:`int`
        Order of the fit.
    breakpoints : `numpy.ndarray`_
        Locations of good breakpoints
    x : :class:`numpy.ndarray`
        Data values, assumed to be monotonically increasing.

    Returns
    -------
    :class:`numpy.ndarray`
        Position of array elements with respect to breakpoints.
    """
    # TODO: Used for testing bspline
#    np.savez_compressed('intrv.npz', nord=nord, breakpoints=breakpoints, x=x)
#    raise ValueError('Entered solution_arrays')
    n = breakpoints.size - nord
    indx = np.zeros(x.size, dtype=int)
    ileft = nord - 1
    for i in range(x.size):
        while x[i] > breakpoints[ileft+1] and ileft < n - 1:
            ileft += 1
        indx[i] = ileft
    return indx

def solution_arrays(nn, npoly, nord, ydata, action, ivar, upper, lower):
    """
    Support function that builds the arrays for Cholesky
    decomposition.

    Args: 
        nn (:obj:`int`):
            Number of good break points.
        npoly (:obj:`int`):
            Polynomial per fit order.
        nord (:obj:`int`):
            Fit order.
        ydata (`numpy.ndarray`_):
            Data to fit.
        action (`numpy.ndarray`_):
            Action matrix. See
            :func:`pypeit.bspline.bspline.bspline.action`. The shape
            of the array is expected to be ``nd`` by ``npoly*nord``.
        ivar (`numpy.ndarray`_):
            Inverse variance in the data to fit.
        upper (`numpy.ndarray`_):
            Vector with the (inclusive) ending indices along the
            second axis of action used to construct the model.
        lower (`numpy.ndarray`_):
            Vector with the starting indices along the second axis of
            action used to construct the model.

    Returns:
        tuple: Returns (1) matrix :math:`A` and (2) vector :math:`b`
        prepared for Cholesky decomposition and used in the solution
        to the equation :math:`Ax=b`.
    """
    # TODO: Used for testing bspline
#    np.savez_compressed('solution_arrays.npz', nn=nn, npoly=npoly, nord=nord, ydata=ydata,
#                        action=action, ivar=ivar, upper=upper, lower=lower)
#    raise ValueError('Entered solution_arrays')
    nfull = nn * npoly
    bw = npoly * nord
    a2 = action * np.sqrt(ivar)[:,None]

    alpha = np.zeros((bw, nfull+bw), dtype=float)
    beta = np.zeros((nfull+bw,), dtype=float)
    bi = np.concatenate([np.arange(i)+(bw-i)*(bw+1) for i in range(bw,0,-1)])
    bo = np.concatenate([np.arange(i)+(bw-i)*bw for i in range(bw,0,-1)])
    upper += 1
    nowidth = np.invert(upper > lower)
    for k in range(nn-nord+1):
        if nowidth[k]:
            continue
        itop = k*npoly
        alpha.T.flat[bo+itop*bw] \
                += np.dot(a2[lower[k]:upper[k],:].T, a2[lower[k]:upper[k],:]).flat[bi]
        beta[itop:min(itop,nfull)+bw] \
                += np.dot(ydata[lower[k]:upper[k]] * np.sqrt(ivar[lower[k]:upper[k]]),
                          a2[lower[k]:upper[k],:])
    upper -= 1
    return alpha, beta


def cholesky_band(l, mininf=0.0):
    """
    Compute Cholesky decomposition of banded matrix.

    This function is pure python.

    Parameters
    ----------
    l : :class:`numpy.ndarray`
        A matrix on which to perform the Cholesky decomposition.
    mininf : :class:`float`, optional
        Entries in the `l` matrix are considered negative if they are less
        than this value (default 0.0).

    Returns
    -------
    :func:`tuple`
        If problems were detected, the first item will be the index or
        indexes where the problem was detected, and the second item will simply
        be the input matrix.  If no problems were detected, the first item
        will be -1, and the second item will be the Cholesky decomposition.
    """
#    # TODO: Used for testing bspline
#    np.savez_compressed('cholesky_band_l.npz', l=l, mininf=mininf)
#    print(l.shape)
#    raise ValueError('Entered band')

    bw, nn = l.shape
    n = nn - bw
    negative = (l[0,:n] <= mininf) | np.invert(np.isfinite(l[0,:n]))
    # JFH changed this below to make it more consistent with IDL version. Not sure
    # why the np.all(np.isfinite(lower)) was added. The code could return an empty
    # list for negative.nonzero() and crash if all elements in lower are NaN.
    # KBW: Added the "or not finite" flags to negative.
    if negative.any():
        nz = negative.nonzero()[0]
        warnings.warn('Found {0} bad entries: {1}'.format(nz.size, nz))
        return nz, l

    lower = l.copy()
    kn = bw - 1
    spot = np.arange(kn, dtype=int) + 1
    bi = np.concatenate([np.arange(i)+(kn-i)*(kn+1) for i in range(kn,0,-1)])
    here = bi[:,None] + (np.arange(n)[None,:] + 1)*bw
    for j in range(n):
        lower[0,j] = np.sqrt(lower[0,j])
        lower[spot,j] /= lower[0,j]
        if not np.all(np.isfinite(lower[spot,j])):
            warnings.warn('NaN found in cholesky_band.')
            return j, l
        hmm = lower[spot,j,None] * lower[None,spot,j]
        lower.T.flat[here[:,j]] -= hmm.flat[bi]
    return -1, lower


def cholesky_solve(a, bb):
    """
    Solve the equation Ax=b where A is a Cholesky-banded matrix.

    This function is pure python.

    Parameters
    ----------
    a : :class:`numpy.ndarray`
        :math:`A` in :math:`A x = b`.
    bb : :class:`numpy.ndarray`
        :math:`b` in :math:`A x = b`.

    Returns
    -------
    :func:`tuple`
        A tuple containing the status and the result of the solution.  The
        status is always -1.
    """
#    # TODO: Used for testing bspline
#    np.savez_compressed('cholesky_solve_abb.npz', a=a, bb=bb)
#    print(a.shape)
#    print(bb.shape)
#    raise ValueError('Entered solve')

    b = bb.copy()
    n = b.shape[0] - a.shape[0]
    kn = a.shape[0] - 1
    spot = np.arange(kn, dtype=int) + 1
    for j in range(n):
        b[j] /= a[0,j]
        b[j+spot] -= b[j]*a[spot,j]
    for j in range(n-1, -1, -1):
        b[j] = (b[j] - np.sum(a[spot,j] * b[j+spot]))/a[0,j]
    return -1, b

