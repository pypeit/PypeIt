# Licensed under a 3-clause BSD style license - see PYDL_LICENSE.rst
# -*- coding: utf-8 -*-
# Also cite https://doi.org/10.5281/zenodo.1095150 when referencing PYDL
import warnings

from IPython import embed

import numpy as np

def cholesky_band(l, mininf=0.0):
    """Compute Cholesky decomposition of banded matrix.

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
    """Solve the equation Ax=b where A is a Cholesky-banded matrix.

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

