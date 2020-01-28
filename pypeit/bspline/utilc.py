# Licensed under a 3-clause BSD style license - see PYDL_LICENSE.rst
# -*- coding: utf-8 -*-
# Also cite https://doi.org/10.5281/zenodo.1095150 when referencing PYDL
import os
import warnings
import ctypes

from IPython import embed

import numpy as np

LIBRARY_PATH = os.path.dirname(__file__)
try:
    _bspline = np.ctypeslib.load_library("_bspline", LIBRARY_PATH)
except Exception:
    raise ImportError('Unable to load bspline C extension.  Try rebuilding pypeit.')

bspline_model_c = _bspline.bspline_model
bspline_model_c.restype = None
bspline_model_c.argtypes = [np.ctypeslib.ndpointer(ctypes.c_double, flags="F_CONTIGUOUS"),
                            np.ctypeslib.ndpointer(ctypes.c_long, flags="C_CONTIGUOUS"),
                            np.ctypeslib.ndpointer(ctypes.c_long, flags="C_CONTIGUOUS"),
                            np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                            ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,
                            np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]

def bspline_model(x, action, lower, upper, coeff, n, nord, npoly):
    # TODO: Can we save some of these objects to self so that we
    # don't have to recreate them?
    # TODO: x is always 1D right?
    yfit = np.zeros(x.size, dtype=x.dtype)
    # TODO: Get rid of this ascontiguousarray call if possible
#    print(action.flags['F_CONTIGUOUS'])
    bspline_model_c(action, lower, upper, coeff.flatten('F'), n, nord, npoly, x.size, yfit)
    return yfit


intrv_c = _bspline.intrv
intrv_c.restype = None
intrv_c.argtypes = [ctypes.c_int, np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                    ctypes.c_int, np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                    ctypes.c_int, np.ctypeslib.ndpointer(ctypes.c_long, flags="C_CONTIGUOUS")]

def intrv(nord, breakpoints, x):
    indx = np.zeros(x.size, dtype=int)
    intrv_c(nord, breakpoints, breakpoints.size, x, x.size, indx)
    return indx


solution_arrays_c = _bspline.solution_arrays
solution_arrays_c.restype = None
solution_arrays_c.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,
                              np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                              np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                              np.ctypeslib.ndpointer(ctypes.c_double, flags="F_CONTIGUOUS"),
                              np.ctypeslib.ndpointer(ctypes.c_long, flags="C_CONTIGUOUS"),
                              np.ctypeslib.ndpointer(ctypes.c_long, flags="C_CONTIGUOUS"),
                              np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                              ctypes.c_int,
                              np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                              ctypes.c_int]

def solution_arrays(nn, npoly, nord, ydata, action, ivar, upper, lower):
    nfull = nn * npoly
    bw = npoly * nord
    # NOTE: Declared as empty because the c code zeros them out
    alpha = np.empty((bw, nfull+bw), dtype=float)
    beta = np.empty((nfull+bw,), dtype=float)
    # NOTE: Beware of the integer types for upper and lower. They must
    # match the argtypes above and in bspline.c explicitly!! np.int32
    # for int and np.int64 for long.
    # NOTE: `action` *must* be stored in fortran-style, column-major contiguous format.
    solution_arrays_c(nn, npoly, nord, ydata.size, ydata, ivar, action,
                      #np.ascontiguousarray(action),
                      upper, lower, alpha, alpha.shape[0], beta, beta.size)
    return alpha, beta


cholesky_band_c = _bspline.cholesky_band
cholesky_band_c.restype = int
cholesky_band_c.argtypes = [np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                            ctypes.c_int, ctypes.c_int]

def cholesky_band(l, mininf=0.0):
    """
    Compute Cholesky decomposition of banded matrix.

    This function is a wrapper for a C function that improves
    performance.

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
    n = np.diff(l.shape)[0]
    negative = (l[0,:n] <= mininf) | np.invert(np.isfinite(l[0,:n]))
    if np.any(negative):
        nz = negative.nonzero()[0]
        warnings.warn('Found {0} bad entries: {1}'.format(nz.size, nz))
        return nz, l
    ll = l.copy()
    err = cholesky_band_c(ll, ll.shape[0], ll.shape[1])
    return err, ll if err == -1 else l


cholesky_solve_c = _bspline.cholesky_solve
cholesky_solve_c.restype = None
cholesky_solve_c.argtypes = [np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                             ctypes.c_int, ctypes.c_int, 
                             np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                             ctypes.c_int]

def cholesky_solve(a, bb):
    r"""
    Solve the equation :math:`Ax=b` where :math:`A` is a
    Cholesky-banded matrix.

    This function is a wrapper for a C function that improves
    performance.

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
    b = bb.copy()
    cholesky_solve_c(a, a.shape[0], a.shape[1], b, b.shape[0])
    return -1, b
