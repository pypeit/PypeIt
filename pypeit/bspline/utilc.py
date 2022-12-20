"""
Implements support methods for
:class:`pypeit.bspline.bspline.bspline`. This module specifically
imports and wrap C functions to improve efficiency.

.. include:: ../include/links.rst

"""

import os
import warnings
import ctypes

from IPython import embed

import numpy as np

# Mimics astropy convention
LIBRARY_PATH = os.path.dirname(__file__)
try:
    _bspline = np.ctypeslib.load_library("_bspline", LIBRARY_PATH)
except Exception:
    raise ImportError('Unable to load bspline C extension.  Try rebuilding pypeit.')

#-----------------------------------------------------------------------
bspline_model_c = _bspline.bspline_model
bspline_model_c.restype = None
bspline_model_c.argtypes = [np.ctypeslib.ndpointer(ctypes.c_double, flags="F_CONTIGUOUS"),
                            np.ctypeslib.ndpointer(ctypes.c_int64, flags="C_CONTIGUOUS"),
                            np.ctypeslib.ndpointer(ctypes.c_int64, flags="C_CONTIGUOUS"),
                            np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                            ctypes.c_int32, ctypes.c_int32, ctypes.c_int32, ctypes.c_int32,
                            np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]

def bspline_model(x, action, lower, upper, coeff, n, nord, npoly):
    """
    Calculate the bspline model.

    This method wraps a C function.

    Args:
        x (`numpy.ndarray`_):
            The independent variable in the fit.
        action (`numpy.ndarray`_):
            Action matrix. See
            :func:`pypeit.bspline.bspline.bspline.action`. The shape
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
        `numpy.ndarray`_: The best fitting bspline model at all
        provided :math:`x`.
    """
    # TODO: Can we save some of these objects to self so that we
    # don't have to recreate them?
    # TODO: x is always 1D right?
    yfit = np.zeros(x.size, dtype=x.dtype)
    upper = np.array(upper, dtype=np.int64)
    lower = np.array(lower, dtype=np.int64)
    # TODO: Get rid of this ascontiguousarray call if possible
#    print(action.flags['F_CONTIGUOUS'])
    bspline_model_c(action, lower, upper, coeff.flatten('F'), n, nord, npoly, x.size, yfit)
    return yfit
#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
intrv_c = _bspline.intrv
intrv_c.restype = None
intrv_c.argtypes = [ctypes.c_int32, np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                    ctypes.c_int32, np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                    ctypes.c_int32, np.ctypeslib.ndpointer(ctypes.c_int64, flags="C_CONTIGUOUS")]

def intrv(nord, breakpoints, x):
    """
    Find the segment between breakpoints which contain each value in
    the array x.

    The minimum breakpoint is nbkptord -1, and the maximum
    is nbkpt - nbkptord - 1.

    This method wraps a C function.

    Parameters
    ----------
    nord : :obj:`int`
        Order of the fit.
    breakpoints : `numpy.ndarray`_
        Locations of good breakpoints
    x : `numpy.ndarray`_
        Data values, assumed to be monotonically increasing.

    Returns
    -------
    indx : `numpy.ndarray`_
        Position of array elements with respect to breakpoints.
    """
    indx = np.zeros(x.size, dtype=np.int64)
    intrv_c(nord, breakpoints, breakpoints.size, x, x.size, indx)
    return indx
#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
solution_arrays_c = _bspline.solution_arrays
solution_arrays_c.restype = None
solution_arrays_c.argtypes = [ctypes.c_int32, ctypes.c_int32, ctypes.c_int32, ctypes.c_int32,
                              np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                              np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                              np.ctypeslib.ndpointer(ctypes.c_double, flags="F_CONTIGUOUS"),
                              np.ctypeslib.ndpointer(ctypes.c_int64, flags="C_CONTIGUOUS"),
                              np.ctypeslib.ndpointer(ctypes.c_int64, flags="C_CONTIGUOUS"),
                              np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                              ctypes.c_int32,
                              np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                              ctypes.c_int32]

def solution_arrays(nn, npoly, nord, ydata, action, ivar, upper, lower):
    """
    Support function that builds the arrays for Cholesky
    decomposition.

    This method wraps a C function.

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
        :obj:`tuple`: Returns (1) matrix :math:`A` and (2) vector :math:`b`
        prepared for Cholesky decomposition and used in the solution
        to the equation :math:`Ax=b`.
    """
    nfull = nn * npoly
    bw = npoly * nord

    alpha = np.zeros((bw, nfull+bw), dtype=float)
    beta = np.zeros((nfull+bw,), dtype=float)
    upper = np.array(upper, dtype=np.int64)
    lower = np.array(lower, dtype=np.int64)

    # need to convert action to fortran-style order and contiguous layout.
    action_t = np.asfortranarray(np.transpose(action))

    # NOTE: Beware of the integer types for upper and lower. They must
    # match the argtypes above and in bspline.c explicitly!! np.int32
    # for int and np.int64 for long.
    # NOTE: `action` *must* be stored in fortran-style, column-major contiguous format.
    solution_arrays_c(nn, npoly, nord, ydata.size, ydata, ivar, action_t,
                      upper, lower, alpha, alpha.shape[0], beta, beta.size)
    return alpha, beta
#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
cholesky_band_c = _bspline.cholesky_band
cholesky_band_c.restype = int
cholesky_band_c.argtypes = [np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                            ctypes.c_int32, ctypes.c_int32]

def cholesky_band(l, mininf=0.0, verbose=False):
    """
    Compute Cholesky decomposition of banded matrix.

    This method wraps a C function.

    Parameters
    ----------
    l : `numpy.ndarray`_
        A matrix on which to perform the Cholesky decomposition.
    mininf : :class:`float`, optional
        Entries in the `l` matrix are considered negative if they are less
        than this value (default 0.0).

    Returns
    -------
    :obj:`tuple`
        If problems were detected, the first item will be the index or
        indexes where the problem was detected, and the second item will simply
        be the input matrix.  If no problems were detected, the first item
        will be -1, and the second item will be the Cholesky decomposition.
    """
    n = np.diff(l.shape)[0]
    negative = (l[0,:n] <= mininf) | np.invert(np.isfinite(l[0,:n]))
    if np.any(negative):
        nz = negative.nonzero()[0]
        if verbose:
            warnings.warn('Found {0} bad entries: {1}'.format(nz.size, nz))
        return nz, l
    ll = l.copy()
    err = cholesky_band_c(ll, ll.shape[0], ll.shape[1])
    return err, ll if err == -1 else l
#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
cholesky_solve_c = _bspline.cholesky_solve
cholesky_solve_c.restype = None
cholesky_solve_c.argtypes = [np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                             ctypes.c_int32, ctypes.c_int32,
                             np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                             ctypes.c_int32]

def cholesky_solve(a, bb):
    r"""
    Solve the equation :math:`Ax=b` where :math:`A` is a
    Cholesky-banded matrix.

    This method wraps a C function.

    Parameters
    ----------
    a : `numpy.ndarray`_
        :math:`A` in :math:`A x = b`.
    bb : `numpy.ndarray`_
        :math:`b` in :math:`A x = b`.

    Returns
    -------
    :obj:`tuple`
        A tuple containing the status and the result of the solution.  The
        status is always -1.
    """
    b = bb.copy()
    cholesky_solve_c(a, a.shape[0], a.shape[1], b, b.shape[0])
    return -1, b
#-----------------------------------------------------------------------
