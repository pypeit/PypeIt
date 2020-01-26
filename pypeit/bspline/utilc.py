# Licensed under a 3-clause BSD style license - see PYDL_LICENSE.rst
# -*- coding: utf-8 -*-
# Also cite https://doi.org/10.5281/zenodo.1095150 when referencing PYDL
import warnings
import ctypes

from IPython import embed

import numpy as np

LIBRARY_PATH = os.path.dirname(__file__)
try:
    _bspline = np.ctypeslib.load_library("_bspline", LIBRARY_PATH)
except Exception:
    raise ImportError('Unable to load bspline C extension.  Try rebuilding pypeit.')

cholesky_solve_c = _bspline.cholesky_solve
cholesky_solve_c.restype = None
cholesky_solve_c.argtypes = [np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                             ctypes.c_int, ctypes.c_int, 
                             np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                             ctypes.c_int]

cholesky_band_c = _bspline.cholesky_band
cholesky_band_c.restype = int
cholesky_band_c.argtypes = [np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                            ctypes.c_int, ctypes.c_int]

def cholesky_solve(a, bb):
    b = bb.copy()
    cholesky_solve_c(a, a.shape[0], a.shape[1], b, b.shape[0])
    return -1, b

def cholesky_band(l, mininf=0.0):
    n = np.diff(l.shape)[0]
    negative = (l[0,:n] <= mininf) | np.invert(np.isfinite(l[0,:n]))
    if np.any(negative):
        nz = negative.nonzero()[0]
        warnings.warn('Found {0} bad entries: {1}'.format(nz.size, nz))
        return nz, l
    ll = l.copy()
    err = cholesky_band_c(ll, ll.shape[0], ll.shape[1])
    return err, ll if err == -1 else l
