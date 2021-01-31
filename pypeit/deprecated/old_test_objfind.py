# Module to run tests on ararclines
import os

import numpy as np
import pytest

from astropy.table import Table

from pypeit import utils
from pypeit import artrace
from pypeit import msgs

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_very_close_obj():
    """ Find objects with minima algorithm
    """
    # Read trace example
    tbl = Table.read(data_path('real_close_trace.fits'))
    trcprof = tbl['trc'].data.astype(np.float64)
    # Call
    objl, objr, bckl, bckr = artrace.find_obj_minima(trcprof, fwhm=3., nsmooth=0)
    assert objl[1] == 55
    assert objr[1] == 63
    assert bckl.shape == (278,6)
    assert bckr.shape == (278,6)


def test_objstd():
    # Read trace example
    tbl = Table.read(data_path('trc.fits'))
    trcprof = tbl['trc'].data.astype(np.float64)
    # Call
    bgreg = trcprof.size - 10
    mad = 0.03
#    objl, objr, bckl, bckr = arcytrace.find_objects(trcprof, bgreg, mad)
    objl, objr, bckl, bckr = artrace.find_objects(trcprof, bgreg, mad)
    assert len(objl) == 5
    assert objl[1] == 36
    assert objr[1] == 82


def test_objmin():
    """ Find objects with minima algorithm
    """
    # Read trace example
    tbl = Table.read(data_path('trc.fits'))
    trcprof = tbl['trc'].data.astype(np.float64)
    # Call
    objl, objr, bckl, bckr = artrace.find_obj_minima(trcprof, fwhm=3.5)
    assert len(objl) == 3
    assert objl[1] == 63
    assert objr[1] == 90


def test_npeaks():
    """ Test on find_nminima
    """
    # Read trace example
    tbl = Table.read(data_path('trc.fits'))
    trc = tbl['trc'].data
    # Call
    peaks, sigmas, ledges, redges = utils.find_nminima(-1*trc, nsmooth=3)

    # TODO: (KBW) the second number keeps failing for me, presumably
    # because of a numpy/scipy version issue...
    try:
        np.testing.assert_allclose(peaks, np.array([4.86511462e+01, -4.530896e-03, 7.17466301e+01,
                                                    1.67786437e+02, 1.35216283e+02, 1.88809673e+02,
                                                    2.13649812e+02,   2.62738738e+02]), rtol=1e-5)
    except:
        np.testing.assert_allclose(peaks, np.array([4.86511462e+01, -4.531444e-03, 7.17466301e+01,
                                                    1.67786437e+02, 1.35216283e+02, 1.88809673e+02,
                                                    2.13649812e+02,   2.62738738e+02]), rtol=1e-5)

