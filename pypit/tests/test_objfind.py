# Module to run tests on ararclines


import os
import numpy as np
import pytest

from astropy.table import Table

from pypit import pyputils
msgs = pyputils.get_dummy_logger()
from pypit import arutils
from pypit import artrace
from pypit import arcytrace


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_objstd():
    # Read trace example
    tbl = Table.read(data_path('trc.fits'))
    trcprof = tbl['trc'].data.astype(np.float64)
    # Call
    bgreg = trcprof.size - 10
    mad = 0.03
    objl, objr, bckl, bckr = arcytrace.find_objects(trcprof, bgreg, mad)
    assert len(objl) == 5


def test_objmin():
    """ Find objects with minima algorithm
    """
    # Read trace example
    tbl = Table.read(data_path('trc.fits'))
    trcprof = tbl['trc'].data.astype(np.float64)
    # Call
    obj = artrace.find_obj_minima(trcprof, fwhm=3.5)
    assert len(obj) == 3


def test_npeaks():
    """ Test on find_nminima
    """
    # Read trace example
    tbl = Table.read(data_path('trc.fits'))
    trc = tbl['trc'].data
    # Call
    peaks, sigmas, ledges, redges = arutils.find_nminima(-1*trc, nsmooth=3)



