# Module to run tests on ararclines
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os

import numpy as np
import pytest

from pypeit import utils
from pypeit import msgs

#def data_path(filename):
#    data_dir = os.path.join(os.path.dirname(__file__), 'files')
#    return os.path.join(data_dir, filename)


def test_func_fit():
    """ Run the parameter setup script
    """
    x = np.pi*np.linspace(0, 1., 100)
    y = np.sin(x)
    # Polynomial
    pcoeff = utils.func_fit(x, y, 'polynomial', 3)
    np.testing.assert_allclose(pcoeff, np.array([ -4.74660344e-02,   1.30745471e+00,
                                                 -4.16175760e-01, 3.08557167e-18]), atol=1e-9)
    # Legendre
    lcoeff = utils.func_fit(x, y, 'legendre', 3)
    np.testing.assert_allclose(lcoeff, np.array([  6.37115652e-01,   6.83317251e-17,
                                                   -6.84581686e-01, -7.59352737e-17]), atol=1e-9)
    # bspline
    bcoeff = utils.func_fit(x, y, 'bspline', 2)
    np.testing.assert_allclose(bcoeff[0], [ 0.        ,  0.        ,  0.        ,  0.31733259,  0.95199777,
        1.58666296,  2.22132814,  3.14159265,  3.14159265,  3.14159265], atol=1e-5)


def test_calc_ivar():
    """ Run the parameter setup script
    """
    x = np.array([-1.0, -0.1, 0.0, 0.1, 1.0])
    res = utils.calc_ivar(x)
    assert np.array_equal(res, np.array([0.0, 0.0, 0.0, 10.0, 1.0]))
    assert np.array_equal(utils.calc_ivar(res), np.array([0.0, 0.0, 0.0, 0.1, 1.0]))
