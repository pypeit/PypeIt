"""
Module to run tests on fitting
"""
import os

from IPython import embed

import pytest

import numpy as np

from pypeit import dataPaths
from pypeit.core import fitting
from pypeit.tests.tstutils import data_output_path


def test_pypeitfit():
    out_file = data_output_path('test_fit.fits')
    if os.path.isfile(out_file):
        os.remove(out_file)
    pypeitFit = fitting.PypeItFit(fitc=np.arange(100).astype(float))
    # Write
    pypeitFit.to_file(out_file)
    # Read
    pypeitFit2 = fitting.PypeItFit.from_file(out_file)
    assert np.array_equal(pypeitFit.fitc, pypeitFit2.fitc)
    pypeitFit2.to_file(out_file, overwrite=True)
    # Finish
    os.remove(out_file)


def test_polynomial():
    """ Run the parameter setup script
    """
    x = np.pi*np.linspace(0, 1., 100)
    y = np.sin(x)
    # Polynomial
    pypeitFit = fitting.PypeItFit(xval=x, yval=y, func='polynomial', order=np.array([3]))
    pypeitFit.fit()
    #pypeitFit = fitting.func_fit(x, y, 'polynomial', 3)
    np.testing.assert_allclose(pypeitFit.fitc, np.array([ -4.74660344e-02,   1.30745471e+00,
                                                  -4.16175760e-01, 3.08557167e-18]), atol=1e-9)
    # Evaluate
    val = pypeitFit.eval(x)
    assert np.isclose(val[0], -0.04746603), 'Bad value'


def test_legendre():
    x = np.pi*np.linspace(0, 1., 100)
    y = np.sin(x)
    # Legendre
    pypeitFit = fitting.PypeItFit(xval=x, yval=y, func='legendre', order=np.array([3]))
    pypeitFit.fit() # = fitting.func_fit(x, y, 'legendre', 3)
    np.testing.assert_allclose(pypeitFit.fitc, np.array([  6.37115652e-01,   6.83317251e-17,
                                                   -6.84581686e-01, -7.59352737e-17]), atol=1e-9)


def test_pypeitfitcollection():

    # These are results generated from the old TraceSet fitter.  This just
    # checks that the results are the same.
    test_data_file = dataPaths.tests.get_file_path('trace_fit_data.npz')
    db = np.load(test_data_file)

    trace_coo = db['trace_coo']
    cen = db['cen']
    trace_fit_ivar = db['trace_fit_ivar']
    xmin = float(db['xmin'])
    xmax = float(db['xmax'])
    function = str(db['function'])
    order = int(db['order'])
    maxdev = float(db['maxdev'])
    maxiter = int(db['maxiter'])
    trace_fit = db['trace_fit']

    fit = fitting.PypeItFitCollection(
        trace_coo, cen, ivar=trace_fit_ivar, func=function, order=order, xmin=xmin, xmax=xmax,
        maxdev=maxdev, maxiter=maxiter
    )

    # NOTE: The array_equal test passed on my machine, but I changed it to
    # allclose, just to be on the safe side.
#    assert np.array_equal(fit.yfit, trace_fit), 'Fit changed!'
    assert np.allclose(fit.yfit, trace_fit), 'Fit changed!'
