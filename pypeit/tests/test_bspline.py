
import os

from IPython import embed

import pytest

import numpy as np

from pypeit import dataPaths
from pypeit import bspline
from pypeit.tests.tstutils import data_output_path
from pypeit.core import fitting


def test_profile_spec():
    """
    Test that bspline_profile (1) is successful and (2) produces the
    same result for a set of data fit spectrally.
    """
    # Files created using `rmtdict` branch (30 Jan 2020)
    files = [dataPaths.tests.get_file_path('gemini_gnirs_32_{0}_spec_fit.npz'.format(slit))
                for slit in [0,1]]
    logrej = 0.5
    spec_samp_fine = 1.2
    for f in files:
        d = np.load(f)
        spec_bspl, spec_gpm_fit, spec_flat_fit, _, exit_status \
                = fitting.bspline_profile(d['spec_coo_data'], d['spec_flat_data'],
                                  d['spec_ivar_data'],
                                  np.ones_like(d['spec_coo_data']), ingpm=d['spec_gpm_data'],
                                  nord=4, upper=logrej, lower=logrej,
                                  kwargs_bspline={'bkspace': spec_samp_fine},
                                  kwargs_reject={'groupbadpix': True, 'maxrej': 5}, quiet=True)
        assert np.allclose(d['spec_flat_fit'], spec_flat_fit), 'Bad spectral bspline result'


def test_io():
    """
    Test that bspline_profile (1) is successful and (2) produces the
    same result for a set of data fit spectrally.
    """
    # Files created using `rmtdict` branch (30 Jan 2020)
    files = [dataPaths.tests.get_file_path('gemini_gnirs_32_{0}_spec_fit.npz'.format(slit))
                for slit in [0,1]]
    logrej = 0.5
    spec_samp_fine = 1.2
    for f in files:
        d = np.load(f)
        spec_bspl, spec_gpm_fit, spec_flat_fit, _, exit_status \
            = fitting.bspline_profile(d['spec_coo_data'], d['spec_flat_data'], d['spec_ivar_data'],
                              np.ones_like(d['spec_coo_data']), ingpm=d['spec_gpm_data'],
                              nord=4, upper=logrej, lower=logrej,
                              kwargs_bspline={'bkspace': spec_samp_fine},
                              kwargs_reject={'groupbadpix': True, 'maxrej': 5}, quiet=True)
    # Write
    ofile = data_output_path('tst_bspline.fits')
    if os.path.isfile(ofile):
        os.remove(ofile)
    spec_bspl.to_file(ofile)
    # Read
    _spec_bspl = bspline.bspline.from_file(ofile)
    # Check that the data read in is the same
    assert np.array_equal(_spec_bspl.breakpoints, spec_bspl.breakpoints), 'Bad read'
    # Evaluate the bsplines and check that the written and read in data
    # provide the same result
    tmp = spec_bspl.value(np.linspace(0., 1., 100))[0]
    _tmp = _spec_bspl.value(np.linspace(0., 1., 100))[0]
    assert np.array_equal(tmp, _tmp), 'Bad bspline evaluate'

    # Test overwrite
    _spec_bspl.to_file(ofile, overwrite=True)

    # None
    bspl = bspline.bspline(None)
    bspl.to_file(ofile, overwrite=True)

    os.remove(ofile)


def test_profile_spat():
    """
    Test that bspline_profile (1) is successful and (2) produces the
    same result for a set of data fit spatially.
    """
    # Files created using `rmtdict` branch (30 Jan 2020)
    files = [dataPaths.tests.get_file_path('gemini_gnirs_32_{0}_spat_fit.npz'.format(slit))
                for slit in [0,1]]
    for f in files:
        d = np.load(f)
        spat_bspl = bspline.bspline(d['spat_coo_data'], nord=4,
                                    bkspace=np.fmax(1.0/d['median_slit_width']/10.0,
                                                    1.2*np.median(np.diff(d['spat_coo_data']))))
        spat_bspl, spat_gpm_fit, spat_flat_fit, _, exit_status \
                = fitting.bspline_profile(d['spat_coo_data'], d['spat_flat_data'],
                                  np.ones_like(d['spat_flat_data']),
                                  np.ones_like(d['spat_flat_data']), nord=4, upper=5.0, lower=5.0,
                                  fullbkpt=spat_bspl.breakpoints, quiet=True)
        assert np.allclose(d['spat_flat_fit'], spat_flat_fit), 'Bad spatial bspline result'


def test_profile_twod():
    """
    Test that bspline_profile (1) is successful and (2) produces the
    same result for a set of data fit two-dimensionally.
    """
    # Files created using `rmtdict` branch (30 Jan 2020)
    files = [dataPaths.tests.get_file_path('gemini_gnirs_32_{0}_twod_fit.npz'.format(slit))
                for slit in [0,1]]
    spec_samp_coarse = 50.0
    twod_sigrej = 4.0
    for f in files:
        d = np.load(f)
        twod_bspl, twod_gpm_fit, twod_flat_fit, _ , exit_status \
                = fitting.bspline_profile(d['twod_spec_coo_data'], d['twod_flat_data'],
                                  d['twod_ivar_data'], d['poly_basis'], ingpm=d['twod_gpm_data'],
                                  nord=4, upper=twod_sigrej, lower=twod_sigrej,
                                  kwargs_bspline={'bkspace': spec_samp_coarse},
                                  kwargs_reject={'groupbadpix': True, 'maxrej': 10}, quiet=True)
        assert np.allclose(d['twod_flat_fit'], twod_flat_fit), 'Bad 2D bspline result'


def test_breakpoints():
    x = np.arange(1000, dtype=float)

    # Should fail because there's insufficient information to construct the breakpoints
    with pytest.raises(ValueError):
        bspline.bspline(x)

    # Provide the breakpoints directly
    #   - Only one break-point is provided
    bspl = bspline.bspline(x, bkpt=np.array([500.]))
    assert bspl.breakpoints.size == 8, 'Incorrect number of breakpoints.'
    assert np.isclose(bspl.breakpoints[bspl.nord-1], 0.), \
            'First breakpoint should be at the first x value'
    assert np.isclose(bspl.breakpoints[-bspl.nord], 999.), \
            'Last breakpoint should be at the last x value'
    #   - Two breakpoints are provided but they're not at the ends of the x vector
    _bspl = bspline.bspline(x, bkpt=np.array([250., 750.]))
    assert np.allclose(_bspl.breakpoints, bspl.breakpoints), 'Breakpoints should be the same'
    #   - Two breakpoints are provided at the end of the x vector
    _bspl = bspline.bspline(x, bkpt=x[[0,-1]])
    assert np.allclose(_bspl.breakpoints, bspl.breakpoints), 'Breakpoints should be the same'
    #   - Increase the spread between the inserted points
    _bspl = bspline.bspline(x, bkpt=x[[0,-1]], bkspread=2.)
    assert np.array_equal(
        np.diff(_bspl.breakpoints)[[0,-1]], 2*np.diff(bspl.breakpoints)[[0,-1]]
    ), 'Incorrect increase in the separation between breakpoints'
    #   - Give a long vector of breakpoints
    bkpt = np.arange(2,x.size,10)
    bspl = bspline.bspline(x, bkpt=bkpt)
    assert bspl.breakpoints.size == bkpt.size + 2*(bspl.nord-1), 'Bad number of breakpoints'
    assert bspl.breakpoints[bspl.nord] == bkpt[1], 'Bad use of input'
    assert bspl.breakpoints[-bspl.nord-1] == bkpt[-2], 'Bad use of input'
    assert bspl.breakpoints[bspl.nord-1] == 0., 'Should move first breakpoint to first value of x'
    assert bspl.breakpoints[-bspl.nord] == 999., 'Should move last breakpoint to last value of x'

    # Provide the breakpoint spacing
    #   - If spacing is too large, just like requesting 1 breakpoint
    bspl = bspline.bspline(x, bkspace=1000.)
    _bspl = bspline.bspline(x, nbkpts=1)
    assert np.array_equal(_bspl.breakpoints, bspl.breakpoints), 'Breakpoints should be the same'
    #   - Set a step of 10
    bspl = bspline.bspline(x, bkspace=10.)
    assert bspl.breakpoints.size == 106, 'Bad number of breakpoints'
    assert np.isclose(bspl.breakpoints[bspl.nord-1], 0.), \
            'First breakpoint should be at the first x value'
    assert np.isclose(bspl.breakpoints[-bspl.nord], 999.), \
            'Last breakpoint should be at the last x value'

    # Provide the number of breakpoints
    _bspl = bspline.bspline(x, nbkpts=100)
    assert np.array_equal(bspl.breakpoints, _bspl.breakpoints), 'Should yield the same breakpoints'
    # Setting the number of breakpoints to 1 is the same as providing on
    # breakpoint directly
    bspl = bspline.bspline(x, bkpt=np.array([500.]))
    _bspl = bspline.bspline(x, nbkpts=1)
    assert np.array_equal(bspl.breakpoints, _bspl.breakpoints), 'Should yield the same breakpoints'

    # Select every Nth x value
    #   - If step is too large, just like requesting 1 breakpoint
    bspl = bspline.bspline(x, everyn=x.size+1)
    _bspl = bspline.bspline(x, nbkpts=1)
    assert np.array_equal(bspl.breakpoints, _bspl.breakpoints), 'Should yield the same breakpoints'
    #   - With a more meaningful value
    bspl = bspline.bspline(x, everyn=10)
    # NOTE: This does NOT result in the same set of breakpoints you would get if
    # you set:
    #   bspl = bspline.bspline(x, bkspace=10.)
    assert bspl.breakpoints.size == 106, 'Bad number of breakpoints'
    assert np.isclose(bspl.breakpoints[bspl.nord-1], 0.), \
            'First breakpoint should be at the first x value'
    assert np.isclose(bspl.breakpoints[-bspl.nord], 999.), \
            'Last breakpoint should be at the last x value'
    assert np.allclose(np.diff(bspl.breakpoints[bspl.nord-1:-bspl.nord-1]), 10.), \
            'Step size is wrong'
    #   - For an irregular grid
    rng = np.random.default_rng(99)
    _x = np.sort(x + 5*rng.normal(size=x.size))
    bspl = bspline.bspline(_x, everyn=10)
    assert np.std(np.diff(bspl.breakpoints) - 10) > 1., 'Breakpoints are not irregularly spaced'

    #   - With a small, floating-point everyn value
    bspl = bspline.bspline(x, everyn=1.5)
    assert np.all(np.diff(bspl.breakpoints) == 1.5), 'Bad float result'
