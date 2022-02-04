
import time
import os
import timeit

from IPython import embed

import numpy as np

from pypeit import bspline
from pypeit.tests.tstutils import bspline_ext_required, data_path
from pypeit.core import fitting

@bspline_ext_required
def test_model_versions_time():
    from pypeit.bspline.utilpy import bspline_model as bspline_model_py
    from pypeit.bspline.utilc import bspline_model as bspline_model_c

    command = "bspline_model(d['x'], d['action'], d['lower'], d['upper'], d['coeff'], \n" \
              "              d['n'], d['nord'], d['npoly'])"
    py_setup = "import numpy as np\n" \
               "from pypeit.tests.tstutils import data_path\n" \
               "from pypeit.bspline.utilpy import bspline_model\n" \
               "d = np.load(data_path('bspline_model.npz'))"
    c_setup = "import numpy as np\n" \
              "from pypeit.tests.tstutils import data_path\n" \
              "from pypeit.bspline.utilc import bspline_model\n" \
              "d = np.load(data_path('bspline_model.npz'))"

    pytime = min(timeit.repeat(stmt=command, setup=py_setup, number=10))
    ctime = min(timeit.repeat(stmt=command, setup=c_setup, number=10))
    assert ctime < 2*pytime, 'C is less efficient!'


@bspline_ext_required
def test_model_versions():
    from pypeit.bspline.utilpy import bspline_model as bspline_model_py
    from pypeit.bspline.utilc import bspline_model as bspline_model_c

    d = np.load(data_path('bspline_model.npz'))

    mod = bspline_model_py(d['x'], d['action'], d['lower'], d['upper'], d['coeff'], d['n'],
                           d['nord'], d['npoly'])
    _mod = bspline_model_c(d['x'], d['action'], d['lower'], d['upper'], d['coeff'], d['n'],
                           d['nord'], d['npoly'])

    assert np.allclose(mod, _mod), 'Differences in index'


@bspline_ext_required
def test_intrv_versions_time():
    command = "intrv(d['nord'], d['breakpoints'], d['x'])"

    py_setup = "import numpy as np\n" \
               "from pypeit.tests.tstutils import data_path\n" \
               "from pypeit.bspline.utilpy import intrv\n" \
               "d = np.load(data_path('intrv.npz'))"
    c_setup = "import numpy as np\n" \
              "from pypeit.tests.tstutils import data_path\n" \
              "from pypeit.bspline.utilc import intrv\n" \
              "d = np.load(data_path('intrv.npz'))"

    pytime = min(timeit.repeat(stmt=command, setup=py_setup, number=10))
    ctime = min(timeit.repeat(stmt=command, setup=c_setup, number=10))
    assert ctime < pytime, 'C is less efficient!'


@bspline_ext_required
def test_intrv_versions():
    from pypeit.bspline.utilpy import intrv as intrv_py
    from pypeit.bspline.utilc import intrv as intrv_c

    d = np.load(data_path('intrv.npz'))

    indx = intrv_py(d['nord'], d['breakpoints'], d['x'])
    _indx = intrv_c(d['nord'], d['breakpoints'], d['x'])
    assert np.allclose(indx, _indx), 'Differences in index'


@bspline_ext_required
def test_solution_array_versions_time():
    command = "solution_arrays(d['nn'], d['npoly'], d['nord'], d['ydata'], d['action'], " \
              "                d['ivar'], d['upper'], d['lower'])"

    py_setup = "import numpy as np\n" \
               "from pypeit.tests.tstutils import data_path\n" \
               "from pypeit.bspline.utilpy import solution_arrays\n" \
               "d = np.load(data_path('solution_arrays.npz'))"
    c_setup = "import numpy as np\n" \
              "from pypeit.tests.tstutils import data_path\n" \
              "from pypeit.bspline.utilc import solution_arrays\n" \
              "d = np.load(data_path('solution_arrays.npz'))"

    pytime = min(timeit.repeat(stmt=command, setup=py_setup, number=10))
    ctime = min(timeit.repeat(stmt=command, setup=c_setup, number=10))
    assert ctime < pytime, 'C is less efficient!'

@bspline_ext_required
def test_solution_array_versions():
    # Import only when the test is performed
    from pypeit.bspline.utilpy import solution_arrays as sol_py
    from pypeit.bspline.utilc import solution_arrays as sol_c

    d = np.load(data_path('solution_arrays.npz'))

    a, b = sol_py(d['nn'], d['npoly'], d['nord'], d['ydata'], d['action'], d['ivar'],
                  d['upper'], d['lower'])
    _a, _b = sol_c(d['nn'], d['npoly'], d['nord'], d['ydata'], d['action'], d['ivar'],
                   d['upper'], d['lower'])

    assert np.allclose(a, _a), 'Differences in alpha'
    assert np.allclose(b, _b), 'Differences in beta'


@bspline_ext_required
def test_cholesky_band_versions_time():
    command = "cholesky_band(d['l'], mininf=d['mininf'])"

    py_setup = "import numpy as np\n" \
               "from pypeit.tests.tstutils import data_path\n" \
               "from pypeit.bspline.utilpy import cholesky_band\n" \
               "d = np.load(data_path('cholesky_band_l.npz'))"
    c_setup = "import numpy as np\n" \
              "from pypeit.tests.tstutils import data_path\n" \
              "from pypeit.bspline.utilc import cholesky_band\n" \
              "d = np.load(data_path('cholesky_band_l.npz'))"

    pytime = min(timeit.repeat(stmt=command, setup=py_setup, number=10))
    ctime = min(timeit.repeat(stmt=command, setup=c_setup, number=10))
    assert ctime < pytime, 'C is less efficient!'


@bspline_ext_required
def test_cholesky_band_versions():
    # Import only when the test is performed
    from pypeit.bspline.utilpy import cholesky_band as cholesky_band_py
    from pypeit.bspline.utilc import cholesky_band as cholesky_band_c

    # Read data
    d = np.load(data_path('cholesky_band_l.npz'))

    # Run python version
    e, l = cholesky_band_py(d['l'], mininf=d['mininf'])
    # Run C version
    e, _l = cholesky_band_c(d['l'], mininf=d['mininf'])

    # Check output arrays
    assert np.allclose(l, _l), 'Differences in cholesky_band'


@bspline_ext_required
def test_cholesky_solve_versions_time():
    command = "cholesky_solve(d['a'], d['bb'])"

    py_setup = "import numpy as np\n" \
               "from pypeit.tests.tstutils import data_path\n" \
               "from pypeit.bspline.utilpy import cholesky_solve\n" \
               "d = np.load(data_path('cholesky_solve_abb.npz'))"
    c_setup = "import numpy as np\n" \
              "from pypeit.tests.tstutils import data_path\n" \
              "from pypeit.bspline.utilc import cholesky_solve\n" \
              "d = np.load(data_path('cholesky_solve_abb.npz'))"

    pytime = min(timeit.repeat(stmt=command, setup=py_setup, number=10))
    ctime = min(timeit.repeat(stmt=command, setup=c_setup, number=10))
    assert ctime < pytime, 'C is less efficient!'


@bspline_ext_required
def test_cholesky_solve_versions():
    # Import only when the test is performed
    from pypeit.bspline.utilpy import cholesky_solve as cholesky_solve_py
    from pypeit.bspline.utilc import cholesky_solve as cholesky_solve_c

    # Read data
    d = np.load(data_path('cholesky_solve_abb.npz'))

    # Run python version
    e, b = cholesky_solve_py(d['a'], d['bb'])
    # Run C version
    e, _b = cholesky_solve_c(d['a'], d['bb'])
    # Check output arrays
    assert np.allclose(b, _b), 'Differences in cholesky_solve'


# NOTE: Used to be in test_pydl.py.
# TODO: Where is the to/from dict functionality used?
def test_bsplinetodict():
    """
    Test for writing a bspline onto a dict (and also reading it out).
    """
    x = np.random.rand(500)

    # Create bspline
    init_bspline = bspline.bspline(x, bkspace=0.01*(np.max(x)-np.min(x)))
    # Write bspline to bspline_dict
    bspline_dict = init_bspline.to_dict()
    # Create bspline from bspline_dict
    bspline_fromdict = bspline.bspline(None, from_dict=bspline_dict)

    assert np.max(np.array(bspline_dict['breakpoints'])-bspline_fromdict.breakpoints) == 0.

def test_profile_spec():
    """
    Test that bspline_profile (1) is successful and (2) produces the
    same result for a set of data fit spectrally.
    """
    # Files created using `rmtdict` branch (30 Jan 2020)
    files = [data_path('gemini_gnirs_32_{0}_spec_fit.npz'.format(slit)) for slit in [0,1]]
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
        assert np.allclose(d['spec_flat_fit'], spec_flat_fit), 'Bad spectral bspline result'

def test_io():
    """
    Test that bspline_profile (1) is successful and (2) produces the
    same result for a set of data fit spectrally.
    """
    # Files created using `rmtdict` branch (30 Jan 2020)
    files = [data_path('gemini_gnirs_32_{0}_spec_fit.npz'.format(slit)) for slit in [0,1]]
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
    ofile = data_path('tst_bspline.fits')
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
    files = [data_path('gemini_gnirs_32_{0}_spat_fit.npz'.format(slit)) for slit in [0,1]]
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
    files = [data_path('gemini_gnirs_32_{0}_twod_fit.npz'.format(slit)) for slit in [0,1]]
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

