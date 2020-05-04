"""
Module to run tests on scripts
"""
import os
import pytest
import shutil

import matplotlib
matplotlib.use('agg')  # For Travis

from astropy.io import fits

from pypeit.scripts import chk_calibs
from pypeit.tests.tstutils import dev_suite_required, cooked_required


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

@dev_suite_required
def test_chk_calibs_not():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/not_alfosc/grism4')
    droot += '/ALD'
    pargs = chk_calibs.parser(['-r', droot, '-s', 'not_alfosc'])
    chk_calibs.main(pargs)

    pytest.set_trace()
    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'not_alfosc*'))
    ext = [f.split('.')[-1] for f in files]
    assert np.all([e in ext for e in expected]), \
        'Did not find all setup file extensions: {0}'.format(expected)

    # Build a PypeIt file
    pargs = setup.parser(['-r', droot, '-s', 'not_alfosc', '-c', 'A', '-d', data_path('')])
    setup.main(pargs)
    pypeit_file = data_path('not_alfosc_A/not_alfosc_A.pypeit')
    pypeIt = pypeit.PypeIt(pypeit_file, calib_only=True)

    # Clean-up
    shutil.rmtree(setup_dir)
    shutil.rmtree(data_path('not_alfosc_A'))

# TODO: Add other instruments!


