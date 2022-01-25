"""
Module to run tests on arsave
"""
import os

from IPython import embed

from pypeit.par.util import make_pypeit_file
from pypeit import pypeitsetup
from pypeit.tests.tstutils import data_path


def test_initialization():
    """ Load input PypeIt file
    """
    # Generate a PYPIT file
    pypit_file = data_path('test.pypeit')
    make_pypeit_file(pypit_file, 'shane_kast_blue', [data_path('b*fits.gz')], setup_mode=True)

    # Perform the setup
    setup = pypeitsetup.PypeItSetup.from_pypeit_file(pypit_file)
    par, spectrograph, fitstbl = setup.run(sort_dir=data_path(''))

    # Test
    assert spectrograph.name == 'shane_kast_blue'
    assert len(fitstbl) == 8

    # Clean-up
    os.remove(data_path('test.calib'))
    os.remove(data_path('test.pypeit'))



