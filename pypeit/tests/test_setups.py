# Module to run tests on scripts
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import glob
import shutil

import numpy as np

import pytest
import yaml

from pypeit.pypmsgs import PypeItError
from pypeit.par.util import parse_pypeit_file
from pypeit.scripts import setup
from pypeit.tests.tstutils import dev_suite_required


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_run_setup():
    """ Test the setup script
    """
    # Remove .setup if needed
    sfiles = glob.glob('*.setups')
    for sfile in sfiles:
        os.remove(sfile)
    #
    droot = data_path('b')
    pargs = setup.parser([droot, 'shane_kast_blue', '-c',
                          '--extension=fits.gz', '--setups_path={:s}'.format(data_path(''))])
    setup.main(pargs)

    setup_file = glob.glob(data_path('setup_files/shane_kast_blue*.setups'))[0]
    # Load
    with open(setup_file, 'r') as infile:
        setup_dict = yaml.load(infile)
    # Test
    assert '01' in setup_dict['A'].keys()
    assert setup_dict['A']['--']['disperser']['name'] == '600/4310'
    # Failures
    pargs2 = setup.parser([droot, 'shane_kast_blu', '-c',
                              '--extension=fits.gz', '--setups_path={:s}'.format(data_path(''))])
    with pytest.raises(ValueError):
        setup.main(pargs2)


def test_setup_made_pypeit_file():
    """ Test the .pypeit file(s) made by pypeit_setup
    This test depends on the one above
    """
    pypeit_file = data_path('shane_kast_blue_setup_A/shane_kast_blue_setup_A.pypeit')
    cfg_lines, data_files, frametype, setups = parse_pypeit_file(pypeit_file)
    print(setups)
    # Test
    assert len(data_files) == 2
    assert frametype['b1.fits.gz'] == 'arc'
    assert setups[0] == 'A'


@dev_suite_required
def test_setup_keck_lris_red():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/Keck_LRIS_red/multi_400_8500_d560')
    droot += '/'
    pargs = setup.parser([droot, 'keck_lris_red'])
    setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'keck_lris_red*'))
    ext = [f.split('.')[-1] for f in files]
    expected = ['lst', 'pypeit', 'setups', 'sorted']
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)


@dev_suite_required
def test_setup_keck_lris_blue():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/Keck_LRIS_blue/multi_600_4000_d560')
    droot += '/'
    pargs = setup.parser([droot, 'keck_lris_blue'])
    setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'keck_lris_blue*'))
    ext = [f.split('.')[-1] for f in files]
    expected = ['lst', 'pypeit', 'setups', 'sorted']
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)

@dev_suite_required
def test_setup_shane_kast_blue():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/Shane_Kast_blue/600_4310_d55')
    droot += '/'
    pargs = setup.parser([droot, 'shane_kast_blue'])
    setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'shane_kast_blue*'))
    ext = [f.split('.')[-1] for f in files]
    expected = ['lst', 'pypeit', 'setups', 'sorted']
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)

@dev_suite_required
def test_setup_shane_kast_red():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/Shane_Kast_red/600_7500_d55')
    droot += '/'
    pargs = setup.parser([droot, 'shane_kast_red'])
    setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'shane_kast_red*'))
    ext = [f.split('.')[-1] for f in files]
    expected = ['lst', 'pypeit', 'setups', 'sorted']
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)

# TODO: We need a test data set for shane_kast_red_ret

@dev_suite_required
def test_setup_keck_deimos():

    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/Keck_DEIMOS/830G_M_8600')
    droot += '/'
    pargs = setup.parser([droot, 'keck_deimos'])
    setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'keck_deimos*'))
    ext = [f.split('.')[-1] for f in files]
    expected = ['lst', 'pypeit', 'setups', 'sorted']
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)

@dev_suite_required
def test_setup_keck_nires():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/Keck_NIRES')
    droot += '/'
    pargs = setup.parser([droot, 'keck_nires'])
    # TODO: There currently is no Echelle PypeIt class, so this should
    # fail
    with pytest.raises(PypeItError):
        setup.main(pargs)

#    cwd = os.getcwd()
#    setup_dir = os.path.join(cwd, 'setup_files')
#    assert os.path.isdir(setup_dir), 'No setup_files directory created'
#
#    files = glob.glob(os.path.join(setup_dir, 'keck_nires*'))
#    ext = [f.split('.')[-1] for f in files]
#    expected = ['lst', 'pypeit', 'setups', 'sorted']
#    assert np.all([e in ext for e in expected]), \
#            'Did not find all setup file extensions: {0}'.format(expected)
#
#    # Clean-up
#    shutil.rmtree(setup_dir)

@dev_suite_required
def test_setup_keck_nirspec():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/Keck_NIRSPEC/NIRSPEC-1')
    droot += '/'
    pargs = setup.parser([droot, 'keck_nirspec'])
    setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'keck_nirspec*'))
    ext = [f.split('.')[-1] for f in files]
    expected = ['lst', 'pypeit', 'setups', 'sorted']
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)

@dev_suite_required
def test_setup_wht_isis_blue():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/WHT_ISIS_blue/long_R300B_d5300')
    droot += '/'
    pargs = setup.parser([droot, 'wht_isis_blue', '--extension', '.fit'])
    setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'wht_isis_blue*'))
    ext = [f.split('.')[-1] for f in files]
    expected = ['lst', 'pypeit', 'setups', 'sorted']
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)

@dev_suite_required
def test_setup_vlt_xshooter_uvb():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/VLT_XSHOOTER/UVB_1x1')
    droot += '/XSHO'
    pargs = setup.parser([droot, 'vlt_xshooter_uvb'])
    # TODO: There currently is no Echelle PypeIt class, so this should
    # fail
    with pytest.raises(PypeItError):
        setup.main(pargs)

#    cwd = os.getcwd()
#    setup_dir = os.path.join(cwd, 'setup_files')
#    assert os.path.isdir(setup_dir), 'No setup_files directory created'
#
#    files = glob.glob(os.path.join(setup_dir, 'vlt_xshooter_uvb*'))
#    ext = [f.split('.')[-1] for f in files]
#    expected = ['lst', 'pypeit', 'setups', 'sorted']
#    assert np.all([e in ext for e in expected]), \
#            'Did not find all setup file extensions: {0}'.format(expected)
#
#    # Clean-up
#    shutil.rmtree(setup_dir)

@dev_suite_required
def test_setup_vlt_xshooter_vis():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/VLT_XSHOOTER/VIS_1x1')
    droot += '/XSHO'
    pargs = setup.parser([droot, 'vlt_xshooter_vis'])
    # TODO: There currently is no Echelle PypeIt class, so this should
    # fail
    with pytest.raises(PypeItError):
        setup.main(pargs)

#    cwd = os.getcwd()
#    setup_dir = os.path.join(cwd, 'setup_files')
#    assert os.path.isdir(setup_dir), 'No setup_files directory created'
#
#    files = glob.glob(os.path.join(setup_dir, 'vlt_xshooter_vis*'))
#    ext = [f.split('.')[-1] for f in files]
#    expected = ['lst', 'pypeit', 'setups', 'sorted']
#    assert np.all([e in ext for e in expected]), \
#            'Did not find all setup file extensions: {0}'.format(expected)
#
#    # Clean-up
#    shutil.rmtree(setup_dir)

@dev_suite_required
def test_setup_vlt_xshooter_nir():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/VLT_XSHOOTER/NIR')
    droot += '/XSHO'
    pargs = setup.parser([droot, 'vlt_xshooter_nir'])
    # TODO: There currently is no Echelle PypeIt class, so this should
    # fail
    with pytest.raises(PypeItError):
        setup.main(pargs)

#    cwd = os.getcwd()
#    setup_dir = os.path.join(cwd, 'setup_files')
#    assert os.path.isdir(setup_dir), 'No setup_files directory created'
#
#    files = glob.glob(os.path.join(setup_dir, 'vlt_xshooter_nir*'))
#    ext = [f.split('.')[-1] for f in files]
#    expected = ['lst', 'pypeit', 'setups', 'sorted']
#    assert np.all([e in ext for e in expected]), \
#            'Did not find all setup file extensions: {0}'.format(expected)
#
#    # Clean-up
#    shutil.rmtree(setup_dir)

if __name__=='__main__':
    test_setup_shane_kast_red()

