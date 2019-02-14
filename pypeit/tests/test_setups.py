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


expected = ['pypeit', 'sorted']

def test_run_setup():
    """ Test the setup script
    """
    # Remove .setup if needed
    sfiles = glob.glob('*.setups')
    for sfile in sfiles:
        os.remove(sfile)
    #
    droot = data_path('b')
    pargs = setup.parser(['-r', droot, '-s', 'shane_kast_blue', '-c=all',
                          '--extension=fits.gz', '--output_path={:s}'.format(data_path(''))])
    setup.main(pargs)

    '''
    setup_file = glob.glob(data_path('setup_files/shane_kast_blue*.setups'))[0]
    # Load
    with open(setup_file, 'r') as infile:
        setup_dict = yaml.load(infile)
    # Test
    assert '01' in setup_dict['A'].keys()
    assert setup_dict['A']['--']['disperser']['name'] == '600/4310'
    '''
    # Failures
    pargs2 = setup.parser(['-r', droot, '-s', 'shane_kast_blu', '-c=all',
                              '--extension=fits.gz', '--output_path={:s}'.format(data_path(''))])
    with pytest.raises(ValueError):
        setup.main(pargs2)


def test_setup_made_pypeit_file():
    """ Test the .pypeit file(s) made by pypeit_setup
    This test depends on the one above
    """
    pypeit_file = data_path('shane_kast_blue_A/shane_kast_blue_A.pypeit')
    cfg_lines, data_files, frametype, usrdata, setups = parse_pypeit_file(pypeit_file)
    # Test
    assert len(data_files) == 2
    assert frametype['b1.fits.gz'] == 'arc'
    assert setups[0] == 'A'


@dev_suite_required
def test_setup_keck_lris_red():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/Keck_LRIS_red/multi_400_8500_d560')
    droot += '/'
    pargs = setup.parser(['-r', droot, '-s', 'keck_lris_red'])
    setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'keck_lris_red*'))
    ext = [f.split('.')[-1] for f in files]
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)
    #os.remove('keck_lris_red.pypeit')


@dev_suite_required
def test_setup_keck_lris_blue():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/Keck_LRIS_blue/multi_600_4000_d560')
    droot += '/'
    pargs = setup.parser(['-r', droot, '-s', 'keck_lris_blue'])
    setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'keck_lris_blue*'))
    ext = [f.split('.')[-1] for f in files]
    #expected = ['lst', 'pypeit', 'setups', 'sorted']
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)
    #os.remove('keck_lris_blue.pypeit')

@dev_suite_required
def test_setup_shane_kast_blue():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/Shane_Kast_blue/600_4310_d55')
    droot += '/'
    pargs = setup.parser(['-r', droot, '-s', 'shane_kast_blue'])
    setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'shane_kast_blue*'))
    ext = [f.split('.')[-1] for f in files]
    #expected = ['lst', 'pypeit', 'setups', 'sorted']
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)
    #os.remove('shane_kast_blue.pypeit')

@dev_suite_required
def test_setup_shane_kast_red():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/Shane_Kast_red/600_7500_d55')
    droot += '/'
    pargs = setup.parser(['-r', droot, '-s', 'shane_kast_red'])
    setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'shane_kast_red*'))
    ext = [f.split('.')[-1] for f in files]
    #expected = ['lst', 'pypeit', 'setups', 'sorted']
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)
    #os.remove('shane_kast_red.pypeit')

# TODO: We need a test data set for shane_kast_red_ret

@dev_suite_required
def test_setup_keck_deimos():

    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/Keck_DEIMOS/830G_M_8600')
    droot += '/'
    pargs = setup.parser(['-r', droot, '-s', 'keck_deimos'])
    setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'keck_deimos*'))
    ext = [f.split('.')[-1] for f in files]
    #expected = ['lst', 'pypeit', 'setups', 'sorted']
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)
    #os.remove('keck_deimos.pypeit')

@dev_suite_required
def test_setup_keck_nires():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/Keck_NIRES/NIRES/')
    droot += '/'
    pargs = setup.parser(['-r', droot, '-s', 'keck_nires'])
    setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'keck_nires*'))
    ext = [f.split('.')[-1] for f in files]
    #expected = ['lst', 'pypeit', 'setups', 'sorted']
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)
    #os.remove('keck_nires.pypeit')

@dev_suite_required
def test_setup_keck_nirspec():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/Keck_NIRSPEC/NIRSPEC-1')
    droot += '/'
    pargs = setup.parser(['-r', droot, '-s', 'keck_nirspec_low'])
    setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'keck_nirspec*'))
    ext = [f.split('.')[-1] for f in files]
    #expected = ['lst', 'pypeit', 'setups', 'sorted']
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)
    #os.remove('keck_nirspec_low.pypeit')

@dev_suite_required
def test_setup_wht_isis_blue():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/WHT_ISIS_blue/long_R300B_d5300')
    droot += '/'
    pargs = setup.parser(['-r', droot, '-s', 'wht_isis_blue', '--extension', '.fit'])
    setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'wht_isis_blue*'))
    ext = [f.split('.')[-1] for f in files]
    #expected = ['lst', 'pypeit', 'setups', 'sorted']
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)
    #os.remove('wht_isis_blue.pypeit')

@dev_suite_required
def test_setup_vlt_xshooter_uvb():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/VLT_XSHOOTER/UVB_1x1')
    droot += '/XSHO'
    pargs = setup.parser(['-r', droot, '-s', 'vlt_xshooter_uvb'])
    setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'vlt_xshooter_uvb*'))
    ext = [f.split('.')[-1] for f in files]
    #expected = ['lst', 'pypeit', 'setups', 'sorted']
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)
    #os.remove('vlt_xshooter_uvb.pypeit')

@dev_suite_required
def test_setup_vlt_xshooter_vis():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/VLT_XSHOOTER/VIS_1x1')
    droot += '/XSHO'
    pargs = setup.parser(['-r', droot, '-s', 'vlt_xshooter_vis'])
    setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'vlt_xshooter_vis*'))
    ext = [f.split('.')[-1] for f in files]
    #expected = ['lst', 'pypeit', 'setups', 'sorted']
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)
    #os.remove('vlt_xshooter_vis.pypeit')

@dev_suite_required
def test_setup_vlt_xshooter_nir():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/VLT_XSHOOTER/NIR')
    droot += '/XSHO'
    pargs = setup.parser(['-r', droot, '-s', 'vlt_xshooter_nir'])
    setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'vlt_xshooter_nir*'))
    ext = [f.split('.')[-1] for f in files]
    #expected = ['lst', 'pypeit', 'setups', 'sorted']
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)
    #os.remove('vlt_xshooter_nir.pypeit')

@dev_suite_required
def test_setup_gemini_gnirs():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/Gemini_GNIRS/GNIRS/')
    droot += '/cN'
    pargs = setup.parser(['-r', droot, '-s', 'gemini_gnirs'])
    setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'gemini_gnirs*'))
    ext = [f.split('.')[-1] for f in files]
    #expected = ['lst', 'pypeit', 'setups', 'sorted']
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)
    #os.remove('gemini_gnirs.pypeit')

# TODO: Add other instruments!

