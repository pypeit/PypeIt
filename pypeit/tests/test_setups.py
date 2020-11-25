"""
Module to run tests on scripts
"""
import os
import glob
import shutil

from IPython import embed

import numpy as np

import pytest
import yaml
from configobj import ConfigObj

from pypeit.pypmsgs import PypeItError
from pypeit.metadata import PypeItMetaData
from pypeit.par import PypeItPar
from pypeit.par.util import parse_pypeit_file
from pypeit.scripts import setup, chk_for_calibs
from pypeit.spectrographs.util import load_spectrograph
from pypeit.tests.tstutils import dev_suite_required
from pypeit import pypeit
from pypeit import pypeitsetup


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def expected_file_extensions():
    return ['pypeit', 'sorted']


def test_run_setup():
    """ Test the setup script
    """
    # Remove .setup if needed
    sfiles = glob.glob('*.setups')
    for sfile in sfiles:
        os.remove(sfile)
    #
    droot = data_path('b')
    pargs = setup.parse_args(['-r', droot, '-s', 'shane_kast_blue', '-c=all',
                              '--extension=fits.gz', '--output_path={:s}'.format(data_path(''))])
    setup.main(pargs)

    #setup_file = glob.glob(data_path('setup_files/shane_kast_blue*.setups'))[0]
    ## Load
    #with open(setup_file, 'r') as infile:
    #    setup_dict = yaml.load(infile)
    ## Test
    #assert '01' in setup_dict['A'].keys()
    #assert setup_dict['A']['--']['disperser']['name'] == '600/4310'
    # Failures
    pargs2 = setup.parse_args(['-r', droot, '-s', 'shane_kast_blu', '-c=all',
                               '--extension=fits.gz', '--output_path={:s}'.format(data_path(''))])
    with pytest.raises(ValueError):
        setup.main(pargs2)
    
    # Cleanup
    shutil.rmtree(data_path('setup_files'))


def test_setup_made_pypeit_file():
    """ Test the .pypeit file(s) made by pypeit_setup
    This test depends on the one above
    """
    pypeit_file = data_path('shane_kast_blue_A/shane_kast_blue_A.pypeit')
    cfg_lines, data_files, frametype, usrdata, setups = parse_pypeit_file(pypeit_file)
    # Test
    assert len(data_files) == 2
    assert sorted(frametype['b1.fits.gz'].split(',')) == ['arc', 'tilt']
    assert setups[0] == 'A'

    # Cleanup
    shutil.rmtree(data_path('shane_kast_blue_A'))


@dev_suite_required
def test_setup_keck_lris_red():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/keck_lris_red/multi_400_8500_d560')
    droot += '/'
    pargs = setup.parse_args(['-r', droot, '-s', 'keck_lris_red'])
    setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'keck_lris_red*'))
    ext = [f.split('.')[-1] for f in files]
    expected = expected_file_extensions()
    assert np.all([e in ext for e in expected]), \
        'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)


@dev_suite_required
def test_setup_keck_lris_red_orig():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/keck_lris_red_orig/long_300_5000')
    droot += '/'
    pargs = setup.parse_args(['-r', droot, '-s', 'keck_lris_red_orig'])
    setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'keck_lris_red_orig*'))
    ext = [f.split('.')[-1] for f in files]
    expected = expected_file_extensions()
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)


@dev_suite_required
def test_setup_keck_lris_blue():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/keck_lris_blue/multi_600_4000_d560')
    droot += '/'
    pargs = setup.parse_args(['-r', droot, '-s', 'keck_lris_blue'])
    setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'keck_lris_blue*'))
    ext = [f.split('.')[-1] for f in files]
    expected = expected_file_extensions()
    assert np.all([e in ext for e in expected]), \
        'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)


@dev_suite_required
def test_setup_keck_lris_blue_orig():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/keck_lris_blue_orig/long_600_4000_d500')
    droot += '/'
    pargs = setup.parse_args(['-r', droot, '-s', 'keck_lris_blue_orig'])
    setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'keck_lris_blue_orig*'))
    ext = [f.split('.')[-1] for f in files]
    expected = expected_file_extensions()
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)


@dev_suite_required
def test_setup_shane_kast_blue():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/shane_kast_blue/600_4310_d55')
    droot += '/'
    pargs = setup.parse_args(['-r', droot, '-s', 'shane_kast_blue'])
    setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'shane_kast_blue*'))
    ext = [f.split('.')[-1] for f in files]
    expected = expected_file_extensions()
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)


@dev_suite_required
def test_setup_shane_kast_red():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/shane_kast_red/600_7500_d55_ret')
    droot += '/'
    pargs = setup.parse_args(['-r', droot, '-s', 'shane_kast_red'])
    setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'shane_kast_red*'))
    ext = [f.split('.')[-1] for f in files]
    expected = expected_file_extensions()
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)

# TODO: We need a test data set for shane_kast_red_ret

@dev_suite_required
def test_setup_keck_deimos():

    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/keck_deimos/830G_M_8600')
    droot += '/'
    pargs = setup.parse_args(['-r', droot, '-s', 'keck_deimos'])
    setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'keck_deimos*'))
    ext = [f.split('.')[-1] for f in files]
    expected = expected_file_extensions()
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)


@dev_suite_required
def test_setup_keck_deimos_multiconfig():

    root = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'keck_deimos')
    files = glob.glob(os.path.join(root, '830G_L_8100', '*fits*'))
    files += glob.glob(os.path.join(root, '830G_L_8400', '*fits*'))

    output_path = os.path.join(os.getcwd(), 'output')
    if os.path.isdir(output_path):
        shutil.rmtree(output_path)
    os.makedirs(output_path)

    ps = pypeitsetup.PypeItSetup(files, spectrograph_name='keck_deimos')
    ps.run(setup_only=True, sort_dir=output_path)
    # Write the automatically generated pypeit data
    pypeit_files = ps.fitstbl.write_pypeit(output_path, cfg_lines=ps.user_cfg,
                                           write_bkg_pairs=True)

    assert len(pypeit_files) == 2, 'Should have created two pypeit files'

    # Test the pypeit files for the correct configuration and
    # calibration group results
    for f, s, c in zip(pypeit_files, ['A', 'B'], ['0', '1']):

        # TODO: All of this front-end stuff, pulled from pypeit.py, should
        # be put into a function.

        # Read the pypeit file
        cfg_lines, data_files, frametype, usrdata, setups = parse_pypeit_file(f, runtime=True)
        # Spectrograph
        cfg = ConfigObj(cfg_lines)
        spectrograph = load_spectrograph(cfg['rdx']['spectrograph'])
        # Configuration-specific parameters
        for idx, row in enumerate(usrdata):
            if 'science' in row['frametype'] or 'standard' in row['frametype']:
                break
        spectrograph_cfg_lines = spectrograph.config_specific_par(data_files[idx]).to_config()
        #  PypeIt parameters
        par = PypeItPar.from_cfg_lines(cfg_lines=spectrograph_cfg_lines, merge_with=cfg_lines)
        #  Metadata
        fitstbl = PypeItMetaData(spectrograph, par, files=data_files, usrdata=usrdata, strict=True)
        fitstbl.finalize_usr_build(frametype, setups[0])

        assert np.all(fitstbl['setup'] == s), 'Setup is wrong'
        assert np.all(fitstbl['calib'] == c), 'Calibration group is wrong'

    # Clean-up
    shutil.rmtree(output_path)


@dev_suite_required
def test_setup_keck_deimos_multiconfig_clean():

    root = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'keck_deimos')
    files = glob.glob(os.path.join(root, '830G_L_8100', '*fits*'))
    files += glob.glob(os.path.join(root, '830G_L_8400', '*fits*'))

    ps = pypeitsetup.PypeItSetup(files, spectrograph_name='keck_deimos')
    ps.build_fitstbl(strict=False)
    ps.get_frame_types(flag_unknown=True)

    # Test that the correct number of configurations are found
    cfgs = ps.fitstbl.unique_configurations()
    assert len(cfgs) == 2, 'Should find 2 configurations'

    # Test that the bias is assigned to the correct configuration
    ps.fitstbl.set_configurations(cfgs)
    biases = np.where(ps.fitstbl.find_frames('bias'))[0]
    assert biases.size == 1, 'Should only be 1 bias'
    assert ps.fitstbl['setup'][biases[0]] == 'B', 'Bias should be in configuration group B'

    # Table should have 25 rows
    assert len(ps.fitstbl) == 25, 'Incorrect number of table rows.'

    # All frames should be from valid configurations
    ps.fitstbl.clean_configurations()
    assert len(ps.fitstbl) == 25, 'Incorrect number of table rows.'

    # Artificially set the amplifier and mode of two frames to be
    # invalid
    ps.fitstbl['amp'][0] = 'SINGLE:A'
    ps.fitstbl['mode'][1] = 'Direct'
    ps.fitstbl.clean_configurations()
    # Those two frames should have been removed
    assert len(ps.fitstbl) == 23, 'Incorrect number of table rows.'


@dev_suite_required
def test_setup_keck_nires():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/keck_nires/NIRES/')
    droot += '/'
    pargs = setup.parse_args(['-r', droot, '-s', 'keck_nires'])
    setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'keck_nires*'))
    ext = [f.split('.')[-1] for f in files]
    expected = expected_file_extensions()
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)


@dev_suite_required
def test_setup_keck_nirspec():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/keck_nirspec/LOW_NIRSPEC-1')
    droot += '/'
    pargs = setup.parse_args(['-r', droot, '-s', 'keck_nirspec_low'])
    setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'keck_nirspec*'))
    ext = [f.split('.')[-1] for f in files]
    expected = expected_file_extensions()
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)


@dev_suite_required
def test_setup_magellan_mage():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/magellan_mage/1x1')
    droot += '/'
    pargs = setup.parse_args(['-r', droot, '-s', 'magellan_mage'])
    setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'magellan_mage*'))
    ext = [f.split('.')[-1] for f in files]
    expected = expected_file_extensions()
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)


@dev_suite_required
def test_setup_wht_isis_blue():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/wht_isis_blue/long_R300B_d5300')
    droot += '/'
    pargs = setup.parse_args(['-r', droot, '-s', 'wht_isis_blue', '--extension', '.fit'])
    setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'wht_isis_blue*'))
    ext = [f.split('.')[-1] for f in files]
    expected = expected_file_extensions()
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)


@dev_suite_required
def test_setup_vlt_xshooter_uvb():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/vlt_xshooter/UVB_1x1')
    droot += '/XSHO'
    pargs = setup.parse_args(['-r', droot, '-s', 'vlt_xshooter_uvb'])
    setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'vlt_xshooter_uvb*'))
    ext = [f.split('.')[-1] for f in files]
    expected = expected_file_extensions()
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)


@dev_suite_required
def test_setup_vlt_xshooter_vis():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/vlt_xshooter/VIS_1x1')
    droot += '/XSHO'
    pargs = setup.parse_args(['-r', droot, '-s', 'vlt_xshooter_vis'])
    setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'vlt_xshooter_vis*'))
    ext = [f.split('.')[-1] for f in files]
    expected = expected_file_extensions()
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)


@dev_suite_required
def test_setup_vlt_xshooter_nir():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/vlt_xshooter/NIR')
    droot += '/XSHO'
    pargs = setup.parse_args(['-r', droot, '-s', 'vlt_xshooter_nir'])
    setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'vlt_xshooter_nir*'))
    ext = [f.split('.')[-1] for f in files]
    expected = expected_file_extensions()
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)


@dev_suite_required
def test_setup_gemini_gnirs():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/gemini_gnirs/32_SB_SXD/')
    droot += '/cN'
    pargs = setup.parse_args(['-r', droot, '-s', 'gemini_gnirs'])
    setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'gemini_gnirs*'))
    ext = [f.split('.')[-1] for f in files]
    expected = expected_file_extensions()
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)


@dev_suite_required
def test_setup_not_alfosc():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/not_alfosc/grism4')
    droot += '/ALD'
    pargs = setup.parse_args(['-r', droot, '-s', 'not_alfosc'])
    setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'not_alfosc*'))
    ext = [f.split('.')[-1] for f in files]
    expected = expected_file_extensions()
    assert np.all([e in ext for e in expected]), \
        'Did not find all setup file extensions: {0}'.format(expected)

    # Build a PypeIt file
    pargs = setup.parse_args(['-r', droot, '-s', 'not_alfosc', '-c', 'A', '-d', data_path('')])
    setup.main(pargs)
    pypeit_file = data_path('not_alfosc_A/not_alfosc_A.pypeit')
    # TODO: Why is this using pypeit.PypeIt and not pypeitsetup.PypeItSetup?
    pypeIt = pypeit.PypeIt(pypeit_file, calib_only=True)

    # Clean-up
    shutil.rmtree(setup_dir)
    shutil.rmtree(data_path('not_alfosc_A'))

@dev_suite_required
def test_setup_vlt_fors2():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/vlt_fors2/300I/')
    droot += '/FORS2'
    pargs = setup.parse_args(['-r', droot, '-s', 'vlt_fors2'])
    setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'vlt_fors2*'))
    ext = [f.split('.')[-1] for f in files]
    expected = expected_file_extensions()
    assert np.all([e in ext for e in expected]), \
        'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)

    # Now chk calib
    pargs = chk_for_calibs.parse_args([droot, '-s', 'vlt_fors2'])
    answers, ps = chk_for_calibs.main(pargs)
    assert answers['pass'][0], 'A must pass!'

# TODO: Add other instruments!


