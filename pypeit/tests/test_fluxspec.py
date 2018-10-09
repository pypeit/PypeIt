# Module to run tests on FluxSpec class
#   Requires files in Development suite (Cooked) and an Environmental variable
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

# TEST_UNICODE_LITERALS

import os

import pytest

from pypeit import fluxspec
from pypeit.scripts import flux_spec
from pypeit.tests.tstutils import dev_suite_required


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


@pytest.fixture
@dev_suite_required
def master_dir():
    # Any test that uses this directory also requires the DevSuite!
#    return data_path('MF_shane_kast_blue') if os.getenv('PYPEIT_DEV') is None \
#            else os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'MF_shane_kast_blue')
    return os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'MF_shane_kast_blue')


# TODO: Not used
@pytest.fixture
@dev_suite_required
def deimos_files():
    return [os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Science',
                         'spec1d_G191B2B_DEIMOS_2017Sep14T152432.fits')]


@pytest.fixture
@dev_suite_required
def kast_blue_files():
    std_file = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Science',
                            'spec1d_Feige66_KASTb_2015May20T041246.96.fits')
    sci_file = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Science',
                            'spec1d_J1217p3905_KASTb_2015May20T045733.56.fits')
    return [std_file, sci_file]


@dev_suite_required
def test_run_from_spec1d(kast_blue_files, master_dir):
    # Instantiate
    std_file, sci_file = kast_blue_files
    FxSpec = fluxspec.FluxSpec(std_spec1d_file=std_file, sci_spec1d_file=sci_file,
                             spectrograph='shane_kast_blue', setup='A_01_aa',
                               master_dir=master_dir)
    assert FxSpec.frametype == 'sensfunc'
    # Find the standard
    std = FxSpec.find_standard()
    #assert std.idx == 'O479-S5009-D01-I0023'
    # Generate the sensitivity function
    sens_dict = FxSpec.generate_sensfunc()
    assert isinstance(sens_dict, dict)
    assert 'FEIGE66' in sens_dict['std_name']
    assert FxSpec.steps[-1] == 'generate_sensfunc'
    # Flux me some science
    FxSpec.flux_science()
    assert 'FLAM' in FxSpec.sci_specobjs[0].optimal.keys()
    # Write
    FxSpec.write_science(data_path('tmp.fits'))
    # Master
    FxSpec.save_master()
    # Load from Master
    sens_dict, _, _ = FxSpec.load_master_frame(force=True)
    assert 'FEIGE66' in sens_dict['std_name']
    

@dev_suite_required
def test_from_sens_func(master_dir):
    """ This test will fail if the previous one does as it need its output
    """
    # TODO: Should change this to a from_sens_file instance.  Most of
    # the class is uninstantiated and methods will fail if you
    # instantiate this way...
    FxSpec3 = fluxspec.FluxSpec(sens_file=os.path.join(master_dir,'MasterSensFunc_A_aa.fits'))
    assert isinstance(FxSpec3.sens_dict, dict)


@dev_suite_required
def test_script(kast_blue_files, deimos_files):
    std_file, sci_file = kast_blue_files
    # Sensitivity function
    pargs = flux_spec.parser(['sensfunc',
                              '--std_file={:s}'.format(std_file),
                              '--instr=shane_kast_blue',
                              '--sensfunc_file={:s}'.format(data_path('tmp.fits'))])
    # Run
    flux_spec.main(pargs)

    # Flux me
    pargs2 = flux_spec.parser(['flux',
                               '--sci_file={:s}'.format(sci_file),
                               '--sensfunc_file={:s}'.format(data_path('tmp.fits')),
                               '--flux_file={:s}'.format(data_path('tmp_fluxed.fits'))])
    flux_spec.main(pargs2)

    # DEIMOS (multi-det)
    #pypeit_flux_spec sensfunc --std_file=spec1d_G191B2B_DEIMOS_2017Sep14T152432.fits  --instr=keck_deimos --sensfunc_file=sens.yaml --multi_det=3,7
    '''
    std_file = deimos_files[0]
    pargs3 = flux_spec.parser(['sensfunc',
                               '--std_file={:s}'.format(std_file),
                               '--instr=keck_deimos',
                               '--sensfunc_file={:s}'.format(data_path('tmp2.yaml')),
                               '--multi_det=3,7'])
    flux_spec.main(pargs3)
    '''


