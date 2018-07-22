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

# These tests are not run on Travis
if os.getenv('PYPIT_DEV') is None:
    skip_test=True
else:
    skip_test=False

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

@pytest.fixture
def deimos_files():
    if not skip_test:
        deimos_std_file = os.path.join(os.getenv('PYPIT_DEV'), 'Cooked', 'Science',
                                                 'spec1d_G191B2B_DEIMOS_2017Sep14T152432.fits')
    else:
        deimos_std_file = None
    return [deimos_std_file]

@pytest.fixture
def kast_blue_files():
    if not skip_test:
        std_file = os.path.join(os.getenv('PYPIT_DEV'), 'Cooked', 'Science',
                                          'spec1d_Feige66_KASTb_2015May20T041246.96.fits')
        sci_file = os.path.join(os.getenv('PYPIT_DEV'), 'Cooked', 'Science',
                                          'spec1d_J1217p3905_KASTb_2015May20T045733.56.fits')
        kast_blue_files = [std_file, sci_file]
    else:
        kast_blue_files = None
    return kast_blue_files

def test_run_from_spec1d(kast_blue_files):
    if skip_test:
        assert True
        return
    # Instantiate
    std_file, sci_file = kast_blue_files
    FxSpec = fluxspec.FluxSpec(std_spec1d_file=std_file, sci_spec1d_file=sci_file,
                               spectrograph='shane_kast_blue', setup='A_01_aa',
                               root_path=data_path('MF'))
    assert FxSpec.frametype == 'sensfunc'
    # Find the standard
    std = FxSpec.find_standard()
    #assert std.idx == 'O479-S5009-D01-I0023'
    # Generate the sensitivity function
    sensfunc = FxSpec.generate_sensfunc()
    assert isinstance(sensfunc, dict)
    assert 'feige66' in sensfunc['std']['file']
    assert FxSpec.steps[-1] == 'generate_sensfunc'
    # Flux me some science
    FxSpec.flux_science()
    assert 'flam' in FxSpec.sci_specobjs[0].optimal.keys()
    # Write
    FxSpec.write_science(data_path('tmp.fits'))
    # Master
    FxSpec.save_master()
    # Load from Master
    sensfunc2, _, _ = FxSpec.load_master_frame(force=True)
    assert 'feige66' in sensfunc2['std']['file']


def test_from_sens_func():
    """ This test will fail if the previous one does as it need its output
    """
    if skip_test:
        assert True
        return
    # TODO: Should change this to a from_sens_file instance.  Most of
    # the class is uninstantiated and methods will fail if you
    # instantiate this way...
    FxSpec3 = fluxspec.FluxSpec(sens_file=data_path('MF_shane_kast_blue/MasterSensFunc_A_aa.yaml'))
    assert isinstance(FxSpec3.sensfunc, dict)


def test_script(kast_blue_files, deimos_files):
    if skip_test:
        assert True
        return
    std_file, sci_file = kast_blue_files
    # Sensitivity function
    pargs = flux_spec.parser(['sensfunc',
                              '--std_file={:s}'.format(std_file),
                              '--instr=shane_kast_blue',
                              '--sensfunc_file={:s}'.format(data_path('tmp.yaml'))])
    # Run
    flux_spec.main(pargs)

    # Flux me
    pargs2 = flux_spec.parser(['flux',
                               '--sci_file={:s}'.format(sci_file),
                               '--sensfunc_file={:s}'.format(data_path('tmp.yaml')),
                               '--flux_file={:s}'.format(data_path('tmp.fits'))])
    flux_spec.main(pargs2)

    # DEIMOS (multi-det)
    #pypeit_flux_spec sensfunc --std_file=spec1d_G191B2B_DEIMOS_2017Sep14T152432.fits  --instr=keck_deimos --sensfunc_file=sens.yaml --multi_det=3,7
    std_file = deimos_files[0]
    pargs3 = flux_spec.parser(['sensfunc',
                               '--std_file={:s}'.format(std_file),
                               '--instr=keck_deimos',
                               '--sensfunc_file={:s}'.format(data_path('tmp2.yaml')),
                               '--multi_det=3,7'])
    flux_spec.main(pargs3)
