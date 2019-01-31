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
from pypeit.spectrographs.util import load_spectrograph


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)



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
                            'spec1d_Feige66_KASTb_2015May20T041246.960.fits')
    sci_file = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Science',
                            'spec1d_J1217p3905_KASTb_2015May20T045733.560.fits')
    return [std_file, sci_file]


@dev_suite_required
def test_gen_sensfunc(kast_blue_files):
    # Get it started
    spectrograph = load_spectrograph('shane_kast_blue')
    par = spectrograph.default_pypeit_par()
    std_file, sci_file = kast_blue_files
    # Instantiate
    FxSpec = fluxspec.FluxSpec(spectrograph, par['fluxcalib'])
    assert FxSpec.frametype == 'sensfunc'
    # Find the standard
    FxSpec.load_objs(std_file, std=True)
    std = FxSpec.find_standard()
    # Generate the sensitivity function
    sens_dict = FxSpec.generate_sensfunc()
    #
    assert isinstance(sens_dict, dict)
    assert 'FEIGE66' in sens_dict['std_name']
    assert FxSpec.steps[-1] == 'generate_sensfunc'
    # Master
    FxSpec.save_sens_dict(FxSpec.sens_dict, outfile=data_path('sensfunc.fits'))


@dev_suite_required
def test_from_sens_func(kast_blue_files ):
    """ This test will fail if the previous one does as it need its output
    """
    spectrograph = load_spectrograph('shane_kast_blue')
    par = spectrograph.default_pypeit_par()
    std_file, sci_file = kast_blue_files
    # Instantiate
    FxSpec = fluxspec.FluxSpec(spectrograph, par['fluxcalib'], sens_file=data_path('sensfunc.fits'))
    assert 'FEIGE66' in FxSpec.sens_dict['std_name']
    # Flux me some science
    FxSpec.flux_science(sci_file)
    assert 'FLAM' in FxSpec.sci_specobjs[0].optimal.keys()
    # Write
    FxSpec.write_science(data_path('tmp.fits'))


@dev_suite_required
def test_script():
    # Sensitivity function
    pargs = flux_spec.parser([data_path('test.flux')])
    # Run
    flux_spec.main(pargs, unit_test=True)

    # Check for output
    assert os.path.isfile('test_sensfunc.fits')


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


