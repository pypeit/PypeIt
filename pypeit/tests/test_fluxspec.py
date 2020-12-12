"""
Module to run tests on SensFunc and FluxCalibrate classes
Requires files in Development suite (Cooked) and an Environmental variable
"""
import os

import pytest

from pypeit import fluxcalibrate
from pypeit import sensfunc
from pypeit.scripts import flux_calib
from pypeit.tests.tstutils import cooked_required
from pypeit.spectrographs.util import load_spectrograph
from pypeit import specobjs


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

## TODO: Not used
#@pytest.fixture
#@dev_suite_required
#def deimos_files():
#    return [os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Science',
#                         'spec1d_G191B2B_DEIMOS_2017Sep14T152432.fits')]

@pytest.fixture
@cooked_required
def kast_blue_files():
    std_file = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Science',
                            'spec1d_b24-Feige66_KASTb_2015May20T041246.960.fits')
    sci_file = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Science',
                            'spec1d_b27-J1217p3905_KASTb_2015May20T045733.560.fits')
    return [std_file, sci_file]


sens_file = data_path('sensfunc.fits')

@cooked_required
def test_gen_sensfunc(kast_blue_files):
    # Get it started
    spectrograph = load_spectrograph('shane_kast_blue')
    par = spectrograph.default_pypeit_par()
    std_file, sci_file = kast_blue_files
    # Instantiate
    sensFunc = sensfunc.UVIS(std_file, sens_file)
    # Test the standard loaded
    assert sensFunc.meta_spec['BINNING'] == '1,1'
    assert sensFunc.meta_spec['TARGET'] == 'Feige 66'

    # Generate the sensitivity function
    sensFunc.run()
    # Test
    assert os.path.basename(sensFunc.meta_table['CAL_FILE'][0]) == 'feige66_002.fits'
    assert 'SENSFUNC' in sensFunc.out_table.keys()
    # Write
    sensFunc.save()


@cooked_required
def test_from_sens_func(kast_blue_files):
    """ This test will fail if the previous one does as it need its output
    """
    spectrograph = load_spectrograph('shane_kast_blue')
    par = spectrograph.default_pypeit_par()
    std_file, sci_file = kast_blue_files
    #
    # Instantiate and run
    outfile = data_path(os.path.basename(sci_file))
    fluxCalibrate = fluxcalibrate.MultiSlitFC([sci_file], [sens_file], par=par['fluxcalib'],
                                              outfiles=[outfile])
    # Test
    sobjs = specobjs.SpecObjs.from_fitsfile(outfile)
    assert 'OPT_FLAM' in sobjs[0].keys()

    os.remove(sens_file)
    os.remove(outfile)

