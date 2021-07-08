"""
Module to run tests on SensFunc and FluxCalibrate classes
Requires files in Development suite (Cooked) and an Environmental variable
"""
import os

import pytest

from IPython import embed

from pypeit import fluxcalibrate
from pypeit import sensfunc
from pypeit.scripts import flux_calib
from pypeit.tests.tstutils import cooked_required, telluric_required
from pypeit.spectrographs.util import load_spectrograph
from pypeit.spectrographs import keck_deimos
from pypeit import specobjs


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


@pytest.fixture
@cooked_required
def kast_blue_files():
    std_file = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Science',
                            'spec1d_b24-Feige66_KASTb_20150520T041246.960.fits')
    sci_file = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Science',
                            'spec1d_b27-J1217p3905_KASTb_20150520T045733.560.fits')
    return [std_file, sci_file]


sens_file = data_path('sensfunc.fits')

@cooked_required
def test_gen_sensfunc(kast_blue_files):

    if os.path.isfile(sens_file):
        os.remove(sens_file)

    # Get it started
    spectrograph = load_spectrograph('shane_kast_blue')
    par = spectrograph.default_pypeit_par()
    std_file, sci_file = kast_blue_files
    # Instantiate
    sensFunc = sensfunc.UVISSensFunc(std_file, sens_file)
    # Test the standard loaded
    assert sensFunc.meta_spec['BINNING'] == '1,1'
    assert sensFunc.meta_spec['TARGET'] == 'Feige 66'

    # Generate the sensitivity function
    sensFunc.run()
    # Test
    assert os.path.basename(sensFunc.std_cal) == 'feige66_002.fits'
    # TODO: @jhennawi, please check this edit
    assert 'SENS_ZEROPOINT' in sensFunc.sens.keys(), 'Bad column names'
    # Write
    sensFunc.to_file(sens_file)

    os.remove(sens_file)


@cooked_required
def test_from_sens_func(kast_blue_files):
    # TODO: Tests need to be self-contained...
    if os.path.isfile(sens_file):
        os.remove(sens_file)

    spectrograph = load_spectrograph('shane_kast_blue')
    par = spectrograph.default_pypeit_par()
    std_file, sci_file = kast_blue_files

    # Build the sensitivity function
    sensFunc = sensfunc.UVISSensFunc(std_file, sens_file)
    sensFunc.run()
    sensFunc.to_file(sens_file)

    # Instantiate and run
    outfile = data_path(os.path.basename(sci_file))
    fluxCalibrate = fluxcalibrate.MultiSlitFC([sci_file], [sens_file], par=par['fluxcalib'],
                                              outfiles=[outfile])
    # Test
    sobjs = specobjs.SpecObjs.from_fitsfile(outfile)
    assert 'OPT_FLAM' in sobjs[0].keys()

    os.remove(sens_file)
    os.remove(outfile)


@telluric_required
def test_wmko_flux_std():

    outfile = data_path('tmp_sens.fits')
    if os.path.isfile(outfile):
        os.remove(outfile)

    # Do it
    wmko_file = data_path('2017may28_d0528_0088.fits')
    spectrograph = load_spectrograph('keck_deimos')

    # Load + write
    spec1dfile = data_path('tmp_spec1d.fits')
    sobjs = keck_deimos.load_wmko_std_spectrum(wmko_file, outfile=spec1dfile)

    # Sensfunc
    #  The following mirrors the main() call of sensfunc.py
    par = spectrograph.default_pypeit_par()
    par['sensfunc']['algorithm'] = "IR"
    par['sensfunc']['multi_spec_det'] = [3,7]

    # Instantiate the relevant class for the requested algorithm
    sensobj = sensfunc.SensFunc.get_instance(spec1dfile, outfile, par['sensfunc'])
    # Generate the sensfunc
    sensobj.run()
    # Write it out to a file
    sensobj.to_file(outfile)

    os.remove(spec1dfile)
    os.remove(outfile)

