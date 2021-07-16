"""
Module to run tests on SensFunc and FluxCalibrate classes
Requires files in Development suite (Cooked) and an Environmental variable
"""
import os

import pytest

from IPython import embed

import numpy as np

from pypeit import fluxcalibrate
from pypeit import sensfunc
from pypeit.par import pypeitpar
from pypeit.scripts import flux_calib
from pypeit.tests.tstutils import cooked_required, telluric_required
from pypeit.spectrographs.util import load_spectrograph
from pypeit.spectrographs import keck_deimos
from pypeit import specobjs, specobj


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

def test_extinction_correction_uvis():
    extinction_correction_tester('UVIS')

def test_extinction_correction_ir():
    extinction_correction_tester('IR')

def extinction_correction_tester(algorithm):
    spec1d_file = data_path('spec1d_test.fits')
    sens_file = data_path('sens_test.fits')

    if os.path.isfile(spec1d_file):
        os.remove(spec1d_file)
    if os.path.isfile(sens_file):
        os.remove(sens_file)

    # make a bogus spectrum that has N_lam = 1
    wave = np.linspace(4000, 6000)
    counts = np.ones_like(wave)
    ivar = np.ones_like(wave)
    sobj = specobj.SpecObj.from_arrays('MultiSlit', wave, counts, ivar)
    sobjs = specobjs.SpecObjs([sobj])

    # choice of PYP_SPEC, DISPNAME and EXPTIME are unimportant here
    # AIRMASS must be > 1
    sobjs.write_to_fits({
        'PYP_SPEC': 'p200_dbsp_blue',
        'DISPNAME': '600/4000',
        'EXPTIME': 1.0,
        'AIRMASS': 1.1}, spec1d_file)

    par = pypeitpar.PypeItPar()

    # set the senfunc algorithm
    par['sensfunc']['algorithm'] = algorithm
    # needed to initiate SensFunc (dummy standard star Feige34)
    par['sensfunc']['star_ra'] = 159.9042
    par['sensfunc']['star_dec'] = 43.1025

    sensobj = sensfunc.SensFunc.get_instance(spec1d_file, sens_file, par['sensfunc'])

    sensobj.wave = np.linspace(3000, 6000, 300).reshape((300, 1))
    sensobj.sens = sensobj.empty_sensfunc_table(*sensobj.wave.T.shape)
    # make the zeropoint such that the sensfunc is flat
    sensobj.zeropoint = 30 - np.log10(sensobj.wave ** 2) / 0.4

    sensobj.to_file(sens_file)

    # now flux our N_lam = 1 specobj
    par['fluxcalib']['extinct_correct'] = None
    fluxCalibrate = fluxcalibrate.MultiSlitFC([spec1d_file], [sens_file], par=par['fluxcalib'])
    # without extinction correction, we should get constant F_lam
    # with extinction correction, the spectrum will be blue

    # make sure that the appropriate default behavior occurred
    sobjs = specobjs.SpecObjs.from_fitsfile(spec1d_file)
    print(sobjs[0].keys())
    if algorithm == 'UVIS':
        assert sobjs[0]['OPT_FLAM'][0] > sobjs[0]['OPT_FLAM'][-1], \
            "UVIS sensfunc was not extinction corrected by default, but should have been"
    elif algorithm == 'IR':
        assert np.isclose(sobjs[0]['OPT_FLAM'][0], sobjs[0]['OPT_FLAM'][-1]), \
            "IR sensfunc was extinction corrected by default, but shouldn't have been"

    # clean up
    os.remove(spec1d_file)
    os.remove(sens_file)


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

