"""
Module to run tests on SensFunc and FluxCalibrate classes
Requires files in PypeIt/pypeit/data
"""
import os

import pytest

from IPython import embed

import numpy as np

from pypeit import fluxcalibrate
from pypeit import sensfunc
from pypeit.par import pypeitpar
from pypeit.tests.tstutils import data_path
from pypeit.spectrographs.util import load_spectrograph
from pypeit.spectrographs import keck_deimos
from pypeit import specobjs, specobj


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


