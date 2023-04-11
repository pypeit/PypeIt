"""
Module to run tests on OneSpec
"""
from pathlib import Path

from IPython import embed

import numpy as np

from pypeit import onespec
from pypeit.tests.tstutils import data_path


def test_init():
    wave = np.linspace(3500,10000,1000)
    flux = np.ones(1000, dtype=float)
    spec = onespec.OneSpec(wave, wave, flux)
    assert np.array_equal(wave, spec.wave), 'Wavelengths munged'
    assert np.array_equal(flux, spec.flux), 'Flux munged'
    assert spec.ivar is None, 'IVAR should not be set'
    assert spec.spectrograph is None, 'Spectrograph should not be set'

    spec = onespec.OneSpec(wave, wave, flux, ivar=2*np.ones_like(flux))
    assert np.allclose(spec.sig, 1/np.sqrt(2)), 'Conversion to sigma is wrong'


def test_io():
    wave = np.linspace(3500,10000,1000)
    flux = np.ones(1000, dtype=float)
    # TODO: PYP_SPEC is required if we want to be able to read the file!
    spec = onespec.OneSpec(wave, wave, flux, PYP_SPEC='shane_kast_blue')
    ofile = Path(data_path('tmp.fits')).resolve()
    spec.to_file(str(ofile), overwrite=True)
    _spec = onespec.OneSpec.from_file(ofile)
    assert np.array_equal(spec.flux, _spec.flux), 'Flux munged'

    ofile.unlink()


