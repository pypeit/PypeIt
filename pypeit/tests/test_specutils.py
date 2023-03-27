"""
Module to test the specutils interface
"""
from pathlib import Path

from IPython import embed

import numpy as np

from astropy.io import fits

from pypeit import onespec
from pypeit import specobj
from pypeit import specobjs
from pypeit.specutils import Spectrum1D, SpectrumList
from pypeit.tests import tstutils

import pytest
specutils_required = pytest.mark.skipif(Spectrum1D is None or SpectrumList is None,
                                        reason='specutils not installed')

@specutils_required
def test_onespec_io():
    rng = np.random.default_rng()
    grid_wave = np.linspace(3500,10000,1000)
    wave = grid_wave + (2*rng.uniform(size=grid_wave.size) - 1)
    flux = np.ones(1000, dtype=float)
    # TODO: PYP_SPEC is required if we want to be able to read the file!
    spec = onespec.OneSpec(wave, grid_wave, flux, PYP_SPEC='shane_kast_blue')
    ofile = Path(tstutils.data_path('tmp.fits')).resolve()
    spec.to_file(str(ofile), overwrite=True)

    _spec = Spectrum1D.read(ofile)
    assert np.array_equal(spec.flux, _spec.flux.data), 'Flux munged'
    assert np.array_equal(spec.wave, _spec.spectral_axis.data), 'Wavelengths munged'

    _spec = Spectrum1D.read(ofile, grid=True)
    assert np.array_equal(spec.wave_grid_mid, _spec.spectral_axis.data), 'Wavelengths munged'

    ofile.unlink()


@specutils_required
def test_spec1d_io():

    ofile = Path(tstutils.data_path('tmp.fits')).resolve()

    spec1 = specobj.SpecObj('MultiSlit', 'DET01', SLITID=0)
    npix_spec = 100

    # NOTE: WAVE, COUNTS, IVAR, and MASK must *all* exist to be able to read the
    # file using specutils because of how the pypeit loader works (i.e., this
    # isn't something intrinsic to specutils).
    spec1['BOX_WAVE'] = np.arange(npix_spec).astype(float)
    spec1['BOX_COUNTS'] = np.ones(npix_spec, dtype=float)
    spec1['BOX_COUNTS_IVAR'] = np.ones(npix_spec, dtype=float)
    spec1['BOX_MASK'] = np.ones(npix_spec, dtype=bool)
    spec1['DETECTOR'] = tstutils.get_kastb_detector()

    allspec = specobjs.SpecObjs([spec1])
    header = fits.PrimaryHDU().header
    header['TST'] = 'TEST'
    allspec.write_to_fits(header, str(ofile), overwrite=True)

    spec = SpectrumList.read(ofile)
    assert len(spec) == 1, 'Should be 1 spectrum'
    assert np.array_equal(spec[0].flux.data, spec1.BOX_COUNTS), 'Bad flux read'
    assert spec[0].meta['name'] == 'SPAT-----SLIT0000-DET01', 'Name changed'
    assert spec[0].meta['extract'] == 'BOX', 'Should have read spectrum as a boxcar extraction'
    assert not spec[0].meta['fluxed'], 'Should have read spectrum as uncalibrated counts'

    spec1['OPT_WAVE'] = spec1['BOX_WAVE'].copy()
    spec1['OPT_COUNTS'] = spec1['BOX_COUNTS'] * 1.1
    spec1['OPT_COUNTS_IVAR'] = spec1['BOX_COUNTS_IVAR'] * 1.1
    spec1['OPT_MASK'] = spec1['BOX_MASK'].copy()

    allspec.write_to_fits(header, str(ofile), overwrite=True)
    spec = SpectrumList.read(ofile)
    assert np.array_equal(spec[0].flux.data, spec1.OPT_COUNTS), 'Bad flux read'
    assert spec[0].meta['extract'] == 'OPT', 'Should have read spectrum as an optimal extraction'

    spec1['OPT_FLAM'] = spec1['OPT_COUNTS'] * 1.1
    spec1['OPT_FLAM_IVAR'] = spec1['OPT_COUNTS_IVAR'] * 1.1

    allspec.write_to_fits(header, str(ofile), overwrite=True)
    spec = SpectrumList.read(ofile)
    assert np.array_equal(spec[0].flux.data, spec1.OPT_FLAM), 'Bad flux read'
    assert spec[0].meta['extract'] == 'OPT', 'Should have read spectrum as an optimal extraction'
    assert spec[0].meta['fluxed'], 'Should have read the flux-calibrated data'

    # Try reading the flux-calibrated boxcar extraction
    spec = SpectrumList.read(ofile, extract='BOX', fluxed=True)
    # TODO: Issue warning when requested data is not available?
    assert spec[0].meta['extract'] == 'BOX', 'Should have read spectrum as a boxcar extraction'
    # BOX_FLAM doesn't exist, so loader resorts to reading the count data
    assert not spec[0].meta['fluxed'], 'Should have read the count data'

    # Add a second spectrum
    spec2 = specobj.SpecObj('MultiSlit', 'DET01', SLITID=1)

    # NOTE: WAVE, COUNTS, IVAR, and MASK must *all* exist to be able to read the
    # file using specutils because of how the pypeit loader works (i.e., this
    # isn't something intrinsic to specutils).
    spec2['BOX_WAVE'] = np.arange(npix_spec).astype(float)
    spec2['BOX_COUNTS'] = 2*np.ones(npix_spec, dtype=float)
    spec2['BOX_COUNTS_IVAR'] = 2*np.ones(npix_spec, dtype=float)
    spec2['BOX_MASK'] = np.ones(npix_spec, dtype=bool)
    spec2['DETECTOR'] = tstutils.get_kastb_detector()

    spec2['OPT_WAVE'] = spec2['BOX_WAVE'].copy()
    spec2['OPT_COUNTS'] = spec2['BOX_COUNTS'] * 1.1
    spec2['OPT_COUNTS_IVAR'] = spec2['BOX_COUNTS_IVAR'] * 1.1
    spec2['OPT_MASK'] = spec2['BOX_MASK'].copy()

    spec2['OPT_FLAM'] = spec2['OPT_COUNTS'] * 1.1
    spec2['OPT_FLAM_IVAR'] = spec2['OPT_COUNTS_IVAR'] * 1.1

    allspec = specobjs.SpecObjs([spec1, spec2])
    allspec.write_to_fits(header, str(ofile), overwrite=True)
    spec = SpectrumList.read(ofile)
    assert len(spec) == 2, 'Should be 2 spectra'
    assert spec[0].meta['name'] == 'SPAT-----SLIT0000-DET01', 'Name changed'
    assert spec[1].meta['name'] == 'SPAT-----SLIT0001-DET01', 'Name changed'
    assert np.array_equal(spec[1].flux.data, spec2.OPT_FLAM), 'Bad flux read'
    assert spec[1].meta['extract'] == 'OPT', 'Should have read spectrum as an optimal extraction'
    assert spec[1].meta['fluxed'], 'Should have read the flux-calibrated data'

    # Try reading the flux-calibrated boxcar extraction
    spec = SpectrumList.read(ofile, extract='BOX', fluxed=True)
    # TODO: Issue warning when requested data is not available?
    assert spec[1].meta['extract'] == 'BOX', 'Should have read spectrum as a boxcar extraction'
    # BOX_FLAM doesn't exist, so loader resorts to reading the count data
    assert not spec[1].meta['fluxed'], 'Should have read the count data'

    ofile.unlink()

