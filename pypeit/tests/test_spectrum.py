
from IPython import embed
import numpy as np
import pytest
from scipy import interpolate
from scipy import special

from pypeit.core import spectrum
from pypeit.pypmsgs import PypeItError


def pixelated_gaussian(x, c=0.0, s=1.0):
    n = np.sqrt(2.)*s
    d = np.asarray(x)-c
    dx = np.mean(np.diff(x))
    return (special.erf((d+dx/2.)/n) - special.erf((d-dx/2.)/n))/2./dx


def get_fake_spec():

    # Wavelengths
    npix = 2001
    wave = np.linspace(4000., 6000., npix, dtype=float)

    # Parameters from some randomly distributed pixelated Gaussians
    rng = np.random.default_rng(99)
    nlines = 10
    a = rng.uniform(size=nlines) * 10 + 0.5
    c = rng.uniform(size=nlines) * npix
    s = np.ones(nlines, dtype=float)
    x = np.arange(npix, dtype=float)

    # Error vector
    err = np.full(npix, 0.2, dtype=float)

    # Build the spectrum
    flux = np.ones(wave.size, dtype=float)
    for i in range(nlines):
        flux += a[i]*pixelated_gaussian(x, c=c[i], s=s[i])
    flux += rng.normal(size=flux.size) * err

    # Mask some random pixels
    gpm = np.ones(flux.shape, dtype=bool)
    gpm[rng.integers(0, high=npix-1, size=50)] = False

    return wave, flux, err, gpm


def test_alloc():

    wave = np.linspace(4000., 6000., 2001, dtype=float)
    flux = np.ones(wave.size, dtype=float)

    spec1D = spectrum.Spectrum(wave, flux)
    assert spec1D.ndim == 1, 'Spectrum should be 1D'

    flux = np.ones((wave.size,10), dtype=float)
    spec2D = spectrum.Spectrum(wave, flux)
    assert spec2D.ndim == 2, 'Spectrum should be 2D'
    assert spec2D.wave.ndim, 'Wavelength should be 1D'
    assert np.array_equal(spec2D.wave, wave), 'Wavelengths changed'


def test_multiply():

    wave = np.linspace(4000., 6000., 2001, dtype=float)
    flux = np.ones(wave.size, dtype=float)

    spec = spectrum.Spectrum(wave, flux)
    spec.multiply(10)
    assert np.all(spec.flux == 10.), 'Scalar multiplication failed'

    # 1D
    spec = spectrum.Spectrum(wave, flux)
    # Should fail when the spectral dimensions do not match
    with pytest.raises(PypeItError):
        spec.multiply(np.arange(10).astype(float))
    # Should fail when the multiplication would increase the dimensionality of the spectrum
    with pytest.raises(PypeItError):
        spec.multiply(np.full((flux.size,2), 2., dtype=float))
    # Should be successful when the array is the right size
    spec.multiply(np.full(flux.size, 2., dtype=float))
    assert np.all(spec.flux == 2.), 'Vector multiplication failed'

    # 2D
    spec = spectrum.Spectrum(wave, np.tile(flux, (2,1)).T)
    a = np.full(flux.size, 2., dtype=float)
    # Should be successful when the array can be broadcast to the right shape
    spec.multiply(a)
    assert np.all(spec.flux == 2.), 'Vector multiplication with array failed'
    # Should be successful when the array shapes match
    spec.multiply(np.tile(a, (2,1)).T)
    assert np.all(spec.flux == 4.), 'Vector multiplication with array failed'
    # Should fault when the arrays are transposed (broadcast fails)
    with pytest.raises(PypeItError):
        spec.multiply(np.tile(a, (2,1)))

    # ND
    spec = spectrum.Spectrum(wave, np.tile(flux, (3,2,1)).T)
    # Should be successful when the array can be broadcast to the right shape
    spec.multiply(a)
    assert np.all(spec.flux == 2.), 'Vector multiplication with array failed'
    spec.multiply(np.tile(a, (2,1)).T)
    assert np.all(spec.flux == 4.), 'Vector multiplication with array failed'
    # Should fault if the order of the shapes is not correct (expand_dims must
    # follow a specific order)
    with pytest.raises(PypeItError):
        spec.multiply(np.tile(a, (3,1)).T)

    # TODO: Check the error and mask propagation
    # TODO: Test multiplying two spectra


def test_resample():

    wave, flux, err, gpm = get_fake_spec()
    ivar = 1/err**2
    new_wave = np.linspace(4220., 5930., 1078, dtype=float)

    # Test without error
    spec = spectrum.Spectrum(wave, flux)
    _spec = spec.resample(new_wave)
    assert np.array_equal(new_wave, _spec.wave), 'Resampled wavelength vector incorrect'
    diff = interpolate.interp1d(spec.wave, spec.flux)(_spec.wave) - _spec.flux
    assert np.std(diff) < 0.06, 'Difference is too large'
    assert _spec.ivar is None, 'Should not include errors'

    # Test with error
    spec = spectrum.Spectrum(wave, flux, ivar=ivar)
    _spec = spec.resample(new_wave)
    assert np.array_equal(new_wave, _spec.wave), 'Resampled wavelength vector incorrect'
    diff = interpolate.interp1d(spec.wave, spec.flux)(_spec.wave) - _spec.flux
    assert np.std(diff) < 0.06, 'Difference is too large'
    assert np.mean(_spec.ivar) > np.mean(spec.ivar), \
        'New spectrum averages pixels so inverse variance should be larger'
    
    # Test 2D
    nspec = 10
    flux2D = np.tile(flux, (nspec,1)).T
    ivar2D = np.tile(ivar, (nspec,1)).T
    spec = spectrum.Spectrum(wave, flux2D, ivar=ivar2D)
    __spec = spec.resample(new_wave)
    assert np.array_equal(_spec.flux, __spec.flux[:,0]), 'Resampling 2D should be same as 1D'

    # Test mask
    new_wave = np.linspace(3500., 6500., flux.size+1000, dtype=float)
    spec = spectrum.Spectrum(wave, flux, ivar=ivar, gpm=gpm)
    _spec = spec.resample(new_wave)
    # Subtle issues with resampling mean that
    #   - the last pixel of the original spectrum is masked in the output
    #     spectrum and
    #   - the flux values are not identical.
    assert np.std(spec.flux[spec.gpm][:-1] - _spec.flux[_spec.gpm]) < 1e-4, 'fluxes too different'
    indx = (_spec.wave < spec.wave[0]) | (_spec.wave > spec.wave[-1])
    assert not np.any(indx & _spec.gpm), 'pixels outside original wavelength range should be masked'
