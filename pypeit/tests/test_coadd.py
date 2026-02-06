
from IPython import embed
import numpy as np
import pytest

from pypeit import PypeItError
from pypeit.core import coadd
from pypeit.core import standard
from pypeit.core import meta

def test_robust_median():
    # Build some test spectra
    ra, dec = meta.convert_radec('00:31:18.49', '-43:36:23')
    std_spec = standard.get_standard_spectrum(ra=ra, dec=dec)
    indx = (std_spec.wave > 3100.) & (std_spec.wave < 9000.)
    wave = std_spec.wave[indx]
    flux_ref = std_spec.flux[indx]
    ivar_ref = (100/flux_ref)**2
    assert np.all(np.isfinite(ivar_ref)), 'Spectrum changed'
    mask_ref = ivar_ref > 0.

    # Add offset and noise
    known_ratio = 2.
    flux = flux_ref / known_ratio
    rng = np.random.default_rng(99)
    sig = np.median(flux) / 10.
    flux += rng.normal(size=flux.size) * sig
    ivar = np.full(flux.size, 1/sig**2, dtype=float)
    mask = ivar > 0.

    ratio = coadd.robust_median_ratio(flux, ivar, flux_ref, ivar_ref, mask=mask, mask_ref=mask_ref)
    assert np.absolute(known_ratio / ratio - 1) < 0.05, 'Median ratio incorrect'

    # Test when mask_ref isn't passed
    ratio = coadd.robust_median_ratio(flux, ivar, flux_ref, ivar_ref, mask=mask)
    assert np.absolute(known_ratio / ratio - 1) < 0.05, 'Median ratio incorrect'

    # Test when mask isn't passed
    ratio = coadd.robust_median_ratio(flux, ivar, flux_ref, ivar_ref, mask_ref=mask_ref)
    assert np.absolute(known_ratio / ratio - 1) < 0.05, 'Median ratio incorrect'

    # Set a high do_not_rescale value
    ratio = coadd.robust_median_ratio(flux, ivar, flux_ref, ivar_ref, snr_do_not_rescale=15.)
    assert ratio == 1., 'S/N should be too low, yielding ratio = 1'
    
    # Set a high min_good
    mask[mask.size//2:] = False
    ratio = coadd.robust_median_ratio(flux, ivar, flux_ref, ivar_ref, mask=mask, min_good=0.6)
    assert ratio == 1., 'Should be too few good pixels, yielding ratio = 1'


def test_median_filt_spec():

    # Read in an example spectrum
    ra, dec = meta.convert_radec('00:31:18.49', '-43:36:23')
    std_spec = standard.get_standard_spectrum(ra=ra, dec=dec)
    indx = (std_spec.wave > 3100.) & (std_spec.wave < 9000.)
    wave = std_spec.wave[indx]
    flux = std_spec.flux[indx]
    ivar = (100/flux)**2
    assert np.all(np.isfinite(ivar)), 'Spectrum changed'
    mask = ivar > 0.

    flux_med, ivar_med = coadd.median_filt_spec(flux, ivar, mask, 31)
    assert np.absolute(np.median(flux - flux_med)) < 1e-10, 'Median of the difference should be 0.'


def test_poly_model_eval():
    wave = np.arange(3000., 9000., 2.)
    # Coefficients for a third order function
    rng = np.random.default_rng(99)
    theta = rng.uniform(size=4) / np.arange(1,5)

    poly_p = coadd.poly_model_eval(theta, 'legendre', 'poly', wave, wave[0], wave[-1])
    assert np.absolute(np.mean(poly_p) - 0.5) < 0.01, 'poly model changed'
    poly_s = coadd.poly_model_eval(theta, 'legendre', 'square', wave, wave[0], wave[-1])
    assert np.absolute(np.mean(poly_s) - 0.3) < 0.01, 'square model changed'
    poly_e = coadd.poly_model_eval(theta, 'legendre', 'exp', wave, wave[0], wave[-1])
    assert np.absolute(np.mean(poly_e) - 1.7) < 0.01, 'exp model changed'
    with pytest.raises(PypeItError):
        poly_r = coadd.poly_model_eval(theta, 'legendre', 'junk', wave, wave[0], wave[-1])


def test_solve_poly_ratio():
    # Build some test spectra
    ra, dec = meta.convert_radec('00:31:18.49', '-43:36:23')
    std_spec = standard.get_standard_spectrum(ra=ra, dec=dec)
    indx = (std_spec.wave > 3100.) & (std_spec.wave < 9000.)
    wave = std_spec.wave[indx]
    flux_ref = std_spec.flux[indx]
    ivar_ref = (100/flux_ref)**2
    assert np.all(np.isfinite(ivar_ref)), 'Spectrum changed'
    mask_ref = ivar_ref > 0.

    # Coefficients for a third order function
    rng = np.random.default_rng(99)
    norder = 3
    theta = rng.uniform(size=norder+1) / np.arange(1,norder+2)

    # "True" ratio
    ratio = coadd.poly_model_eval(theta, 'legendre', 'square', wave, wave[0], wave[-1])

    # Add polynomial and noise
    flux = flux_ref / ratio
    rng = np.random.default_rng(99)
    sig = np.median(flux) / 10.
    flux += rng.normal(size=flux.size) * sig
    ivar = np.full(flux.size, 1/sig**2, dtype=float)
    mask = ivar > 0.

    # Check that it fails when norder=0
    with pytest.raises(PypeItError):
        coadd.solve_poly_ratio(wave, flux, ivar, flux_ref, ivar_ref, 0)

    # Check that it fails when model is not recognized
    with pytest.raises(PypeItError):
        coadd.solve_poly_ratio(wave, flux, ivar, flux_ref, ivar_ref, norder, model='junk')

    # Perform the fit
    scale_min = 0.05
    _ratio, fit_tuple, flux_rescale, ivar_rescale, outmask = coadd.solve_poly_ratio(
        wave, flux, ivar, flux_ref, ivar_ref, norder, mask=mask, mask_ref=mask_ref,
        scale_min=scale_min
    )
    assert np.all(np.absolute(fit_tuple[0]/theta - 1) < 0.15), \
        'Coefficients should match input to within 15%.'
    assert np.all(_ratio >= scale_min), 'scale_min not imposed correctly'

    scale_min = 1e-10
    _ratio, fit_tuple, flux_rescale, ivar_rescale, outmask = coadd.solve_poly_ratio(
        wave, flux, ivar, flux_ref, ivar_ref, norder, mask=mask, mask_ref=mask_ref,
        scale_min=scale_min
    )
    assert abs(np.median(flux_rescale / flux_ref - 1)) < 2e-4, 'Fit not as good as it should be'
