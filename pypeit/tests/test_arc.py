"""
Module to run tests on ararclines
"""
import pytest
import numpy as np

from pypeit.core import arc
from pypeit.core import skyspec
from pypeit.core.wavecal import wvutils
from pypeit.core import tracewave
from pypeit.par.pypeitpar import PypeItPar, WavelengthSolutionPar


def test_detect_lines():
    # Using Paranal night sky as an 'arc'
    arx_sky = skyspec.load_sky_spectrum('paranal_sky.fits')
    arx_amp_true, arx_amp, arx_cent, arx_wid, arx_centerr, arx_w, arx_yprep, _ \
            = arc.detect_lines(arx_sky.flux)
    # The default rejects lines narrower than half the median width
    assert 3250 < len(arx_w) < 3350
    assert np.all(arx_wid[arx_w] > 0.5 * np.median(arx_wid[arx_w]) * 0.9)
    # Without the floor (previous behaviour)
    old = arc.detect_lines(arx_sky.flux, min_frac_med_width=0.)
    assert len(old[5]) > 3700


def _synthetic_arc(nspec=4096, spacing=20., fwhm=4.0, noise=1.0, seed=42):
    """
    Isolated Gaussian lines on a regular grid with random offsets, amplitudes,
    and widths, plus Gaussian noise.
    """
    rng = np.random.default_rng(seed)
    x = np.arange(nspec, dtype=float)
    cen = np.arange(20., nspec - 20., spacing)
    cen += rng.uniform(-2., 2., cen.size)
    amp = 10.**rng.uniform(1.5, 3.5, cen.size)
    sig = fwhm / 2.355 * rng.uniform(0.9, 1.1, cen.size)
    y = np.sum(amp[:,None] * np.exp(-0.5 * ((x[None,:] - cen[:,None]) / sig[:,None])**2),
               axis=0)
    y += rng.normal(scale=noise, size=nspec)
    return x, y, amp, cen, sig


def test_fit_arcspec_recovers_truth():
    x, y, amp, cen, sig = _synthetic_arc(nspec=1024)
    pixt = np.round(cen).astype(int)
    a, c, w, ce = arc.fit_arcspec(x, y, pixt, 5)
    good = c > 0
    assert good.sum() > 0.95 * cen.size
    # The faintest lines have S/N ~ 30
    assert np.median(np.absolute(c[good] - cen[good])) < 0.01
    assert np.max(np.absolute(c[good] - cen[good])) < 0.1
    assert np.median(np.absolute(w[good] / sig[good] - 1.)) < 0.02
    assert np.max(np.absolute(w[good] / sig[good] - 1.)) < 0.1
    assert np.max(np.absolute(a[good] / amp[good] - 1.)) < 0.1
    assert np.all(ce[good] > 0)


def test_fit_arcspec_edge_cases():
    nspec = 64
    x = np.arange(nspec, dtype=float)
    y = np.zeros(nspec)
    y[32] = 100.
    y[:8] = -5.                                    # all-negative window
    pixt = np.array([0, 2, 32, 63])
    a, c, w, ce = arc.fit_arcspec(x, y, pixt, 5)
    for v in [a, c, w, ce]:
        assert not np.any(np.isnan(v))             # no NaN leaks to the caller
    # Failed fits are flagged with the -999 sentinel
    assert np.all(c[[0, 1, 3]] == -999.)

    # Lines at both ends of the spectrum, with truncated windows
    y = 100. * np.exp(-0.5 * ((x - 0.4) / 1.7)**2) \
            + 100. * np.exp(-0.5 * ((x - (nspec - 1.2)) / 1.7)**2)
    pixt = np.array([0, 1, nspec - 2, nspec - 1])
    a, c, w, ce = arc.fit_arcspec(x, y, pixt, 5)
    assert np.allclose(c, [0.4, 0.4, nspec - 1.2, nspec - 1.2], atol=1e-4)
    assert np.allclose(w, 1.7, rtol=1e-4)

    # A NaN in the window does not raise (curve_fit would); the line is
    # flagged as a failed fit
    y = 100. * np.exp(-0.5 * ((x - 30.3) / 1.7)**2)
    y[29] = np.nan
    a, c, w, ce = arc.fit_arcspec(x, y, np.array([30, 50]), 5)
    assert c[0] == -999.
    assert c[1] != -999.

    # No lines
    a, c, w, ce = arc.fit_arcspec(x, y, np.array([], dtype=int), 5)
    assert all(v.size == 0 for v in [a, c, w, ce])


def _arc_with_spikes(nspec=2048, seed=3):
    """A synthetic arc (sigma ~1.7 px) plus narrow spikes."""
    x, y, amp, cen, sig = _synthetic_arc(nspec=nspec, spacing=40., seed=seed)
    # Spikes (sigma = 0.4 px) half-way between lines, bright enough to be
    # detected
    spikes = 0.5 * (cen[1:] + cen[:-1])[::3]
    y += np.sum(300. * np.exp(-0.5 * ((x[None,:] - spikes[:,None]) / 0.4)**2), axis=0)
    return y, cen, spikes


def _found(cent, good, pos, tol=1.):
    return np.array([np.any(np.absolute(cent[good] - p) < tol) for p in pos])


def test_detect_lines_rejects_narrow_features():
    """By default, features narrower than half the median width are rejected."""
    y, cen, spikes = _arc_with_spikes()
    _, _, cent, wid, _, good, _, _ = arc.detect_lines(y, fwhm=4.)
    assert np.all(_found(cent, good, cen))                 # every arc line is kept
    assert not np.any(_found(cent, good, spikes))           # no spike survives
    assert np.all(wid[good] > 0.5 * np.median(wid[good]))

    # Without the floor, most spikes are accepted (curve_fit fails on a few)
    _, _, cent, wid, _, good, _, _ = arc.detect_lines(y, fwhm=4., min_frac_med_width=0.)
    assert np.sum(_found(cent, good, spikes)) > 0.5 * spikes.size


def test_detect_lines_floor_needs_enough_lines():
    """The median-width floor is not applied to spectra with very few lines."""
    x = np.arange(256, dtype=float)
    rng = np.random.default_rng(1)
    # Three broad lines and one narrow one
    y = sum(a * np.exp(-0.5 * ((x - c) / s)**2) for a, c, s in
            [(500., 40., 3.), (400., 100., 3.), (600., 160., 3.), (500., 220., 1.)])
    y += rng.normal(scale=1., size=x.size)
    _, _, cent, wid, _, good, _, _ = arc.detect_lines(y, fwhm=7., cont_subtract=False,
                                                      input_thresh=50.)
    assert good.size < 5
    assert np.any(np.absolute(cent[good] - 220.) < 1.)


def test_detect_lines_floor_explicit_value():
    """An explicit min_frac_med_width sets the floor."""
    y, cen, spikes = _arc_with_spikes()
    _, _, cent, wid, _, good, _, _ = arc.detect_lines(y, fwhm=4., min_frac_med_width=0.)
    med = np.median(wid[good])
    for frac in [0.3, 0.5, 0.9]:
        _, _, c, w, _, g, _, _ = arc.detect_lines(y, fwhm=4., min_frac_med_width=frac)
        assert np.all(w[g] > frac * med)
        assert np.array_equal(g, good[w[good] > frac * med])


def test_get_xcorr_arc_no_floor(monkeypatch):
    """The synthetic cross-correlation arc is built without the width floor."""
    y, cen, spikes = _arc_with_spikes()
    seen = []
    _detect_lines = arc.detect_lines
    def _spy(*a, **k):
        seen.append(k.get('min_frac_med_width'))
        return _detect_lines(*a, **k)
    monkeypatch.setattr(arc, 'detect_lines', _spy)
    wvutils.get_xcorr_arc(y, fwhm=4.)
    assert seen == [0.]

def test_med_width_pars():
    """The median-width floor is set by WavelengthSolutionPar."""
    par = WavelengthSolutionPar()
    assert par['min_frac_med_width'] == 0.5
    assert par['min_nlines_med_width'] == 5
    assert wvutils.med_width_kwargs(par) == dict(min_frac_med_width=0.5, min_nlines_med_width=5)
    assert wvutils.med_width_kwargs(None) == {}
    cfg = ['[calibrations]', '    [[wavelengths]]', '        min_frac_med_width = 0.',
           '        min_nlines_med_width = 10']
    par = PypeItPar.from_cfg_lines(merge_with=(cfg,))['calibrations']['wavelengths']
    assert par['min_frac_med_width'] == 0.
    assert par['min_nlines_med_width'] == 10
    with pytest.raises(ValueError):
        WavelengthSolutionPar(min_frac_med_width=-1.)
    with pytest.raises(ValueError):
        WavelengthSolutionPar(min_nlines_med_width=0)


def test_tilts_find_lines_med_width():
    """The tilt line detection applies the median-width floor it is given."""
    y, cen, spikes = _arc_with_spikes()
    slit_cen = np.full(y.size, 50.)
    nspike = []
    for frac in [0.5, 0.]:
        lines_spec, _, good = tracewave.tilts_find_lines(y, slit_cen, tracethresh=5., fwhm=4.,
                                                         min_frac_med_width=frac)
        nspike.append(np.sum(_found(lines_spec, good, spikes)))
    assert nspike[0] == 0
    assert nspike[1] > 0

# Many more functions in pypeit.core.arc that need tests!
