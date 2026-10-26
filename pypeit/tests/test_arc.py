"""
Module to run tests on ararclines
"""
import pytest
import numpy as np

from pypeit.core import arc
from pypeit.core import fitting
from pypeit.core import skyspec
from pypeit.par.pypeitpar import PypeItPar, WavelengthSolutionPar
from pypeit.wavecalib import BuildWaveCalib


def test_detect_lines():
    # Using Paranal night sky as an 'arc'
    arx_sky = skyspec.load_sky_spectrum('paranal_sky.fits')
    arx_amp_true, arx_amp, arx_cent, arx_wid, arx_centerr, arx_w, arx_yprep, _ \
            = arc.detect_lines(arx_sky.flux)
    assert (len(arx_w) > 3275)


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


def _legacy_fit(monkeypatch, *args):
    """Run fit_arcspec with the per-line curve_fit implementation."""
    with monkeypatch.context() as m:
        m.setattr(arc, 'USE_VECTORIZED_ARCFIT', False)
        return arc.fit_arcspec(*args)


def test_fit_arcspec_recovers_truth():
    x, y, amp, cen, sig = _synthetic_arc()
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


def test_fit_arcspec_matches_curve_fit(monkeypatch):
    """The vectorized fit must agree with the legacy per-line curve_fit."""
    x, y, _, cen, _ = _synthetic_arc()
    pixt = np.round(cen).astype(int)
    new = arc.fit_arcspec(x, y, pixt, 5)
    old = _legacy_fit(monkeypatch, x, y, pixt, 5)
    good = (new[1] > 0) & (old[1] > 0)
    assert good.sum() > 0.95 * cen.size
    assert np.max(np.absolute(new[1][good] - old[1][good])) < 0.02   # centroid, pixels
    assert np.max(np.absolute(new[2][good] / old[2][good] - 1.)) < 0.02
    assert np.max(np.absolute(new[0][good] / old[0][good] - 1.)) < 0.02
    # Same definition of the centroid variance as curve_fit
    assert np.max(np.absolute(new[3][good] / old[3][good] - 1.)) < 0.02


def test_fit_arcspec_analytic_only(monkeypatch):
    """The pure analytic log-parabola fit is close to curve_fit on clean data."""
    x, y, _, cen, _ = _synthetic_arc()
    pixt = np.round(cen).astype(int)
    old = _legacy_fit(monkeypatch, x, y, pixt, 5)
    monkeypatch.setattr(arc, 'ARCFIT_GN_ITER', 0)
    monkeypatch.setattr(arc, 'ARCFIT_RESCUE_ITER', 0)
    new = arc.fit_arcspec(x, y, pixt, 5)
    good = (new[1] > 0) & (old[1] > 0)
    assert good.sum() > 0.95 * cen.size
    assert np.median(np.absolute(new[1][good] - old[1][good])) < 0.005
    assert np.max(np.absolute(new[1][good] - old[1][good])) < 0.05
    assert np.max(np.absolute(new[2][good] / old[2][good] - 1.)) < 0.1


def test_fit_arcspec_kill_switch(monkeypatch):
    """With the module switch off, fit_arcspec is the per-line fit."""
    x, y, _, cen, _ = _synthetic_arc(nspec=512)
    pixt = np.round(cen).astype(int)
    old = _legacy_fit(monkeypatch, x, y, pixt, 5)
    # fitp = 5 gives a fit half-width of 3 pixels
    ref = [np.full(pixt.size, -999.) for _ in range(4)]
    arc._fit_arcspec_curvefit(x, y, pixt, 3, y.size, np.arange(pixt.size), *ref)
    for o, r in zip(old, ref):
        assert np.array_equal(o, r)


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
    assert c[2] == pytest.approx(32., abs=0.5)

    # Lines at both ends of the spectrum, with truncated windows
    y = 100. * np.exp(-0.5 * ((x - 0.4) / 1.7)**2) \
            + 100. * np.exp(-0.5 * ((x - (nspec - 1.2)) / 1.7)**2)
    pixt = np.array([0, 1, nspec - 2, nspec - 1])
    a, c, w, ce = arc.fit_arcspec(x, y, pixt, 5)
    assert np.allclose(c, [0.4, 0.4, nspec - 1.2, nspec - 1.2], atol=1e-4)
    assert np.allclose(w, 1.7, rtol=1e-4)

    # A NaN in the window does not raise (curve_fit would), and is ignored
    y = 100. * np.exp(-0.5 * ((x - 30.3) / 1.7)**2)
    y[29] = np.nan
    a, c, w, ce = arc.fit_arcspec(x, y, np.array([30]), 5)
    assert c[0] == pytest.approx(30.3, abs=1e-4)

    # No lines
    a, c, w, ce = arc.fit_arcspec(x, y, np.array([], dtype=int), 5)
    assert all(v.size == 0 for v in [a, c, w, ce])


def test_fit_arcspec_negative_wings(monkeypatch):
    """Negative pixels in the wings are used, as in the per-line fit."""
    nspec = 64
    x = np.arange(nspec, dtype=float)
    y = 100. * np.exp(-0.5 * ((x - 30.3) / 1.7)**2) - 30.
    # The outermost pixels of the 7-pixel window are negative
    assert np.sum(y[27:34] < 0) == 2
    pixt = np.array([30])
    def _boom(*a, **k):
        raise AssertionError('fit_gauss should not be called')
    # The model cannot match the negative pedestal, so the iterations converge
    # only slowly; test the vectorized fit itself, without handing the line to
    # curve_fit
    with monkeypatch.context() as m:
        m.setattr(fitting, 'fit_gauss', _boom)
        m.setattr(arc, 'ARCFIT_MATCH_CURVEFIT', False)
        new = arc.fit_arcspec(x, y, pixt, 5)
    old = _legacy_fit(monkeypatch, x, y, pixt, 5)
    for n, o in zip(new[:3], old[:3]):
        assert n[0] == pytest.approx(o[0], rel=1e-6)


def test_fit_arcspec_narrow_features(monkeypatch):
    """Features much narrower than an arc line are left to curve_fit."""
    x, y, _, cen, _ = _synthetic_arc(nspec=512)
    # A hot pixel between two lines
    hot = int(np.round(cen[5] + cen[6])) // 2
    y[hot] += 500.
    pixt = np.append(np.round(cen).astype(int), hot)
    ncalls = [0]
    _fit_gauss = fitting.fit_gauss
    def _count(*a, **k):
        ncalls[0] += 1
        return _fit_gauss(*a, **k)
    monkeypatch.setattr(fitting, 'fit_gauss', _count)
    new = arc.fit_arcspec(x, y, pixt, 5)
    assert ncalls[0] >= 1
    old = _legacy_fit(monkeypatch, x, y, pixt, 5)
    # The hot pixel gets exactly the curve_fit result
    for n, o in zip(new, old):
        assert n[-1] == o[-1]
    # Without the hand-off, the vectorized fit returns a very narrow line
    ncalls[0] = 0
    monkeypatch.setattr(arc, 'ARCFIT_MATCH_CURVEFIT', False)
    fast = arc.fit_arcspec(x, y, pixt, 5)
    assert ncalls[0] == 0
    assert 0 < fast[2][-1] < arc.ARCFIT_MIN_WIDTH * 3


def test_detect_lines_matches_curve_fit(monkeypatch):
    """With the default hand-off, the accepted lines are those of curve_fit."""
    flux = skyspec.load_sky_spectrum('paranal_sky.fits').flux
    new = arc.detect_lines(flux)
    with monkeypatch.context() as m:
        m.setattr(arc, 'USE_VECTORIZED_ARCFIT', False)
        old = arc.detect_lines(flux)
    assert np.array_equal(new[5], old[5])
    assert np.max(np.absolute(new[2][new[5]] - old[2][old[5]])) < 0.02


def test_fit_arcspec_no_curve_fit_on_clean_data(monkeypatch):
    """On clean data the vectorized path must not call curve_fit at all."""
    x, y, _, cen, _ = _synthetic_arc()
    def _boom(*a, **k):
        raise AssertionError('fit_gauss should not be called')
    monkeypatch.setattr(fitting, 'fit_gauss', _boom)
    a, c, w, ce = arc.fit_arcspec(x, y, np.round(cen).astype(int), 5)
    assert np.all(c > 0)

def test_arcfit_options():
    """The context manager sets the option temporarily, even on error."""
    assert arc.ARCFIT_MATCH_CURVEFIT
    with arc.arcfit_options(match_curvefit=False):
        assert not arc.ARCFIT_MATCH_CURVEFIT
    assert arc.ARCFIT_MATCH_CURVEFIT
    with pytest.raises(ValueError):
        with arc.arcfit_options(match_curvefit=False):
            raise ValueError('test')
    assert arc.ARCFIT_MATCH_CURVEFIT
    with arc.arcfit_options():
        assert arc.ARCFIT_MATCH_CURVEFIT


def test_fast_line_fit_par():
    """The wavelength-solution parameter parses and defaults to the exact fit."""
    assert WavelengthSolutionPar()['fast_line_fit'] is False
    cfg = ['[calibrations]', '    [[wavelengths]]', '        fast_line_fit = True']
    par = PypeItPar.from_cfg_lines(merge_with=(cfg,))
    assert par['calibrations']['wavelengths']['fast_line_fit'] is True


def test_wavecalib_applies_fast_line_fit(monkeypatch):
    """BuildWaveCalib.run applies fast_line_fit to the line fits it performs."""
    seen = []
    def _extract_arcs(self):
        seen.append(arc.ARCFIT_MATCH_CURVEFIT)
        raise RuntimeError('stop')
    monkeypatch.setattr(BuildWaveCalib, 'extract_arcs', _extract_arcs)
    for fast in [True, False]:
        wvc = object.__new__(BuildWaveCalib)
        wvc.par = WavelengthSolutionPar(fast_line_fit=fast)
        with pytest.raises(RuntimeError):
            wvc.run()
        # Restored afterwards
        assert arc.ARCFIT_MATCH_CURVEFIT
    assert seen == [False, True]

# Many more functions in pypeit.core.arc that need tests!
