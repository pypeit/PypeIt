import numpy as np
import pytest

from astropy.table import Table

from pypeit import PypeItError
from pypeit.core.wavecal import autoid
from pypeit.core.wavecal import wvutils
from pypeit.par import pypeitpar


def test_reidentify_xcorr_only_uses_best_archive_and_forces_linear(monkeypatch):
    nspec = 50
    spec = np.ones(nspec)
    spec_arxiv = np.column_stack((np.full(nspec, 2.0), np.full(nspec, 3.0)))
    pixel = np.arange(nspec, dtype=float)
    wave_arxiv = np.column_stack((
        5000.0 + 2.0*pixel,
        6000.0 + 3.0*pixel + 0.01*pixel**2
    ))

    def fake_arc_lines_from_spec(input_spec, **kwargs):
        empty = np.array([], dtype=float)
        return empty, empty, empty, np.array([], dtype=int), input_spec.copy()

    calls = []

    def fake_xcorr_shift_stretch(input_spec, archive_spec, **kwargs):
        calls.append(kwargs['stretch_func'])
        if archive_spec[0] == 2.0:
            return 1, 4.0, 0.99, 0.0, 0.8, 4.0, 0.8
        return 1, -2.0, 1.01, 0.0, 0.9, -2.0, 0.9

    monkeypatch.setattr(wvutils, 'arc_lines_from_spec', fake_arc_lines_from_spec)
    monkeypatch.setattr(wvutils, 'xcorr_shift_stretch', fake_xcorr_shift_stretch)

    _, _, patt_dict = autoid.reidentify(
        spec, spec_arxiv, wave_arxiv, Table({'wave': [5000.0]}), 1,
        xcorr_only=True, stretch_func='quadratic')

    expected = wvutils.shift_and_stretch(
        wave_arxiv[:, 1], -2.0, 1.01, 0.0, stretch_func='linear')

    assert calls == ['linear', 'linear']
    assert patt_dict['acceptable']
    assert patt_dict['xcorr_only']
    assert patt_dict['ibest'] == 1
    assert patt_dict['xcorr_cc'] == 0.9
    assert np.array_equal(patt_dict['xcorr_valid'], expected > 1.0)
    assert np.allclose(patt_dict['xcorr_wave'], expected)


def test_fit_xcorr_wave_solution():
    nspec = 100
    pixel = np.arange(nspec, dtype=float)
    wave = 5000.0 + 2.0*pixel + 0.002*pixel**2
    patt_dict = {
        'xcorr_wave': wave,
        'xcorr_valid': np.ones(nspec, dtype=bool),
        'xcorr_shift': 3.5
    }

    wave_fit = autoid.fit_xcorr_wave_solution(
        np.ones(nspec), patt_dict, func='legendre', order=2)

    assert wave_fit is not None
    assert wave_fit.pypeitfit is not None
    assert wave_fit.shift == 3.5
    assert wave_fit.rms < 1e-8
    assert np.allclose(wave_fit.pypeitfit.eval(pixel/(nspec-1)), wave)


def test_echelle_rejects_xcorr_only():
    par = pypeitpar.WavelengthSolutionPar(echelle=True)
    par['xcorr_only'] = True

    with pytest.raises(PypeItError, match='not available for echelle'):
        autoid.echelle_wvcalib(
            np.ones((10, 1)), np.array([1]), np.ones((10, 1)),
            np.ones((10, 1)), [], par)


def test_full_template_xcorr_only(monkeypatch):
    nspec = 20
    spec = np.ones((nspec, 1))
    wave = 5000.0 + 2.0*np.arange(nspec)
    template_dict = {
        'wave': wave,
        'spec': np.ones(nspec),
        'bin': 1,
        'order': None,
        'lines_pix': [np.array([5.0, 10.0])],
        'lines_wav': [wave[[5, 10]]],
        'lines_fit_ord': [1]
    }
    par = pypeitpar.WavelengthSolutionPar(
        method='full_template', xcorr_only=True, fwhm_fromlines=False)

    monkeypatch.setattr(
        autoid.waveio, 'load_line_lists',
        lambda *args, **kwargs: (Table({'wave': [5000.0]}), None, None))
    monkeypatch.setattr(
        wvutils, 'arc_lines_from_spec',
        lambda input_spec, **kwargs: (
            np.array([]), np.array([]), np.array([]), np.array([], dtype=int),
            input_spec.copy()))
    monkeypatch.setattr(wvutils, 'xcorr_shift', lambda *args, **kwargs: (0.0, 0.9))

    calls = {}

    def fake_reidentify(input_spec, archive_spec, archive_wave, *args, **kwargs):
        calls['reidentify'] = kwargs
        return np.array([]), input_spec, {
            'acceptable': True,
            'xcorr_wave': archive_wave,
            'xcorr_valid': archive_wave > 1.0,
            'xcorr_shift': 0.0
        }

    expected = {'sentinel': True}

    def fake_fit(input_spec, patt_dict, **kwargs):
        calls['fit'] = kwargs
        return expected

    monkeypatch.setattr(autoid, 'reidentify', fake_reidentify)
    monkeypatch.setattr(autoid, 'fit_xcorr_wave_solution', fake_fit)

    result, order = autoid.full_template(
        spec, [], par, np.array([0]), 1, 1, measured_fwhms=np.array([4.0]),
        template_dict=template_dict)

    assert result['0'] == expected
    assert order is None
    assert calls['reidentify']['xcorr_only']
    assert calls['reidentify']['stretch_func'] == 'linear'


def test_full_template_rejects_echelle_xcorr_only():
    par = pypeitpar.WavelengthSolutionPar(echelle=True)
    par['xcorr_only'] = True

    with pytest.raises(PypeItError, match='not available for echelle'):
        autoid.full_template(
            np.ones((10, 1)), [], par, np.array([0]), 1, 1,
            measured_fwhms=np.array([4.0]))
