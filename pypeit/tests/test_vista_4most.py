"""Unit tests for the experimental 4MOST commissioning adapter."""

from types import SimpleNamespace

import numpy as np
import pytest

from pypeit import PypeItError
from pypeit.spectrographs.vista_4most import VISTA4MOSTLRSBRedSpectrograph


def _spectrum(fiber, name, column, n=400):
    return SimpleNamespace(
        MASKDEF_ID=fiber, MASKDEF_OBJNAME=name, TRACE_SPAT=np.full(n, float(column)),
        OPT_WAVE=np.linspace(8440., 8859., n), BOX_COUNTS=np.full(n, 1500.),
        OPT_COUNTS=np.full(n, 1400.), BOX_MASK=np.ones(n, dtype=bool),
        OPT_MASK=np.ones(n, dtype=bool), OPT_COUNTS_SIG=np.full(n, 3.),
        OPT_COUNTS_SIG_DET=np.full(n, 2.), OPT_COUNTS_SKY=np.full(n, 80.),
        OPT_COUNTS_IVAR=np.full(n, 1/9), OPT_COUNTS_NIVAR=np.full(n, 1/4)
    )


def test_aperture_calibration_selects_only_bright_isolated_science(monkeypatch):
    monkeypatch.setenv('PYPEIT_4MOST_APERTURE_SCALE', '1.01')
    spec = VISTA4MOSTLRSBRedSpectrograph()
    isolated = _spectrum(1, 'Object', 8)
    bright_left = _spectrum(2, 'Object', 20)
    bright_right = _spectrum(3, 'Object', 26)
    faint = _spectrum(4, 'Object', 40)
    sky = _spectrum(5, 'SKY', 52)
    image = np.zeros((400, 60))
    image[:, [8, 20, 26, 40, 52]] = [1600., 1800., 1500., 200., 1000.]

    spec.calibrate_fiber_aperture([isolated, bright_left, bright_right, faint, sky], image)

    factor = 1.01 * 1500. / 1400.
    np.testing.assert_allclose(isolated.OPT_COUNTS, 1515.)
    np.testing.assert_allclose(isolated.OPT_COUNTS_IVAR, (1/9) / factor**2)
    np.testing.assert_allclose(isolated.OPT_COUNTS_NIVAR, (1/4) / factor**2)
    np.testing.assert_allclose(isolated.OPT_COUNTS_SIG, 3. * factor)
    np.testing.assert_allclose(isolated.OPT_COUNTS_SIG_DET, 2. * factor)
    np.testing.assert_allclose(isolated.OPT_COUNTS_SKY, 80. * factor)
    np.testing.assert_array_equal(isolated.BOX_COUNTS, np.full(400, 1500.))
    np.testing.assert_array_equal(isolated.OPT_MASK, np.ones(400, dtype=bool))
    np.testing.assert_array_equal(isolated.OPT_WAVE, np.linspace(8440., 8859., 400))
    for obj in (bright_left, bright_right, faint, sky):
        np.testing.assert_array_equal(obj.OPT_COUNTS, np.full(400, 1400.))


def test_aperture_calibration_disabled_without_opt_in(monkeypatch):
    monkeypatch.delenv('PYPEIT_4MOST_APERTURE_SCALE', raising=False)
    spec = VISTA4MOSTLRSBRedSpectrograph()
    obj = _spectrum(1, 'Object', 8)
    spec.calibrate_fiber_aperture([obj], np.full((400, 20), 1000.))
    np.testing.assert_array_equal(obj.OPT_COUNTS, np.full(400, 1400.))


@pytest.mark.parametrize('scale', ('0', '-1', 'nan', 'inf', 'invalid'))
def test_aperture_calibration_rejects_invalid_scale(monkeypatch, scale):
    monkeypatch.setenv('PYPEIT_4MOST_APERTURE_SCALE', scale)
    with pytest.raises(PypeItError, match='positive number'):
        VISTA4MOSTLRSBRedSpectrograph().aperture_scale
