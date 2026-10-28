"""Regression tests for SciPy signal calls with FITS byte-order arrays."""

import numpy as np

from pypeit.core import telluric
from pypeit.core.wavecal import wvutils


def test_xcorr_shift_normalizes_byte_order(monkeypatch):
    pixels = np.arange(128)
    spectrum = np.exp(-0.5 * ((pixels - 64) / 4)**2).astype(np.float32)
    expected = wvutils.xcorr_shift(spectrum, spectrum, do_xcorr_arc=False,
                                   cont_subtract=False)

    correlate = wvutils.scipy.signal.correlate

    def check_correlate(first, second, **kwargs):
        assert first.dtype.isnative and second.dtype.isnative
        return correlate(first, second, **kwargs)

    monkeypatch.setattr(wvutils.scipy.signal, 'correlate', check_correlate)
    swapped = spectrum.astype(spectrum.dtype.newbyteorder('S'))
    result = wvutils.xcorr_shift(swapped, swapped, do_xcorr_arc=False,
                                cont_subtract=False)
    np.testing.assert_allclose(result, expected)


def test_conv_telluric_normalizes_byte_order(monkeypatch):
    spectrum = np.exp(-0.5 * ((np.arange(128) - 64) / 8)**2)
    expected = telluric.conv_telluric(spectrum, 1e-5, 20000)

    convolve = telluric.scipy.signal.convolve

    def check_convolve(first, second, **kwargs):
        assert first.dtype.isnative and second.dtype.isnative
        return convolve(first, second, **kwargs)

    monkeypatch.setattr(telluric.scipy.signal, 'convolve', check_convolve)
    swapped = spectrum.astype(spectrum.dtype.newbyteorder('S'))
    result = telluric.conv_telluric(swapped, 1e-5, 20000)
    np.testing.assert_allclose(result, expected)
