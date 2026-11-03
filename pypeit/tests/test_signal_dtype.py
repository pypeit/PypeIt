"""Regression tests for SciPy signal calls with FITS byte-order arrays."""

import numpy as np
from scipy import signal

from pypeit.core import telluric
from pypeit.core.wavecal import wvutils


def test_xcorr_shift_normalizes_byte_order(monkeypatch):
    pixels = np.arange(128)
    spectrum = np.exp(-0.5 * ((pixels - 64) / 4)**2).astype(np.float32)
    expected = wvutils.xcorr_shift(spectrum, spectrum, do_xcorr_arc=False,
                                   cont_subtract=False)

    correlate = signal.correlate

    def check_correlate(first, second, **kwargs):
        assert first.dtype.isnative and second.dtype.isnative
        return correlate(first, second, **kwargs)

    monkeypatch.setattr(signal, 'correlate', check_correlate)
    swapped = spectrum.astype(spectrum.dtype.newbyteorder('S'))
    result = wvutils.xcorr_shift(swapped, swapped, do_xcorr_arc=False,
                                cont_subtract=False)
    np.testing.assert_allclose(result, expected)



def test_xcorr_shift_stretch_normalizes_byte_order(monkeypatch):
    pixels = np.arange(256)
    spectrum = np.zeros(pixels.size, dtype=float)
    for cen, amp in zip([40.3, 97.8, 151.2, 203.6], [1.0, 3.0, 2.0, 1.5]):
        spectrum += amp * np.exp(-0.5 * ((pixels - cen) / 1.5)**2)
    spectrum = (spectrum + 1.).astype(np.float32)
    expected = wvutils.xcorr_shift_stretch(spectrum, spectrum, fwhm=3.5,
                                           stretch_func='linear')

    correlate = signal.correlate

    def check_correlate(first, second, **kwargs):
        assert first.dtype.isnative and second.dtype.isnative
        return correlate(first, second, **kwargs)

    monkeypatch.setattr(signal, 'correlate', check_correlate)
    swapped = spectrum.astype(spectrum.dtype.newbyteorder('S'))
    result = wvutils.xcorr_shift_stretch(swapped, swapped, fwhm=3.5,
                                         stretch_func='linear')
    assert result[0] == 1
    np.testing.assert_allclose(result[1:], expected[1:])

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
