"""
Tests for Binospec overscan subtraction and nonlinearity correction.
"""
import numpy as np
import pytest

from pypeit.core.procimg import clean_overscan_vector
from pypeit.spectrographs.mmt_binospec import (binospec_read_amp,
                                                MMTBINOSPECSpectrograph)


class TestCleanOverscanVector:
    """Tests for clean_overscan_vector."""

    def test_clean_no_outliers(self):
        """A smooth vector should be returned unchanged."""
        rng = np.random.default_rng(42)
        vec = 1000.0 + rng.normal(scale=1.0, size=100)
        cleaned = clean_overscan_vector(vec)
        np.testing.assert_array_equal(cleaned, vec)

    def test_clean_single_outlier(self):
        """A single large outlier should be interpolated over."""
        vec = np.full(100, 1000.0)
        vec[50] = 2000.0  # outlier: deviation = 1000 >> nsig*rdnoise = 4.0
        cleaned = clean_overscan_vector(vec)
        # Outlier should be replaced with interpolated value (~1000)
        assert abs(cleaned[50] - 1000.0) < 1.0
        # Non-outlier pixels should be unchanged
        assert cleaned[0] == 1000.0
        assert cleaned[99] == 1000.0

    def test_clean_edge_outlier(self):
        """Outliers at edges should be handled (extrapolation)."""
        vec = np.full(50, 500.0)
        vec[0] = 1500.0  # outlier at left edge
        vec[49] = 1500.0  # outlier at right edge
        cleaned = clean_overscan_vector(vec)
        assert abs(cleaned[0] - 500.0) < 1.0
        assert abs(cleaned[49] - 500.0) < 1.0

    def test_clean_constant_vector(self):
        """A constant vector with strict threshold should be unchanged."""
        vec = np.full(20, 500.0)
        # Even with nsig=0, constant vector has zero deviation
        # from its median-filtered version, so nothing is rejected
        cleaned = clean_overscan_vector(vec, nsig=0.0)
        np.testing.assert_array_equal(cleaned, vec)

    def test_clean_custom_params(self):
        """Custom window and nsig should be respected."""
        vec = np.full(20, 100.0)
        vec[10] = 200.0
        # With nsig=100, threshold = 100*4 = 400, so 200-100=100 < 400
        # Outlier should NOT be cleaned
        cleaned = clean_overscan_vector(vec, w=5, nsig=100.0)
        assert cleaned[10] == 200.0


class TestBinospecReadAmpOverscan:
    """Tests for overscan subtraction in binospec_read_amp."""

    def _make_fake_amp_hdu(self, bias_level=1000.0, signal=500.0):
        """Create a fake Binospec amplifier HDU with known overscan.

        Layout matches real data: NAXIS1=2114, NAXIS2=2072,
        DATASEC=[51:2098,1:2056].  The image is stored in standard
        FITS orientation (not transposed).

        The data section is filled with bias_level + signal.
        The overscan regions contain only bias_level.
        """
        from astropy.io import fits as pyfits

        nx, ny = 2114, 2072
        img = np.full((ny, nx), bias_level, dtype=np.float32)
        # Data section: cols 50:2098, rows 0:2056 (0-indexed)
        img[0:2056, 50:2098] += signal

        hdr = pyfits.Header()
        hdr['DATASEC'] = '[51:2098,1:2056]'
        hdr['DETSEC'] = '[1:2048,1:2056]'
        hdr['NAXIS1'] = nx
        hdr['NAXIS2'] = ny

        hdu_primary = pyfits.PrimaryHDU()
        hdu_ext = pyfits.ImageHDU(data=img, header=hdr)
        hdulist = pyfits.HDUList([hdu_primary, hdu_ext])
        return hdulist

    def test_bias_subtracted(self):
        """Overscan subtraction should remove the bias level."""
        bias = 1000.0
        signal = 500.0
        hdulist = self._make_fake_amp_hdu(bias_level=bias, signal=signal)
        data, overscan, datasec, biassec = binospec_read_amp(hdulist, 1)

        # After overscan subtraction, data section should be close to
        # signal only (bias removed)
        med_data = np.median(data)
        assert abs(med_data - signal) < 5.0, \
            f"Bias not removed: median={med_data}, expected ~{signal}"

    def test_output_shape(self):
        """Output data should have datasec dimensions (2048 x 2056)."""
        hdulist = self._make_fake_amp_hdu()
        data, overscan, datasec, biassec = binospec_read_amp(hdulist, 1)
        # Note: binospec_read_amp transposes the image, so shape is
        # (x, y) = (2048, 2056) after cropping to datasec
        assert data.shape == (2048, 2056), f"Unexpected shape: {data.shape}"

    def test_zero_fake_overscan(self):
        """Returned overscan should be all zeros (fake)."""
        hdulist = self._make_fake_amp_hdu()
        data, overscan, datasec, biassec = binospec_read_amp(hdulist, 1)
        assert np.all(overscan == 0), "Overscan should be fake zeros"

    def test_row_dependent_bias_removed(self):
        """Row-dependent bias structure should be removed by overscan."""
        from astropy.io import fits as pyfits

        nx, ny = 2114, 2072
        # Create a bias pattern that varies along FITS rows (axis 0)
        row_bias = np.linspace(990, 1010, ny).astype(np.float32)
        img = np.broadcast_to(row_bias[:, None], (ny, nx)).copy()
        # Add signal to data section
        img[0:2056, 50:2098] += 500.0

        hdr = pyfits.Header()
        hdr['DATASEC'] = '[51:2098,1:2056]'
        hdr['DETSEC'] = '[1:2048,1:2056]'
        hdr['NAXIS1'] = nx
        hdr['NAXIS2'] = ny
        hdu_primary = pyfits.PrimaryHDU()
        hdu_ext = pyfits.ImageHDU(data=img, header=hdr)
        hdulist = pyfits.HDUList([hdu_primary, hdu_ext])

        data, _, _, _ = binospec_read_amp(hdulist, 1)

        # After overscan subtraction, the row-dependent bias pattern
        # should be mostly removed. Check that the row-wise variation
        # in the data is much less than the original 20 ADU range.
        col_medians = np.median(data, axis=0)  # median along x for each y
        assert np.ptp(col_medians) < 5.0, \
            f"Row-dependent bias not removed: range={np.ptp(col_medians)}"


class TestNonlinearityCorrection:
    """Tests for per-amplifier nonlinearity correction."""

    def test_coefficients_shape(self):
        """Nonlinearity coefficients should be 8x5 (8 amps, degree 4)."""
        coeffs = MMTBINOSPECSpectrograph.nonlinearity_coeffs
        assert coeffs.shape == (8, 5)

    def test_coefficients_zero_constant(self):
        """All constant terms should be zero."""
        coeffs = MMTBINOSPECSpectrograph.nonlinearity_coeffs
        np.testing.assert_array_equal(coeffs[:, 0], 0.0)

    def test_coefficients_near_unity_linear(self):
        """Linear terms should be close to 1.0 (small correction)."""
        coeffs = MMTBINOSPECSpectrograph.nonlinearity_coeffs
        assert np.all(np.abs(coeffs[:, 1] - 1.0) < 0.01)

    def test_correction_applied_in_read_amp(self):
        """binospec_read_amp should apply nonlinearity correction."""
        from astropy.io import fits as pyfits

        nx, ny = 2114, 2072
        bias_level = 1000.0
        signal = 10000.0
        img = np.full((ny, nx), bias_level, dtype=np.float32)
        # Data section gets bias + signal; overscan has only bias
        img[0:2056, 50:2098] += signal

        hdr = pyfits.Header()
        hdr['DATASEC'] = '[51:2098,1:2056]'
        hdr['DETSEC'] = '[1:2048,1:2056]'
        hdr['NAXIS1'] = nx
        hdr['NAXIS2'] = ny
        hdu_primary = pyfits.PrimaryHDU()
        hdu_ext = pyfits.ImageHDU(data=img, header=hdr)
        hdulist = pyfits.HDUList([hdu_primary, hdu_ext])

        data, _, _, _ = binospec_read_amp(hdulist, 1)

        # After overscan subtraction, data ~ signal. Nonlinearity
        # correction then maps signal -> polyval(signal, coeffs).
        coeffs = MMTBINOSPECSpectrograph.nonlinearity_coeffs[0]
        expected = np.polynomial.polynomial.polyval(signal, coeffs)
        med_data = np.median(data)
        assert abs(med_data - expected) < 5.0, \
            f"Nonlinearity not applied: got {med_data}, expected ~{expected}"

    def test_correction_is_small(self):
        """At typical science levels (~1000 ADU), correction < 1%."""
        coeffs = MMTBINOSPECSpectrograph.nonlinearity_coeffs
        test_counts = 1000.0
        for i in range(8):
            corrected = np.polynomial.polynomial.polyval(test_counts, coeffs[i])
            ratio = corrected / test_counts
            assert 0.99 < ratio < 1.01, \
                f"Amp {i+1}: correction too large: {ratio}"
