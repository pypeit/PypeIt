# Binospec Overscan & Nonlinearity Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development
> (if subagents available) or superpowers:executing-plans to implement this plan.
> Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Update `binospec_read_amp()` to match the IDL pipeline's overscan
subtraction (resistant_mean + cleaning, both overscan regions, Y-first) and add
per-amplifier polynomial nonlinearity correction.

**Architecture:** All changes are in `pypeit/spectrographs/mmt_binospec.py`.
A new module-level helper function (`clean_overscan_vector`)
plus a rewrite of the overscan block in the existing
`binospec_read_amp()` function, and a new class attribute for nonlinearity
coefficients on `MMTBINOSPECSpectrograph`.

**Tech Stack:** numpy, scipy.ndimage (median_filter), astropy.stats
(sigma_clipped_stats), astropy.io.fits

**Spec:** `claude_docs/specs/2026-03-16-binospec-overscan-nonlinearity-design.md`

**IDL reference:** `~/MMT/binospec/pipeline/bino_mosaic.pro` (functions
`clean_overscan_vector` and `remove_overscan`, lines 1-70; `bino_mosaic`
procedure lines 115-126)

---

## File Map

- **Modify:** `pypeit/spectrographs/mmt_binospec.py`
  - Add import: `from astropy.stats import sigma_clipped_stats`
  - Add import: `from scipy.ndimage import median_filter`
  - Add module-level function: `clean_overscan_vector()` (before
    `binospec_read_amp`)
  - Rewrite overscan logic in `binospec_read_amp()` (lines 1080-1094)
  - Add nonlinearity correction in `binospec_read_amp()` (after overscan,
    before return)
  - Add class attribute `nonlinearity_coeffs` on `MMTBINOSPECSpectrograph`
- **Create:** `pypeit/tests/test_binospec_overscan.py` (unit tests)

---

## Task 1: Add `clean_overscan_vector()` with Tests

**Files:**
- Modify: `pypeit/spectrographs/mmt_binospec.py` (add function before
  `binospec_read_amp`, around line 1040)
- Create: `pypeit/tests/test_binospec_overscan.py`

- [ ] **Step 1: Write the test file with tests for `clean_overscan_vector`**

Create `pypeit/tests/test_binospec_overscan.py`:

```python
"""
Tests for Binospec overscan subtraction and nonlinearity correction.
"""
import numpy as np
import pytest

from pypeit.spectrographs.mmt_binospec import clean_overscan_vector


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

    def test_clean_all_rejected(self):
        """If all pixels are rejected, return unmodified vector."""
        # With rdnoise=4.0 and nsig=1.0, threshold is 4.0.
        # Make a vector where every pixel differs from its local median
        # by more than 4.0 -- practically impossible with smooth data,
        # so we test the edge case by using a very noisy vector with
        # nsig=0 (everything rejected).
        vec = np.array([100.0, 200.0, 300.0, 400.0, 500.0])
        # nsig=0 means threshold=0, so any deviation rejects
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
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /Users/tim/MMT/pypeit && python -m pytest pypeit/tests/test_binospec_overscan.py -v`
Expected: ImportError — `clean_overscan_vector` does not exist yet.

- [ ] **Step 3: Add imports and implement `clean_overscan_vector`**

Add to the import block at the top of `mmt_binospec.py` (after existing
imports, around line 14):

```python
from scipy.ndimage import median_filter
```

Add the function before `binospec_read_amp` (insert before line 1041):

```python
def clean_overscan_vector(overscan, w=9, nsig=1.0, rdnoise=4.0):
    """
    Clean a 1D overscan vector by median-filtering and interpolating
    over outliers.

    Replicates the IDL ``clean_overscan_vector`` function from
    ``bino_mosaic.pro``.

    Parameters
    ----------
    overscan : `numpy.ndarray`_
        1D overscan vector to clean.
    w : :obj:`int`, optional
        Window size for median filtering. Must be >= 3. Default is 9.
    nsig : :obj:`float`, optional
        Sigma threshold for outlier rejection. Pixels deviating from
        the median-filtered vector by more than ``nsig * rdnoise`` are
        replaced by interpolation. Default is 1.0.
    rdnoise : :obj:`float`, optional
        Read noise in ADU, used to set the outlier threshold.
        Default is 4.0.

    Returns
    -------
    clean : `numpy.ndarray`_
        Cleaned overscan vector with outliers interpolated over.
    """
    w = max(w, 3)
    m_overscan = median_filter(overscan, size=w, mode='reflect')
    bad = np.abs(overscan - m_overscan) > rdnoise * nsig
    good = ~bad
    if not np.any(bad):
        return overscan.copy()
    if not np.any(good):
        return overscan.copy()
    clean = overscan.copy()
    good_idx = np.where(good)[0]
    bad_idx = np.where(bad)[0]
    clean[bad_idx] = np.interp(bad_idx, good_idx, overscan[good_idx])
    return clean
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd /Users/tim/MMT/pypeit && python -m pytest pypeit/tests/test_binospec_overscan.py::TestCleanOverscanVector -v`
Expected: All 5 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add pypeit/spectrographs/mmt_binospec.py pypeit/tests/test_binospec_overscan.py
git commit -m "Add clean_overscan_vector helper for Binospec overscan cleaning"
```

---

## Task 2: Rewrite Overscan Subtraction in `binospec_read_amp()`

**Files:**
- Modify: `pypeit/spectrographs/mmt_binospec.py:1080-1099` (the overscan
  block inside `binospec_read_amp`)
- Modify: `pypeit/tests/test_binospec_overscan.py` (add overscan tests)

- [ ] **Step 1: Write tests for the new overscan subtraction**

Add to `pypeit/tests/test_binospec_overscan.py`:

```python
from pypeit.spectrographs.mmt_binospec import binospec_read_amp


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
```

- [ ] **Step 2: Run tests to verify the new tests fail**

Run: `cd /Users/tim/MMT/pypeit && python -m pytest pypeit/tests/test_binospec_overscan.py::TestBinospecReadAmpOverscan -v`
Expected: Some tests fail (the current code uses simple median and wrong
axis order, and doesn't use both X overscan regions).

- [ ] **Step 3: Add `astropy.stats` import and rewrite overscan in `binospec_read_amp`**

Add to the import block at top of `mmt_binospec.py`:

```python
from astropy.stats import sigma_clipped_stats
```

Replace lines 1080-1099 of `binospec_read_amp()` (the overscan subtraction
block, from the `# NOTE:` comment through the `overscan = np.zeros_like`
line) with:

```python
    # Overscan subtraction following IDL pipeline (bino_mosaic.pro):
    # Y-axis first, then X-axis. Uses sigma-clipped mean (resistant_mean)
    # with outlier cleaning, matching IDL defaults (clean_w=9, clean_nsig=1.0).

    # Y-axis overscan: postscan rows after datasec
    if y2 < nyt:
        overscan_y = temp[:, y2:nyt]
        overscan_vec, _, _ = sigma_clipped_stats(overscan_y, sigma=3.0, axis=1)
        overscan_vec = clean_overscan_vector(overscan_vec, w=9, nsig=1.0)
        temp = temp - overscan_vec[:, None]

    # X-axis overscan: prescan + postscan columns
    overscan_x_regions = []
    if x1 > 1:
        overscan_x_regions.append(temp[0:x1-1, :])
    if x2 < nxt:
        overscan_x_regions.append(temp[x2:nxt, :])
    if len(overscan_x_regions) > 0:
        overscan_x = np.concatenate(overscan_x_regions, axis=0)
        overscan_x_vec, _, _ = sigma_clipped_stats(overscan_x, sigma=3.0, axis=0)
        overscan_x_vec = clean_overscan_vector(overscan_x_vec, w=9, nsig=1.0)
        temp = temp - overscan_x_vec[None, :]

    # Crop to datasec
    data = temp[x1-1:x2, y1-1:y2]

    # Fake overscan for PypeIt's general pipeline (effectively a no-op)
    biassec = f'[0:{x1-1},{y1-1}:{y2}]'
    xos1, xos2, yos1, yos2 = chain.from_iterable(parse.load_sections(biassec, fmt_iraf=False))
    overscan = np.zeros_like(temp[xos1:xos2, yos1:yos2])
```

Note on axis mapping: `binospec_read_amp` transposes the raw FITS image at
line 1067 (`temp = hdu[ext].data.transpose()`), so after transpose:
- axis 0 = X (columns in FITS = NAXIS1 direction, 2114 pixels)
- axis 1 = Y (rows in FITS = NAXIS2 direction, 2072 pixels)

The IDL `remove_overscan` operates on `image[x, y]` with `datasec =
[x1, x2, y1, y2]` (1-indexed). After PypeIt's transpose, `temp[x, y]` has
the same orientation.

- Y overscan: `temp[:, y2:nyt]` = all X pixels, postscan Y rows.
  `sigma_clipped_stats(..., axis=1)` collapses along Y → 1D vector of
  length nxt. Subtracted as `temp - vec[:, None]` (column vector broadcast).
  This matches IDL line 44: `image_sub = image_sub - (overscan_vec #
  (fltarr(s_im[2])+1))`.

- X overscan: `temp[0:x1-1, :]` (prescan) + `temp[x2:nxt, :]` (postscan),
  concatenated along axis 0. `sigma_clipped_stats(..., axis=0)` collapses
  along X → 1D vector of length nyt. Subtracted as `temp - vec[None, :]`
  (row vector broadcast). This matches IDL line 65:
  `image_sub = image_sub - ((fltarr(s_im[1])+1) # overscan_x_vec)`.

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd /Users/tim/MMT/pypeit && python -m pytest pypeit/tests/test_binospec_overscan.py -v`
Expected: All tests PASS (both `TestCleanOverscanVector` and
`TestBinospecReadAmpOverscan`).

- [ ] **Step 5: Commit**

```bash
git add pypeit/spectrographs/mmt_binospec.py pypeit/tests/test_binospec_overscan.py
git commit -m "Rewrite Binospec overscan subtraction to match IDL pipeline"
```

---

## Task 3: Add Nonlinearity Correction

**Files:**
- Modify: `pypeit/spectrographs/mmt_binospec.py` (add class attribute +
  apply in `binospec_read_amp`)
- Modify: `pypeit/tests/test_binospec_overscan.py` (add nonlinearity tests)

- [ ] **Step 1: Write tests for nonlinearity correction**

Add to `pypeit/tests/test_binospec_overscan.py`:

```python
from pypeit.spectrographs.mmt_binospec import MMTBINOSPECSpectrograph


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
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /Users/tim/MMT/pypeit && python -m pytest pypeit/tests/test_binospec_overscan.py::TestNonlinearityCorrection -v`
Expected: AttributeError — `nonlinearity_coeffs` does not exist yet.

- [ ] **Step 3: Add nonlinearity coefficients and apply in `binospec_read_amp`**

Add the class attribute to `MMTBINOSPECSpectrograph` (after `name =
'mmt_binospec'`, around line 35):

```python
    # Per-amplifier nonlinearity correction coefficients from IDL pipeline
    # calibration (scicam_bino_sep2017.fits). Applied as:
    #   C_corr = polyval(C_raw, coeffs[amp_index])
    # Amps 1-4 = detector A (DET01), amps 5-8 = detector B (DET02).
    nonlinearity_coeffs = np.array([
        [0.00000000e+00, 1.00400089e+00, -1.39235362e-06,
         8.31711824e-12, -1.20653479e-17],
        [0.00000000e+00, 1.00361458e+00, -1.29223833e-06,
         6.93723177e-12, -9.67406255e-18],
        [0.00000000e+00, 1.00269542e+00, -9.29361806e-07,
         5.97902827e-12, -2.30257302e-17],
        [0.00000000e+00, 1.00339616e+00, -8.47134521e-07,
         7.92441693e-12, -4.46542834e-17],
        [0.00000000e+00, 1.00727205e+00, -1.69093388e-06,
         2.07225055e-11, -1.62655178e-16],
        [0.00000000e+00, 1.00858745e+00, -2.35668901e-06,
         2.40641019e-11, -1.50286358e-16],
        [0.00000000e+00, 1.00728526e+00, -1.80779473e-06,
         1.73427719e-11, -1.01685780e-16],
        [0.00000000e+00, 1.00845168e+00, -2.02050567e-06,
         2.97587091e-11, -2.65508521e-16],
    ])
```

Add nonlinearity correction in `binospec_read_amp()`, after the `data =
temp[x1-1:x2, y1-1:y2]` line and before the `biassec` line:

```python
    # Apply per-amplifier nonlinearity correction (IDL: poly(im_cur, c_poly))
    data = np.polynomial.polynomial.polyval(
        data, MMTBINOSPECSpectrograph.nonlinearity_coeffs[ext - 1])
```

- [ ] **Step 4: Run all tests to verify they pass**

Run: `cd /Users/tim/MMT/pypeit && python -m pytest pypeit/tests/test_binospec_overscan.py -v`
Expected: All tests PASS.

- [ ] **Step 5: Commit**

```bash
git add pypeit/spectrographs/mmt_binospec.py pypeit/tests/test_binospec_overscan.py
git commit -m "Add per-amplifier nonlinearity correction for Binospec"
```

---

## Task 4: Verify Against Existing Test Suite

**Files:** None modified — verification only.

- [ ] **Step 1: Run the full PypeIt test suite**

Run: `cd /Users/tim/MMT/pypeit && python -m pytest pypeit/tests/ -x -q --timeout=120`
Expected: All tests pass. No shared infrastructure was modified, so only
Binospec-specific tests are new.

- [ ] **Step 2: Verify imports work correctly**

Run: `cd /Users/tim/MMT/pypeit && python -c "from pypeit.spectrographs.mmt_binospec import MMTBINOSPECSpectrograph, MMTBINOSPECIFUSpectrograph, binospec_read_amp, clean_overscan_vector; print('All imports OK'); print('IFU inherits nonlinearity_coeffs:', hasattr(MMTBINOSPECIFUSpectrograph, 'nonlinearity_coeffs'))"`
Expected: `All imports OK` and `IFU inherits nonlinearity_coeffs: True`

- [ ] **Step 3: Commit (if any fixes were needed)**

Only if test failures required fixes in earlier code.
