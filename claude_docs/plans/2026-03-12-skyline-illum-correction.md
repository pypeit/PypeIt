# Sky-Line-Based Fiber Illumination Correction

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Measure sky emission line fluxes across all fibers and derive a wavelength-dependent throughput correction that equalizes sky/science fiber response before joint sky subtraction.

**Architecture:** Add a `compute_skyline_illum()` method on `MMTBINOSPECIFUSpectrograph` that measures bright sky line fluxes per fiber, normalizes to the median, and interpolates to build a 2D correction image. A no-op hook on the base `Spectrograph` class (`skyline_illum_correct()`) is called from `joint_skysub()` in `find_objects.py` just before the iterative B-spline fit loop. The hook only modifies `sciimg` (a working copy); variance correction is applied once via `apply_relative_scale()` after the iterative loop converges. This avoids double-applying the variance correction.

**Tech Stack:** numpy, scipy.interpolate (interp1d), PypeIt internals (find_objects, spectrograph, flat)

---

## File Map

| File | Action | Responsibility |
|------|--------|---------------|
| `pypeit/spectrographs/spectrograph.py` | Modify | Add no-op `skyline_illum_correct()` hook |
| `pypeit/spectrographs/mmt_binospec.py` | Modify | Implement `compute_skyline_illum()` and override `skyline_illum_correct()` |
| `pypeit/find_objects.py` | Modify | Call `skyline_illum_correct()` hook in `joint_skysub()` |
| `pypeit/tests/test_skyline_illum.py` | Create | Unit tests for the correction algorithm |

---

## Task 1: Add base-class hook `skyline_illum_correct()` to Spectrograph

**Files:**
- Modify: `pypeit/spectrographs/spectrograph.py` (near `modify_pixelflat`, ~line 1383)

- [ ] **Step 1: Add the no-op hook method**

Add after `modify_pixelflat()` (around line 1383):

```python
def skyline_illum_correct(self, sciimg, waveimg, slits, slitmask):
    """
    Spectrograph-specific sky-line-based illumination correction.

    Measures bright sky emission line fluxes across fibers and builds
    a wavelength-dependent throughput correction.  The default
    implementation does nothing.

    This method should only modify ``sciimg`` in place (dividing by
    the correction).  Variance propagation is handled separately by
    the caller via :meth:`~pypeit.find_objects.FindObjects.apply_relative_scale`.

    Args:
        sciimg (`numpy.ndarray`_):
            2D flat-fielded science image (nspec, nspat).  Modified
            in place (divided by the correction).
        waveimg (`numpy.ndarray`_):
            Wavelength image in Angstroms.
        slits (:class:`~pypeit.slittrace.SlitTraceSet`):
            Slit traces.
        slitmask (`numpy.ndarray`_):
            2D image mapping pixels to slit ``spat_id``.

    Returns:
        `numpy.ndarray`_: 2D correction image that was applied
        (1.0 everywhere if no correction).
    """
    return np.ones_like(sciimg)
```

- [ ] **Step 2: Verify import exists**

`numpy` is already imported in `spectrograph.py`. No new imports needed.

- [ ] **Step 3: Commit**

```bash
git add pypeit/spectrographs/spectrograph.py
git commit -m "Add no-op skyline_illum_correct hook to Spectrograph base class"
```

---

## Task 2: Implement `compute_skyline_illum()` on Binospec IFU

**Files:**
- Modify: `pypeit/spectrographs/mmt_binospec.py` (after `modify_pixelflat`, ~line 1577)

- [ ] **Step 1: Add the sky line list as a class attribute**

Add to `MMTBINOSPECIFUSpectrograph` class attributes (after `sky_fiber_ids`, ~line 1132):

```python
# Bright sky emission lines for throughput correction (Angstroms).
# These must be isolated enough to measure reliably with
# continuum sidebands.  Source: Binospec IDL pipeline
# (Chilingarian et al. 2025, arXiv:2501.01528).
# Note: 4358.335 (Hg I) omitted -- unreliable at dark sites.
skyline_list_ang = np.array([
    5577.34, 6300.304, 6863.951, 7340.881,
    7993.327, 8465.353, 8885.843, 9502.808
])
```

- [ ] **Step 2: Add the `compute_skyline_illum()` method**

This is the core algorithm. It accepts `spat_ids` (numpy array) rather than
a `SlitTraceSet` so that unit tests can call it directly with mock data.
Add after `modify_pixelflat()`:

```python
def compute_skyline_illum(self, sciimg, waveimg, slitmask, spat_ids):
    """
    Compute per-fiber throughput correction from sky emission lines.

    For each sky line within the wavelength range, extracts boxcar
    flux per fiber, subtracts local continuum, and normalizes by
    the median across fibers.  Interpolates between lines to build
    a wavelength-dependent correction per fiber.

    Fibers with no valid line measurements (e.g. dead fibers) are
    left uncorrected (correction = 1.0).

    Args:
        sciimg (`numpy.ndarray`_):
            2D flat-fielded science image.
        waveimg (`numpy.ndarray`_):
            Wavelength image in Angstroms.
        slitmask (`numpy.ndarray`_):
            2D slit image (pixel -> spat_id).
        spat_ids (`numpy.ndarray`_):
            Array of spatial IDs for each fiber.

    Returns:
        `numpy.ndarray`_: 2D correction image (same shape as
        sciimg).  Values > 1 for fibers brighter than median.
    """
    from pypeit import log
    from scipy.interpolate import interp1d

    nfibers = len(spat_ids)

    # Determine wavelength range from the data
    valid = waveimg > 0
    if not np.any(valid):
        log.warning("No valid wavelength data; skipping skyline "
                    "illumination correction")
        return np.ones_like(sciimg)
    wmin = waveimg[valid].min()
    wmax = waveimg[valid].max()

    # Filter sky lines to those within the wavelength range with
    # enough margin for continuum estimation (50 Ang each side)
    margin = 50.0
    usable = ((self.skyline_list_ang > wmin + margin)
              & (self.skyline_list_ang < wmax - margin))
    sky_lines = self.skyline_list_ang[usable]
    if len(sky_lines) == 0:
        log.warning("No sky lines in wavelength range "
                    f"[{wmin:.0f}, {wmax:.0f}] Ang; skipping "
                    "skyline illumination correction")
        return np.ones_like(sciimg)
    log.info(f"Measuring {len(sky_lines)} sky lines for fiber "
             f"throughput correction")

    # Extract 1D boxcar spectrum per fiber
    nspec = sciimg.shape[0]
    fiber_flux = np.zeros((nfibers, nspec))
    fiber_wave = np.zeros((nfibers, nspec))
    for i, spat_id in enumerate(spat_ids):
        slit_pix = slitmask == spat_id
        for j in range(nspec):
            row_pix = slit_pix[j, :]
            if np.any(row_pix):
                fiber_flux[i, j] = np.sum(sciimg[j, row_pix])
                fiber_wave[i, j] = np.mean(waveimg[j, row_pix])

    # Measure each sky line in each fiber
    line_window = 4.0  # Angstroms half-width for line flux
    cont_inner = 8.0   # Angstroms from line center to start of
                        # continuum window
    cont_outer = 20.0  # Angstroms from line center to end of
                        # continuum window
    min_valid_fibers = max(10, nfibers // 10)

    # Shape: (n_lines, n_fibers)
    line_ratios = np.full((len(sky_lines), nfibers), np.nan)

    for k, wl in enumerate(sky_lines):
        fiber_line_flux = np.zeros(nfibers)
        for i in range(nfibers):
            wave_i = fiber_wave[i]
            flux_i = fiber_flux[i]
            good = wave_i > 0

            if not np.any(good):
                continue

            # Line window
            in_line = good & (np.abs(wave_i - wl) <= line_window)
            # Continuum windows (blue and red sidebands)
            in_cont = good & (np.abs(wave_i - wl) >= cont_inner) \
                           & (np.abs(wave_i - wl) <= cont_outer)

            if np.sum(in_line) < 2 or np.sum(in_cont) < 3:
                continue

            cont_level = np.median(flux_i[in_cont])
            line_sum = np.sum(flux_i[in_line] - cont_level)
            fiber_line_flux[i] = line_sum

        # Normalize by median across fibers (exclude zeros/negatives)
        valid_flux = fiber_line_flux > 0
        if np.sum(valid_flux) < min_valid_fibers:
            log.warning(f"Sky line {wl:.1f} Ang: too few valid "
                        f"fibers ({np.sum(valid_flux)}), skipping")
            continue
        med_flux = np.median(fiber_line_flux[valid_flux])
        ratios = np.where(valid_flux,
                          fiber_line_flux / med_flux, np.nan)
        line_ratios[k] = ratios
        n_valid = np.sum(valid_flux)
        rmin = np.nanmin(ratios)
        rmax = np.nanmax(ratios)
        log.info(f"  {wl:.1f} Ang: {n_valid} fibers, "
                 f"range {rmin:.3f} - {rmax:.3f}")

    # Check we have at least one usable line
    usable_lines = ~np.all(np.isnan(line_ratios), axis=1)
    if not np.any(usable_lines):
        log.warning("No sky lines measured successfully; skipping "
                    "skyline illumination correction")
        return np.ones_like(sciimg)
    sky_lines = sky_lines[usable_lines]
    line_ratios = line_ratios[usable_lines]

    # Build wavelength-dependent correction per fiber.
    # For fibers with missing measurements at some lines, use
    # nearest valid value via fill_value extrapolation.
    corr_2d = np.ones_like(sciimg)
    for i, spat_id in enumerate(spat_ids):
        slit_pix = slitmask == spat_id
        if not np.any(slit_pix):
            continue

        # Get this fiber's ratios across lines
        ratios_i = line_ratios[:, i]
        good_lines = ~np.isnan(ratios_i)

        if not np.any(good_lines):
            continue

        if np.sum(good_lines) == 1:
            # Single line: constant correction
            corr_val = ratios_i[good_lines][0]
            corr_2d[slit_pix] = corr_val
        else:
            # Interpolate between lines
            interp_func = interp1d(
                sky_lines[good_lines], ratios_i[good_lines],
                kind='linear', bounds_error=False,
                fill_value=(ratios_i[good_lines][0],
                            ratios_i[good_lines][-1]))
            # Get wavelength at each pixel in this fiber
            wave_pix = waveimg[slit_pix]
            corr_2d[slit_pix] = interp_func(wave_pix)

    # Safety: clip extreme corrections
    corr_2d = np.clip(corr_2d, 0.3, 3.0)

    return corr_2d
```

- [ ] **Step 3: Add the `skyline_illum_correct()` override**

This override only modifies `sciimg` in place. It does NOT touch ivar —
variance propagation is handled by the caller via `apply_relative_scale()`.

```python
def skyline_illum_correct(self, sciimg, waveimg, slits, slitmask):
    """
    Apply sky-line-based illumination correction.

    Corrects for throughput differences between sky fibers (bare
    fibers) and science fibers (lenslet-fed) that the dome-flat-based
    illumination correction cannot capture.  Only modifies ``sciimg``
    in place; variance is handled by the caller.

    See :meth:`compute_skyline_illum` for the algorithm.
    """
    from pypeit import log

    corr = self.compute_skyline_illum(sciimg, waveimg, slitmask,
                                      slits.spat_id)
    if np.allclose(corr, 1.0):
        return corr

    # Apply: divide science image only (variance handled by caller)
    good = corr > 0.1
    sciimg[good] /= corr[good]
    log.info("Applied sky-line illumination correction "
             f"(range {corr[good].min():.3f} - "
             f"{corr[good].max():.3f})")
    return corr
```

- [ ] **Step 4: Commit**

```bash
git add pypeit/spectrographs/mmt_binospec.py
git commit -m "Implement sky-line-based fiber illumination correction for Binospec IFU

Measures bright sky emission line fluxes in each fiber, normalizes
by the median across fibers, and interpolates to build a wavelength-
dependent throughput correction. This accounts for the different
optical paths of sky fibers (bare fibers) vs science fibers (lenslet
array), which the dome-flat-based f_illum cannot capture."
```

---

## Task 3: Call the hook from `joint_skysub()`

**Files:**
- Modify: `pypeit/find_objects.py` (in `joint_skysub`, ~line 1099)

- [ ] **Step 1: Insert the hook call before the iterative loop**

Insert after line 1099 (`sciimg = skysub.convolve_skymodel(...)`) and before
line 1101 (`numiter = 4`).  At this point `sciimg` is a copy returned by
`convolve_skymodel`, so in-place modification is safe.  The hook does NOT
touch `model_ivar` — variance is corrected once later via
`apply_relative_scale`.

```python
        # Apply spectrograph-specific sky-line illumination correction.
        # For fiber-fed spectrographs (e.g. Binospec IFU), this equalizes
        # throughput differences between sky and science fibers that are
        # not captured by the dome-flat-based illumination correction.
        # Only sciimg (a copy) is modified here; variance propagation to
        # self.sciImg is deferred to apply_relative_scale after the loop.
        # NOTE: This interacts with the illum_profile_spectral_poly
        # iteration below — both corrections converge together, which is
        # the desired behavior.
        skyline_corr = self.spectrograph.skyline_illum_correct(
            sciimg, self.waveimg, self.slits, slitmask)
```

- [ ] **Step 2: Apply the correction to the original science image after the loop**

After line 1161 (`self.apply_relative_scale(scaleImg)`), add:

```python
        # Apply the skyline illumination correction to the original
        # science image (and propagate to variance) so the final sky
        # recalculation uses corrected data.
        if not np.allclose(skyline_corr, 1.0):
            self.apply_relative_scale(skyline_corr)
```

- [ ] **Step 3: Commit**

```bash
git add pypeit/find_objects.py
git commit -m "Call skyline_illum_correct hook in joint_skysub before sky fitting"
```

---

## Task 4: Unit tests

**Files:**
- Create: `pypeit/tests/test_skyline_illum.py`

Tests call the production `compute_skyline_illum()` method directly, passing
`spat_ids` as a numpy array (the method accepts this without needing a full
`SlitTraceSet`).

- [ ] **Step 1: Write tests**

```python
"""Tests for sky-line-based fiber illumination correction."""
import numpy as np
import pytest

from pypeit.spectrographs.mmt_binospec import MMTBINOSPECIFUSpectrograph


def _make_mock_data(nfibers=40, nspec=200, wmin=5000.0, wmax=9000.0,
                    sky_fiber_boost=1.3, n_sky_fibers=4):
    """Create synthetic 2D data with sky lines and fiber throughput
    variation.

    Returns sciimg, waveimg, slitmask, spat_ids, sky_indices.
    """
    nspat_per_fiber = 3
    nspat = nfibers * nspat_per_fiber
    sciimg = np.zeros((nspec, nspat))
    waveimg = np.zeros((nspec, nspat))
    slitmask = np.zeros((nspec, nspat), dtype=int)

    wave_arr = np.linspace(wmin, wmax, nspec)
    spat_ids = np.arange(1, nfibers + 1)
    sky_indices = np.arange(n_sky_fibers)

    # Throughput varies smoothly across fibers, sky fibers are brighter
    throughput = np.ones(nfibers)
    throughput[sky_indices] = sky_fiber_boost
    # Add some variation to science fibers
    rng = np.random.default_rng(42)
    throughput[n_sky_fibers:] += 0.05 * rng.standard_normal(
        nfibers - n_sky_fibers)

    # Sky continuum + emission lines
    sky_lines_ang = np.array([5577.34, 6300.304, 6863.951, 7340.881,
                              7993.327, 8465.353])
    line_sigma = 2.0  # pixels

    for i in range(nfibers):
        cols = slice(i * nspat_per_fiber, (i + 1) * nspat_per_fiber)
        slitmask[:, cols] = spat_ids[i]
        waveimg[:, cols] = wave_arr[:, None]

        # Sky spectrum: continuum + lines
        sky_spec = np.full(nspec, 100.0)  # continuum
        for wl in sky_lines_ang:
            pix = np.argmin(np.abs(wave_arr - wl))
            sky_spec += 500.0 * np.exp(
                -0.5 * ((np.arange(nspec) - pix) / line_sigma) ** 2)

        # Apply throughput
        sky_spec *= throughput[i]
        sciimg[:, cols] = sky_spec[:, None]

    return sciimg, waveimg, slitmask, spat_ids, sky_indices


@pytest.fixture
def spectrograph():
    return MMTBINOSPECIFUSpectrograph()


class TestComputeSkylineIllum:
    """Test the core correction algorithm."""

    def test_uniform_throughput_gives_unit_correction(self, spectrograph):
        """If all fibers have same throughput, correction should be ~1."""
        sciimg, waveimg, slitmask, spat_ids, _ = _make_mock_data(
            sky_fiber_boost=1.0)
        corr = spectrograph.compute_skyline_illum(
            sciimg, waveimg, slitmask, spat_ids)
        assert np.allclose(corr[slitmask > 0], 1.0, atol=0.05)

    def test_sky_fibers_get_higher_correction(self, spectrograph):
        """Sky fibers with 30% more throughput should get corr > 1."""
        sciimg, waveimg, slitmask, spat_ids, sky_idx = _make_mock_data(
            sky_fiber_boost=1.3)
        corr = spectrograph.compute_skyline_illum(
            sciimg, waveimg, slitmask, spat_ids)
        # Sky fiber pixels should have correction > 1
        for i in sky_idx:
            sky_pix = slitmask == spat_ids[i]
            assert np.median(corr[sky_pix]) > 1.1

    def test_correction_normalizes_sky_level(self, spectrograph):
        """After applying correction, sky line flux should be
        equalized across fibers."""
        sciimg, waveimg, slitmask, spat_ids, sky_idx = _make_mock_data(
            sky_fiber_boost=1.3)
        corr = spectrograph.compute_skyline_illum(
            sciimg, waveimg, slitmask, spat_ids)
        corrected = sciimg / np.where(corr > 0.1, corr, 1.0)

        # Measure flux in a sky line region for sky vs science fibers
        wave_arr = waveimg[:, 1]
        line_pix = np.abs(wave_arr - 5577.34) < 10
        sky_flux = []
        sci_flux = []
        nspat_per = 3
        for i in range(len(spat_ids)):
            cols = slice(i * nspat_per, (i + 1) * nspat_per)
            f = np.sum(corrected[line_pix, cols])
            if i in sky_idx:
                sky_flux.append(f)
            else:
                sci_flux.append(f)
        # After correction, sky and science fibers should have similar
        # flux (within 5%)
        ratio = np.mean(sky_flux) / np.mean(sci_flux)
        assert abs(ratio - 1.0) < 0.05

    def test_no_lines_in_range(self, spectrograph):
        """If no sky lines fall in wavelength range, return ones."""
        sciimg, waveimg, slitmask, spat_ids, _ = _make_mock_data(
            wmin=3000.0, wmax=4000.0)
        corr = spectrograph.compute_skyline_illum(
            sciimg, waveimg, slitmask, spat_ids)
        assert np.allclose(corr, 1.0)

    def test_extreme_corrections_clipped(self, spectrograph):
        """Corrections should be clipped to [0.3, 3.0]."""
        sciimg, waveimg, slitmask, spat_ids, _ = _make_mock_data(
            sky_fiber_boost=5.0)
        corr = spectrograph.compute_skyline_illum(
            sciimg, waveimg, slitmask, spat_ids)
        assert np.all(corr >= 0.3)
        assert np.all(corr <= 3.0)
```

- [ ] **Step 2: Run the tests**

```bash
pytest pypeit/tests/test_skyline_illum.py -v
```

Expected: All 5 tests pass.

- [ ] **Step 3: Commit**

```bash
git add pypeit/tests/test_skyline_illum.py
git commit -m "Add unit tests for sky-line-based fiber illumination correction"
```

---

## Task 5: Integration test with real data

This task verifies the full pipeline works with actual Binospec IFU data.

- [ ] **Step 1: Delete cached spec2d files to force re-reduction**

```bash
rm ~/MMT/bino_ifu/JADES_1031022/Science/spec2d_*.fits
```

(Calibrations can stay — the pixel flat already has f_illum baked in.)

- [ ] **Step 2: Run the pipeline**

```bash
cd ~/MMT/bino_ifu/JADES_1031022
run_pypeit jades_1031022.pypeit -o
```

- [ ] **Step 3: Inspect results**

Check the log output for:
- `Measuring N sky lines for fiber throughput correction`
- Per-line range printout (expect sky fibers to have ratios > 1)
- `Applied sky-line illumination correction (range X - Y)`

Inspect spec2d files:
- SKYMODEL extension should no longer be over-subtracted
- Residuals (science - sky) should be flat across fiber types

- [ ] **Step 4: Commit any adjustments**

If parameters need tuning (line_window, cont margins, clipping), update and commit.
