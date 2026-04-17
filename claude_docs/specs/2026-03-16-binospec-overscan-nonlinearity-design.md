# Binospec Overscan Subtraction and Nonlinearity Correction

**Date**: 2026-03-16
**Status**: Approved
**Scope**: `pypeit/spectrographs/mmt_binospec.py` only

## Problem

PypeIt's Binospec overscan subtraction (`binospec_read_amp()`) differs from
the IDL pipeline (`remove_overscan` in `bino_mosaic.pro`) in several ways that
produce measurable differences in processed images:

1. Only uses one overscan region per axis (prescan for X, postscan for Y)
   instead of combining both sides
2. Uses simple `np.median` instead of sigma-clipped mean with outlier cleaning
3. Subtracts X-axis first, then Y-axis (IDL does Y first, then X)
4. No nonlinearity correction (IDL applies a per-amplifier polynomial)

## Binospec Detector Geometry

Each of the 8 amplifier extensions has identical layout (from FITS headers):

- **Raw image**: 2114 x 2072 pixels
- **DATASEC**: `[51:2098, 1:2056]` (2048 x 2056 data pixels)
- **X prescan**: columns 1-50 (50 pixels)
- **X postscan**: columns 2099-2114 (16 pixels)
- **Y prescan**: none
- **Y postscan**: rows 2057-2072 (16 pixels)

## Design

All changes are in `pypeit/spectrographs/mmt_binospec.py`.

### 1. Add `clean_overscan_vector()` Module-Level Helper

Replicates the IDL `clean_overscan_vector` function. Defined as a module-level
function (consistent with `binospec_read_amp` which is also module-level):

1. Median-filter the 1D overscan vector (window=9)
2. Identify outliers: pixels deviating from the median-filtered version by
   more than `nsig * rdnoise` (nsig=1.0, rdnoise=4.0)
3. Interpolate over rejected pixels using `np.interp` from good pixels.
   Edge pixels with no good neighbors are extrapolated (nearest good value).
   If all pixels are rejected, return the unmodified vector.

Parameters match IDL defaults: `w=9` (clean window), `nsig=1.0`,
`rdnoise=4.0`.

### 2. Rewrite Overscan Subtraction in `binospec_read_amp()`

Replace the current overscan logic with the IDL `remove_overscan` approach:

**Y-axis overscan first**:
1. Collect postscan rows (rows after datasec in Y). For current geometry:
   rows 2057-2072 (16 rows). No prescan exists in Y.
2. Collapse to 1D vector using sigma-clipped mean at 3 sigma, implemented
   via `astropy.stats.sigma_clipped_stats(data, sigma=3.0, axis=...)`,
   taking the mean. This matches IDL's `resistant_mean` with 3-sigma
   iterative clipping.
3. Clean the vector with `clean_overscan_vector()`.
4. Subtract from the full image.

**X-axis overscan second**:
1. Collect prescan columns (cols 1-50) and postscan columns (cols 2099-2114),
   combining both regions into a single array.
2. Collapse to 1D vector using sigma-clipped mean (same as above).
3. Clean the vector.
4. Subtract from the (already Y-corrected) image.

**Crop** to datasec and return.

The function continues to return a fake zero overscan array so that PypeIt's
general `subtract_overscan()` step is effectively a no-op. The `datasec` and
`biassec` return values remain unchanged, as `get_rawimage()` uses them to
populate `rawdatasec_img` and `oscansec_img`.

### 3. Add Nonlinearity Correction Coefficients

Store the per-amplifier polynomial coefficients from the IDL calibration file
(`scicam_bino_sep2017.fits`, measured 2017-09) as a class attribute on
`MMTBINOSPECSpectrograph` (inherited by `MMTBINOSPECIFUSpectrograph`).
These coefficients are stable detector properties not expected to change
unless the detectors are replaced. The coefficients are an 8x5 array
(8 amplifiers, degree-4 polynomial):

```python
nonlinearity_coeffs = np.array([
    [0.00000000e+00, 1.00400089e+00, -1.39235362e-06, 8.31711824e-12, -1.20653479e-17],
    [0.00000000e+00, 1.00361458e+00, -1.29223833e-06, 6.93723177e-12, -9.67406255e-18],
    [0.00000000e+00, 1.00269542e+00, -9.29361806e-07, 5.97902827e-12, -2.30257302e-17],
    [0.00000000e+00, 1.00339616e+00, -8.47134521e-07, 7.92441693e-12, -4.46542834e-17],
    [0.00000000e+00, 1.00727205e+00, -1.69093388e-06, 2.07225055e-11, -1.62655178e-16],
    [0.00000000e+00, 1.00858745e+00, -2.35668901e-06, 2.40641019e-11, -1.50286358e-16],
    [0.00000000e+00, 1.00728526e+00, -1.80779473e-06, 1.73427719e-11, -1.01685780e-16],
    [0.00000000e+00, 1.00845168e+00, -2.02050567e-06, 2.97587091e-11, -2.65508521e-16],
])
```

The polynomial form is `C_corr = c[0] + c[1]*C + c[2]*C^2 + c[3]*C^3 + c[4]*C^4`
(constant term is zero, linear term ~1.004-1.009, with small higher-order
corrections). Compatible with `np.polynomial.polynomial.polyval`.

### 4. Apply Nonlinearity in `binospec_read_amp()`

After overscan subtraction and before returning data:

```python
data = np.polynomial.polynomial.polyval(data, nonlinearity_coeffs[ext - 1])
```

The amplifier extension index (1-8) is already passed to the function and
selects the correct coefficients.

**Order**: overscan subtraction -> nonlinearity correction -> (gain applied
later by PypeIt's general pipeline). This matches the IDL pipeline order
confirmed in `bino_mosaic.pro` lines 123-125: `remove_overscan` then
`poly(im_cur, nonlin_str)` then `im_cur * gain`.

Extensions 1-4 correspond to detector A (DET01), extensions 5-8 to
detector B (DET02). The nonlinearity coefficients array is ordered to match
(index 0 = ext 1, index 7 = ext 8).

### 5. No Other Pipeline Changes

- `get_rawimage()`: No changes needed.
- `default_pypeit_par()`: `correct_nonlinear` remains `None` (nonlinearity
  handled pre-assembly).
- `rawimage.py`, `procimg.py`, `pypeitpar.py`: No changes.
- PypeIt's general `apply_gain` step continues to apply gain after assembly,
  maintaining the correct order.

## Testing

- Compare PypeIt-processed images against IDL pipeline output for the same
  raw data to verify the overscan and nonlinearity corrections match.
- Existing unit tests in `pypeit/tests/` should continue to pass (no shared
  infrastructure is modified).
- Dev suite Binospec IFU tests should be re-run to verify end-to-end results.

## References

- IDL `remove_overscan` and `clean_overscan_vector`: `bino_mosaic.pro`
- IDL nonlinearity calibration: `calib_Bino/detector/nonlinearity/scicam_bino_sep2017.fits`
- PypeIt overscan infrastructure: `pypeit/images/rawimage.py`, `pypeit/core/procimg.py`
