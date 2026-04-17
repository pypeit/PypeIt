# Binospec Bad Pixel Mask Update — Design Spec

## Problem

PypeIt's Binospec BPM masks only a handful of hard-coded bad columns (2 for
det=1, 4 for det=2), while the IDL pipeline (`bino_mosaic.pro`) applies three
layers of masking totaling ~12,500 individual bad pixels plus bad columns and
detector trap regions. This causes artifacts in PypeIt reductions that are
absent from IDL-processed data.

## Solution

Generate pre-built compressed FITS BPM images from the IDL pipeline's masking
data, ship them as static calibration files, and load them in `bpm()`.

## Coordinate Mapping

The IDL mosaic and PypeIt BPM use different coordinate systems. Empirically
verified mapping:

- **Side A (det=1)**: `x_pyp = (x_idl + 2048) % 4096`, `y_pyp = y_idl`
- **Side B (det=2)**: `x_pyp = x_idl`, `y_pyp = 4111 - y_idl`

Both produce 4096×4112 images.

## IDL Masking Layers (All Three Ported)

### Layer 1: Calibration File (`badpix_binospec.fits`)

Individual bad pixel coordinates stored as (X, Y) tables:
- Extension 1 (Side A): 4,141 bad pixels
- Extension 2 (Side B): 8,335 bad pixels

Coordinates are in IDL mosaic space; must be transformed to PypeIt BPM space.

### Layer 2: Hard-Coded Bad Columns (`bino_mosaic.pro` lines 255-281)

**Side A:**
- Top half (y≥2056): 13 columns — `[1446, 1621, 1780, 1781, 1786, 1834, 1941, 1942, 1944, 1946, 1964, 2912, 3311]`
- Bottom half (y<2056): 2 columns — `[3236, 3391]`

**Side B:**
- Top half (y≥2056): 7 columns — `[113, 210, 2099, 3336, 3337, 4056, 4057]`

### Layer 3: Detector Trap Regions (Side A only, lines 260-277)

Rectangular regions at specific (x, y) ranges:
```
(480, 1635:1783), (366, 1670:1711), (402, 1480:1576), (403, 1450:1480),
(404:405, 1345:1624), (406:409, 1148:1472), (410, 1037:1150),
(411, 925:1011), (412, 790:906), (413, 718:790), (414, 617:722),
(415, 545:616), (1964, 323:470), (2748, 349:460), (3037:3038, 0:260)
```

Diagonal defect (lines 275-277):
```
yr1bad = 763 + range(1433 - 763 + 1)
xr1bad = int(2801.0 + (2795 - 2801) * (yr1bad - 763) / (1433. - 763.))
```

## File Storage

Two compressed FITS files in `pypeit/data/static_calibs/mmt_binospec/`:
- `bpm_binospec_det1.fits.gz` — Side A BPM (4096×4112, int8)
- `bpm_binospec_det2.fits.gz` — Side B BPM (4096×4112, int8)

These are shipped with the package (no remote download), following the
X-Shooter BPM precedent. Expected size: ~30-50KB each after gzip compression.

## BPM Generation

A one-time conversion script (`claude_docs/scripts/generate_binospec_bpm.py`)
reads the IDL calibration data, applies all three masking layers, transforms
coordinates to PypeIt space, and writes the FITS files. This script is not
shipped with PypeIt — only the output BPM files are.

## Code Changes

### `mmt_binospec.py` — `MMTBINOSPECSpectrograph.bpm()`

Replace the existing hard-coded pixel logic with:

```python
def bpm(self, filename, det, shape=None, msbias=None):
    bpm_img = super().bpm(filename, det, shape=shape, msbias=msbias)

    bpm_file = dataPaths.static_calibs.get_file_path(
        f'mmt_binospec/bpm_binospec_det{det}.fits.gz')
    static_bpm = fits.getdata(bpm_file)
    bpm_img |= static_bpm

    return bpm_img
```

Requires adding `from pypeit import dataPaths` to the imports in
`mmt_binospec.py`.

### Binning

Binospec IFU data is always unbinned (1×1). The static BPM files are at native
resolution (4096×4112). No binning support is implemented. If binned modes are
ever used, the BPM loading should be extended with block-max rebinning.

### FITS Shape Convention

The generation script must write the BPM array with NumPy shape `(4096, 4112)`
so that `fits.getdata()` returns the same shape as `super().bpm()`. The test
`test_bpm_shape` verifies the NumPy array shape after loading, not the FITS
NAXIS values.

### Dropped PypeIt-Only BPM Entries

The existing hard-coded values (det=1: x=2111, 2447; det=2: x=3011) that are
not in the IDL mask are dropped — they do not correspond to any known bad
pixels in the IDL calibration and appear to have been tentative entries (the
code has `# TODO: Fix this` on the det=1 block and `#ToDo: Need to double
check` on det=2). Values that overlap with IDL (det=2: 3336, 3337, 4056) are
already covered by the static files.

## Testing

New file `pypeit/tests/test_binospec_bpm.py` with:

1. **test_bpm_file_exists** — Both static BPM files exist in `static_calibs`
2. **test_bpm_shape** — Each file is (4096, 4112) int8
3. **test_bpm_nonzero_count** — Masked pixel count is in expected range
4. **test_bpm_applied** — `bpm()` incorporates the static mask (called with
   `shape` parameter to avoid requiring raw data files)

## Existing PypeIt BPM Entries

The current hard-coded entries will be removed:

| Detector | PypeIt Column | In IDL? | Action |
|----------|--------------|---------|--------|
| det=1 | 2447 (top) | No | Drop |
| det=1 | 2111 (top) | No | Drop |
| det=2 | 3336 (bot) | Yes | Covered by static file |
| det=2 | 3337 (bot) | Yes | Covered by static file |
| det=2 | 4056 (bot) | Yes | Covered by static file |
| det=2 | 3011 (top) | No | Drop |
