# Binospec Bad Pixel Mask Update — Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace Binospec's incomplete hard-coded BPM with comprehensive static BPM images ported from the IDL pipeline.

**Architecture:** A one-time generation script converts IDL bad pixel data (calibration file + hard-coded columns + detector trap regions) to PypeIt coordinates and writes compressed FITS BPM images. These are stored in `pypeit/data/static_calibs/mmt_binospec/` and loaded by a rewritten `bpm()` method.

**Tech Stack:** numpy, astropy.io.fits, pypeit.dataPaths

**Spec:** `claude_docs/specs/2026-03-16-binospec-bpm-design.md`

---

## File Structure

| File | Action | Purpose |
|------|--------|---------|
| `claude_docs/scripts/generate_binospec_bpm.py` | Create | One-time script to generate BPM FITS files from IDL data |
| `pypeit/data/static_calibs/mmt_binospec/bpm_binospec_det1.fits.gz` | Create | Static BPM for Side A (det=1) |
| `pypeit/data/static_calibs/mmt_binospec/bpm_binospec_det2.fits.gz` | Create | Static BPM for Side B (det=2) |
| `pypeit/spectrographs/mmt_binospec.py` | Modify | Rewrite `bpm()` to load static files, add `dataPaths` import |
| `pypeit/tests/test_binospec_bpm.py` | Create | Tests for BPM files and loading |

---

### Task 1: Write and run the BPM generation script

**Files:**
- Create: `claude_docs/scripts/generate_binospec_bpm.py`
- Create: `pypeit/data/static_calibs/mmt_binospec/bpm_binospec_det1.fits.gz`
- Create: `pypeit/data/static_calibs/mmt_binospec/bpm_binospec_det2.fits.gz`

- [ ] **Step 1: Create the generation script**

The script reads the IDL calibration data, applies all three masking layers,
transforms coordinates to PypeIt space, and writes compressed FITS files.

IDL source file: `~/MMT/binospec/pipeline/calib_Bino/detector/badpix_binospec.fits`

Coordinate mappings (empirically verified):
- Side A (det=1): `x_pyp = (x_idl + 2048) % 4096`, `y_pyp = y_idl`
- Side B (det=2): `x_pyp = x_idl`, `y_pyp = 4111 - y_idl`

```python
#!/usr/bin/env python
"""
Generate Binospec BPM FITS files from IDL pipeline calibration data.

Reads:
  ~/MMT/binospec/pipeline/calib_Bino/detector/badpix_binospec.fits

Writes:
  pypeit/data/static_calibs/mmt_binospec/bpm_binospec_det1.fits.gz
  pypeit/data/static_calibs/mmt_binospec/bpm_binospec_det2.fits.gz

Run once from the repo root:
  python claude_docs/scripts/generate_binospec_bpm.py
"""
import os
import numpy as np
from astropy.io import fits
from pathlib import Path

# Paths
IDL_BADPIX = os.path.expanduser(
    '~/MMT/binospec/pipeline/calib_Bino/detector/badpix_binospec.fits')
OUTDIR = Path('pypeit/data/static_calibs/mmt_binospec')

NX, NY = 4096, 4112  # Detector dimensions


def build_side_a_bpm():
    """Build BPM for Side A (det=1) in PypeIt coordinates."""
    bpm = np.zeros((NX, NY), dtype=np.int8)

    # --- Layer 1: Individual bad pixels from calibration file ---
    hdu = fits.open(IDL_BADPIX)
    bpix = hdu[1].data  # Side A: 4141 pixels
    x_idl, y_idl = bpix['X'], bpix['Y']
    x_pyp = (x_idl + 2048) % 4096
    y_pyp = y_idl
    bpm[x_pyp, y_pyp] = 1
    hdu.close()

    # --- Layer 2: Hard-coded bad columns (bino_mosaic.pro lines 255-258) ---
    # IDL: bx1top=[1447,1622,...]-1  →  0-indexed columns, top half (y>=2056)
    bx1top = np.array([1447, 1622, 1781, 1782, 1787, 1835,
                       1942, 1943, 1945, 1947, 1965, 2913, 3312]) - 1
    for x_idl in bx1top:
        x_pyp = (x_idl + 2048) % 4096
        bpm[x_pyp, 2056:] = 1

    # IDL: bx1bot=[3237,3392]-1  →  bottom half (y<2056)
    bx1bot = np.array([3237, 3392]) - 1
    for x_idl in bx1bot:
        x_pyp = (x_idl + 2048) % 4096
        bpm[x_pyp, :2056] = 1

    # --- Layer 3: Detector trap regions (Side A only, lines 260-277) ---
    # All coordinates are IDL mosaic space; transform x only.
    trap_regions = [
        # (x_slice, y_slice) in IDL coordinates
        (slice(480, 481), slice(1635, 1784)),
        (slice(366, 367), slice(1670, 1712)),
        (slice(402, 403), slice(1480, 1577)),
        (slice(403, 404), slice(1450, 1481)),
        (slice(404, 406), slice(1345, 1625)),
        (slice(406, 410), slice(1148, 1473)),
        (slice(410, 411), slice(1037, 1151)),
        (slice(411, 412), slice(925, 1012)),
        (slice(412, 413), slice(790, 907)),
        (slice(413, 414), slice(718, 791)),
        (slice(414, 415), slice(617, 723)),
        (slice(415, 416), slice(545, 617)),
        (slice(1964, 1965), slice(323, 471)),
        (slice(2748, 2749), slice(349, 461)),
        (slice(3037, 3039), slice(0, 261)),
    ]
    for x_sl, y_sl in trap_regions:
        x_indices = np.arange(x_sl.start, x_sl.stop)
        y_indices = np.arange(y_sl.start, y_sl.stop)
        x_pyp = (x_indices + 2048) % 4096
        bpm[np.ix_(x_pyp, y_indices)] = 1

    # Diagonal defect (lines 275-277)
    yr1bad = np.arange(763, 1434)
    xr1bad = (2801.0 + (2795 - 2801) * (yr1bad - 763) / (1433. - 763.)).astype(int)
    x_pyp = (xr1bad + 2048) % 4096
    bpm[x_pyp, yr1bad] = 1

    return bpm


def build_side_b_bpm():
    """Build BPM for Side B (det=2) in PypeIt coordinates."""
    bpm = np.zeros((NX, NY), dtype=np.int8)

    # --- Layer 1: Individual bad pixels from calibration file ---
    hdu = fits.open(IDL_BADPIX)
    bpix = hdu[2].data  # Side B: 8335 pixels
    x_idl, y_idl = bpix['X'], bpix['Y']
    x_pyp = x_idl
    y_pyp = 4111 - y_idl
    bpm[x_pyp, y_pyp] = 1
    hdu.close()

    # --- Layer 2: Hard-coded bad columns (bino_mosaic.pro line 279-280) ---
    # IDL: bx2top=[114,211,2100,3337,3338,4057,4058]-1  →  top half (y>=2056)
    bx2top = np.array([114, 211, 2100, 3337, 3338, 4057, 4058]) - 1
    for x_idl in bx2top:
        # IDL top half y>=2056 → PypeIt y = 4111-y → y <= 2055 (bottom half)
        bpm[x_idl, :2056] = 1

    # No Layer 3 for Side B (no trap regions in IDL code)

    return bpm


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    for det, builder in [(1, build_side_a_bpm), (2, build_side_b_bpm)]:
        bpm = builder()
        n_bad = np.sum(bpm > 0)
        outfile = OUTDIR / f'bpm_binospec_det{det}.fits.gz'
        hdu = fits.PrimaryHDU(data=bpm)
        hdu.header['DETECTOR'] = (det, 'PypeIt detector number')
        hdu.header['SIDE'] = ('A' if det == 1 else 'B', 'Binospec detector side')
        hdu.header['ORIGIN'] = ('IDL pipeline', 'Source of bad pixel data')
        hdu.writeto(outfile, overwrite=True)
        print(f'det={det}: {n_bad} bad pixels -> {outfile} '
              f'({outfile.stat().st_size / 1024:.1f} KB)')


if __name__ == '__main__':
    main()
```

- [ ] **Step 2: Run the generation script**

Run: `python claude_docs/scripts/generate_binospec_bpm.py`

Expected output: Two `.fits.gz` files created with bad pixel counts and file sizes printed.

- [ ] **Step 3: Verify the generated files**

```bash
python -c "
from astropy.io import fits
import numpy as np
for det in [1, 2]:
    f = f'pypeit/data/static_calibs/mmt_binospec/bpm_binospec_det{det}.fits.gz'
    bpm = fits.getdata(f)
    print(f'det={det}: shape={bpm.shape}, dtype={bpm.dtype}, '
          f'n_bad={np.sum(bpm > 0)}')
"
```

Expected: Both files have shape `(4096, 4112)`, dtype `int8`, and thousands of bad pixels.

- [ ] **Step 4: Commit**

```bash
git add claude_docs/scripts/generate_binospec_bpm.py \
        pypeit/data/static_calibs/mmt_binospec/bpm_binospec_det1.fits.gz \
        pypeit/data/static_calibs/mmt_binospec/bpm_binospec_det2.fits.gz
git commit -m "Add static BPM files for Binospec from IDL pipeline calibration"
```

---

### Task 2: Write BPM tests

**Files:**
- Create: `pypeit/tests/test_binospec_bpm.py`

- [ ] **Step 1: Write the test file**

```python
"""Tests for Binospec static bad pixel masks."""
import numpy as np
import pytest
from astropy.io import fits

from pypeit import dataPaths
from pypeit.spectrographs.mmt_binospec import MMTBINOSPECSpectrograph


class TestBinospecBPMFiles:
    """Tests for the static BPM FITS files."""

    @pytest.mark.parametrize('det', [1, 2])
    def test_bpm_file_exists(self, det):
        """Both static BPM files exist in static_calibs."""
        bpm_path = dataPaths.static_calibs.get_file_path(
            f'mmt_binospec/bpm_binospec_det{det}.fits.gz')
        assert bpm_path.exists(), f'BPM file not found: {bpm_path}'

    @pytest.mark.parametrize('det', [1, 2])
    def test_bpm_shape(self, det):
        """Each BPM file has the correct NumPy shape and dtype."""
        bpm_path = dataPaths.static_calibs.get_file_path(
            f'mmt_binospec/bpm_binospec_det{det}.fits.gz')
        bpm = fits.getdata(bpm_path)
        assert bpm.shape == (4096, 4112)
        assert bpm.dtype == np.int8

    @pytest.mark.parametrize('det,min_bad,max_bad', [
        (1, 30000, 80000),
        (2, 15000, 50000),
    ])
    def test_bpm_nonzero_count(self, det, min_bad, max_bad):
        """Each BPM has a reasonable number of masked pixels."""
        bpm_path = dataPaths.static_calibs.get_file_path(
            f'mmt_binospec/bpm_binospec_det{det}.fits.gz')
        bpm = fits.getdata(bpm_path)
        n_bad = np.sum(bpm > 0)
        assert min_bad < n_bad < max_bad, \
            f'det={det}: {n_bad} bad pixels outside expected range [{min_bad}, {max_bad}]'


class TestBinospecBPMApplied:
    """Test that bpm() loads and applies the static mask."""

    @pytest.mark.parametrize('det', [1, 2])
    def test_bpm_applied(self, det):
        """bpm() returns a mask that includes the static BPM pixels."""
        spec = MMTBINOSPECSpectrograph()
        # Use shape parameter to avoid needing a raw data file
        bpm_img = spec.bpm(filename=None, det=det, shape=(4096, 4112))

        # Load the static file directly for comparison
        bpm_path = dataPaths.static_calibs.get_file_path(
            f'mmt_binospec/bpm_binospec_det{det}.fits.gz')
        static_bpm = fits.getdata(bpm_path)

        # Every pixel marked in the static file should be marked in bpm_img
        assert np.all(bpm_img[static_bpm > 0] > 0)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest pypeit/tests/test_binospec_bpm.py -v`

Expected: `test_bpm_file_exists` and `test_bpm_shape` pass (files exist from Task 1).
`test_bpm_applied` fails because `bpm()` still uses hard-coded logic and doesn't
load the static files.

- [ ] **Step 3: Commit**

```bash
git add pypeit/tests/test_binospec_bpm.py
git commit -m "Add tests for Binospec static bad pixel masks"
```

---

### Task 3: Rewrite bpm() to load static files

**Files:**
- Modify: `pypeit/spectrographs/mmt_binospec.py:22-32` (imports)
- Modify: `pypeit/spectrographs/mmt_binospec.py:341-406` (`bpm()` method)

- [ ] **Step 1: Add the dataPaths import**

Add `from pypeit import dataPaths` to the imports in `mmt_binospec.py`
(after line 22, alongside the other `from pypeit import ...` lines):

```python
from pypeit import dataPaths
```

- [ ] **Step 2: Rewrite the bpm() method**

Replace the entire `bpm()` method (lines 341-406) with:

```python
    def bpm(self, filename, det, shape=None, msbias=None):
        """
        Generate a default bad-pixel mask.

        Loads a pre-built static BPM from the IDL pipeline calibration data
        (``badpix_binospec.fits`` + hard-coded bad columns and detector trap
        regions from ``bino_mosaic.pro``).

        Even though they are both optional, either the precise shape for
        the image (``shape``) or an example file that can be read to get
        the shape (``filename`` using :func:`get_image_shape`) *must* be
        provided.

        Args:
            filename (:obj:`str` or None):
                An example file to use to get the image shape.
            det (:obj:`int`):
                1-indexed detector number to use when getting the image
                shape from the example file.
            shape (tuple, optional):
                Processed image shape
                Required if filename is None
                Ignored if filename is not None
            msbias (`numpy.ndarray`_, optional):
                Processed bias frame used to identify bad pixels

        Returns:
            `numpy.ndarray`_: An integer array with a masked value set
            to 1 and an unmasked value set to 0.  All values are set to
            0.
        """
        # Call the base-class method to generate the empty bpm
        bpm_img = super().bpm(filename, det, shape=shape, msbias=msbias)

        # Load and apply the static BPM from IDL pipeline calibration
        bpm_file = dataPaths.static_calibs.get_file_path(
            f'mmt_binospec/bpm_binospec_det{det}.fits.gz')
        static_bpm = fits.getdata(bpm_file)
        bpm_img |= static_bpm

        return bpm_img
```

- [ ] **Step 3: Run the BPM tests**

Run: `pytest pypeit/tests/test_binospec_bpm.py -v`

Expected: All tests pass.

- [ ] **Step 4: Run the full Binospec test suite**

Run: `pytest pypeit/tests/test_binospec_overscan.py pypeit/tests/test_binospec_bpm.py -v`

Expected: All tests pass.

- [ ] **Step 5: Update release notes**

Add a bullet under the "Instrument-specific Updates" section in
`doc/releases/2.1.0dev.rst`, after the existing Binospec overscan bullet:

```rst
- Updated the MMT Binospec bad pixel mask to match the IDL pipeline, adding
  ~12,500 individual bad pixels, additional bad columns, and detector trap
  region masking from static calibration files.
```

- [ ] **Step 6: Commit**

```bash
git add pypeit/spectrographs/mmt_binospec.py \
        doc/releases/2.1.0dev.rst
git commit -m "Rewrite Binospec bpm() to load static BPM from IDL calibration"
```
