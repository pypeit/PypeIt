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
    with fits.open(IDL_BADPIX) as hdu:
        bpix = hdu[1].data  # Side A: 4141 pixels
        x_idl, y_idl = bpix['X'], bpix['Y']
        x_pyp = (x_idl + 2048) % 4096
        y_pyp = y_idl
        bpm[x_pyp, y_pyp] = 1

    # --- Layer 2: Hard-coded bad columns (bino_mosaic.pro lines 255-258) ---
    # IDL: bx1top=[1447,1622,...]-1  ->  0-indexed columns, top half (y>=2056)
    bx1top = np.array([1447, 1622, 1781, 1782, 1787, 1835,
                       1942, 1943, 1945, 1947, 1965, 2913, 3312]) - 1
    for x_idl in bx1top:
        x_pyp = (x_idl + 2048) % 4096
        bpm[x_pyp, 2056:] = 1

    # IDL: bx1bot=[3237,3392]-1  ->  bottom half (y<2056)
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
    with fits.open(IDL_BADPIX) as hdu:
        bpix = hdu[2].data  # Side B: 8335 pixels
        x_idl, y_idl = bpix['X'], bpix['Y']
        x_pyp = x_idl
        y_pyp = 4111 - y_idl
        bpm[x_pyp, y_pyp] = 1

    # --- Layer 2: Hard-coded bad columns (bino_mosaic.pro line 279-280) ---
    # IDL: bx2top=[114,211,2100,3337,3338,4057,4058]-1  ->  top half (y>=2056)
    bx2top = np.array([114, 211, 2100, 3337, 3338, 4057, 4058]) - 1
    for x_idl in bx2top:
        # IDL top half y>=2056 -> PypeIt y = 4111-y -> y <= 2055 (bottom half)
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
