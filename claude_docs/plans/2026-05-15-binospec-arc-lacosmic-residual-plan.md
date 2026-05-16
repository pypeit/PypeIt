# Binospec arc/tilt CR cleanup via row-median residual — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the Hough-binning body of `MMTBINOSPECIFUSpectrograph.clean_calibration_image` with row-median subtraction followed by `procimg.lacosmic` on the residual (with `remove_compact_obj=False`, `objlim=0.0`). Delete the now-unused `_long_diagonal_arc_cr_mask` and its test, replacing the latter with an end-to-end test of `clean_calibration_image`.

**Architecture:** The detector pipeline keeps both the base `Spectrograph.clean_calibration_image` hook and its call sites in `pypeit/calibrations.py`. Only the Binospec-IFU override changes: instead of Hough-binning across slope candidates, the new body subtracts a row-local median (which removes near-horizontal arc/tilt line signal), then runs LA Cosmic on the residual without the compact-object discriminator. That turns LA Cosmic into a high-pass + threshold detector, which catches the body of an extended trail — the exact failure mode of Option 3.

**Tech Stack:** Python, NumPy, SciPy `ndimage`, PypeIt (`procimg.lacosmic`, `utils.inverse`, `MMTBINOSPECIFUSpectrograph`, `PypeItImage`).

**Spec reference:** `claude_docs/specs/2026-05-15-binospec-arc-lacosmic-residual-design.md`
**Supersedes the Option 3 plan at:** `claude_docs/plans/2026-05-15-binospec-arc-lacosmic-plan.md`

---

## Branch hygiene note

Before starting, run `git status` and note: the working tree already has uncommitted edits in `pypeit/spectrographs/mmt_binospec.py`, `pypeit/tests/test_fiber_block_extraction.py`, and two `doc/spectrographs/*.rst` files (tilt-order tuning and documentation work, unrelated to this plan). When committing Task 2 / Task 3 work, **stage by file path and review the staged hunks** (`git diff --cached`) before each commit — confirm the staged hunks match this plan and *include* whatever uncommitted edits to those files are already in place (the user has not asked to revert them). If a hunk is unfamiliar and not from this plan, stop and ask.

---

## File Structure

| File | Status | Responsibility |
| --- | --- | --- |
| `claude_docs/scratch/validate_arc_lacosmic_residual.py` | create (uncommitted) | Param sweep on residual-based detector; not added to git |
| `pypeit/spectrographs/mmt_binospec.py` | modify | New `clean_calibration_image` body + module-level CR constants; delete `_long_diagonal_arc_cr_mask` |
| `pypeit/tests/test_fiber_block_extraction.py` | modify | Replace `test_long_diagonal_arc_cr_mask` with `test_clean_calibration_image_residual_cr` |

Task 1 is the empirical validation gate. Tasks 2–3 are the TDD-style code change. Task 4 is the full-suite check. Task 5 is the optional end-to-end re-run on the JADES dataset.

---

## Task 1: Empirical validation of `(sigclip, maxiter, grow)` on the row-median residual

**Purpose:** Pick the LA Cosmic parameters that satisfy the spec §3.3 pass criteria on the residual. Validation output gates Task 2. If no combination passes, escalate per spec §3.4.

**Files:**
- Create (uncommitted): `claude_docs/scratch/validate_arc_lacosmic_residual.py`

- [ ] **Step 1: Create the scratch directory if missing**

```bash
mkdir -p /Users/tim/MMT/pypeit/claude_docs/scratch
```

- [ ] **Step 2: Write the validation script**

Create `/Users/tim/MMT/pypeit/claude_docs/scratch/validate_arc_lacosmic_residual.py` with:

```python
"""Validate (sigclip, maxiter, grow) for lacosmic on row-median residual.

Not part of the repository. Run from the pypeit project root.
"""
from pathlib import Path

import numpy as np
from scipy import ndimage

from pypeit.spectrographs.mmt_binospec import MMTBINOSPECIFUSpectrograph
from pypeit.images.buildimage import buildimage_fromlist
from pypeit.core import procimg
from pypeit import utils

RAW = Path("/Users/tim/MMT/bino_ifu/JADES_1031022/2025.0524/sci_img_2025.0524.054037.fits")
MEDIAN_WIDTH = 51
SIGCLIP_GRID = [4.0, 5.0, 6.0, 7.0, 8.0, 10.0]
MAXITER_GRID = [1, 2, 3]
GROW_GRID = [1.5, 2.0, 3.0]
ARC_LINE_ROW_PERCENTILE = 90.0

PASS_TRAIL_COVERAGE = 0.90
PASS_ARC_LINE_MASKED_FRACTION = 0.001
PASS_DET02_TOTAL_FRACTION = 0.005
PASS_DET01_TOTAL_FRACTION = 0.0005
PASS_TILT_EXTRA_FRACTION = 0.001


def build_calib_image(spec, det, frametype):
    par = spec.default_pypeit_par()
    frame_par = par["calibrations"][f"{frametype}frame"]
    frame_par["process"]["mask_cr"] = False
    bpm = spec.bpm(str(RAW), det)
    if bpm.ndim == 3:
        bpm = bpm[0]
    return buildimage_fromlist(spec, det, frame_par, [str(RAW)], bpm=bpm)


def _detector_value(detector, key):
    if isinstance(detector, list):
        detector = detector[0]
    try:
        return float(detector[key])
    except (TypeError, KeyError):
        return float(getattr(detector, key))


def arc_line_rows(image: np.ndarray) -> np.ndarray:
    row_medians = np.nanmedian(image, axis=1)
    threshold = np.nanpercentile(row_medians, ARC_LINE_ROW_PERCENTILE)
    return np.broadcast_to((row_medians >= threshold)[:, None], image.shape).copy()


def main() -> None:
    spec = MMTBINOSPECIFUSpectrograph()
    for det in (1, 2):
        for frametype in ("arc", "tilt"):
            print(f"\n=== DET{det:02d} {frametype} ===")
            calib = build_calib_image(spec, det, frametype)
            image = calib.image
            if image.ndim == 3:
                image = image[0]
            replacement = ndimage.median_filter(image, size=(1, MEDIAN_WIDTH), mode="nearest")
            resid = image - replacement
            var = utils.inverse(calib.ivar)
            if var.ndim == 3:
                var = var[0]

            if det == 2 and frametype == "arc":
                ref_mask, _ = spec._long_diagonal_arc_cr_mask(image)
                print(f"Reference (_long_diagonal_arc_cr_mask): {ref_mask.sum()} pixels")
            else:
                ref_mask = np.zeros(image.shape, dtype=bool)

            bpm = None
            if hasattr(calib, "fullmask") and calib.fullmask is not None:
                try:
                    bpm = calib.fullmask.flagged("BPM")
                    if bpm.ndim == 3:
                        bpm = bpm[0]
                except Exception:
                    bpm = None

            arc_rows = arc_line_rows(image)
            saturation = _detector_value(calib.detector, "saturation")
            nonlinear = _detector_value(calib.detector, "nonlinear")

            header = f"{'sigclip':>8} {'maxiter':>8} {'grow':>5} {'cov':>6} {'arc_frac':>9} {'total':>7}"
            print(header)
            for sigclip in SIGCLIP_GRID:
                for maxiter in MAXITER_GRID:
                    for grow in GROW_GRID:
                        crmask = procimg.lacosmic(
                            resid, saturation=saturation, nonlinear=nonlinear,
                            bpm=bpm, varframe=var,
                            sigclip=sigclip, sigfrac=0.3,
                            objlim=0.0, remove_compact_obj=False,
                            maxiter=maxiter, grow=grow,
                        )
                        ref_n = int(ref_mask.sum())
                        cov = float((ref_mask & crmask).sum()) / ref_n if ref_n else float("nan")
                        arc_pixels_total = int(arc_rows.sum())
                        arc_frac = (
                            int((crmask & arc_rows).sum()) / arc_pixels_total
                            if arc_pixels_total else float("nan")
                        )
                        total = float(crmask.sum()) / crmask.size
                        marker = ""
                        if (
                            det == 2 and frametype == "arc"
                            and cov >= PASS_TRAIL_COVERAGE
                            and arc_frac <= PASS_ARC_LINE_MASKED_FRACTION
                            and total <= PASS_DET02_TOTAL_FRACTION
                        ):
                            marker = "  <-- PASS"
                        if det == 1 and total <= PASS_DET01_TOTAL_FRACTION:
                            marker += "  [DET01 ok]"
                        if frametype == "tilt" and total <= PASS_TILT_EXTRA_FRACTION:
                            marker += "  [tilt ok]"
                        print(
                            f"{sigclip:>8.1f} {maxiter:>8d} {grow:>5.1f} "
                            f"{cov:>6.3f} {arc_frac:>9.5f} {total:>7.5f}{marker}"
                        )


if __name__ == "__main__":
    main()
```

- [ ] **Step 3: Run the validation script**

```bash
cd /Users/tim/MMT/pypeit
python claude_docs/scratch/validate_arc_lacosmic_residual.py 2>&1 \
    | grep -vE "^\s*\[INFO\]" \
    | tee claude_docs/scratch/validate_arc_lacosmic_residual_output.txt
```

The `grep -v` filters out LA Cosmic's per-iteration `[INFO]` log spam so the output table is readable. Expected runtime: 5–15 minutes on a developer laptop (4 detector/frametype combos × 54 lacosmic calls).

- [ ] **Step 4: Pick the parameters**

From the output, choose the least aggressive `(sigclip, maxiter, grow)` that satisfies all four pass criteria for the relevant rows:

1. **DET02 / arc:** `<-- PASS` marker is present.
2. **DET01 / arc:** `[DET01 ok]` marker is present.
3. **DET02 / tilt:** `[tilt ok]` marker is present.
4. **DET01 / tilt:** `[tilt ok]` *and* `[DET01 ok]` markers are both present.

Tie-breaking (least aggressive first): highest `sigclip` → lowest `grow` → lowest `maxiter`.

Fill the chosen values into this plan inline by replacing the block below before proceeding:

```
CHOSEN_SIGCLIP = 10.0
CHOSEN_MAXITER = 2
CHOSEN_GROW    = 2.0
```

**Note:** the spec's strict pass criteria (arc_frac ≤ 0.001, DET01 total ≤ 0.0005, tilt total ≤ 0.001) were not reachable with `objlim=0, remove_compact_obj=False` at any setting in the grid; the residual is not pure Gaussian noise (sub-σ leakage on arc rows) and DET01 also contains many real compact CRs/hot pixels that the detector legitimately flags. The criteria were loosened per user decision: replacement with the row-median is benign on arc-line and tilt-line pixels, so a higher `total` mask rate is acceptable. Selection rationale: highest `sigclip` (=10) that still yields DET02 arc coverage ≥ 0.95 (achieved cov=0.979 at maxiter=2, grow=2.0); maxiter=2 over maxiter=1 picks up the trail edges; grow=2.0 keeps mask compact. DET01 total at this setting is 0.0024 (0.24% of pixels), tilt totals are the same.

- [ ] **Step 5: Failure case**

- If **only the tilt criteria** (criteria 3 and/or 4) fail for every parameter combo, set `CHOSEN_*` based on the arc-only criteria and proceed to Task 2 with the implementation restricted to `frametype == "arc"` (spec §3.4 fallback). Note this restriction in the commit message and add a one-line comment in the new `clean_calibration_image` body.
- If **DET02 arc** itself fails (criterion 1), stop. Report BLOCKED and ask the user before continuing — this is the case the spec said requires rethinking.

- [ ] **Step 6: Commit the plan with the chosen values filled in**

```bash
cd /Users/tim/MMT/pypeit
git add claude_docs/plans/2026-05-15-binospec-arc-lacosmic-residual-plan.md
git diff --cached --stat
git commit -m "docs: record validated lacosmic-on-residual params for Binospec arc/tilt"
```

Verify with `git status` that `claude_docs/scratch/` is untracked (not staged).

---

## Task 2: TDD rewrite of `clean_calibration_image`

**Files:**
- Modify: `pypeit/spectrographs/mmt_binospec.py` (around lines 1364–1501)
- Modify: `pypeit/tests/test_fiber_block_extraction.py` (replace `test_long_diagonal_arc_cr_mask`)

In every code/test block below, substitute `<CHOSEN_SIGCLIP>`, `<CHOSEN_MAXITER>`, `<CHOSEN_GROW>` with the literal numeric values recorded in Task 1 Step 4.

- [ ] **Step 1: Write the new failing test**

In `pypeit/tests/test_fiber_block_extraction.py`, locate `test_long_diagonal_arc_cr_mask` (around lines 226–245) and replace it with:

```python
def test_clean_calibration_image_residual_cr():
    """clean_calibration_image must catch a synthetic diagonal CR
    without masking synthetic arc lines."""
    from types import SimpleNamespace

    rng = np.random.default_rng(7)
    image = 10.0 + rng.normal(scale=1.0, size=(120, 160))
    image[35, :] += 500.0   # synthetic arc line at row 35
    image[85, :] += 300.0   # synthetic arc line at row 85
    xs = np.arange(20, 135)
    ys = np.rint(15 + 0.48 * xs).astype(int)
    image[ys, xs] += 1000.0  # synthetic diagonal CR trail

    arc_line_pixels = image[[35, 85]].copy()

    bpm_2d = np.zeros(image.shape, dtype=bool)
    fullmask = SimpleNamespace(flagged=lambda flag: bpm_2d.copy())
    calib_image = SimpleNamespace(
        image=image,
        ivar=None,
        fullmask=fullmask,
        update_mask_cr=lambda mask: None,
    )

    spec.clean_calibration_image(calib_image, "arc", det=2)

    # Trail pixels are replaced by the row-median (~10), no longer at +1000.
    assert np.nanmedian(calib_image.image[ys, xs]) < 50.0
    # Arc-line rows are untouched.
    np.testing.assert_allclose(calib_image.image[35, :], arc_line_pixels[0], atol=2.0)
    np.testing.assert_allclose(calib_image.image[85, :], arc_line_pixels[1], atol=2.0)
```

The `SimpleNamespace` stub keeps the test pure module-level — no `PypeItImage` construction required and no DataContainer plumbing.

- [ ] **Step 2: Run the test and verify it fails**

```bash
cd /Users/tim/MMT/pypeit
pytest pypeit/tests/test_fiber_block_extraction.py::test_clean_calibration_image_residual_cr -v
```

Expected: FAIL. With the *existing* `_long_diagonal_arc_cr_mask`-based `clean_calibration_image`, the test will fail because:
- The existing detector has `det != 2` skip behavior keyed on DET02-only — the test passes `det=2` so this part is fine.
- BUT the existing detector uses `min_length=500`, while our synthetic trail is only ~120 pixels long (xs from 20 to 135) — so the existing detector will return early without flagging anything.
- Result: cleaned trail pixels still have +1000 added, `np.nanmedian` will be >>50.

If the test happens to PASS against the existing code, stop and double-check the synthetic geometry — the test must fail before the new implementation goes in.

- [ ] **Step 3: Replace `clean_calibration_image` body and delete `_long_diagonal_arc_cr_mask`**

In `pypeit/spectrographs/mmt_binospec.py`:

First, near the top of the `MMTBINOSPECIFUSpectrograph` class (just below the existing class-level constants, e.g., right after the `sky_fiber_indices_0based` block around line 1140 — pick a position adjacent to other CR-related constants if any exist), add:

```python
    # Cosmic-ray cleanup constants for arc/tilt frames.  See
    # claude_docs/specs/2026-05-15-binospec-arc-lacosmic-residual-design.md.
    # Row-median subtraction removes the near-horizontal arc/tilt line
    # signal so that the residual contains only noise plus CRs; LA Cosmic
    # then runs without compact-object filtering as a simple high-pass
    # threshold detector that catches trail bodies (which Laplacian S/N
    # alone misses on raw arc/tilt images).
    _ARC_CR_MEDIAN_WIDTH = 51
    _ARC_CR_SIGCLIP = <CHOSEN_SIGCLIP>
    _ARC_CR_MAXITER = <CHOSEN_MAXITER>
    _ARC_CR_GROW = <CHOSEN_GROW>
```

Then replace the body of `clean_calibration_image` (currently lines 1364–1385) AND delete the entire `_long_diagonal_arc_cr_mask` method that follows (currently lines 1387–1501). The replacement:

```python
    def clean_calibration_image(self, calib_image, frametype, det):
        """
        Remove bright cosmic-ray artifacts from Binospec IFU arc/tilt images.

        Arc and tilt frames are taken as single exposures per observing block
        and cannot rely on median-combine CR rejection.  Standard LA Cosmic on
        the raw frame either misses the body of long thin trails (Laplacian
        S/N is low for non-edge pixels along an extended feature) or masks
        arc lines themselves.  Subtracting a row-local median first removes
        the (near-)horizontal arc/tilt line signal, leaving a residual on
        which sigma-clipped LA Cosmic cleanly identifies CRs.
        """
        if frametype not in ("arc", "tilt"):
            return calib_image

        from scipy import ndimage

        image = calib_image.image
        replacement = ndimage.median_filter(
            image, size=(1, self._ARC_CR_MEDIAN_WIDTH), mode="nearest")
        resid = image - replacement

        if getattr(calib_image, "ivar", None) is not None:
            var = utils.inverse(calib_image.ivar)
        else:
            var = np.maximum(np.abs(image), 1.0)

        bpm = None
        fullmask = getattr(calib_image, "fullmask", None)
        if fullmask is not None:
            try:
                bpm = fullmask.flagged("BPM")
            except Exception:
                bpm = None

        crmask = procimg.lacosmic(
            resid, varframe=var, bpm=bpm,
            sigclip=self._ARC_CR_SIGCLIP, sigfrac=0.3,
            objlim=0.0, remove_compact_obj=False,
            maxiter=self._ARC_CR_MAXITER, grow=self._ARC_CR_GROW,
        )
        if not np.any(crmask):
            return calib_image

        calib_image.image[crmask] = replacement[crmask]
        if hasattr(calib_image, "update_mask_cr"):
            calib_image.update_mask_cr(crmask)
        log.info(
            f"Cleaned {int(crmask.sum())} CR pixels from Binospec IFU "
            f"{frametype} image on DET{det:02d} "
            f"(row-median residual + LA Cosmic)."
        )
        return calib_image
```

(If Task 1 §5 set the implementation to be arc-only, change the frametype check to `if frametype != "arc": return calib_image` and add a one-line `# Tilt frames not handled here — see plan Task 1 fallback.` comment above it.)

Verify required imports exist at the top of `mmt_binospec.py`:
- `from pypeit.core import procimg`
- `from pypeit import utils`
- `from pypeit import log`
- `import numpy as np`

If any of these is missing, add it in the standard import block at the top of the file.

- [ ] **Step 4: Run the test and verify it passes**

```bash
cd /Users/tim/MMT/pypeit
pytest pypeit/tests/test_fiber_block_extraction.py::test_clean_calibration_image_residual_cr -v
```

Expected: PASS.

- [ ] **Step 5: Run the full file**

```bash
cd /Users/tim/MMT/pypeit
pytest pypeit/tests/test_fiber_block_extraction.py -v
```

Expected: all tests PASS — including any existing tests untouched by this task. If a previously-passing test now fails, investigate before proceeding.

- [ ] **Step 6: Commit**

```bash
cd /Users/tim/MMT/pypeit
git add pypeit/spectrographs/mmt_binospec.py pypeit/tests/test_fiber_block_extraction.py
git diff --cached --stat
```

Review `git diff --cached` and confirm the staged hunks cover (a) the new constants block, (b) the rewritten `clean_calibration_image`, (c) the deletion of `_long_diagonal_arc_cr_mask`, (d) the replacement of the test function — plus any pre-existing uncommitted edits to these files that were on the branch already (those carry over). Then:

```bash
git commit -m "replace Binospec arc/tilt CR detector with row-median + LA Cosmic"
```

---

## Task 3: Full unit-test suite check

- [ ] **Step 1: Run the unit-test directory**

```bash
cd /Users/tim/MMT/pypeit
pytest pypeit/tests/ -x -q
```

Expected: all PASS. `-x` stops on first failure so any regression is immediately visible.

- [ ] **Step 2: If anything fails, investigate**

The most likely culprit for a non-Binospec test failure would be an unintended import-time change in `mmt_binospec.py`. Check `git diff HEAD~1 pypeit/spectrographs/mmt_binospec.py` for anything outside the planned hunks.

- [ ] **Step 3: No commit needed if the suite passes (no code changed in this task).**

---

## Task 4: Documentation update

**Files:**
- Modify: `doc/spectrographs/mmt_binospec.rst` *or* `doc/spectrographs/mmt_binospec_pipeline_comparison.rst` — wherever the existing CR-handling notes live

- [ ] **Step 1: Locate any existing notes**

```bash
cd /Users/tim/MMT/pypeit
grep -n "long diagonal\|long.diagonal\|clean_calibration_image\|cosmic" \
    doc/spectrographs/mmt_binospec.rst \
    doc/spectrographs/mmt_binospec_pipeline_comparison.rst 2>&1
```

- [ ] **Step 2: Rewrite the notes**

If any match exists from Step 1, edit it to describe the new behavior: arc/tilt frames go through `clean_calibration_image`, which subtracts a row-local median (width 51 px) and runs LA Cosmic on the residual with `sigclip=<CHOSEN_SIGCLIP>`, `maxiter=<CHOSEN_MAXITER>`, `grow=<CHOSEN_GROW>`, `objlim=0`, `remove_compact_obj=False`. Briefly mention that this catches the trail body that standard LA Cosmic structurally cannot see on raw arc images. If no existing notes match, skip this step.

- [ ] **Step 3: Add a changelog entry if a current dev release file exists**

```bash
ls /Users/tim/MMT/pypeit/doc/releases/ | tail -3
```

If a current `*dev.rst` file exists, add under "Spectrograph Changes" or similar:

```rst
- Binospec IFU: ``clean_calibration_image`` now uses row-median
  residual + LA Cosmic to clean cosmic-ray trails on arc and tilt
  frames; replaces the previous Hough-binning detector.
```

If no dev release file exists, skip.

- [ ] **Step 4: Commit (only if there were doc changes)**

```bash
cd /Users/tim/MMT/pypeit
git add doc/
git diff --cached --stat
git commit -m "docs: note Binospec arc/tilt CR cleanup via lacosmic-on-residual"
```

---

## Task 5: End-to-end re-run on JADES dataset (optional but recommended)

This is a manual verification against real data. Does not change repository contents.

- [ ] **Step 1: Re-run the Binospec IFU calibration on the JADES setup**

```bash
cd /Users/tim/MMT/bino_ifu/JADES_1031022/mmt_binospec_ifu_A
run_pypeit mmt_binospec_ifu_A.pypeit -o -r calib
```

(Adjust the `-r` flag if your installed `run_pypeit` uses a different option to limit the run to calibrations.)

- [ ] **Step 2: Compare the new arc/tilt images against the previously-committed ones**

```bash
cd /Users/tim/MMT/bino_ifu/JADES_1031022/mmt_binospec_ifu_A
python -c "
from astropy.io import fits
import numpy as np
for name in ['Arc_A_0_DET02.fits', 'Tiltimg_A_0_DET02.fits']:
    new = fits.getdata(name, ext=1)
    old = fits.getdata('../Calibrations/' + name, ext=1)
    diff = new - old
    print(f'{name}: max |diff| = {np.nanmax(np.abs(diff)):.1f}, '
          f'rms diff = {np.nanstd(diff):.2f}, '
          f'masked pixels (new) = {(np.isnan(new) | (new == 0)).sum()}')
"
```

Expected: differences at the CR trail location (both methods masked it but the new code may pick a slightly larger or smaller pixel set). RMS diff elsewhere near noise. Adjust paths if the reference Calibrations location differs.

- [ ] **Step 3: Inspect the tilt-tracing log**

```bash
grep -i "cosmic\|tilt.*track\|wave.*fit" mmt_binospec_ifu_A.log | tail -40
```

Expected: any CR-related messages now come from PypeIt's `procimg.lacosmic` and the new `clean_calibration_image` log line. No new tilt-tracing warnings versus the previous run.

- [ ] **Step 4: If a regression appears**

If trail removal looks worse than before, or arc lines appear masked in `Arc_A_0_DET02.fits`, revisit the chosen `(sigclip, maxiter, grow)` from Task 1. Either widen the grid, tighten one of the pass thresholds, or accept the validation as a lower-confidence pass and loop back to update Task 2's constants (one focused fixup commit on top is fine).

---

## Self-Review Notes

- Spec §1 motivation + scope reminder — context only; no implementation step needed.
- Spec §2.1 implementation — covered by Task 2 Step 3.
- Spec §2.2 deletions (`_long_diagonal_arc_cr_mask`) — covered by Task 2 Step 3.
- Spec §2.3 kept (base hook, call sites) — by *not* removing them; Task 3 confirms no test regressions on `pypeit/spectrographs/spectrograph.py` or `pypeit/calibrations.py`.
- Spec §3 empirical validation — Task 1.
- Spec §3.3 four pass criteria — Task 1 Step 4 enumerates all four; Task 1 Step 5 handles each failure mode.
- Spec §3.4 fallback (arc-only) — Task 1 Step 5 and Task 2 Step 3 final paragraph.
- Spec §4 test replacement — Task 2 Steps 1, 2, 4.
- Spec §5 risks — manual end-to-end check in Task 5 surfaces real-data regressions (tilt residual leakage, dense CR clusters); the `objlim=0.0` rationale is captured in the constants comment in Task 2 Step 3.
- Spec §6 out of scope — none of the tasks touch science frames, MOS arcs, or other spectrographs.
