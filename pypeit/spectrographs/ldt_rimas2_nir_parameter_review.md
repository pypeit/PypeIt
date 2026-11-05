# LDT/RIMAS2 near-IR parameter review

This is a read-only comparison of the active and commented-out parameter
overrides in `ldt_rimas2.py` as of 2026-08-24.  The closest useful comparison
spectrographs are the long-slit NIR instruments Gemini/Flamingos,
MMT/MMIRS, LBT/LUCI, and the TripleSpec variants (ARC, P200, and SOAR).
The RIMAS grism mode is also compared with cross-dispersed TripleSpec,
GNIRS, and FIRE where the echelle-specific settings are relevant.

The document is not a prescription to make every setting match another
instrument: slit geometry, detector behavior, and the available calibration
frames must take precedence over numerical similarity.

## Shared RIMAS defaults

| Area | Active RIMAS2 setting | Commented alternative | Near-IR precedent | Selection suggestion |
|---|---|---|---|---|
| Calibration darks | Use a dark for arcs, tilts, and standards; do not CR-mask darks | None | Normal for IR detectors; the TripleSpec modes use dark processing | Keep. This matches the thermal/dark-current character of NIR data. Verify a same-exposure dark exists; do not turn CR masking on for darks without checking hot-pixel masks. |
| Wavelength lamps | `OH_GNIRS` by default | None | Flamingos, LUCI, and TripleSpec similarly use OH/sky-line lists | Keep for HK/grism sky calibration. The VPH-YJ Xe override is appropriate when the science spectrum lacks usable OH coverage. |
| Slit edges | `edge_thresh=50`, `trace_thresh=10`, `fit_order=2`, `fit_min_spec_length=0.5`, `left_right_pca=True` | VPH notes propose `edge_thresh=20`, `trace_thresh=20`, Sobel enhancement 3, and several geometry-specific limits | Flamingos: trace 10, edge 200, min length 0.4; TripleSpec: trace 5, min length 0.3 | Retain the active shared values unless edge QA shows missed/false edges. Use the existing VPH-specific overrides rather than reactivating the generic commented block. |
| Arc/tilt combination | `clip=False`, `combine="mean"` | Tilt `clip=False`, `combine="mean"` is commented in the VPH block | No common universal override; NIR arc/sky frames are instrument-specific | The active and commented values agree. Keep if a small number of frames are clean; otherwise test sigma clipping because mean-without-clipping preserves transient defects. |
| Arc/tilt continuum | `subtract_continuum=False` for VPH | Shared comments propose `subtract_continuum=True` and `tilts.rm_continuum=True` | The grism mode actively uses `rm_continuum=True`; other NIR definitions vary | For VPH, keep continuum subtraction off unless broad lamp/sky structure demonstrably compromises centroiding. For grism/echelle, retain the active continuum removal because order tracing and dense OH fitting benefit more. |
| Flat-field sampling | `slit_illum_finecorr=False`, `spec_samp_coarse=3`, `spat_samp=2`, `tweak_slits=False` | Historical VPH block proposes `pixelflat_min_wave=3000`, fine illumination off, `spec_samp_fine=30`, tweak off | LUCI modifies slit tweaking; the other NIR definitions generally leave sampling near defaults | Keep the active coarse sampling only if flat QA is smooth and stable. `spec_samp_fine=30` is much more aggressive than normal and should remain disabled until it improves a demonstrated artifact. |
| Object finding/extraction | Seeing-derived `find_fwhm`, `trace_npoly=3`, one standard, up to five science objects, `boxcar_radius=1.92"`, no 2D-model mask | VPH comments repeat these values (and label them against generic defaults) | TripleSpec: one standard, two science objects, fixed 0.75–2.0" boxcar radii; LUCI: one standard | Keep the active physically motivated seeing conversion. Lower `maxnumber_sci` only for narrow-slit/echelle modes where multiple detections are known to be spurious. |
| Sky rejection | `sky_sigrej=4` | Same value repeated in VPH comments | Flamingos uses 5; generic default is 3 | Start at 4. Raise to 5 for bright/structured residuals only after confirming it does not retain cosmic rays; lower only if genuine sky pixels are being rejected. |
| Spectral flexure | `spec_method="boxcar"`, `spec_maxshift=30` (binning-scaled) | Comments repeat this choice | TripleSpec, Flamingos, and LUCI commonly use `skip`; RIMAS has large allowed shifts | Keep `boxcar` if RIMAS shows measurable exposure-to-exposure shifts. Move to `skip` only after flexure QA shows a stable wavelength solution; it is faster and avoids fitting noise when no correction is needed. |

## Sensitivity function and telluric fitting

| Area | Active RIMAS2 setting | Commented alternative | Near-IR precedent | Selection suggestion |
|---|---|---|---|---|
| Algorithm and polynomial | `algorithm="IR"`, `polyorder=8` | No active UVIS alternative remains | Flamingos, MMIRS, TripleSpec, and NIRES use IR; their orders are commonly 8 (GNIRS uses 6) | Keep IR/8. It is the prevailing NIR configuration and jointly models telluric absorption. Lower the order only if QA shows polynomial ringing. |
| PCA grid and rejection passes | `TellPCA_3000_26000_R10000.fits`, `maxiter=2`, `lower=upper=3` | Earlier historical R25000/five-pass setup is no longer active | Flamingos, MMIRS, LUCI, TripleSpec, MOSFIRE, and NIRES use R10000; MOSFIRE/NIRES explicitly use two passes | Keep. This matches low-resolution RIMAS data and avoids needless high-resolution convolution and repeated differential-evolution fits. |
| Sensfunc shift bounds | VPH and grism use `(-5, 5)` pixels | No active alternative; old RIMAS/X-Shooter-style settings used wider ranges | The TelluricPar default is `(-5, 5)`; high-resolution NIRSPEC/X-Shooter use roughly `(-8, 8)` or `(-10, 10)` | Keep the default-width search unless flexure QA repeatedly reaches a bound. Broaden only for demonstrated calibration shifts. |
| Sensfunc resolution guess | VPH30: 30; VPH300: 300; grism: 4000 | VPH formerly set inactive UVIS resolution values (YJ VPH30 400, HK VPH30 40, VPH300 800) | The IR fit expects `resln_guess`; similar definitions usually rely on sampling rather than setting it | VPH300 = 300 is supported by the HK WaveCalib measurement (~315 at the central wavelength). Recheck VPH30 and grism with measured arc FWHM before treating 30 and 4000 as established values; they are the highest-priority remaining measurements. |
| Sensfunc extrapolation | `extrap_blu=extrap_red=0.4` | None | PypeIt default is 0.1; most NIR definitions do not override it | Reduce toward 0.1 if science/standard wavelength coverage agrees. This does not materially accelerate fitting, but limits less-reliable extrapolated calibration. Retain 0.4 only when the standard is routinely narrower than science coverage. |
| Telluric (post-flux) grism bounds | `pix_shift_bounds=(-10,10)`, `resln_frac_bounds=(0.4,2.0)` | None | Wider than SensFunc/IR defaults, as can be appropriate for echelle orders | Keep separate from sensfunc bounds. Tighten only after fitting representative data shows the solution is well inside narrower ranges. |

## VPH wavelength and trace calibration

| Area | Active RIMAS2 setting | Commented alternative | Near-IR precedent | Selection suggestion |
|---|---|---|---|---|
| VPH parent wavelength fit | `rms_thresh_frac_fwhm=0.4`, `sigdetect=5`, `n_first=n_final=1` | `reidentify`, `cc_thresh=0.6`, GNIRS archive; `fwhm=3`, `nsnippet=1` | Flamingos/LUCI use fractional RMS thresholds 0.07–0.1 and final order 4; TripleSpec uses reidentification with archives | The 0.4 RMS threshold is deliberately permissive and should be retained only while calibration QA requires it. Once stable templates exist, test tightening toward 0.1–0.2. Use the active configuration-specific archives rather than the generic GNIRS archive. |
| VPH YJ calibration | Xe, `holy-grail`, `ldt_nihts.fits`, `fwhm_fromlines=True`, one snippet, S/N threshold 50 | Same ingredients appear in a commented historical block | LUCI uses holy-grail with OH; Flamingos uses arc/OH solutions | Keep the active YJ configuration. Xe and fwhm-from-lines are appropriate when the spectral region lacks dense sky lines. Reduce the S/N threshold only if valid faint standards/science traces are missed. |
| VPH HK calibration | `full_template`, `fwhm=5`, one snippet, `sigdetect=1`, OH archive; VPH300 raises `sigdetect` to 5 | VPH30 comments propose no-local-sky/no-poly | LUCI uses full templates for specific low-resolution settings; TripleSpec uses reidentification archives | Keep the mode-specific archive strategy. The general `sigdetect=1` is permissive; retain VPH300's 5 unless line-identification QA shows too few lines. Enable `no_local_sky`/`no_poly` only for a clear failure of the local sky model. |
| VPH polynomial orders | VPH30 and VPH300 use `n_first=3`, `n_final=5`; VPH30 splits the range (`nsnippet=2`) | Parent comments list `n_first/n_final` alternatives only indirectly | TripleSpec commonly uses final order 3–4 | Keep the higher order/split VPH30 setting for its wide range. For VPH300, reduce only if residuals show overfit structure; otherwise its sparse/curved solution can justify order 5. |
| VPH tilt fitting | Parent uses `tracethresh=10`, `sig_neigh=5`, `nfwhm_neigh=2`; HK uses spatial/spectral orders 4/5 | Historical comments propose 3/4 or the same 4/5 | TripleSpec commonly sets trace threshold 10; echelle orders vary with format | Keep active mode-specific orders. Select lower orders only if tilt QA shows edge oscillations or sparse-line instability. |

## Grism/echelle calibration

| Area | Active RIMAS2 setting | Commented alternative | Near-IR precedent | Selection suggestion |
|---|---|---|---|---|
| Echelle wavelength solution | `holy-grail`, echelle enabled, FWHM 4 then split-arm override 3, final order 4, 5x5 echelle coefficients, 3-sigma rejection | None | TripleSpec/SOAR uses reidentify plus echelle coefficients 4x6 and 3-sigma rejection | Keep the current explicit echelle configuration while the holy-grail solution is needed. Compare against a template/reidentify path only after a stable RIMAS grism archive exists. |
| Resolution from lines | `fwhm_fromlines=False` | None | Many NIR modes use a measured fixed FWHM; YJ VPH enables line measurement | Keep disabled as documented in the source: measured line FWHM currently degrades the grism wavelength solution. Re-enable only with a regression dataset showing better residuals. |
| Grism slit/tilt treatment | `det_min_spec_length=0.10`, nearest synchronization, continuum removal, `tracethresh=25`, tilt orders 3/4 then split-arm 4/5 | None | TripleSpec uses trace threshold 5–10 and min length ~0.3; its format differs substantially | Retain the RIMAS-specific values because the source notes noisy YJ tracing. Use order QA, not numerical matching, to decide whether to lower the trace threshold or polynomial orders. |
| Grism extraction/sky | `bspline_spacing=0.8`, `global_sky_std=False`, `model_full_slit=True`, standard CR masking off | None | ARC/SOAR TripleSpec use `bspline_spacing=0.8` and full-slit modeling; LUCI differs by camera | Keep; this is strongly aligned with cross-dispersed NIR practice. Change `global_sky_std` only if the standard's local background is visibly unstable. |

## Recommended next measurements before changing code

1. Measure arc-line FWHM and dispersion for VPH30 and grism, as was done for
   VPH300, and use those values to validate `IR.resln_guess`.
2. Inspect wavelength-fit RMS in units of the measured FWHM.  If VPH results
   are consistently much better than 0.4 FWHM, tighten
   `rms_thresh_frac_fwhm` incrementally.
3. Compare flexure QA with `boxcar` and `skip` on the same RIMAS sequence;
   retain the correction only if it changes sky-line residuals materially.
4. Examine flat/illumination QA before changing the coarse sampling or
   reactivating the historical `spec_samp_fine=30` setting.
5. Use sensfunc QA to decide whether 40% extrapolation is actually required;
   do not infer that need from fitting runtime.
