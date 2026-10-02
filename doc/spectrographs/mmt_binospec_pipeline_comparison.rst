
.. include:: ../include/links.rst

.. _mmt_binospec_pipeline_comparison:

**********************************************
Binospec Pipeline Comparison: IDL vs PypeIt
**********************************************

Comparison of the Binospec IDL reduction pipeline and PypeIt for all
Binospec spectroscopic modes: long-slit, multi-object (MOS), and
integral field unit (IFU).  Both pipelines share many algorithmic
approaches — cross-correlation for fiber identification, B-spline sky
fitting, grating-dependent calibration parameters — while differing in
error propagation, extraction methods, and scattered light treatment.

.. contents:: Table of Contents
   :depth: 2
   :local:


Supported Modes
===============

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Mode
     - PypeIt
     - IDL Pipeline
   * - Long-slit
     - ``MultiSlit`` pypeline (``mmt_binospec``)
     - Standard slit processing
   * - Multi-object (MOS)
     - ``MultiSlit`` pypeline (``mmt_binospec``)
     - Standard slit processing (``*_ms.pro``)
   * - IFU
     - ``Fiber`` pypeline (``mmt_binospec_ifu``)
     - Fiber-specific processing (``bino_ifu_*.pro``)

Both pipelines use the same detector readout, nonlinearity correction,
and basic image processing for all modes.  The processing paths diverge
at tracing (slits vs fibers), sky subtraction, and extraction.


Image Processing
================

Shared Steps
++++++++++++

Both pipelines perform the same fundamental image processing:

- **Bias subtraction** from overscan regions (8 amplifiers per detector)
- **Gain correction** to convert ADU to electrons (per-amplifier gains)
- **Nonlinearity correction** using degree-4 polynomial coefficients
  measured from detector characterization data
- **Flat fielding** from internal lamp flats
- **Cosmic ray rejection** (PypeIt uses L.A.Cosmic; IDL uses custom
  sigma-clipping)
- **Bad pixel masking** from static masks and runtime detection

Scattered Light
+++++++++++++++

Scattered light from inter-fiber and inter-slit gaps is a significant
contaminant for Binospec, particularly in IFU mode where fibers are
closely packed.

**IDL pipeline** (``bino_ifu_model_sc_light.pro``):

- Samples scattered light in inter-fiber gaps on a coarse grid
  (128 x 128 bins)
- Uses ``resistant_mean()`` at 3-sigma in each gap region
- Interpolates between gap measurements with spline fitting,
  working in quadrants to avoid edge effects
- Applied to IFU data only

**PypeIt** (``pypeit/images/rawimage.py``):

- Three methods controlled by the ``scattlight`` parameter:

  - ``method='model'``: Uses a pre-fit parametric scattered light model
    from a dedicated calibration frame
  - ``method='frame'``: Fits a scattered light model to each science
    frame independently using inter-slit/inter-fiber gap pixels.
    Falls back to model or archival parameters if the fit fails
  - ``method='archive'``: Uses archival model parameters for the
    grating/binning configuration

- Optional fine correction step for residual structure
- Binospec IFU uses ``method='frame'`` for science data (no separate
  scattered light calibration frame required)

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Aspect
     - PypeIt
     - IDL Pipeline
   * - Approach
     - Parametric model fit to gap pixels
     - Grid sampling + spline interpolation in gaps
   * - IFU handling
     - Per-frame model fit (``method='frame'``)
     - Per-frame gap sampling
   * - MOS handling
     - Same framework (``method='model'`` or ``method='frame'``)
     - Not applied (inter-slit gaps are wider)
   * - Fallback
     - Archival or calibration-frame model parameters
     - None


Error Propagation
=================

PypeIt: Full Variance Model
++++++++++++++++++++++++++++

PypeIt tracks errors through every processing step using inverse variance
(``ivar``) arrays.  The variance model
(:func:`~pypeit.core.procimg.variance_model`) is:

.. math::

   V = s^2 \left[\max(0,C) + \frac{D \cdot t}{3600} + V_\mathrm{rn}
   + V_\mathrm{proc}\right] + \varepsilon^2 \max(0,c)^2

where *s* is the inverse sensitivity (1/gain), *C* is counts, *D* is dark
current, *V_rn* is read noise variance, *V_proc* is processing variance,
and *eps* is the noise floor (fractional error).

Key properties:

- **Base variance** (detector-level: read noise + dark + processing) is
  tracked separately from signal-dependent Poisson noise
- ``ivar`` arrays propagated through every step: image processing, flat
  fielding, sky subtraction, and extraction
- Optimal extraction weights computed as ``ivar * profile**2``
- Two output ivar columns: ``OPT_COUNTS_IVAR`` (full variance including
  object shot noise) and ``OPT_COUNTS_NIVAR`` (sky + readnoise only,
  for S/N estimation)
- Model variance iteratively updated during local sky subtraction when
  ``model_noise=True`` (MOS/long-slit mode)
- Relevant code: ``pypeit/core/procimg.py``, ``pypeit/core/skysub.py``,
  ``pypeit/core/extract.py``

IDL Pipeline: Incomplete Error Tracking
++++++++++++++++++++++++++++++++++++++++

The IDL pipeline does not propagate errors through most processing steps:

- **No error tracking** through fiber tracing, flat fielding, wavelength
  calibration, or spectral linearization
- **Sky subtraction** uses inverse-variance weighting internally:

  .. code-block:: none

     isky = gain^2 / (gain * |sky| + rdnoise^2)

  but the resulting sky model errors are not propagated to downstream
  products

- **Wavelength calibration** stores ``wl_s_err`` (RMS per fiber) but
  never propagates it
- **IFU fiber extraction** via bounded least-squares
  (``bounded_least_squares.pro``) solves for fiber fluxes with
  positivity constraints but does not compute covariance
- **Only at cube building** (``bino_ifu_cube.pro``) are errors finally
  estimated:

  .. code-block:: none

     err = sqrt(flux * flat + rdnoise^2) / flat

  This is a pure Poisson + read noise assumption that ignores
  accumulated errors from all prior processing steps

- No covariance tracking between adjacent fibers (relevant since fiber
  profiles overlap)

Comparison Table
++++++++++++++++

.. list-table::
   :header-rows: 1
   :widths: 25 40 35

   * - Processing Step
     - PypeIt
     - IDL Pipeline
   * - Raw to processed image
     - Full variance model (Poisson+RN+dark+proc)
     - None
   * - Flat fielding
     - ``ivar`` updated for flat noise
     - No error map
   * - Sky subtraction
     - B-spline weighted by ``ivar``; IFU uses extracted fiber-flat and
       static illumination corrections before fitting
     - Weighted internally; no error output
   * - Extraction
     - Optimal with model variance; outputs ``ivar``
     - MOS: optimal with ``ivar``; IFU: bounded LS, no error bars
   * - Wavelength cal
     - Propagated through resampling
     - Stores ``wl_s_err`` (unused)
   * - Cube building
     - Full ``ivar`` propagation via ``coadd3d``
     - Poisson+RN assumption only

Implication
+++++++++++

PypeIt's variance model provides proper error propagation across all
modes.  This is particularly important for:

- Faint emission-line science where reliable S/N estimates are critical
- Combining exposures with different conditions (proper inverse-variance
  weighting)
- Flagging unreliable spaxels from dead/weak fibers (IFU mode)


Slit and Fiber Tracing
======================

MOS / Long-slit
++++++++++++++++

Both pipelines trace slit edges from flat field exposures.  PypeIt uses
Sobel-filter edge detection with polynomial fits; the IDL pipeline uses
peak detection with iterative profile fitting.

IFU Fiber Tracing
+++++++++++++++++

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Aspect
     - PypeIt
     - IDL Pipeline
   * - Trace model
     - Block-level edge detection (21 blocks per side, each containing 8-20 fibers)
     - Gaussian-Hermite profile per fiber (h3-h6)
   * - Detection method
     - Sobel filter at high threshold (100σ) detects block boundaries at ~70px inter-block gaps
     - Peak detection + iterative profile fitting
   * - Profile model
     - Empirical from flat field for extraction
     - 8-parameter Gaussian-Hermite or Moffat
   * - Dead fiber handling
     - Excluded during matching (distance set to infinity)
     - Interpolation from reference catalog

For PypeIt Binospec reductions, arc and tilt images on both detectors
pass through an instrument-specific cleanup that subtracts a row-local
median and runs LA Cosmic on the residual.  This removes extended
cosmic-ray trails (including the body of long trails that the standard
3x3 Laplacian misses) before tilt fitting; such features can otherwise
drive spurious curvature in the wavelength image.  The cleanup is tuned
via the standard ``[calibrations][arcframe/tiltframe][process]`` blocks
(see ``cr_median_width``, ``sigclip``, ``lamaxiter``, ``grow``, etc.).


Fiber Identification (IFU)
==========================

Both pipelines use cross-correlation against a reference fiber profile
to assign physical fiber IDs to detected traces.  The reference profile
(``fiber_ref_profile.fits``) contains expected pixel positions and
Gaussian-Hermite profile parameters for each fiber, obtained from a
high-quality flat field observation.  The algorithm is adapted from the
IDL pipeline's ``bino_ifu_fiber_id.pro``.

Shared Algorithm
++++++++++++++++

1. Build a synthetic spatial profile from detected fiber positions
   (sum of Gaussians at each detected position)
2. Cross-correlate the synthetic profile against the reference profile
   in 5 overlapping segments with cosine apodization, each covering
   a portion of the detector while skipping 150 pixels at each edge
3. Fit a linear polynomial to the segment offsets to capture
   position-dependent shifts from flexure, scale, or distortion
4. Match individual fibers using a distance threshold after applying
   the per-position polynomial shift

PypeIt Enhancements
+++++++++++++++++++

PypeIt extends the IDL algorithm with:

- **Two-pass iterative matching**: A tight threshold (1.7 px) is
  applied first.  The polynomial is then refit using only matched pairs,
  and a second pass with a relaxed threshold (3.5 px) catches
  physically displaced fibers (e.g., those adjacent to dead fibers
  on side B)
- **Sub-pixel accuracy**: Float-valued slit center positions from the
  midpoint of traced edges are passed to the matching algorithm,
  avoiding the ~0.5 px rounding error from integer ``spat_id`` values
- **Sky fiber identification by name**: Sky fibers are identified by
  their ``FIB_NAME`` prefix (``SKY*``) from the reference profile,
  rather than hardcoded index lists.  The ``FIB_TYPE`` field in the
  reference file is unreliable (all fibers are marked as ``SKY``)

Matching Results
++++++++++++++++

- Side A (DET01): 360/360 fibers matched (100%), including 40 sky fibers
- Side B (DET02): 352/356 fibers matched (98.9%), including 40 sky fibers.
  4 unmatched fibers (B26, B78, B160, B200) are each adjacent to a dead
  fiber and physically displaced by 5-7 pixels — consistent with the IDL
  pipeline, which also fails to match these fibers


Sky Subtraction
===============

MOS / Long-slit
++++++++++++++++

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Aspect
     - PypeIt
     - IDL Pipeline
   * - Method
     - B-spline (1D spectral + polynomial spatial)
     - Kelson-style B-spline (2D/3D in wavelength + focal plane coords)
   * - Spatial model
     - Polynomial basis within slit
     - 3D spline over focal plane for TARGET slits; 2D for BOX slits
   * - Variance weighting
     - Full ``ivar`` from variance model
     - Approximate: ``gain^2 / (gain * |sky| + rdnoise^2)``
   * - Knot spacing
     - Configurable per grating
     - 1.05/0.50/0.35 Angstrom for 270/600/1000 gpm
   * - Object handling
     - Joint sky + object B-spline fit (local sky subtraction)
     - Separate sky and object steps
   * - Extended slits
     - Full slit modeled jointly
     - Two-stage fit (blue/red halves separately)

Both pipelines use B-spline fitting for the sky model; the IDL pipeline
additionally models spatial variation across the focal plane for MOS
observations using 3D splines.

IFU
++++

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Aspect
     - PypeIt
     - IDL Pipeline
   * - Strategy
     - All fibers extracted and throughput-corrected; a single 2D
       B-spline (joint over sky and science fibers by default) builds a
       per-pixel sky model subtracted before extraction
     - Dedicated sky fibers with ``resistant_mean`` or B-spline
   * - Sky fibers
     - 40 per side, identified by fiber name (``SKY*``)
     - 40 per side, at outermost ring of each sub-bundle
   * - Knot spacing
     - Grating-dependent (1.05/0.50/0.35 Angstrom)
     - Grating-dependent (same values)
   * - Sky line correction
     - Per-fiber wavelength refinement (2-pass); LSF/PSF kernel not applied
     - PSF difference kernel applied (``sigdiff > 0.2`` px)
   * - Flexure
     - Spectral flexure disabled (active flexure control)
     - Cross-correlation offset in wavelength zero-point

Both pipelines use the same 40 dedicated sky fibers per side (8 fibers
at each of 5 radial positions on the outermost ring of each sub-bundle)
and the same grating-dependent B-spline knot spacing.

In PypeIt, the IFU sky subtraction is handled by
:class:`~pypeit.find_objects.FiberFindObjects`, which inherits the joint
sky fitting framework from :class:`~pypeit.find_objects.SlicerIFUFindObjects`.
The sky subtraction flow is:

1. One spectrum is created for each identified fiber, and each fiber is
   boxcar pre-extracted from the un-subtracted image to provide the
   aperture sums and wavelengths that feed the sky fit.
2. Each fiber's extracted flat spectrum is collapsed to a single
   ``fiber_throughput`` scalar (not a wavelength-dependent division), and
   the aperture sums are divided by it onto a common surface.
3. For Binospec IFU, the static ``fiber_illumination.fits`` vector is
   multiplied into the per-fiber scalar.  The DET02 vector is mirrored
   before the fiber-ID lookup to account for the side-B
   detector/extraction-axis flip.
4. A single 2D B-spline (wavelength plus a low-order spatial coordinate)
   is fit to the throughput-normalized fibers — by default jointly over
   both the dedicated sky fibers and the science fibers.  The fit is run
   in two passes with a per-fiber wavelength refinement in between.
5. Each fiber's predicted sky is scaled back to detector counts and
   distributed across its aperture using the empirical flat profile,
   forming a per-pixel 2D ``SKYMODEL`` that is subtracted from ``SCIIMG``
   before extraction.

The ``SKYMODEL`` frame is therefore the actual model used for
subtraction, in the same detector-pixel units as ``SCIIMG`` — not a
post-hoc reconstruction.  The wavelength-dependent spectral response is
left for flux calibration rather than removed by the flat.


Extraction
==========

MOS / Long-slit
++++++++++++++++

Both pipelines use optimal extraction for slit spectroscopy:

- **PypeIt**: Horne (1986) optimal extraction with joint local sky +
  object fitting, iterative model variance, and full ``ivar`` output
- **IDL pipeline**: Standard optimal extraction with ``ivar`` weighting

IFU
++++

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Aspect
     - PypeIt
     - IDL Pipeline
   * - Method
     - Per-fiber boxcar + optimal (Horne 1986) extraction with
       flat-derived profiles
     - Bounded least-squares (positivity constraint)
   * - Sky subtraction
     - Global 1D fiber sky model used directly (no local sky)
     - Sky subtracted during linearization
   * - Cross-talk
     - Mitigated by flat-derived profiles, but not solved
       simultaneously across neighboring fibers
     - Simultaneous multi-fiber solve per block
   * - Profile
     - Empirical from flat field, normalized per spectral row
     - Gaussian-Hermite from flat field
   * - Error output
     - Full ``ivar`` per spectral pixel
     - None

PypeIt's IFU extraction (:class:`~pypeit.extraction.FiberExtract`)
performs both boxcar and Horne (1986) optimal extraction for each fiber.
The global fiber sky model is used directly — there is no local sky
subtraction step.  Spatial profiles are built empirically from the flat
field within each fiber aperture and normalized to unit sum per spectral
row.  If fewer than half the fibers have valid empirical profiles, the
code falls back to simpler extraction behavior.

The IDL pipeline solves for fiber fluxes simultaneously in blocks using
bounded least-squares with a positivity constraint.  This naturally
handles cross-talk between adjacent fibers whose Gaussian-Hermite
profiles overlap, at the cost of not providing per-pixel error estimates.


Wavelength Calibration
======================

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Aspect
     - PypeIt
     - IDL Pipeline
   * - Method
     - Full template matching + reidentification
     - Per-block Legendre polynomial (degree 3)
   * - Spatial model
     - 2D fit (spectral + spatial)
     - Independent per-fiber solution
   * - Refinement
     - Flexure correction from sky lines (MOS); disabled for IFU
     - Sky line cross-correlation adjustment
   * - Arc lamps
     - HeI, NeI, ArI, ArII
     - Same (HeNe + Ar)

For IFU mode, PypeIt now runs wavelength calibration per block-slit (42
total: 21 blocks per detector × 2 detectors) rather than per individual
fiber (720 previously).  Each block-slit contains 8–20 fibers that share
a common wavelength solution, providing a significant runtime improvement
while retaining the accuracy of the 2D spatial fit within each block.


Cube Building (IFU)
===================

Both pipelines produce 3D datacubes from extracted IFU spectra by
mapping fiber positions to sky coordinates and interpolating onto a
regular spatial grid.

- **PypeIt** (``pypeit/scripts/binospec_ifu_cube.py``): Combines both
  detectors (640 science fibers total, 320 per side) using fiber sky
  positions from the IFU layout file.  Reads already-extracted fiber
  spectra from spec1d files.  Full ``ivar`` propagation through the
  cube-building process
- **IDL pipeline** (``bino_ifu_cube.pro``): Similar fiber-to-sky
  mapping.  Errors estimated at cube-building time using a Poisson +
  read noise assumption


.. rubric:: Summary

The IDL and PypeIt pipelines share many core algorithmic approaches for
Binospec data reduction — B-spline sky fitting, cross-correlation fiber
identification, grating-dependent calibration parameters, and dedicated
sky fiber subtraction for IFU mode.  They differ primarily in:

- **Error propagation**: PypeIt provides full variance tracking through
  every processing step; the IDL pipeline estimates errors only at the
  final cube-building stage
- **Extraction**: PypeIt uses Horne (1986) optimal extraction with
  empirical profiles for IFU mode; the IDL pipeline uses simultaneous
  bounded least-squares which more explicitly models inter-fiber
  cross-talk
- **Scattered light**: PypeIt fits a parametric model per-frame;
  the IDL pipeline samples inter-fiber gaps on a coarse grid
- **Fiber matching**: Both use the same cross-correlation approach;
  PypeIt adds two-pass iterative matching for displaced fibers and
  sub-pixel position accuracy
- **Throughput mapping**: PypeIt applies the extracted fiber flat plus
  the static IDL ``fiber_illumination.fits`` vector, including the
  side-B mirror mapping needed for DET02
- **Sky line sharpness**: The IDL pipeline applies a PSF difference
  kernel correction not yet implemented in PypeIt
