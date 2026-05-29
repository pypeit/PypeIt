
.. include:: ../include/links.rst

.. _mmt_binospec:

************
MMT Binospec
************

Overview
========

This file summarizes several instrument specific
items for the MMTO's Binospec spectrograph.

PypeIt supports Binospec in three modes:

- **Multi-Object Spectroscopy (MOS)**: ``mmt_binospec``
- **Longslit**: ``mmt_binospec`` (auto-detected from the ``MASK`` header keyword)
- **IFU**: ``mmt_binospec_ifu`` (fiber-fed integral field unit)

Wavelength Calibration
++++++++++++++++++++++

Templates were created in 2023 that cover the usable range for each of the 3
gratings; the 1000 l/mm template was extended in 2026 to support red central
wavelengths.  All wavelengths are vacuum.  The current template coverage is:

==============  ==========================  ==============
Grating         Template file               Range (Å)
==============  ==========================  ==============
270 l/mm        ``mmt_binospec_270.fits``   3825 – 9216
600 l/mm        ``mmt_binospec_600.fits``   4042 – 10047
1000 l/mm       ``mmt_binospec_1000.fits``  3757 – 7211
==============  ==========================  ==============

For the 1000 l/mm grating, central wavelengths from roughly 4500 Å through
6450 Å are covered with full overlap on both detector sides.  Settings near
the blue end of the grating range (e.g. central wavelength ~4000 Å) extend
slightly below the template's blue limit but should still be usable, since
``full_template`` only requires sufficient overlap to anchor the
cross-correlation; better test data are needed to confirm this corner.

Bad pixel mask
++++++++++++++

The static bad pixel mask is loaded from FITS files derived from the IDL
pipeline calibration data (``badpix_binospec.fits`` plus the hard-coded bad
columns and detector trap regions defined in ``bino_mosaic.pro``).  This
adds roughly 12,500 individual bad pixels per detector, along with bad
columns and detector trap region masking, replacing the small set of
hard-coded bad columns previously identified from 2019 flat and bias
observations.  The mask files are distributed with the package as
``static_calibs/mmt_binospec/bpm_binospec_det{1,2}.fits.gz``.

Cosmic-ray cleanup for arc and tilt frames
++++++++++++++++++++++++++++++++++++++++++

For Binospec arc and tilt calibrations on both detectors, PypeIt applies
an instrument-specific cosmic-ray cleanup before the wavelength-tilt fit.
The cleanup subtracts a row-local median (default 31-pixel window along
the dispersion axis) and runs LA Cosmic on the residual.  This catches
the body of extended cosmic-ray trails, which standard LA Cosmic's 3x3
Laplacian structurally cannot see on raw arc images.  Without this
cleanup, long trails can be mistaken for tilted arc lines and produce
distorted wavelength images.

The cleanup is tuned via the standard ``[calibrations][arcframe][process]``
and ``[calibrations][tiltframe][process]`` blocks.  Defaults for Binospec
are ``cr_median_width=31``, ``sigclip=10``, ``sigfrac=0.3``, ``objlim=0``,
``lamaxiter=2``, ``grow=2.0``, ``rmcompact=False``.  These can be
overridden in a user ``.pypeit`` file; setting ``cr_median_width=0``
disables the cleanup.

IFU Mode
========

The Binospec IFU is a fiber-fed integral field unit with a hexagonal
lenslet array feeding approximately 360 fibers per side into the
spectrograph.  Each side has 40 dedicated sky fibers located at the
outermost ring of each hexagonal sub-bundle, which are used for sky
subtraction.  See the `Binospec IFU instrument page
<https://www.mmto.org/instrument-suite/binospec/binospec-ifu-information/>`__
for more details on the hardware and observing modes.

IFU data are identified automatically from the FITS header keyword
``MASK = 'IFU'`` and reduced using the ``mmt_binospec_ifu`` spectrograph
class with the ``Fiber`` pypeline.

Unlike the ``SlicerIFU`` pypeline used for slicer-based IFUs, the
``Fiber`` pypeline treats each fiber as a distinct object and performs
1D spectral extraction as part of the standard pipeline run.  This
produces both spec2d and spec1d output files.

Key IFU parameters set by default:

- Joint sky fitting using both the dedicated sky fibers and the science
  fibers (``joint_fit_use_sci = True``); see `Tuning the joint sky fit`_
- Grating-dependent B-spline spacing for sky subtraction:
  1.05 Angstrom (270 gpm), 0.5 Angstrom (600 gpm), 0.35 Angstrom
  (1000 gpm)
- Fiber edge detection tuned for densely-packed traces
  (``edge_thresh = 5``)
- Slit edge tweaking using the gradient method
- Spectral flexure correction disabled (Binospec has active flexure
  control)

For a detailed comparison of the PypeIt and IDL IFU pipeline approaches,
see :ref:`mmt_binospec_pipeline_comparison`.

Reducing IFU data
+++++++++++++++++

IFU data are reduced using the standard ``run_pypeit`` workflow.  PypeIt
will automatically detect IFU frames from the ``MASK = 'IFU'`` header
keyword.  Both detectors (``DET01`` for side A and ``DET02`` for side B)
are processed and written to separate extensions within a single spec2d
output file.

.. code-block:: bash

   pypeit_setup -r /path/to/raw -s mmt_binospec_ifu -c all
   run_pypeit mmt_binospec_ifu_A/mmt_binospec_ifu_A.pypeit

The pipeline produces spec1d files containing one extracted spectrum per
fiber.  Each spectrum is identified by its instrument fiber name (e.g.,
``A1``, ``B1``, ``SKY6-1``) via cross-correlation against a reference
profile.  Both boxcar (``BOX``) and optimal Horne (1986) (``OPT``)
extractions are performed.  The spec1d files can be inspected with
``pypeit_show_1dspec``.

.. note::

   The pipeline wavelength calibration is done independently for each
   fiber, which can be time-consuming (~360 fibers per detector).
   A typical reduction with both detectors takes several hours.

Fiber throughput correction and sky subtraction
++++++++++++++++++++++++++++++++++++++++++++++++

Fiber-to-fiber throughput is corrected with a single scalar per fiber
rather than a wavelength-dependent flat division.  The wavelength-dependent
spectral response (lamp × system) is *not* removed here; it is left for
flux calibration downstream.  The reduction proceeds as follows:

1. The ``FiberFlatImages`` calibration is loaded and each fiber's flat
   spectrum is collapsed to a single ``fiber_throughput`` scalar — the
   median over the central wavelength region, divided by the typical
   science-fiber median so science fibers sit near unity and sky fibers
   scale up by their larger on-sky footprint.  The full extracted flat
   spectra (``normflat``) are retained for diagnostics only.

2. When the spectrograph supplies a static ``fiber_illumination.fits``
   vector (Binospec does), it is multiplied into the per-fiber scalar as
   a refinement of the flat-derived throughput.  The side-B (``DET02``)
   vector is mirrored before the fiber-ID lookup, because the IDL
   calibration vector is stored in the extracted side-B fiber-axis order
   whereas PypeIt's DET02 reference profile is in detector order.

3. The combined scalar is attached to each fiber's ``SpecObj`` so that
   the sky-model fit and the post-extraction correction apply the same
   value.

Sky subtraction uses a single 2D B-spline rather than a per-fiber 1D
model.  Every fiber is boxcar pre-extracted from the un-subtracted image,
the aperture sums are divided by each fiber's throughput scalar onto a
common surface, and one B-spline in ``(wavelength, normalized spatial
position)`` is fit — by default jointly over the sky *and* science fibers
(see `Tuning the joint sky fit`_).  Each fiber's predicted sky spectrum is
scaled back to detector counts and distributed across the fiber aperture
using the empirical spatial profile from the flat, producing a genuine
per-pixel 2D ``SKYMODEL``.  Standard extraction then runs on
``sciimg - skymodel``, after which the per-fiber throughput scalar is
divided out of the extracted ``BOX`` and ``OPT`` spectra.

To suppress dipole residuals on bright OH lines, the B-spline fit is run
in two passes: after the first pass each fiber's wavelength solution is
shifted by the small offset that best matches the joint model, and the
fit is repeated.

The IDL pipeline additionally matches the (typically wider) science-fiber
line-spread function by broadening the sky model per fiber.  This LSF
matching step is not currently applied by PypeIt.

Tuning the joint sky fit
++++++++++++++++++++++++

The Binospec IFU sky model is a single 2D B-spline.  By default it is fit
jointly to the dedicated sky fibers *and* all science fibers, relying on
the standard iterative sigma rejection to down-weight pixels that contain
real source flux.  Pulling the science fibers into the fit gives it many
more samples per wavelength and generally produces a smoother, better
constrained sky model.  Two parameters in the ``[reduce][skysub]`` block
control which fibers contribute:

``joint_fit_use_sci`` (:obj:`bool`, default ``True``)
   When ``True``, the global sky B-spline is fit jointly to the dedicated
   sky fibers and the science fibers.  When ``False``, only the dedicated
   sky fibers (those whose ``MASKDEF_OBJNAME`` starts with ``SKY``) are
   used — the legacy behaviour.  Set this to ``False`` when most science
   fibers carry significant source flux (e.g. a target that fills much of
   the field), so that source light cannot leak into the sky model.

``sci_exclude_radius`` (:obj:`float`, default ``0.0``, in arcsec)
   When greater than zero *and* ``joint_fit_use_sci = True``, any science
   fiber whose on-sky position lies within this radius of the IFU
   geometric center is dropped from the fit.  This is the middle ground
   between the two ``joint_fit_use_sci`` extremes: keep the bulk of the
   science fibers in the joint fit, but exclude the central fibers that a
   bright, extended target would dominate.  The default of ``0.0`` keeps
   all science fibers.

To disable the joint fit and use only the dedicated sky fibers, add the
following to your ``.pypeit`` file:

.. code-block:: ini

    [reduce]
        [[skysub]]
            joint_fit_use_sci = False

To keep the joint fit but exclude science fibers within 3 arcsec of the
field center (e.g. for a bright extended target):

.. code-block:: ini

    [reduce]
        [[skysub]]
            joint_fit_use_sci = True
            sci_exclude_radius = 3.0

Producing datacubes
+++++++++++++++++++

Because the Binospec IFU is fiber-fed rather than slicer-based, the
general-purpose ``pypeit_coadd_datacube`` script (designed for slicer
IFUs like KCWI) does not produce correct results.  Instead, use the
dedicated ``pypeit_binospec_ifu_cube`` script that reads
extracted 1D fiber spectra from spec1d files and
builds a datacube.

By default, optimal (``OPT``) extraction columns are used; pass
``--boxcar`` to use boxcar (``BOX``) extraction columns instead.  Sky
subtraction is already applied by the pipeline, so no additional sky
subtraction is performed.

The script then:

1. Identifies sky vs. science fibers from fiber metadata
2. Resamples all fiber spectra onto a common linear wavelength grid
3. Combines both detectors (up to 640 science fibers total)
4. Maps each fiber to its on-sky position using the IFU layout
   calibration file (``bino_IFU_sky_layout.fits``)
5. Interpolates the irregularly-spaced fiber positions onto a regular
   spatial grid using ``scipy.interpolate.griddata``

.. note::

   The two Binospec detectors produce mirror-image spectra.  The script
   automatically accounts for this by reversing the fiber-to-sky mapping
   for side B so that the two halves of the IFU tile correctly.

.. note::

   Fiber-to-fiber throughput correction is handled entirely by the
   pipeline (the per-fiber throughput scalar, optionally refined by the
   static fiber illumination vector, divided out of the extracted
   spectra).  No additional illumination correction is needed before
   building datacubes.  Spectral response (flux calibration) from
   standard star observations is not yet implemented.

Basic usage
^^^^^^^^^^^

To build datacubes from spec1d files:

.. code-block:: bash

   pypeit_binospec_ifu_cube spec1d_*.fits

Each input file produces a separate output datacube named
``cube_sci_img_*.fits``.  Each cube combines both detectors and
contains ``FLUX`` and ``VAR`` extensions.

Command-line options
^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   pypeit_binospec_ifu_cube spec1d_*.fits [options]

``--output FILENAME``
   Output FITS filename.  Only valid when processing a single input
   file.  Default is auto-generated from the input filename
   (e.g., ``cube_sci_img_*.fits``).

``--spatial_scale SCALE``
   Output spatial pixel scale in arcsec.  Default is 0.27, which
   matches the IDL pipeline (``scl = 0.269461``).

``--boxcar``
   Use boxcar (``BOX``) extraction columns from the spec1d file
   instead of the default optimal (``OPT``) columns.

``--method METHOD``
   Spatial interpolation method: ``nearest``, ``linear`` (default), or
   ``cubic``.  The ``linear`` method is recommended for most use cases.

Output format
^^^^^^^^^^^^^

The output FITS file contains:

- **Extension 0** (``PRIMARY``): Empty primary HDU
- **Extension 1** (``FLUX``): 3D datacube with axes
  (NAXIS1=RA, NAXIS2=DEC, NAXIS3=wavelength) and a full 3-axis WCS
  (``RA---TAN``, ``DEC--TAN``, ``WAVE``)
- **Extension 2** (``VAR``): Variance datacube with the same shape and
  WCS

The cube can be viewed directly with tools like ds9 or QFitsView.
Typical output dimensions for the default spatial scale are approximately
63 x 47 spatial pixels, with the number of wavelength pixels depending
on the grating and wavelength coverage.

.. note::

   The datacube construction can be sped up significantly by using an
   accelerated BLAS library.  On macOS (13.3+), switching to Apple's
   ``newaccelerate`` backend yielded a ~3x speedup.  See the
   `conda-forge BLAS documentation
   <https://conda-forge.org/docs/maintainer/knowledge_base/#switching-blas-implementation>`__
   for instructions on selecting a BLAS implementation with conda, e.g.:

   .. code-block:: bash

      conda install "libblas=*=*newaccelerate"

Interactive 1D fiber extraction
++++++++++++++++++++++++++++++++

The ``pypeit_binospec_ifu_extract`` script opens an interactive matplotlib
GUI showing the IFU as a hexagonal grid of fibers, colored by per-fiber
integrated flux.  Users select fibers via rectangle, circle, or
individual-click modes and write the combined extracted 1D spectrum as a
:class:`~pypeit.onespec.OneSpec` FITS file readable by
``pypeit_show_1dspec``, ``pypeit_coadd_1dspec``, and other 1D-spectrum
tools.

Basic usage
^^^^^^^^^^^

.. code-block:: bash

   pypeit_binospec_ifu_extract spec1d_obj.fits

By default the script uses optimal (``OPT_*``) extraction columns; pass
``--boxcar`` to use boxcar (``BOX_*``) instead.  The output filename
defaults to the input basename with ``spec1d_`` replaced by
``extract1d_`` (e.g. ``spec1d_obj.fits`` → ``extract1d_obj.fits``).

Command-line options
^^^^^^^^^^^^^^^^^^^^

``-o, --output FILENAME``
   Output OneSpec FITS filename.  Default is auto-generated from the
   input filename.

``--boxcar``
   Use boxcar (``BOX_*``) extraction columns from the spec1d file
   instead of the default optimal (``OPT_*``) columns.

GUI controls
^^^^^^^^^^^^

- **Rectangle / Circle**: drag in the fiber view to select a spatial
  region.  Selected fibers are highlighted in cyan as you drag.
- **Click Fibers**: click individual hexagons to toggle them in or out
  of the selection.  Combine with rectangle/circle selections — both
  contribute to the final extraction.
- **Wavelength slider**: restrict the per-fiber integrated flux used to
  colour the hex grid.
- **Extract Spectrum**: resample all selected fibers onto a common
  wavelength grid and sum, propagating inverse variance.
- **Reset Display**: clear the selection and restore the full-range
  wavelength view.
- **Save Spectrum**: write the extracted spectrum to the output OneSpec
  FITS file.

Sky-line masking
^^^^^^^^^^^^^^^^

Per-fiber integrated flux used to colour the hex grid masks out the
brightest optical sky lines (5577 Å, Na I D 5890/5896 Å, [OI]
6300/6363 Å, plus the OH band centers around 7340–8867 Å) within
±10 Å.  This prevents residual sky-subtraction artifacts from biasing
the visualization.  The masking applies *only* to the hex display:
extracted spectra are summed over the full wavelength range.
