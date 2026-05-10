
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

- Joint sky fitting across all fibers using dedicated sky fibers
  (``joint_fit = True``)
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
``SCI1-1``, ``SKY6-1``) via cross-correlation against a reference
profile.  Both boxcar (``BOX``) and optimal Horne (1986) (``OPT``)
extractions are performed.  The spec1d files can be inspected with
``pypeit_show_1dspec``.

.. note::

   The pipeline wavelength calibration is done independently for each
   fiber, which can be time-consuming (~360 fibers per detector).
   A typical reduction with both detectors takes several hours.

Fiber illumination correction
+++++++++++++++++++++++++++++

The pipeline applies two levels of fiber-to-fiber throughput correction
during reduction:

1. **Dome-flat illumination correction**: The relative throughput from
   flat field observations is baked into the pixel flat during
   flat-fielding.

2. **Sky-line illumination correction**: During joint sky subtraction,
   the pipeline measures bright sky emission lines across all fibers
   to derive a per-fiber throughput correction.  This accounts for the
   different optical paths of sky fibers (bare fibers) versus science
   fibers (lenslet-fed), which cannot be fully captured by dome flat
   illumination alone.

Both corrections are applied automatically by the pipeline and require
no additional post-processing steps.

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
   pipeline (dome-flat illumination + sky-line correction during sky
   subtraction).  No additional illumination correction is needed
   before building datacubes.  Spectral response (flux calibration)
   from standard star observations is not yet implemented.

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

