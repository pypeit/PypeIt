.. include:: ../include/links.rst

.. _jwst_nirspec:

************
JWST NIRSpec
************

Overview
========

This file summarizes the instrument-specific settings and the reduction
workflow for the JWST/NIRSpec spectrograph (PypeIt name ``jwst_nirspec``).

PypeIt can reduce NIRSpec fixed-slit (FS, ``EXP_TYPE = NRS_FIXEDSLIT``) and
multi-object (MSA, ``EXP_TYPE = NRS_MSASPEC``) spectroscopy (see
:ref:`jwst_nirspec_tested` for the modes that have been tested so far).
Instead of processing the raw ramps directly,
PypeIt starts from the products of the first two stages of the `jwst`_ pipeline
(``calwebb_detector1`` and ``calwebb_spec2``):

    - the ramp fitting, the detector-level corrections, and the WCS assignment
      are performed by the JWST pipeline;
    - the wavelength calibration, the slit positions, and the flat-field,
      path-loss, bar-shadow, and photometric corrections are taken from the
      JWST calibrated (``_cal.fits``) and interpolated flat-field
      (``_interpolatedflat.fits``) files;
    - PypeIt then performs the sky subtraction (or the nod subtraction), the
      object finding, the optimal extraction, and the 1D and 2D coadding of
      the exposures.

The data are reduced **one slit at a time**: each slit (or MSA slitlet) is
reduced independently, using the cutout of that slit produced by the JWST
pipeline to define the region of the full-frame detector image(s) to reduce.

The workflow is:

    #. :ref:`jwst_nirspec_preprocess`: run the JWST pipeline with
       ``pypeit_preprocess_jwst``.
    #. :ref:`jwst_nirspec_setup`: create the :ref:`pypeit_file` with
       :ref:`pypeit_setup` and select the slits to reduce.
    #. :ref:`jwst_nirspec_run`: reduce the data with :ref:`run-pypeit`.
    #. :ref:`jwst_nirspec_coadd`: coadd the exposures with
       :ref:`pypeit_coadd_2dspec <coadd2d>` and :ref:`pypeit_coadd_1dspec
       <coadd1d>`.

.. _jwst_nirspec_tested:

Tested modes
------------

The following modes are included in the PypeIt development suite and have
been tested:

.. list-table::
   :header-rows: 1
   :widths: 10 18 12 12 12 36

   * - Mode
     - Aperture (``decker``)
     - Disperser
     - Filter
     - Program
     - Notes
   * - FS
     - ``NRS_S200A2_SLIT``
     - G395M
     - F290LP
     - 1967
     - 3-point nod; the S200A2 slit falls on NRS1 only (``DET01``).
   * - MSA
     - ``NRS_FULL_MSA``
     - G395M
     - F290LP
     - 3543
     - 3-shutter nod; the selected slitlet falls on both NRS1 and NRS2
       and is reduced as a mosaic (``MSC01``).
   * - MSA
     - ``NRS_FULL_MSA``
     - PRISM
     - CLEAR
     - 1287
     - 3-shutter nod (two sets of three exposures); the selected slitlets
       fall on NRS1 only (``DET01``).

Other gratings and filters are expected to work in the same way, but they have
not been tested yet.  Please report any problem on the PypeIt `GitHub
<https://github.com/pypeit/PypeIt/issues>`__ issue tracker.

Installation
============

Two additional packages are needed to reduce NIRSpec data:

    - The `jwst`_ pipeline package is required to run
      ``pypeit_preprocess_jwst``.  It is **not** used by any of the other
      PypeIt scripts.

    - The `stdatamodels`_ package is **strongly recommended** for the PypeIt
      reduction.  PypeIt uses it to read the GWCS of each slit stored in the
      JWST ``_cal.fits`` files, from which the along-slit platescale of the slit
      is computed (see :ref:`jwst_nirspec_platescale`).  PypeIt does not import
      it directly, but it needs the ``asdf`` converters that ``stdatamodels``
      provides for the JWST-specific transforms.  If ``stdatamodels`` is not
      installed, PypeIt issues a warning and uses the nominal NIRSpec
      platescale for all slits.

Both are installed with:

.. code-block:: console

    pip install "pypeit[jwst]"

If you have already preprocessed your data (or you preprocess them in a
separate environment), you only need ``stdatamodels``:

.. code-block:: console

    pip install stdatamodels

See :ref:`optional-dependencies`.

.. _jwst_nirspec_preprocess:

Preprocessing with the JWST pipeline
====================================

The raw, uncalibrated JWST files (``*_uncal.fits``) can be downloaded from the
`MAST archive <https://mast.stsci.edu/>`__.  The ``pypeit_preprocess_jwst``
script runs the JWST ``Detector1Pipeline`` and ``Spec2Pipeline`` on them, with
the parameters needed to produce the files used by PypeIt.  The JWST pipeline
requires the reference files from the `JWST Calibration Reference Data System
(CRDS) <https://jwst-crds.stsci.edu/>`__; the local CRDS cache is set either
with the ``CRDS_PATH`` environment variable or with the ``--crds_path``
option.  For example:

.. code-block:: console

    pypeit_preprocess_jwst raw_data -o jwst_output --crds_path ~/crds_cache

processes all the ``*_uncal.fits`` files in the ``raw_data`` directory and
writes the results to the ``jwst_output`` directory.  Use ``--list`` to see the
files that would be processed, and the ``--pid``, ``--obs``, ``--vis``, ``--vg``,
``--s``, ``--a``, and ``--exp`` options to select a subset of the files based on
the JWST file naming convention.

Compared to the default JWST pipeline, ``pypeit_preprocess_jwst``:

    - applies the ``clean_flicker_noise`` step in stage 1;
    - skips the background subtraction, the master background, the 2D
      resampling, and the 1D extraction in stage 2 (these are performed by
      PypeIt);
    - always applies the bar-shadow correction, so that the illumination
      along the MSA slitlets is flat;
    - saves the ``_assign_wcs.fits``, ``_extract_2d.fits``, and
      ``_interpolatedflat.fits`` intermediate products.

The source type used by the ``srctype`` step (which selects, e.g., the
path-loss correction) is ``EXTENDED`` by default, and can be changed with
``--source_type POINT``.  All other parameters of the JWST pipelines can be
changed by providing JSON files with ``--det1_config`` and ``--spec2_config``;
use ``--save_config`` to write the default parameters to JSON files that can be
edited.  The parameters used in each run are always saved in the output
directory.

For each exposure and detector, PypeIt uses three of the files in the output
directory:

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - JWST file
     - PypeIt ``frametype``
     - Used for
   * - ``*_assign_wcs.fits``
     - ``science``
     - Full-frame count-rate image to be reduced.
   * - ``*_cal.fits``
     - ``arc,tilt``
     - Slit cutouts, wavelength maps, path-loss and bar-shadow corrections,
       photometric conversion, and slit GWCS.
   * - ``*_interpolatedflat.fits``
     - ``pixelflat,trace``
     - Flat field.

The full script usage is:

.. include:: ../help/pypeit_preprocess_jwst.rst

.. _jwst_nirspec_setup:

PypeIt File
===========

Run :ref:`pypeit_setup` on the directory with the output of the JWST pipeline:

.. code-block:: console

    pypeit_setup -r jwst_output -s jwst_nirspec -c all

The frame types are assigned based on the suffix of the file names (see the
table above).  The other files in the directory (e.g., ``*_rate.fits`` or
``*_extract_2d.fits``) are not classified and are commented out in the
:ref:`pypeit_file`.  The instrument configurations are defined by the grating
(``dispname``), the filter (``filter1``), the aperture (``decker``), and the
target.

The NRS1 and NRS2 files of the same exposure are grouped together: they are
assigned to the same calibration group (``calib``), and the science files get
the same ``comb_id``.  Each exposure is calibrated using its own ``_cal.fits``
and ``_interpolatedflat.fits`` files.  Here is an example of the
:ref:`data_block` for an FS observation with a 3-point nod:

.. code-block:: console

    data read
     path /path/to/jwst_output
                                                filename |       frametype | ... |     dithpat | dithpos | dithoff | calib | comb_id | bkg_id
                 jw01967012001_03102_00001_nrs1_cal.fits |        arc,tilt | ... |        None |    None |     0.0 |     1 |      -1 |     -1
                 jw01967012001_03102_00001_nrs2_cal.fits |        arc,tilt | ... |        None |    None |     0.0 |     1 |      -1 |     -1
    jw01967012001_03102_00001_nrs1_interpolatedflat.fits | pixelflat,trace | ... |        None |    None |     0.0 |     1 |      -1 |     -1
    jw01967012001_03102_00001_nrs2_interpolatedflat.fits | pixelflat,trace | ... |        None |    None |     0.0 |     1 |      -1 |     -1
          jw01967012001_03102_00001_nrs1_assign_wcs.fits |         science | ... | 3-POINT-NOD |       1 |  -0.919 |     1 |       1 |     -1
          jw01967012001_03102_00001_nrs2_assign_wcs.fits |         science | ... | 3-POINT-NOD |       1 |  -0.919 |     1 |       1 |     -1
          ...
    data end

The ``dithoff`` column gives the offset of the exposure along the slit (in
arcsec).  It is the negative of the ``YOFFSET`` header keyword (the offset in
the NIRSpec ideal coordinate frame), so that it follows the PypeIt sign
convention.  This value is only used as a fallback by ``pypeit_coadd_2dspec``;
see :ref:`jwst_nirspec_dithoff`.

Selecting the slits
-------------------

By default, PypeIt reduces all the slits found in the ``_cal.fits`` files.
To reduce only some of them, set the ``maskIDs`` parameter in the
:ref:`parameter_block` of the :ref:`pypeit_file`.  The slits can be selected by
slit name (``SLTNAME``, e.g., ``S200A2`` for a fixed slit, or the slitlet
number for the MSA), by source ID (``SOURCEID``), or by source alias
(``SRCALIAS``):

.. code-block:: ini

    [rdx]
        spectrograph = jwst_nirspec
        maskIDs = 20013731, 183348

The detector(s) used for each slit are determined automatically: a slit that
falls only on NRS1 is reduced as ``DET01``, a slit that falls only on NRS2 as
``DET02``, and a slit that falls on both detectors is reduced as a mosaic of
the two (``MSC01``).  The ``detnum`` parameter is not used.

Sky and nod subtraction
-----------------------

By default (``bkg_id = -1``), the sky is modeled and subtracted by PypeIt.
For nodded observations, you can instead subtract the other nod positions by
setting the ``comb_id`` and ``bkg_id`` columns (see :ref:`a-b_differencing`).
For example, for an MSA 3-shutter nod:

.. code-block:: console

                                          filename | frametype | ... | dithpos | calib | comb_id | bkg_id
    jw03543001001_07101_00002_nrs1_assign_wcs.fits |   science | ... |       1 |     1 |       1 |    3,5
    jw03543001001_07101_00002_nrs2_assign_wcs.fits |   science | ... |       1 |     1 |       1 |    3,5
    jw03543001001_07101_00003_nrs1_assign_wcs.fits |   science | ... |       2 |     2 |       3 |    1,5
    jw03543001001_07101_00003_nrs2_assign_wcs.fits |   science | ... |       2 |     2 |       3 |    1,5
    jw03543001001_07101_00004_nrs1_assign_wcs.fits |   science | ... |       3 |     3 |       5 |    1,3
    jw03543001001_07101_00004_nrs2_assign_wcs.fits |   science | ... |       3 |     3 |       5 |    1,3

In this case, the global sky subtraction is usually skipped:

.. code-block:: ini

    [reduce]
        [[findobj]]
            skip_skysub = True

.. _jwst_nirspec_run:

Running PypeIt
==============

Run :ref:`run-pypeit` as usual:

.. code-block:: console

    run_pypeit jwst_nirspec_A.pypeit

For each calibration group, PypeIt loops over the selected slits.  The
calibrations of each slit are built from the JWST products of the exposure
(see :class:`~pypeit.calibrations.NIRSpecSlitCalibrations`):

    - the slit cutout (the ``SLTSTRT1``, ``SLTSIZE1``, ``SLTSTRT2``, and
      ``SLTSIZE2`` header keywords of the ``_cal.fits`` file) defines the
      region of the full-frame science image that is reduced;
    - the wavelength calibration and the tilts are derived from the
      ``WAVELENGTH`` extension of the ``_cal.fits`` file;
    - the bad-pixel mask flags the pixels without a valid wavelength (in the
      ``_cal.fits`` file) or flat-field value (in the
      ``_interpolatedflat.fits`` file), and the slit edges are traced from
      this mask;
    - the flat field combines the interpolated flat, the path-loss
      (``PATHLOSS_PS`` or ``PATHLOSS_UN``, depending on the source type) and
      bar-shadow (``BARSHADOW``) corrections, and the photometric conversion
      (``PHOTMJSR``) from the ``_cal.fits`` file.  The reduced spectra are
      therefore already flux calibrated, and no sensitivity function is needed.

The variance of the science image is computed from the ``VAR_RNOISE`` and
``VAR_POISSON`` extensions, and pixels flagged as ``DO_NOT_USE``, ``HOT``, or
``DEAD`` in the ``DQ`` extension are masked.

The name of the slit is added to all the output files.  For example, for slit
``S200A2`` the calibration files are named ``Slits_A_1_DET01_slitS200A2.fits.gz``,
``Flat_A_1_DET01_slitS200A2.fits``, etc., and the reduced files are named:

.. code-block:: console

    spec2d_jw01967012001_03102_00001-J2255+0251_NIRSPEC_20221028T062125.115_slitS200A2.fits
    spec1d_jw01967012001_03102_00001-J2255+0251_NIRSPEC_20221028T062125.115_slitS200A2.fits

.. _jwst_nirspec_coadd:

Coadding
========

2D coadd
--------

Use :ref:`pypeit_setup_coadd2d` to create a :ref:`coadd2d_file` for each slit,
selecting the spec2d files with ``--slitname`` and the detector (or mosaic)
on which the slit was reduced with ``--det``:

.. code-block:: console

    pypeit_setup_coadd2d -f jwst_nirspec_A.pypeit --slitname S200A2 --det 1

and then run :ref:`pypeit_coadd_2dspec <coadd2d>` on it.  The offsets between
the exposures can be measured from the brightest object (``offsets = auto``,
the default), or computed from the dither offsets recorded in the headers:

.. code-block:: ini

    [coadd2d]
        offsets = header

In the latter case, the dither offsets recorded in the headers of the spec2d
files (see :ref:`jwst_nirspec_dithoff`) are converted to pixels using the
platescale of the slit in each exposure (see :ref:`jwst_nirspec_platescale`).

.. _jwst_nirspec_dithoff:

Dither offsets
++++++++++++++

The ``YOFFSET`` header keyword (and therefore the ``dithoff`` column of the
:ref:`pypeit_file`) is the offset of the telescope pointing from the
reference point of the aperture (``APERNAME``).  For the fixed slits, this is
the center of the slit.  For the MSA (``APERNAME = NRS_FULL_MSA``), however,
it is a fixed point of the whole MSA, and ``YOFFSET`` includes the pointing of
the visit, which changes for different visits or MSA configurations.  Within a
visit, the differences between the ``YOFFSET`` values of the exposures are
the nod offsets, but the ``YOFFSET`` values of exposures taken in different
visits cannot be compared.

Instead, for each slit and exposure, PypeIt measures the position of the
source along the slit, relative to the center of the slit, using the source
position recorded by the JWST pipeline (``SRCYPOS``) and the slit GWCS in the
``_cal.fits`` file of the same exposure (see
:func:`~pypeit.spectrographs.jwst_nirspec.slit_source_offset_from_gwcs`).  The
negative of this offset is written to the ``DITHOFF`` keyword of the spec2d
and spec1d files, and it is used by ``pypeit_coadd_2dspec`` when ``offsets =
header``.  Because it is measured from the center of the slit in every
exposure, it can be used to coadd exposures taken in different visits, even
if the source fell in different slitlets.  For example, for an MSA source
observed in two visits with a 3-shutter nod, the offsets computed from
``YOFFSET`` differ from the actual ones by one nod step (:math:`\sim`\ 5
pixels) for the exposures of the second visit, while the offsets computed by
PypeIt agree with the actual ones to within a few hundredths of a pixel.

The JWST pipeline records the position of the source only for the slits that
include a planned source.  For the other slits (i.e., the fixed slits other
than the one targeted by the observation, ``FXD_SLIT``, and the MSA background
slitlets), and if the position of the source cannot be computed (e.g., because
``stdatamodels`` is not installed), PypeIt writes the ``dithoff`` value of the
:ref:`pypeit_file` (i.e., the negative of ``YOFFSET``), which should only be
used to coadd exposures of the same visit.

.. warning::

    The definition of the dither offset of NIRSpec data changed in PypeIt
    version 2.1.0 (both its sign and its reference point).  Spec2d files
    created with previous versions should be recreated before coadding them
    with ``offsets = header``.

1D coadd
--------

The 1D spectra of the individual exposures can be coadded with
:ref:`pypeit_coadd_1dspec <coadd1d>`.  Because the spectra are already flux
calibrated, set ``flux_value = False`` in the :ref:`coadd1d_file`:

.. code-block:: ini

    [coadd1d]
        coaddfile = J2255+0251_1dcoadd.fits
        flux_value = False

Detector parameters
===================

.. _jwst_nirspec_platescale:

Platescale
----------

The platescale of NIRSpec is not uniform.  Along the slit, the nominal value
is 0.105 arcsec/pixel, while in the dispersion direction it is 0.103
arcsec/pixel (`Jakobsen et al. 2022, A&A, 661, A80
<https://ui.adsabs.harvard.edu/abs/2022A%26A...661A..80J>`__).  However, the
actual along-slit platescale varies by a few percent across the field of
view: it is :math:`\sim`\ 0.1053 arcsec/pixel for the fixed slits, and it
ranges from :math:`\sim`\ 0.103 to :math:`\sim`\ 0.108 arcsec/pixel across the
MSA, depending on the quadrant and the shutter row.

PypeIt therefore computes the along-slit platescale of each slit from the slit
GWCS stored in the ``_cal.fits`` files (see
:func:`~pypeit.spectrographs.jwst_nirspec.slit_platescale_from_gwcs`), and uses
it for that slit instead of the nominal value.  For a slit that falls on both
detectors, the average of the values measured on the two detectors is used.
The value is recorded in the detector parameters of the spec2d files and of
the calibration files (e.g., ``Edges``), and is used, e.g., by the object
finding and by ``pypeit_coadd_2dspec`` when ``offsets = header``.

.. warning::

    Reading the slit GWCS requires the `stdatamodels`_ package.  If it is not
    installed, PypeIt issues a warning and uses the nominal platescale (0.105
    arcsec/pixel) for all slits, which can be off by up to :math:`\sim`\ 3%
    for some MSA slitlets.

Detector gap
------------

When a slit falls on both detectors, the NRS1 and NRS2 cutouts are combined
into a single mosaic image, with a gap of 146 pixels along the dispersion
direction between the two detectors.  This is the separation between the full
2048-pixel arrays given by the CRDS focal-plane-assembly (FPA) reference file,
and corresponds to the :math:`\sim`\ 2.8 mm gap between the light-sensitive
pixels quoted by Jakobsen et al. (2022) once the four reference columns on each
side are excluded.
