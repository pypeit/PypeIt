
.. include:: ../include/links.rst

.. _keck_osiris:

***********
Keck OSIRIS
***********

Overview
========

This page summarizes the instrument-specific items for the spectrograph of
`Keck OSIRIS <https://www2.keck.hawaii.edu/inst/osiris/>`__ (Larkin et al.
2006), the Keck I adaptive-optics near-infrared *lenslet* integral-field
spectrograph (``keck_osiris``).  For a step-by-step reduction of an example
night see the :doc:`../tutorials/keck_osiris_howto` tutorial.

The OSIRIS lenslet array produces up to ~3000 interleaved spectra on the
2048 x 2048 H2RG detector, only two pixels apart, so the individual lenslet
spectra cannot be isolated by slit tracing.  As in the Keck OSIRIS data
reduction pipeline (DRP; Lyke et al. 2017; Lockhart et al. 2019), PypeIt
separates them with the *rectification matrices* measured by Keck (a model
of the detector footprint of every lenslet), and then runs its ``Fiber``
pypeline on the rectified frame, in which every lenslet is one "fiber":

#. Each raw frame is read (``SCI`` extension, ADU/s/coadd, converted to ADU),
   corrected for the column offsets measured on the unilluminated reference
   rows, cleaned of cosmic rays with L.A.Cosmic, and *rectified*: a
   least-squares solve with the rectification matrix turns the raw image
   into a *lenslet image* with one spectral row per detector column and
   one spatial pixel per lenslet.
#. The lenslet image is cut into *block-slits*, one per Keck lenslet column
   (19 for the broadband filters, 51 for the narrowband filters), each
   holding the 64 lenslets of that column.
#. Wavelengths come from the Keck per-lenslet wavelength coefficients
   (vacuum), refined on the OH sky lines of every science frame.  The tilts
   are derived from the same coefficients; no arc or flat frames are used.
#. Background (sky) frames are subtracted frame-by-frame (A-B), the residual
   OH lines are scaled with the Davies (2007) method, and every lenslet is
   extracted to a spec1d spectrum.
#. :ref:`pypeit_fiber_cube <keck_osiris_cube>` turns the spec1d files into
   datacubes, or sums the lenslets of a point source into a single spectrum
   for the sensitivity function and the telluric correction.

The rectification, the lenslet geometry and the Davies OH scaling follow the
OSIRIS DRP (https://github.com/Keck-DataReductionPipelines/OsirisDRP); we
gratefully acknowledge its authors.

Supported data
==============

- Only data taken with the H2RG spectrograph detector, i.e. on or after
  2016-01-01 (MJD 57388), are supported.  Earlier (Hawaii-2) frames are
  flagged by ``pypeit_setup`` and refused by ``run_pypeit`` with a clear
  error; use the Keck DRP for them.
- Tested modes: Kbb/0.050" (2024, the development-suite science-mode setup),
  Hn3/0.100" (2016) and Kn3/0.035" (2024).  The other broadband (Zbb, Jbb,
  Hbb, Kbb, Kcb) and narrowband filters and the 0.020" and 0.035"/0.050"/
  0.100" scales use the same code paths and only need the matching
  rectification matrix.
- Headers of all eras are read, including the truncated 2016 headers
  without ``ISSKY``, ``PA_SPEC``, ``DTMP7``, ``INSTRUME`` or weather cards
  (defaults are used, with a warning where it matters).  A header without
  ``SFILTER`` or ``MJD-OBS`` cannot be reduced and gives an error naming
  the missing card.

Data requirements
=================

OSIRIS observers do not take arcs or flats.  A reduction needs:

- the science frames and their **sky frames** (nod-off frames, ``ISSKY =
  1``) or a nodded partner (A-B pairs within the field);
- an A0V **telluric standard** (optional; same filter and scale);
- **darks** only for frames without a background frame (see
  :ref:`keck_osiris_darks`);
- the **rectification matrix** of the filter, scale and epoch (see
  :ref:`keck_osiris_matrices`).

Raw files are available from the `Keck Observatory Archive
<https://koa.ipac.caltech.edu>`__ (``lev0`` files; public data can be
downloaded without a login).

.. _keck_osiris_matrices:

Rectification matrices
======================

Keck builds one rectification matrix per (filter, scale) and *validity
epoch* from calibration scans taken whenever the instrument changes.  The
matrices are public and listed by epoch at
http://tkserver.keck.hawaii.edu/osiris/ (plain HTTP only; the server does
not support modern TLS, so use ``http://``).  The files (152 MB each,
30-90 MB gzipped) are named ``sYYMMDD_cNNN___infl_<filter>_<scale>.fits.gz``
and live in ``data/<epoch>/`` sub-directories, e.g.:

.. code-block:: bash

    mkdir -p osiris_matrices
    curl -o osiris_matrices/s240524_c006___infl_Kbb_050.fits.gz \
        http://tkserver.keck.hawaii.edu/osiris/data/230701-250310/s240524_c006___infl_Kbb_050.fits.gz

Put the matrices (gzipped or not) in one directory and point PypeIt to it
in the :ref:`pypeit_file`; environment variables are expanded:

.. code-block:: ini

    [calibrations]
        osiris_rectmat_dir = /path/to/osiris_matrices

If the parameter is not set, the directory of the raw files and its
``calib`` sub-directory are searched (with a warning).

**Epoch rule.**  The file date of a matrix is the date of the calibration
scan, which can be *after* the observations it serves (e.g.
``s240524_c006___infl_Kbb_050`` serves 2024-02-03).  PypeIt therefore
selects, among the matrices of the filter and scale of the frame,

#. the one whose file date lies inside the validity epoch containing the
   observation (``MJD-OBS``; the epochs are shipped with PypeIt as
   ``rectmat_epochs.ecsv``, copied from the tkserver listing; an
   observation that falls in a gap between epochs uses the nearest epoch,
   with a warning);
#. otherwise the latest matrix dated on or before the observation, with a
   warning naming the epoch to download.

If no matrix qualifies, the error lists the matrices and the other files
found in the directory (names must follow the Keck convention), the epoch
that applies, and its tkserver sub-directory.  The selected matrix is
reported in the log (``Rectifying with the lstsq solver and ...``).

Setups and frame typing
=======================

A configuration (setup) is defined by the filter (``dispname``, header
``SFILTER``) and the lenslet scale (``decker``, header ``SSCALE``).
``pypeit_setup`` types the frames as follows:

- ``dark``: ``SFILTER = 'Drk'``;
- ``standard``: the target name carries a spectral type (e.g.
  ``HD101060(A0V)``) or is a bright-star catalogue name (HD, HIP, HR, ...)
  with an exposure time <= 60 s;
- ``science`` (and ``arc``, ``tilt``, ``trace``): every other on-sky frame;
  frames with ``ISSKY = 1`` are typed ``science,sky``.

No lamp frames exist: the ``arc``, ``tilt`` and ``trace`` types only tell
PypeIt which frame provides the header for the fixed calibrations (block
edges, Keck tilts and wavelengths), and its rectified image is kept as the
``Arc`` calibration frame (used by the OH-refinement fallback below).  One
long exposure per setup is enough; assigning these types to a single A/B
pair avoids rectifying and combining many frames for nothing.

.. _keck_osiris_bkg:

The ``bkg_id`` convention
-------------------------

Sky subtraction is done by frame differencing (see
:ref:`a-b_differencing`).  Give every frame its own ``comb_id`` and set the
``bkg_id`` of each science (or standard) frame to the ``comb_id`` of its
sky frame.  Sky frames that should not be reduced on their own are typed
``sky`` only.  For example (columns trimmed):

.. code-block:: console

                     filename |              frametype | ... | issky | calib | comb_id | bkg_id
    OS.20240203.42705.45.fits | arc,science,trace,tilt | ... |     0 |     0 |      12 |     13
    OS.20240203.43383.79.fits |     arc,sky,tilt,trace | ... |     1 |     0 |      13 |     -1
    OS.20240203.44064.64.fits |                science | ... |     0 |     0 |      14 |     13
    OS.20240203.44757.73.fits |                science | ... |     0 |     0 |      15 |     16
    OS.20240203.45436.78.fits |                    sky | ... |     1 |     0 |      16 |     -1

Nodded standard-star pairs (both positions on the field) can use each
other as background: ``bkg_id`` of the first is the ``comb_id`` of the
second and vice versa.

.. _keck_osiris_darks:

Darks and ITIME grouping
------------------------

No dark is subtracted from A-B pairs: the frame difference removes the dark
current and the detector glow (``use_darkimage = False`` by default).  For
frames without a background frame, turn the dark subtraction on and give
the darks of the same exposure time (PypeIt does not rescale darks) the
calibration group of those frames.  Darks are never rectified: they are
subtracted from the raw frame before its rectification.  The development
suite groups the darks of the 2024-02-03 night by ITIME in calibration
groups 1-3 (600, 1.5 and 8.85 s).

Processing details
==================

Rectification
-------------

The default solver (``[calibrations] osiris_rect_method = lstsq``) is an
exact block-tridiagonal weighted least-squares solve with flux-independent
weights (read noise plus a smoothed background variance), whose variance is
propagated exactly; ``drp`` is a port of the iterative DRP solver
(``spatrectif_000.c``).  Both reproduce the DRP reference cube of the
development suite (median spaxel difference 0.7-0.8 %).  Before the solve, the per-column median of
the unilluminated rows 2032-2044 is subtracted and cosmic rays are flagged
with L.A.Cosmic on each raw frame (``osiris_rect_crmask``); flagged pixels
are excluded from the solve.  A rectification takes ~10-15 s per frame.

Hits that L.A.Cosmic misses on the single raw frames survive as narrow
positive (A) or negative (B) spikes in the A-B spectra.  A second pass on
the sky-subtracted lenslet spectra flags spikes above
``[reduce][skysub] fiber_cr_nsig`` (8 sigma) that exceed three times the
residual of the neighbouring lenslets at the same wavelength (real lines
are shared by the neighbours).  Standards are not searched.

Wavelengths and tilts
---------------------

``[calibrations][wavelengths] method = fixed``: the Keck per-lenslet
wavelength coefficients of the epoch of the observation (shipped with
PypeIt) give the vacuum wavelength of every lenslet; the spec1d files carry
``AIRORVAC = 'vac'``.  Within a block the lenslets are staggered by ~1.9
pixels per lenslet row; the tilts of each block are a 2D Legendre fit to
the Keck solutions (0.02 px robust rms) and every lenslet keeps its exact
Keck solution in the spec1d.

The Keck solutions are refined on the OH lines of every science frame
(``[calibrations][wavelengths] fiber_sky_refine = zeropoint``; line list
``OH_OSIRIS`` in K, ``OH_NIRES`` otherwise): isolated lines are centroided
in every lenslet of the un-subtracted frame and the median offset is
applied per lenslet (on the 2024 epoch the Keck solution is off by +1.5 A,
+0.54 px).  Lenslets with fewer than three lines get the frame median; short
frames with few measurable lenslets use the shift of the median-stacked
spectrum.  **Frames without any measurable sky line** (e.g. 1.5-s standard
stars) adopt the per-lenslet shifts measured on the ``Arc`` frame of their
calibration group -- a long science or sky exposure of the same night and
configuration -- with a warning.  A QA plot of the shifts is written to
``QA/PNGs/*_fiber_sky_refine.png``.

The DRP grating-temperature term (``DTMP7``) moves the OH lines away from
their vacuum wavelengths by +0.4 A on the 2016 and 2024 test data and is off
by default (``[calibrations] osiris_grating_temp = False``).

OH-line scaling
---------------

The OH lines vary between the A and B exposures.  With
``[reduce][skysub] osiris_oh_scaling = global`` (the default) the residual
lines of every A-B frame are corrected as in Davies (2007) and the DRP
(``scaledskysub_000.pro``): the lines of the sky frame are grouped by
vibrational transition and one scale factor per group, pooled over all
lenslets (``block``: per block; ``fiber``: per lenslet), is fit between the
A-B and the sky spectra; the rescaled lines (not the continuum) form the
sky model of the A-B frame (spec2d ``SKYMODEL``).  Of the significant
factors, the combination that most lowers the median integrated residual of
the isolated OH lines (relative to their sky flux) is kept, and none if the
median does not decrease; this keeps short exposures with weak lines from
being degraded.  ``osiris_oh_min_exptime`` (seconds)
skips the scaling for shorter frames altogether.  Standard stars are not
scaled.

Parameters
----------

The OSIRIS defaults are set in
:func:`~pypeit.spectrographs.keck_osiris.KeckOSIRISSpectrograph.default_pypeit_par`;
see :ref:`parameters` for their descriptions.  The OSIRIS-specific and
fiber-pypeline parameters are:

=====================================================  ===============  ==========================================
Parameter                                              OSIRIS default   Purpose
=====================================================  ===============  ==========================================
``[calibrations] osiris_rectmat_dir``                  None             directory of the Keck matrices
``[calibrations] osiris_rect_method``                  ``lstsq``        rectification solver (``lstsq``, ``drp``)
``[calibrations] osiris_rect_niter``, ``_relax``       40, 1.0          ``drp`` solver iterations / relaxation
``[calibrations] osiris_rect_crmask``                  ``lacosmic``     raw-frame cosmic rays (``drp2016``, ``none``)
``[calibrations] osiris_grating_temp``                 False            DRP ``DTMP7`` expansion term
``[calibrations][wavelengths] method``                 ``fixed``        Keck wavelength coefficients
``[calibrations][wavelengths] fiber_sky_refine``       ``zeropoint``    per-lenslet OH refinement (``block``)
``[reduce][skysub] osiris_oh_scaling``                 ``global``       Davies OH scaling pooling
``[reduce][skysub] osiris_oh_min_exptime``             None             no OH scaling below this exposure time
``[reduce][skysub] fiber_cr_nsig``                     8.0              second cosmic-ray pass on A-B spectra
``[sensfunc] polyorder``, ``[telluric] polyorder``     7                instrument response (filter ripples)
=====================================================  ===============  ==========================================

Outputs
=======

The spec2d file holds the rectified lenslet image (one column per lenslet)
and the spec1d file one spectrum per lenslet (1216 for the broadband
filters, up to 3072 for the narrowband filters), with
``MASKDEF_ID = row * 51 + col`` and ``MASKDEF_OBJNAME = L<col>_<row>`` of
the Keck lenslet grid.  For A-B frames the spec1d sky columns are zero
(the sky was removed by the frame difference); the Davies model is in the
spec2d ``SKYMODEL``.

.. _keck_osiris_cube:

Datacubes and point sources
===========================

``pypeit_fiber_cube`` builds a datacube (wavelength, lenslet row, lenslet
column) from one or more spec1d files, with the plate scale from ``SSCALE``
and a WCS referenced to the lenslet-grid origin (the DRP reference spaxel)
and rotated by the position angle of the frame:

.. code-block:: bash

    pypeit_fiber_cube Science/spec1d_OS.20240203.42705.45-NGC3245_OSIRIS_20240203T115145.453.fits

With ``--point_source R`` it also sums the lenslets within ``R`` arcsec of
the brightest source into a single-object spec1d file (``*_ps.fits``), the
input of the sensitivity-function and telluric tools.  One cube is built
per exposure; co-adding cubes is not yet supported.

.. _pypeit_fiber_cube:

pypeit_fiber_cube
-----------------

The script is generic for the ``Fiber`` pypeline (``pypeit_binospec_ifu_cube``
is an alias for MMT Binospec); the spectrograph supplies the plate scale and
WCS through its ``fiber_cube_setup`` hook.

.. include:: ../help/pypeit_fiber_cube.rst

Telluric correction
===================

Use an A0V star of the same filter and scale, reduced as an A-B pair:

.. code-block:: bash

    # Point-source spectrum of the standard (no cube)
    pypeit_fiber_cube Science/spec1d_OS.20240203.59407.32-HD101060\(A0V\)_OSIRIS_20240203T163007.324.fits \
        --point_source 0.4 --no_cube
    # Sensitivity function (IR algorithm)
    pypeit_sensfunc Science/spec1d_OS.20240203.59407.32-HD101060\(A0V\)_OSIRIS_20240203T163007.324_ps.fits \
        -s keck_osiris_kbb_050.sens
    # Telluric correction with the stellar model
    pypeit_tellfit Science/spec1d_OS.20240203.59407.32-HD101060\(A0V\)_OSIRIS_20240203T163007.324_ps.fits \
        -t keck_osiris_kbb_050.tell

with the input files

.. code-block:: ini

    # keck_osiris_kbb_050.sens
    [sensfunc]
        algorithm = IR
        star_type = A0
        star_mag = 8.78

    # keck_osiris_kbb_050.tell
    [telluric]
        objmodel = star
        star_type = A0
        star_mag = 8.78

The default telluric grid is ``TellPCA_3000_26000_R10000.fits`` and the
instrument-response polynomial has order 7 (the broadband filters have
few-percent transmission ripples).  On the development-suite star the
corrected spectrum is flat to 1.5 % rms in 50-A bins.  See
:ref:`fluxing` and :ref:`telluric_correction`.

Limitations
===========

- Pre-2016 (Hawaii-2) data are not supported.
- The rectification matrix of each (filter, scale, epoch) must be
  downloaded by the user (above); PypeIt does not fetch it.
- Lamp flats are not used: lenslet-to-lenslet throughput variations are
  those of the Keck matrices (which are normalised per lenslet).
- One cube per exposure; dithered exposures must be registered before
  being combined (``pypeit_fiber_cube`` does not co-add yet).  Maps built
  by stacking spec1d files in lenslet coordinates mix the dither positions.
- Emission-line maps measured in a fixed wavelength window can show sharp,
  straight edges where the line velocity leaves the window (e.g. the
  rotating nuclear gas of NGC 3245 in the development-suite data, +-230
  km/s); integrate over the velocity range of the line.  This is not an
  instrumental artefact: the development-suite analysis found no detector
  or rectification feature at the corresponding lenslet rows.
- The ``drp`` solver and the ``drp2016`` cosmic-ray rule are kept for
  comparison with the Keck DRP; the defaults are recommended.
- The 2026 TRICK-dichroic epoch (``260604-YYMMDDTRICK`` on the tkserver
  listing) has not been tested.

References
==========

- Larkin, J. et al. 2006, SPIE 6269, 62691A (OSIRIS)
- Lyke, J. et al. 2017, Astrophysics Source Code Library, ascl:1710.021
  (OSIRIS DRP)
- Lockhart, K. et al. 2019, AJ 157, 75 (the OSIRIS DRP; characterisation
  and improvements)
- Davies, R. I. 2007, MNRAS 375, 1099 (OH-line scaling)

