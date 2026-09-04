.. include:: ../include/links.rst

.. _mmt_mmirs:

*********
MMT MMIRS
*********

Overview
========

This file summarizes several instrument specific
items for the MMT/MMIRS spectrograph.


up-the-ramp fitting
+++++++++++++++++++

MMIRS raw frames store the non-destructive reads of the HAWAII-2 detector as
separate FITS extensions.  For frames with at least 5 reads, PypeIt performs
up-the-ramp fitting with likelihood-based jump (cosmic-ray) detection using
the algorithm of `Brandt (2024) <https://arxiv.org/abs/2404.01326>`__,
replacing the correlated double sampling (first minus last read) used
previously.  Frames with fewer reads still use correlated double sampling
(fewer than 5 reads do not sample the ramp well enough to fit reliably).
Each read is reference-pixel corrected before fitting.  This implementation
was inspired by the prototype at
`mmt-mmirs-up-the-ramp-pypeit
<https://github.com/zhechenghu/mmt-mmirs-up-the-ramp-pypeit>`__.

The per-read noise needed by the fit is not the instantaneous detector read
noise but an **effective** noise: the read noise plus the shot noise from dark
current (and any flux) that accumulates over the ramp.  It therefore *grows
with exposure time* -- measured on MMIRS darks, roughly 5.7 e- for a 3-read
(3 s) ramp, 8.3 e- for 8 reads (10 s), and 9.4 e- for 69 reads (300 s) -- so it
must be measured at the science exposure time.

PypeIt calibrates it from the dark frames listed in the :ref:`pypeit_file`
(include darks in your raw-data directory when running :ref:`pypeit_setup` to
enable this).  Only darks whose exposure time matches the science frame and
which have at least 10 reads are used; each is calibrated independently by
rescaling the ramp-fit chi-squared, and the results are combined as an
**inverse-variance weighted mean** (each dark weighted by the bootstrap
uncertainty on its own calibrated noise).  With no matching dark, a frame with
at least 10 reads self-calibrates from itself; a frame with fewer reads is too
poorly constrained and instead uses a starting-guess effective noise of 9 e-.

The calibrated value is floored at the instantaneous read noise from the header
``RDNOISE`` card (a derived noise below the physical read noise is unphysical);
``RDNOISE`` is also used as the correlated-double-sampling read noise for
short-ramp frames.  The effective read noise of the fitted image,
``sigma * sqrt(12 (N-1) / (N (N+1)))`` for ``N`` reads, is propagated to the
detector parameters.

An explicit dark supplied to ``pypeit_fit_ramp --dark`` is used as given,
without the exposure-time match.

The per-pixel fit is independent across the detector, so it is performed in
parallel: the rows are split into small blocks that are fit concurrently in a
thread pool (NumPy releases the GIL during the element-wise arithmetic that
dominates the fit).  With the default of ``min(6, os.cpu_count())`` threads,
expect roughly 30-40 seconds and ~5 GB of memory to fit a 69-read cube of
2048x2048 frames (about 3x faster than the serial fit); cubes with few reads
take seconds.  See `Preprocessed ramp images`_ below for how PypeIt avoids
repeating this cost each time a science frame is loaded.

Performance tuning
^^^^^^^^^^^^^^^^^^^

Two :ref:`parameters` in the ``[rdx]`` block control the parallel fit:

.. code-block:: ini

    [rdx]
        spectrograph = mmt_mmirs
        ramp_fit_cores = 12        # worker threads (None -> min(6, os.cpu_count()); 1 disables)
        ramp_fit_chunk_rows = 16   # detector rows fit per call (expert knob)

The values that ship as defaults were tuned on a laptop (Apple M1, 16 GB).
Because the fit is **memory-bandwidth bound** rather than compute bound
(it is dominated by element-wise array arithmetic, not linear algebra), the
speedup plateaus at roughly ``ramp_fit_cores = 6`` on such a machine and can
*regress* beyond it as the memory bus saturates.  Raising ``ramp_fit_cores``
therefore pays off mainly on a workstation with more memory bandwidth (more
memory channels), where the bandwidth ceiling is higher.  Peak memory scales
roughly as ``ramp_fit_cores * ramp_fit_chunk_rows``, so increase the core
count only if you have both the cores and the RAM to spare.

``ramp_fit_chunk_rows`` is an expert knob: small blocks keep each thread's
working set resident in cache, and 16 rows was empirically near-optimal and
largely machine-independent.  Larger blocks raise peak memory and usually
reduce throughput.

Preprocessed ramp images
^^^^^^^^^^^^^^^^^^^^^^^^

Fitting a many-read ramp takes minutes, and PypeIt loads each science frame
several times during a reduction.  To avoid re-fitting, the fitted 2D
count-rate image (units of e-/s) is written to a ``RampFit`` directory
inside the reduction directory (alongside ``Calibrations``, ``Science``,
and ``QA``) the first time each cube is loaded, and reused on subsequent
loads.  A preprocessed image is re-fit automatically if the raw cube's
modification time changes; if the ``RampFit`` directory is not writable,
the fit proceeds in memory with a warning and nothing is cached.

The fitted images can also be created (and inspected) ahead of a reduction
with:

.. code-block:: bash

    pypeit_fit_ramp mmt_mmirs raw/*.fits

which writes into the ``RampFit`` directory under ``--odir`` (default: the
current directory, so run it from the reduction directory) and accepts
``--sig`` to force the single-read noise, ``--dark`` to calibrate it from a
dark cube, and ``--force`` to re-fit existing outputs.  Preprocessed files
carry the fit parameters in header cards (``RAMPSIG``, ``RAMPRON``,
``NGROUPS``) and preserve all raw metadata, so :ref:`pypeit_setup` can be
run directly on a ``RampFit`` directory if preferred.  ``run_pypeit`` never
requires the manual step.

Changing the noise-calibration source between runs (e.g. adding dark frames
to the raw-data directory, or forcing a different ``--sig``) does not
invalidate an existing preprocessed image, since freshness is only judged
against the raw cube's modification time; use ``pypeit_fit_ramp --force``
to re-fit with the new calibration.

Multislit observations
++++++++++++++++++++++

PypeIt supports both MMIRS longslit and multi-object (MOS) slitmask
observations. For MOS masks, PypeIt ingests the MMIRS ``.msk`` mask-design
file so that slit tracing, alignment-box registration, and target
identities all follow the mask design; see the mask-definition section below.

Supported grism/filter modes
++++++++++++++++++++++++++++

MMIRS offers several grism/filter combinations for longslit and slitmask
spectroscopy; see
the observatory's `MMIRS recommendations page
<https://www.mmto.org/mmirs-recommendations/>`__ for the full list and for
guidance on choosing among them.  For each mode below PypeIt ships an archival
wavelength solution (``reid_arxiv``) and automatically selects the
``full_template`` wavelength-calibration method (see
:func:`~pypeit.spectrographs.mmt_mmirs.MMTMMIRSSpectrograph.config_specific_par`):

=========  =========  ===============================
Grism      Filter     Approx. coverage
=========  =========  ===============================
J          zJ         0.95-1.50 um
H3000      H          1.44-1.86 um
K3000      Kspec      1.90-2.51 um
HK         HK3        1.14-2.51 um
HK         zJ         0.94-1.29 um (2nd order; below)
=========  =========  ===============================

All of the observatory's recommended spectroscopic modes, plus the ``HK``/``zJ``
setup, have been reduced end-to-end and verified against real data in the PypeIt
development suite.

HK Grism and zJ filter
++++++++++++++++++++++

This is an unusual setup that people might not use frequently.  The HK grism's
first order is the HK band itself (roughly 1.14-2.51 um), so it cannot reach the
z/J region in first order; that coverage is only available in **second order**,
where the dispersion is also doubled.  PypeIt reduces the second-order spectrum,
so the wavelength axis it produces (roughly 0.94-1.29 um) is the true
second-order wavelength.  On that axis, anything appearing below 9000 angstrom is
first-order leakage at double the indicated wavelength, not real z/J flux (recall
that second-order light at wavelength :math:`\lambda` lands at the same detector
position as first-order light at :math:`2\lambda`).  The major advantage of this
setup is that the second-order spectrum provides the bluest wavelength coverage,
which is important for some science cases.

A-B nod background subtraction
++++++++++++++++++++++++++++++

MMIRS science frames are commonly taken as an A-B nod sequence (dithering the
target along the slit between exposures), for both longslit and multi-object
(MOS) masks. PypeIt derives the along-slit dither offset of each frame from the
telescope pointing relative to the catalog target (the MMIRS header has no
dither keyword) and records it in the ``dithoff`` column, along with ``dithpos``
(A/B/A'/B' labels), ``dithpat`` (e.g. ``ABA'B'``), and ``frameno``.

During :ref:`pypeit_setup`, PypeIt automatically fills the ``comb_id`` and
``bkg_id`` columns so that each science frame is background-subtracted using the
temporally-adjacent frame at the opposite nod position. Pairing is done
separately for each target within an instrument configuration (so a standard
star, or a second target sharing the same longslit setup, is never paired
against an unrelated object). Each frame keeps a unique ``comb_id`` (frames are
not 2D-coadded across the small sub-dither), so PypeIt produces one
``spec2d``/``spec1d`` per nod pair. See :ref:`a-b_differencing` for how to edit
these columns by hand.

A sequence whose dither offsets span less than ~1 arcsec is treated as a stare
and left unpaired (``bkg_id = -1``).

Combine the resulting 1D spectra with :ref:`pypeit_coadd_1dspec`.

Mask definition (MOS)
+++++++++++++++++++++

For multi-object (MOS) masks, PypeIt can ingest the MMIRS ``.msk`` mask-design
file (the xfitmask output delivered with the raw data) so that every slit,
``SpecObj``, and coadd carries the real target identity. During
:ref:`pypeit_setup`, if a file named ``<MOSID>.msk`` (matching the ``APERTURE``
header keyword and the ``decker`` metadata column) is found in the raw-data
directory, PypeIt automatically:

- enables slitmask-design tracing
  (``[calibrations][slitedges] use_maskdesign = True``,
  ``maskdesign_filename`` pointing at the ``.msk``);
- stamps ``MASKDEF_ID``, ``RA``, ``DEC``, and ``MASKDEF_OBJNAME`` (the catalog
  target name) on each extracted object
  (``[reduce][slitmask] assign_obj = True``);
- flags the alignment-box (``BOX``) slits as alignment rather than science, and
  uses them to register the mask-to-detector offset
  (``use_alignbox = True``); and
- force-extracts designed targets that were not auto-detected, at their
  predicted positions (``extract_missing_objs = True``).

The mask ``y`` coordinate (mm) maps linearly, and anti-aligned, to the detector
spatial pixel (``spat = C - scale * y``, with
``scale = (1/arc2mm)/platescale/bin_spat`` px/mm); the mask ``x`` coordinate is
the dispersion axis and does not enter the spatial edge prediction. The
absolute registration offset is fit by PypeIt's
:func:`~pypeit.edgetrace.EdgeTraceSet.maskdesign_matching` from the detected
slits, so only the relative scale must be correct.

The ``SpecObj`` name stays ``SPATxxxx-SLITyyyy-DET`` (PypeIt's internal object
identifier, as for every multislit spectrograph); the catalog identity lives in
the ``MASKDEF_OBJNAME``, ``MASKDEF_ID``, ``RA``, and ``DEC`` attributes. Those
are printed alongside the name by ``pypeit_show_1dspec <spec1d> --list`` and in
the ``.txt`` object summary written next to each ``spec1d`` file.

To coadd across nod pairs by target, use :ref:`collate1d`, which groups objects
by sky coordinate. For MMIRS the ``[collate1d] outfile_from`` parameter defaults
to ``maskdef_objname``, so each coadded file is named after the catalog target
(e.g. ``qso_172239.955+655201.69_MMIRS_<date>.fits``) rather than its sky
coordinate; set ``outfile_from = coord`` to restore coordinate-based names.

.. note::

    A coadded 1D spectrum (``OneSpec``) can contain masked bins where no input
    exposure contributed a good pixel (e.g. a persistent bad pixel or an OH
    residual that lands at the same wavelength in every flexure-controlled nod
    exposure). These bins are correctly flagged in the good-pixel mask (``flux``,
    ``ivar``, and ``mask`` are all zero there), but their stored wavelength is
    left at 0 rather than the grid wavelength. A viewer that plots the raw
    ``wave`` array *without* honoring the mask (such as ``lt_xspec``) therefore
    draws a spurious spike back toward ``wave = 0``. The data are fine — this is
    purely a display artifact. Load the file through `specutils`_ instead, which
    applies the mask and returns a strictly monotonic wavelength axis::

        from specutils import Spectrum
        import pypeit.specutils   # registers the PypeIt loaders
        spec = Spectrum.read('coadd/qso_172239.955+655201.69_MMIRS_<date>.fits')
        # spec.spectral_axis (Angstrom) and spec.flux have the masked bins dropped

When no ``.msk`` file is present, PypeIt falls back to generic edge tracing with
``SPATxxxx-SLITyyyy`` object names; group the resulting 1D spectra for coaddition
by slit id in the :ref:`pypeit_coadd_1dspec` input file.

