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
separate FITS extensions.  For frames with at least 3 reads, PypeIt performs
up-the-ramp fitting with likelihood-based jump (cosmic-ray) detection using
the algorithm of `Brandt (2024) <https://arxiv.org/abs/2404.01326>`__,
replacing the correlated double sampling (first minus last read) used
previously.  Frames with 2 reads still use correlated double sampling.  Each
read is reference-pixel corrected before fitting.  This implementation was
inspired by the prototype at
`mmt-mmirs-up-the-ramp-pypeit
<https://github.com/zhechenghu/mmt-mmirs-up-the-ramp-pypeit>`__.

The single-read noise needed by the fit is calibrated from a dark frame
listed in the :ref:`pypeit_file` when one with at least 10 reads is
available (include darks in your raw-data directory when running
:ref:`pypeit_setup` to enable this); otherwise it is self-calibrated from
each frame by rescaling the fit chi-squared.  Both start from the measured
per-read noise of 11 e- (at gain 1) reported on the `MMIRS instrument
statistics page
<https://lweb.cfa.harvard.edu/mmti/mmirs/instrstats.html>`__.  The effective read noise of
the fitted image, ``sigma * sqrt(12 (N-1) / (N (N+1)))`` for ``N`` reads, is
propagated to the detector parameters.

Expect roughly 2 minutes of processing and ~5 GB of memory to fit a 69-read
frame; flats with few reads take seconds.  See `Preprocessed ramp images`_
below for how PypeIt avoids repeating this cost each time a science frame
is loaded.

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

    pypeit_mmirs_ramp raw/*.fits

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
against the raw cube's modification time; use ``pypeit_mmirs_ramp --force``
to re-fit with the new calibration.

Multislit observations
++++++++++++++++++++++

The PypeIt developers only tested MMIRS longslit observations, but,
in principle, PypeIt should also work for multi-slit data.

HK Grism and zJ filter
++++++++++++++++++++++

This is an unusual setup that people might not use frequently.
This setup has both first and second order light. Wavelengths
greater than 9000 angstrom are in second order at the wavelength shown.
Wavelengths shown below 9000 are in first order at double the
wavelength shown. The major advantage of this setup is that
it provides the bluest wavelength coverage by the second
order spectrum and therefore are important for some science
cases. Currently, PypeIt only reduces the second order spectrum
of this setup.

