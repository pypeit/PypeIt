
.. include:: ../include/links.rst

.. _keck_osiris_howto:

=================
Keck OSIRIS HOWTO
=================

Overview
========

This tutorial reduces one night of `Keck OSIRIS
<https://www2.keck.hawaii.edu/inst/osiris/>`__ lenslet integral-field data
with PypeIt: the nucleus of NGC 3245 observed in Kbb at 0.050"/lenslet on
2024-02-03, with nodded sky frames, the A0V telluric standard HD101060 and
darks.  The same files form the ``keck_osiris/Kbb_050`` setup of the
:ref:`dev-suite`.  Read the :ref:`keck_osiris` page first for the
instrument-specific concepts (rectification matrices, block-slits, the
``bkg_id`` convention).

The reduction takes ~40 minutes for the 50 frames reduced here on a laptop
(~15 s per frame for the rectification, the rest for the extraction of
~1200 lenslet spectra per frame).

Data
====

The raw files come from the `Keck Observatory Archive
<https://koa.ipac.caltech.edu>`__: the 600-s NGC 3245 object frames
(``ISSKY = 0``) and their sky frames (``ISSKY = 1``), the HD101060 frames and
the darks of the night.  Download the Kbb/0.050" rectification matrix of
the validity epoch of the night (2023-07-01 to 2025-03-10) from the Keck
listing:

.. code-block:: bash

    mkdir -p Kbb_050/calib
    curl -o Kbb_050/calib/s240524_c006___infl_Kbb_050.fits.gz \
        http://tkserver.keck.hawaii.edu/osiris/data/230701-250310/s240524_c006___infl_Kbb_050.fits.gz

The matrix is dated after the observations; PypeIt selects it because its
date lies inside the validity epoch of 2024-02-03 (see
:ref:`keck_osiris_matrices`).

Setup
=====

Run :ref:`pypeit_setup` on the raw directory:

.. code-block:: bash

    pypeit_setup -r /path/to/Kbb_050/ -s keck_osiris -c all

This writes one :ref:`pypeit_file` per configuration (filter and lenslet
scale; here ``keck_osiris_A`` for Kbb/0.050").  The darks, taken at the
0.020" scale with ``SFILTER = Drk``, are typed ``dark``; the standards are
recognised from their target names; frames with ``ISSKY = 1`` are typed
``science,sky``.

Edit the pypeit file
====================

Four edits turn the automatically generated file into the one used for the
reduction (the development-suite file, ``keck_osiris_kbb_050.pypeit``, can
be used as a template):

#. **Matrix directory**:

   .. code-block:: ini

       [calibrations]
           osiris_rectmat_dir = /path/to/Kbb_050/calib

#. **Sky frames**: type the ``ISSKY = 1`` frames ``sky`` only, so they serve
   as backgrounds and are not reduced on their own.

#. **Pairs**: give each frame its own ``comb_id`` and set the ``bkg_id`` of
   every object frame to the ``comb_id`` of the nearest sky frame
   (:ref:`keck_osiris_bkg`).  Nodded standard-star pairs use each other.

#. **Calibration frame**: keep the ``arc,tilt,trace`` types on one long A/B
   pair only (here ``OS.20240203.42705.45`` and its sky
   ``OS.20240203.43383.79``).  OSIRIS has no lamp calibrations; these frames
   only provide the header (filter, date) for the fixed calibrations and the
   OH-line reference for very short exposures.

The darks can stay in the file (calibration groups 1-3 by exposure time);
they are not used for A-B pairs.

Run the reduction
=================

.. code-block:: bash

    run_pypeit keck_osiris_A.pypeit -o

The log reports, for every frame, the rectification (``Rectifying with the
lstsq solver and s240524_c006___infl_Kbb_050.fits.gz``), the per-lenslet OH
refinement (``Sky-line refinement (zeropoint) ... median shift +0.54 px``),
the Davies OH scaling of A-B frames (``OH scaling (global) group 8-6: r - 1
= ...``) and the residual cosmic rays flagged on the A-B spectra.  The
1.5-s HD89239 frames have no measurable sky line; they adopt the shifts of
the arc frame, with a warning.

Inspect the outputs
===================

- ``Science/spec2d_*.fits``: the rectified lenslet image (one column per
  lenslet, 19 block-slits of 64 lenslets), the Davies sky model and the
  wavelength image; view them with ``pypeit_show_2dspec``.
- ``Science/spec1d_*.fits``: one spectrum per lenslet (``MASKDEF_ID = row
  * 51 + col`` of the Keck grid; ``MASKDEF_OBJNAME = L<col>_<row>``),
  vacuum wavelengths.  ``pypeit_show_1dspec --list`` lists the spectra;
  ``--obj`` selects one by name.
- ``QA/PNGs/*_fiber_sky_refine.png``: the per-lenslet OH shifts.

Datacube
========

.. code-block:: bash

    pypeit_fiber_cube Science/spec1d_OS.20240203.42705.45-NGC3245_OSIRIS_20240203T115145.453.fits

writes ``cube_OS.20240203.42705.45-NGC3245_OSIRIS_20240203T115145.453.fits``,
a (wavelength, 66, 21) cube with 0.05" spaxels, 2.8-A steps and a WCS
referenced to the lenslet-grid origin.  See :ref:`keck_osiris_cube`.

Telluric correction
===================

Extract the point-source spectrum of the A0V star, compute the sensitivity
function and fit the telluric model:

.. code-block:: bash

    pypeit_fiber_cube "Science/spec1d_OS.20240203.59407.32-HD101060(A0V)_OSIRIS_20240203T163007.324.fits" \
        --point_source 0.4 --no_cube
    pypeit_sensfunc "Science/spec1d_OS.20240203.59407.32-HD101060(A0V)_OSIRIS_20240203T163007.324_ps.fits" \
        -s keck_osiris_kbb_050.sens
    pypeit_tellfit "Science/spec1d_OS.20240203.59407.32-HD101060(A0V)_OSIRIS_20240203T163007.324_ps.fits" \
        -t keck_osiris_kbb_050.tell

The input files are listed on the :ref:`keck_osiris` page (IR algorithm,
``star_type = A0``, ``star_mag = 8.78``, ``objmodel = star``).  The
corrected spectrum of the standard is flat to ~1.5 % rms in 50-A bins.

