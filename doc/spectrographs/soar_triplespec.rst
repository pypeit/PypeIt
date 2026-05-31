***************
SOAR TripleSpec
***************

Overview
========

This file summarizes several instrument specific settings that are related to
the SOAR/TripleSpec spectrograph (PypeIt name ``soar_tspec``).

TripleSpec is a fixed-format, cross-dispersed near-infrared spectrograph.  PypeIt
reduces it with the Echelle pipeline, covering five orders (order numbers 7
through 3) spanning roughly 0.8--2.47 microns (:math:`\sim`\ 8000--24700
Angstrom, i.e. the *z/J/H/K* range).  Because the format is fixed, the order
positions, order numbers, and approximate spectral ranges are hard-coded in the
spectrograph class.


PypeIt File
===========

Here is some advice on how to setup your :ref:`pypeit_file`.  First, run:

.. code-block:: console

    pypeit_setup -r absolute_path -s soar_tspec -b -c A

where ``-b`` indicates that the data use sky subtraction and the ``calib``,
``comb_id``, and ``bkg_id`` columns are added to the :ref:`data_block`.  See
:ref:`pypeit_setup` and :doc:`../A-B_differencing` for the syntax used for the
data in these columns and how PypeIt uses them.

Here is an example of the :ref:`data_block` of the PypeIt file:

.. code-block:: console

    # Data block
    data read
     path ./
                                    filename |       frametype | ... |  dispname |  decker | binning | exptime | calib | comb_id | bkg_id
             SPEC_DFLAT_04-02-2025_0088.fits |            dark | ... |      spec | default |     1,1 |     2.0 |     0 |      -1 |     -1
             SPEC_DFLAT_04-02-2025_0042.fits | pixelflat,trace | ... |      spec | default |     1,1 |     2.0 |     0 |      -1 |     -1
            SPEC_2025cy_05-02-2025_0116.fits | arc,tilt,science | ... |     spec | default |     1,1 |   300.0 |     0 |       1 |     -1
            SPEC_2025cy_05-02-2025_0117.fits | arc,tilt,science | ... |     spec | default |     1,1 |   300.0 |     0 |       2 |     -1
    SPEC_2025cy_05-02-2025_standard0120.fits |        standard | ... | spec | default |     1,1 |    9.25 |     0 |       3 |     -1
    data end

``frametype`` is automatically assigned to each frame using the values of
various header keywords (see :meth:`~pypeit.spectrographs.soar_tspec.SOARTSPECSpectrograph.check_frame_type`):
the dome flats taken with the lamp **on** are typed ``pixelflat,trace``, the
dome flats taken with the lamp **off** are typed ``dark``, and the science
exposures are typed ``science``.


Calibrations
============

Wavelength Calibration
----------------------

TripleSpec has no dedicated arc lamp that is useful across the full near-IR
range; instead, the dense forest of **OH sky-emission lines** in the science
frames is used as the wavelength reference.  You must therefore add the ``arc``
and ``tilt`` frame types to the science frames in the PypeIt file, as shown in
the data block above, and comment out (or omit) any short dedicated arc
exposures.  The long (e.g. 300 s) science exposures contain far more usable OH
lines than a short arc and give a much better solution, especially in the
reddest order.

PypeIt ships a native reidentification archive,
``soar_triplespec.fits``, that covers all five orders.  The default parameters
(set in :meth:`~pypeit.spectrographs.soar_tspec.SOARTSPECSpectrograph.default_pypeit_par`)
use it via the ``reidentify`` method with the ``OH_NIRES`` line list, which
carries OH lines all the way out to :math:`\sim`\ 24980 Angstrom and so anchors
the reddest order (order 3, *K* band) with real lines.  No user intervention is
normally required.

.. warning::

    If you change the frame typing (for example, switching which frames are
    used as the ``arc``), you must **delete the affected files in
    ``Calibrations/``** before re-running.  ``run_pypeit -o`` overwrites only the
    *science* products, not the cached calibration frames, so a stale
    ``Arc_*`` / ``WaveCalib_*`` will silently be reused.

Flat Fielding
-------------

Dome flats are taken in pairs: lamp **on** and lamp **off**.  PypeIt types the
lamp-on flats as ``pixelflat,trace`` and the lamp-off flats as ``dark``; the
dark is subtracted to remove the thermal/dome background before the field flat
and order edges are constructed.  Bias subtraction and overscan correction are
turned off for this detector.


Additional Reading
==================

Here are additional docs relevant to reducing SOAR/TripleSpec data:

- :doc:`../tutorials/soar_triplespec_howto` -- a full worked example reduction.
- :doc:`../calibrations/wave_calib`
- :doc:`../calibrations/flat`
- :doc:`../A-B_differencing`
