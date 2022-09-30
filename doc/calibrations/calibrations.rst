
.. _calibrations:

============
Calibrations
============

Overview
========

This doc summarizes the main calibrations performed by PypeIt.  As per standard
spectroscopic reduction, there are a series of calibrations files generated to
correct the detector response and perform wavelength calibration.

Although fluxing is certainly a calibration step, this is *not* performed as a
standard part of PypeIt's main data reduction procedure, :ref:`run-pypeit`.
Instead it is treated as a post-processing procedure; see :doc:`../fluxing`.

.. _calibration-groups:

Calibration Groups
==================

After :ref:`pypeit_setup` has organized your files into distinct instrument
configurations, the simplest approach (and PypeIt's default approach) is to use
all calibrations for all science frames; however, that isn't always the best
approach.  For example, you may want to associate different science frames with
calibrations taken either in the evening or the morning.  To ensure the best
calibrations are used with each science frame, you could split the data between,
e.g., two :ref:`PypeIt Files<pypeit_file>`, or you can use PypeIt's calibration
group designations.

To assign specific calibration frames with each science frame, you need to
include the ``-b`` option when running :ref:`pypeit_setup` to add the ``calib``,
``comb_id`` and ``bkg_id`` columns to the :ref:`data_block` of the
:ref:`pypeit_file`, although you only need the ``calib`` column and can leave
the other two at their default values (-1).

The ``calib`` column is used to pair science frames with calibration frames.
Here's an example:
    
.. code-block:: console

                     filename |                 frametype | ... | calib | comb_id | bkg_id
    DE.20170425.09554.fits.gz |                  arc,tilt | ... |   all |      -1 |     -1
    DE.20170425.09632.fits.gz | pixelflat,illumflat,trace | ... |   1,2 |      -1 |     -1
    DE.20170425.09722.fits.gz | pixelflat,illumflat,trace | ... |   1,2 |      -1 |     -1
    DE.20170425.09803.fits.gz | pixelflat,illumflat,trace | ... |     3 |      -1 |     -1
    DE.20170425.50487.fits.gz |                   science | ... |     1 |      -1 |     -1
    DE.20170425.51771.fits.gz |                   science | ... |     2 |      -1 |     -1
    DE.20170425.53065.fits.gz |                   science | ... |     3 |      -1 |     -1

Here, the ``arc,tilt`` frame is used for the calibration of every science frame,
so the value in its ``calib`` entry can be set to ``all`` (``1,2,3`` also
works). The first two ``illumflat,pixelflat,trace`` frames will be combined and
are used for the calibration of the first two science frames, while the third is
used for the calibration of the third science frame.

.. note::
    
    Importantly, because all of the ``comb_id`` and ``bkg_id`` values are -1,
    none of the science frames will be combined.  To perform a simple
    combination of science frames (i.e., combine the pixel values in the frames
    after :ref:`image_proc` **without** accounting for any shifts in the
    wavelength calibration or slit traces), see :ref:`2d_combine`.

Calibration Steps
=================

The primary calibration procedures are, in the order they're performed:

.. toctree::
   :maxdepth: 1

   ../image_proc
   bias_dark
   slit_tracing
   flexure
   wave_calib
   Slit Alignment (IFU only) <alignment>
   flat_fielding


Products
========

The main products of calibrations are :doc:`masters` which
are placed in the Masters/ folder.  Here are the full set
that may be created (not all are required; depends on the
instrument):

.. toctree::
   :maxdepth: 1

   masters
   master_align
   master_arc
   master_bias
   master_dark
   master_edges
   master_slits
   master_flat
   master_tilt
   master_tilts
   master_wvcalib

Modifications
=============

Here are the global modifications one may make
for calibrations:

* Add/Suppress bias/dark frame generation. See :doc:`bias_dark`
* Add/Suppress flexure correction.  See :doc:`flexure`
* Add/Suppress aspects of flat fielding.  See :doc:`flat_fielding`
* Associate different calibration frames to different science frames. See :ref:`2d_combine_calibs`

