
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
Instead, it is treated as a post-processing procedure; see :doc:`../fluxing`.

----

.. _calibration-groups:

Calibration Groups
==================

After :ref:`pypeit_setup` has organized your files into distinct instrument
configurations, the simplest approach (and PypeIt's default approach) is to use
all calibrations for all science frames; however, that isn't always the best
approach.  For example, you may want to associate different science frames with
calibrations taken either in the evening or the morning.  To ensure the best
calibrations are used with each science frame, you could split the data between
two :ref:`PypeIt Files<pypeit_file>` or you can use PypeIt's *calibration group*
designations.

To assign specific calibration frames to each science frame, you need to
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

.. _calibrations-calibfile:

Output calib File
-----------------

When PypeIt runs successful, a file with a ``.calib`` extension is generated
that provides a summary of the calibration groups used during the reduction.
The file lists all the source file of each frame type ordered by their
calibration group.  It looks roughly like this:

.. code-block:: console

    A:
    --:
        SETUP DETAILS
    1:
        align: []
        arc:
        - /arc/file/1
        - /arc/file/2
        bias: []
        ...

You may generate a similar (and perhaps more readable) file using the
:ref:`pypeit-parse-calib-id` script.

----

Calibration Steps
=================

The primary calibration procedures are, in the order they're performed:

.. toctree::
   :maxdepth: 1

   image_proc
   slit_tracing
   flexure
   wave_calib
   Slit Alignment (IFU only) <alignment>
   flat_fielding

Follow the above links for a detailed description of each step, and how the
default calibration procedures may be **modified** to suit your needs.

----

Products
========

As the calibrations are completed, PypeIt will save the results to files in the
``Masters/`` folder (see :ref:`outputs-dir`).  Unless otherwise indicated, all
output files are in FITS format.

Saving the results of each calibration step to a file allows:

 - the user to inspect the calibrations, and

 - for a quicker re-reduction (i.e. these steps can be skipped) including in
   cases where the code crashed, things were fixed, and then one re-runs the
   script.

Below is the full list of possible master frame produced by PypeIt.  For any
given run, the files actually produced will depend on the :doc:`spectrograph
<../spectrographs/spectrographs>` and the files listed in the
:ref:`pypeit_file`.

.. toctree::
   :maxdepth: 1

   Alignment Image (MasterAlignment; IFU only) <master_align>
   Processed arc spectral image (MasterArc) <master_arc>
   Processed bias image (MasterBias) <master_bias>
   Processed dark image (MasterDark) <master_dark>
   Images used to trace slit edges (MasterEdges) <master_edges>
   Consolidated slit traces (MasterSlits) <master_slits>
   Normalized flat field images (MasterFlat) <master_flat>
   Image used to trace wavelengths within each slit (MasterTiltimg) <master_tilt>
   Mapping of pixels to constant wavelength (MasterTilts) <master_tilts>
   Solution of 1D wavelength calibration (MasterWaveCalib) <master_wvcalib>

.. _master-naming:

Master Frame Naming
-------------------

All reduced calibration frames are given the file name prefix ``Master`` (e.g.,
``MasterArc``).  They are also assigned a unique identifier that is a combination of:

    - the instrument configuration (setup) identifier (e.g., ``A``),

    - the calibration group **bit** identifier (e.g., ``7``), and

    - the detector or mosaic identifier (e.g., ``DET01`` or ``MSC01``).

Importantly, note that the calibration group **bit** identifier is used and not
a specific calibration group ID number.  This is because master frames can
belong to multiple calibration groups.  For example, if there are 3 calibration
groups (groups 1, 2, and 3) and a set of bias frames are used for all three
groups, the resulting master bias frame will be called
``MasterBias_A_7_DET01.fits`` for detector 1 in setup A.  The ``7`` comes from
the bitwise calculation :math:`2^{1-1} + 2^{2-1} + 2^{3-1} = 7`.

Although convenient from a development perspective, this abstraction of the
master file names can make it difficult to associate each master frame with its
source input frames.  To help with this, see the :ref:`calibrations-calibfile`.


