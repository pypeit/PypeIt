
.. include:: ../include/links.rst

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

To assign specific calibration frames to each science frame, you must edit the
``calib`` column in the :ref:`data_block` of the :ref:`pypeit_file`, which is 
used to pair science frames with calibration frames.  Here's an example:
    
.. code-block:: console

                     filename |                 frametype | ... | calib |
    DE.20170425.09554.fits.gz |                  arc,tilt | ... |   all |
    DE.20170425.09632.fits.gz | pixelflat,illumflat,trace | ... |   1,2 |
    DE.20170425.09722.fits.gz | pixelflat,illumflat,trace | ... |   1,2 |
    DE.20170425.09803.fits.gz | pixelflat,illumflat,trace | ... |     3 |
    DE.20170425.50487.fits.gz |                   science | ... |     1 |
    DE.20170425.51771.fits.gz |                   science | ... |     2 |
    DE.20170425.53065.fits.gz |                   science | ... |     3 |

Here, the ``arc,tilt`` frame is used for the calibration of every science frame,
so the value in its ``calib`` entry can be set to ``all`` (``1,2,3`` also
works). The first two ``illumflat,pixelflat,trace`` frames will be combined and
are used for the calibration of the first two science frames, while the third is
used for the calibration of the third science frame.

.. note::

    - Calibration group numbers **must** be integers and unique for each group;
      however, they do not need to follow any particular sequence.  

    - There is currently a limit of no more than 63 calibration groups.

    - Calibration frames can be assigned to multiple calibration groups, but
      science frames can only be assigned to *one* calibration group.

.. _calibrations-calibfile:

".calib" File
-------------

A PypeIt ``.calib`` file is written both when running :ref:`pypeit_setup` and
:ref:`run-pypeit`, a yaml file that provides a summary of the calibration groups
used during the reduction.  The file lists all the source file of each frame
type ordered by their calibration group and looks roughly like this:

.. code-block:: yaml

    # Auto-generated calibration association file using PypeIt version:  1.12.2.dev137+g95b15695a
    # UTC 2023-03-28T17:09:20.359
    # NOTE: DET01 is a placeholder for the reduced detectors/mosaics
    A:
      --:
        dichroic: d55
        dispname: 600/4310
      0:
        arc:
          proc:
          - /rdx/path/shane_kast_blue_A/Calibrations/Arc_A_0_DET01.fits
          - /rdx/path/shane_kast_blue_A/Calibrations/WaveCalib_A_0_DET01.fits
          raw:
          - /raw/path/b1.fits.gz
        bias:
          proc:
          - /rdx/path/shane_kast_blue_A/Calibrations/Bias_A_0_DET01.fits
          raw:
          - /raw/path/b14.fits.gz
          ...
        science:
        - /raw/path/b27.fits.gz
        - /raw/path/b28.fits.gz
        standard:
        - /raw/path/b24.fits.gz
        ...

The outermost element gives the setup/configuration identifier, ``A`` in this
case.  This is immediately followed by the instrument configuration details.
Then, for each calibration group, the file provides the associated raw science
and standard files, the raw calibration files, and the processed calibration
files the code *expects* to create.  Whether or not the processed calibration
files exist after executing :ref:`run-pypeit` will depend on the success of the
run and any relevant user-based parameters.

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
   scattlight

Follow the above links for a detailed description of each step, and how the
default calibration procedures may be **modified** to suit your needs.

----

Products
========

As the calibrations are completed, PypeIt will save the results to files in the
``Calibrations/`` folder (see :ref:`outputs-dir`).  Unless otherwise indicated, all
output files are in FITS format.

Saving the results of each calibration step to a file allows:

 - the user to inspect the calibrations, and

 - for a quicker re-reduction (i.e. these steps can be skipped) including in
   cases where the code crashed, things were fixed, and then one re-runs the
   script.

Below is the full list of possible calibration frames produced by PypeIt.  For any
given run, the files actually produced will depend on the :doc:`spectrograph
<../spectrographs/spectrographs>` and the files listed in the
:ref:`pypeit_file`.

.. toctree::
   :maxdepth: 1

   Alignment Image (Alignment; IFU only) <align>
   Processed arc spectral image (Arc) <arc>
   Processed bias image (Bias) <bias>
   Processed dark image (Dark) <dark>
   Images used to trace slit edges (Edges) <edges>
   Consolidated slit traces (Slits) <slits>
   Normalized flat field images (Flat) <flat>
   Image used to trace wavelengths within each slit (Tiltimg) <tilt>
   Mapping of pixels to constant wavelength (Tilts) <tilts>
   Solution of 1D wavelength calibration (WaveCalib) <wvcalib>

.. _calib-naming:

Calibration Frame Naming
------------------------

All reduced calibration frames are named according to their primary calibration
type (e.g., ``Arc``).  They are also assigned a unique identifier that is a
combination of:

    #. the instrument configuration (setup) identifier (e.g., ``A``),

    #. a compressed list of associated calibration groups (e.g., ``1+2`` or ``all``), and

    #. the detector or mosaic identifier (e.g., ``DET01`` or ``MSC01``).

For the second component, sequential numbers are reduced to a range; e.g.,
``'0-1-2-3-4'`` becomes ``'0+4'`` and ``'3-5-6-10-11-12-15-18-19'`` becomes
``'3-5+6-10+12-15-18+19'``.

.. warning::

    If you have a lot of calibration groups in your pypeit file, you may end up
    with very long file names!  This may cause a fault when the file name is
    included in the header of the output fits files.  If using the calibration
    group ``all`` doesn't solve the problem or isn't possible given your
    application, please `Submit an issue`_.


