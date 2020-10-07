=============
Master Frames
=============

Overview
========

As calibrations are generated,
by default PypeIt will write these data
to the hard-drive in the :doc:`masters` folder named Masters/.
This is to:

 - Allow the user to inspect the calibrations.
 - Allow for a quicker re-reduction (i.e. these steps can be skipped) including in cases where
   the code crashed, things were fixed and one re-runs the script.

Here is a listing of all the of MasterFrame files that
can be generated.  The ones that are made depend on the
:doc:`spectrographs` and the calibrations files provided.

================= ========= ===========================================
Type              Format    Description
================= ========= ===========================================
MasterAlign       FITS      Alignment image (IFU)
MasterArc         FITS      Processed arc spectral image
MasterBias        FITS      Processed bias image
MasterEdges       FITS      Several images describing the slit traces
MasterFlat        FITS      Normalized flat field image
MasterPinhole     FITS      Pinhole flat image
MasterSlits       FITS      Reduced output of slit tracing
MasterTiltimg     FITS      Mapping of pixel to constant wavelength
MasterTilts       FITS      Image used to trace wavelengths in the slits
MasterWave        FITS      Wavelength image (in air and Angstroms)
MasterWaveCalib   FITS      Solution of 1D wavelength calibration
================= ========= ===========================================


.. _master-naming:

Masters Naming
==============

The naming convention for MasterFrames is a bit obscure.
Here is an example of one::

    MasterEdges_A_1_01.fits

Here is how we break it down:

  - The prefix is as described in the Table above
  - The **A** specifies the spectrograph :doc:`setup`
  - The **1** specifies the detector number (one-based indexing)
  - The **01** specifies a bit-wise description of the `calib`
