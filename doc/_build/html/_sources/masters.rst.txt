.. highlight:: rest

*************
Master Frames
*************

By default, PypeIt will write a series of
calibration files per setup to the hard-drive in the
MastersFrames folder named MF_instrument (e.g.
MF_lris_blue).  This is to allow the user
to inspect the calibations.  It is also to allow
a quicker re-reduction (i.e. these steps can be
skipped) including in cases where the code crashed,
things were fixed and one re-runs the script.

Here is a listing of the typical set of MasterFrame files
(there are differences between AMRLSD and ARMED):

================= ========= ===========================================
Type              Format    Description
================= ========= ===========================================
MasterArc         2d image  Processed arc spectral image
MasterBadPix      2d image  2d image of bad pixels
MasterBias        2d image  Processed bias image
MasterFlatField   2d image  Normalized flat field image
MasterPinhole     2d image  Pinhole flat image
MasterSensFunc    YAML      Standard star sensitivity function
MasterSlitProfile 2d image  Pinhole flat image
MasterTilts       2d image  Mapping of pixel to constant wavelength
MasterTrace       2d images Several images describing the slit traces
MasterWave        2d image  Wavelength image (in air and Angstroms)
MasterWaveCalib   JSON      Solution of 1D wavelength calibration
================= ========= ===========================================


Reusing Masters
===============

There are 2 standard modes to direct the code to use any
existing MasterFrame files from the hard-drive when reducing.
These are described separately below (ReUse and Force).

Note that by default, the code reuses any MasterFrames already in memory,
i.e. those produced during the course of the reductions.

ReUse
+++++

The softer approach (recommended) is to ReUse any existing
MasterFrames in the folder and to generate new calibration
files when there is no MasterFrame.  Of course, this requires
the raw calibration files exist and have been properly
identified by the code.

Command Line
------------

When executing run_pypeit, call with -m or --use_masters, e.g.::

    run_pypeit pypeit_file.pypeit -m

PypeIt file
----------

Alternatively, you can add `reduce masters reuse True` to your
PypeIt file.

Force
+++++

There may be cases where you wish to use only MasterFrame files
that have been previously generated, e.g. on a previous night.
In this case, take the following steps:

1. Generate the properly named MF directory (e.g. MF_lris_red).
2. Copy all the necessary MasterFrame files into this directory
3. Add `reduce masters force True` to the PypeIt file
4. Add `reduce masters setup SETUP_NAME` to the PypeIt file, where
SETUP_NAME needs to match the setup used in the MasterFrames, e.g.
C_01_aa or A_02_ab.  Note that the detector number is not used
but should be 2 digits.
5. Run

You will note that this can only be performed on a single, specific
setup.  Also, the code will crash out if any required MasterFrame file is missing.
Finally, note that all calibration files are ignored (and QA files are
generally not recreated), including standard stars.
Only science frames are processed.


Internal
========

Here are some notes on the internal workings of MasterFrames.

Reduce dict
-----------

settings.argflag holds a dict of basic info on the MasterFrames.

====== ===== ============================================
Key    Type  Description
====== ===== ============================================
file   str   Points to the .setup file
loaded ??    ??
reuse  bool  Flag to specify whether to use MasterFrames
setup  str   Name of setup, e.g. '01'
====== ===== ============================================
