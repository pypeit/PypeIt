.. highlight:: rest

*************
Master Frames
*************

By default, PYPIT will write a series of
calibration files per setup to the hard-drive in the
MastersFrames folder.  This is to allow the user
to inspect the calibations.  It is also to allow
a quicker re-reduction (i.e. these steps can be
skipped).

Here is a listing of the standard MasterFrame files:

=============== ======== ===========================================
Type            Format   Description
=============== ======== ===========================================
MasterArc       2d image Processed arc spectral image
MasterBadPix    2d image 2d image of bad pixels
MasterBias      2d image Processed bias image
MasterFlatField 2d image Normalized flat field image
MasterTilts     2d image Mapping of pixel to constant wavelength
MasterWaveCalib JSON     Solution of 1D wavelength calibration
MasterWave      2d image Wavelength image (in air and Angstroms)
=============== ======== ===========================================


**Warning:**  If the code has exited prematurely, one or
more of these frames can be essentially empty.

Reusing Masters
===============

There are 2 standard ways to direct the code to reuse any
existing MasterFrame files from the hard-drive when reducing.

By default, the code reuses any MasterFrames already in memory,
i.e. those produced during the course of the reductions.

Command Line
------------

When executing run_pypit, call with -m or --use_masters, e.g.::

    run_pypit pypit_file.pypit -m

PYPIT file
----------

Alternatively, you can add `reduce masters reuse True` to your
PYPIT file.

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
