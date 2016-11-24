.. highlight:: rest

*****************
Calibration Check
*****************

Overview
========

We *strongly recommend* that one perform a calibration
check with the .pypit file before proceeding to run the
reduction.  This simply verfies that the number of desired
calibration files exist.  It does **not** check the
sanctity of the files nor process the calibrations in any manner.

Procedure
=========

The procedure is simple.  Add the following line to your
.pypit file::

    run calcheck True

You must also verify that your .pypit file does **not**
include this line::

    run setup True   # Cannot be set for calcheck or reduction

Either set to False, comment it out, or remove it altogether.

You may then run PYPIT, e.g.::

    run_pypit kast_blue_setup_01_02.pypit

The code will exist with error if there are insufficient calibration
frames.  Otherwise, it will exit after organizing the files and
will produce a new .group file for your inspection.


Settings
========

The primary settings you need to specify at this stage are:

#.  The number of calibration files required of each frametype

#.  Over-ride any frametype designations, as necessary.

For the second issue, see :ref:`modifying_frametype`.

For the first, add a series of lines (or edit the existing ones)
in the :ref:`spect_block` of the .pypit file.
One line per calibration frametype desired.
Here is a standard block for LRISb::

     pixelflat number 5
     arc number 1
     slitflat number 5
     bias number 10
     standard number 1

When a positive, non-zero value is used, the code will require
that there be that many calibration frames for each science
frame reduced.  And, PYPIT will restrict to precisely that many
calibration files.

If you wish to you *at least* that many frames, specify with a
negative integer value, e.g.::

     pixelflat number 5
     arc number 1
     slitflat number -5
     bias number -5
     standard number -1


