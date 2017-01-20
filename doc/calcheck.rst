.. highlight:: rest

*****************
Calibration Check
*****************

Overview
========

We *strongly recommend* that one perform a calibration
check with the .pypit file before proceeding to run the
reduction.  This verifies that the number of desired
calibration files exist and allows the user to examine
how the code will group calibration files with science
frames.  It does **not** check the
sanctity of the files nor process the calibrations in any manner.

Procedure
=========

The procedure is simple.  Add the following line to your
.pypit file::

    run calcheck True

You must also verify that your .pypit file does **not**
include this line::

    run setup True   # Cannot be set for calcheck or full reduction

Either set 'run setup' to False, comment it out, or remove it altogether.

You may then run PYPIT, e.g.::

    run_pypit kast_blue_setup_A.pypit

The code will exit with error if there are insufficient calibration
frames.  Otherwise, it will exit after organizing the files and
will produce a new .group file for your inspection.

You should confirm that the correct number of science and
exposure standard files have been identified.


Settings
========

PYPIT identifies calibration files that are closest in time to every
individual science frame.
You can place an upper limit on the time window that PYPIT uses to search
for calibrations but setting the keyword::

     fits calwin 12.0

which will search for calibrations that were taken within +/-12 hours
from a science frame.  See docs on :ref:`calwin` for a further
discussion.

The primary settings you need to specify at this stage are:

#.  The number of calibration files required of each frametype

#.  Over-ride any frametype designations, as necessary.

#.  Modify the method(s) for bias subtraction, flat fielding etc.

For the second issue, see :ref:`modifying_frametype`.

For the first issue see below:

.. _calwin:

calwin
------

When associating calibration files to a given science frame,
PYPIT will restrict to data within a window in time.  This
is specified by a `calwin` parameter which has a default
value of 12 (hours) for most instruments.  One can turn
off this restriction by setting the value to 0 in
the :ref:`spect_block`::

    fits calwin 0

This is the default for :ref:`LRISb` and
may become the default for all instruments.

Calib number
------------

The user can specify and/or over-ride defaults
for the number of calibration frames required
by adding a series of lines (or edit the existing ones)
in the :ref:`spect_block` of the .pypit file.
One line per calibration frametype, as desired.
Here is a block one might use for LRISb::

    # Spect
    spect read
     arc number 1
     trace number 5
     bias number 10
     standard number 1
     pixelflat number 3
    spect end

When a positive, non-zero value is used, the code will require
that there be that many calibration frames for each science
frame reduced.  And, PYPIT will restrict to precisely that many
calibration files.

If you wish to use *at least* an input number of frames (and
more if they exist), then specify the calibration number
with a negative integer value, e.g.::

     pixelflat number 5
     arc number 1
     trace number -5
     bias number -5
     standard number -1


