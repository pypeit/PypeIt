*****************
Calibration Check
*****************

Overview
========

We *strongly recommend* that one perform a calibration
check with the .pypeit file before proceeding to run the
reduction.  This verifies that the number of desired
calibration files exist and allows the user to examine
how the code will group calibration files with science
frames.  It does **not** check the
sanctity of the files nor process the calibrations in any manner.

.. _run-calcheck:

calcheck
========

The procedure is simple.  Add the following line to your
.pypeit file::

    run calcheck True

You must also verify that your .pypeit file does **not**
include this line::

    run setup True   # Cannot be set for calcheck or full reduction

Either set 'run setup' to False, comment it out, or remove it altogether.

You may then run PypeIt, e.g.::

    run_pypeit kast_blue_setup_A.pypeit

The code will exit with error if there are insufficient calibration
frames.  Otherwise, it will exit after organizing the files and
will produce a new .group file for your inspection.

You should confirm that the correct number of science and
exposure standard files have been identified.


