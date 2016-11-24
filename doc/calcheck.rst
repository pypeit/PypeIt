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


