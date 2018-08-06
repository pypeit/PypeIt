.. highlight:: rest

====================
Calibration Overview
====================

This document gives an overview of the calibration
steps of PypeIt.  There is additional documentation
for several of the more complex steps.

Sequence of Events
==================

Here is the sequence of events:

========== ============= ===========================================
Step       Products      Description
========== ============= ===========================================
datasec    datasec image 2D image describing the detector pixels for analysis
========== ============= ===========================================

datasec
=======

In this step, PypeIt parses the instrument settings file
(or user input) to establish the region on each detector
for analysis.  The overscan section is also established
and is included in the datasec if one exists.

A 2D image defining the datasec pixels
is generated and stored internally (in _datasec).

Standard
--------

The standard approach to defining the datasec is to set these
in the instrument settings file.  It is necessary to generate
one set per amplifier as each of these may have distinct
properties (e.g. readnoise, gain).

Here are the values for Kast blue::

    det01 numamplifiers 2                 # Number of amplifiers
    det01 datasec01 [:,0:1024]
    det01 oscansec01 [:,2049:2080]
    det01 datasec02 [:,1024:2048]
    det01 oscansec02 [:,2080:2111]
    det01 gain 1.2,1.2                    # Inverse gain (e-/ADU)
    det01 ronoise 3.7,3.7                 # Read-out noise (e-)

LRIS
----

The FITS file writing for LRIS is sufficiently complex that the
datasec definition (and file loading)
is guided by a custom method: arlris.read_lris()
