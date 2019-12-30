.. highlight:: rest

************
VLT XShooter
************


Overview
========

This file summarizes several instrument specific
settings that are related to the VLT/XShooter spectrograph.


Calibrations
============

It is common for the VLT to observe calibration frames,
especially arc images, at the native (i.e. unbinned)
resolution.  The code will handle this with the exception
that you cannot bias subtract (only overscan).

To tell the code to skip bias subtraction, do::

    [calibrations]
        [[arcframe]]
            [[[process]]]
                bias = skip

This may become the default eventually.

