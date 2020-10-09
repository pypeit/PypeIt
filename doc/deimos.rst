***********
Keck DEIMOS
***********

Overview
========

This file summarizes several instrument specific
settings that are related to the Keck/DEIMOS spectrograph.


Deviations
==========

Here are the main deviations from the default settings
for DEIMOS
(see :func:`~pypeit.spectrographs.keck_deimos.KeckDEIMOSSpectrograph.default_pypeit_par`).

These are tuned to the standard calibration
set taken with DEIMOS.

Calibrations
============

Edge Tracing
------------

It has been reported that the default `edge_thresh` of 50
for DEIMOS is too high for some setups.  If some of your
'fainter' slits on the blue side of the spectrum are missing,
try::

    [calibrations]
      [[slitedges]]
         edge_thresh = 10

It is possible, however, that our new implementation of using
the slitmask design file has alleviated this issue.

Flat Fielding
-------------

When using the *LVMslitC* mask, it is common for the
widest slits to have saturated flat fields.  If so, the
code will exit during flat fielding. You can skip over them
as described in :ref:`flat_fielding:Saturated Slits`.


Fluxing
-------

If you use the LVMslitC (common), avoid placing your standard
star in the right-most slit as you are likely to collide with
a bad column.

Additional Reading
==================

Here are additional docs related to Keck/DEIMOS:

.. toctree::
   :maxdepth: 1

   dev/deimosframes
