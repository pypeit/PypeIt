=========
Keck LRIS
=========


Overview
========

This file summarizes several instrument specific
settings that are related to the Keck/LRIS spectrograph.

Common Items
============

Flexure
+++++++

There is substantial flexure in the LRIS instrument and
the default settings attemps to characterize both the spectral
and spatial effects.

See the :doc:`flexure` notes if you wish
to turn either of these off.

.. _LRISb:

keck_lris_blue
==============

LRISb Default Settings
++++++++++++++++++++++

See :ref:`pypeit_par:KECK LRISb (``keck_lris_blue``)` for
a listing of modifications to the default settings.

Taking Calibrations for LRISb
+++++++++++++++++++++++++++++

Arcs
----

We recommend that you turn on *all* of the standard
arc lamps,  including those slated for the red side.

These are::

    Ne,Ar,Cd,Kr,Xe,Zn,Hg

The archived solutions expect all of these lamps.

Pixel Flat
----------

It is recommend to correct for pixel-to-pixel variations using a slitless
flat.  If you did not take such calibration frames or cannot process them,
you may wish to use an archival.  Request it from the developers.

And then set the following in your :doc:`pypeit_file`::

    [calibrations]
      [[flatfield]]
           frame = path_to_the_file/PYPEIT_LRISb_pixflat_B600_2x2_17sep2009.fits.gz

WARNING: Internal flats may be too bright and need to be tested.

Trace Flat
----------

We strongly recommend on-the-sky trace flats through full instrument
setup.  Aim for 1000 counts per pixel above bias.
These are best achieved by taking twilight flats within 15 minutes
of sunset/sunrise.

WARNING: Internal/dome flats are likely to be too faint in the
very blue.


.. _keck-lris-red:

keck_lris_red
=============

Taking Calibrations for LRISr
=============================

LRISr Default Settings
++++++++++++++++++++++

Here are the deviations from the default settings
for LRISr::

    settings trace slits sigdetect 50.0   # Good for relatively bright dome flats
    settings trace slits pca params [3,2,1,0]

Known issues
++++++++++++

Multi-slit
----------

The code may identify a 'ghost' slit in empty detector real
estate if your mask does not fill most of the field.  Be prepared
to ignore it.
