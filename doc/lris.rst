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

See :ref:`pypeit_par:Instrument-Specific Default Configuration` for
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
you may wish to use an archival.
`This link <https://drive.google.com/drive/folders/1YmDgCgXrsRbkuH_Pc_MLShWVdSrMkoFP?usp=sharing>`_
has the existing ones staged by the PypeIt team.

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

LRISb Slit Edges
----------------

When observing in long-slit mode, PypeIt might set the slit incorrectly
for detector 2.  This may occur if the counts from the flat field
are too low (e.g., using internal flats rather than twilight
flats with a higher signal in the blue).
Therefore, if you use internal flats, be careful to inspect the
slits defined by PypeIt as described in :doc:`master_edges`.

If the defined slit(s) does not cover the portion of
the illuminated detector where your source falls, you
can manually define the slit position as described
in :ref:`slit_tracing:Missing A Slit`.


Here is an example for the PypeIt file::

    [calibrations]
       [[slitedges]]
         add_slits = 2:788:10:650
         sync_predict = nearest

This will force a slit onto the detector for reduction.

Multi-slit
----------

The code may identify a 'ghost' slit in empty detector real
estate if your mask does not fill most of the field.  Be prepared
to ignore it.
