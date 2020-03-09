.. _master-edges:

===========
MasterEdges
===========

Overview
========

This file describes the data model for the MasterEdges object.
This is written to disk as a multi-extension FITS file prefixed by
MasterEdges in the Masters/ folder.

It contains the core information required to describe the edges
of each slit on a given detector (or orders for echelle).

.. _pypeit-chk-edges:

pypeit_chk_edges
================

There are currently 2 options for viewing the slit edges on the image
used to construct them.  Each uses the :ref:`pypeit-chk-edges`
script.  Use `pypeit_chk_edges -h` to view its current usage.

ginga
-----

This is the default mode and it requires you to first launch
a `ginga` viewer in remote control mode::

    ginga --modules=RC

You should see `RC` in a tab on the upper right of the ginga viewer.

You may then run the script::

    pypeit_chk_edges Masters/MasterEdges_A_1_01.fits.gz


Here is an zoom-in screen shot for the :ref:`keck-lris-red` spectrograph.

.. image:: figures/slit_edges_ginga.png

A few things to note from this good-performing example:

 - The slits run vertically on the image
 - The left/right edge of each slit identified is a green/red line
 - There were 13 slits identified (0 indexing)
 - The brightest `slit` is an alignment box and was discarded by the code

matplotlib
----------

To avoid ginga, use the `--mpl` flag::

    pypeit_chk_edges Masters/MasterEdges_A_1_01.fits.gz --mpl

The color scheme is distinct and the labeling
now includes -1 or +1 for left/right edges.

SlitTrace
=========

TODO:
Describe that object here (or link to its own doc??)

Current EdgeTrace Data Model
============================

Internally, the image is held in EdgeTrace DataContainer.
The datamodel written to disk is:

TO APPEAR HERE

