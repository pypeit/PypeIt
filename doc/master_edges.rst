.. _master-edges:

===========
MasterEdges
===========

Overview
========

This file describes the MasterEdges object.
It contains the core information required to describe the edges
of each slit on a given detector (or orders for echelle).

See below for the `Current EdgeTrace Data Model`_.
This is written to disk as a multi-extension FITS file prefixed by
MasterEdges in the Masters/ folder.
See the :ref:`master-naming` docs for more.

We also describe how to view the slit edges
using `pypeit_chk_edges`_.

Viewing
=======

The preferred way to view the slit edges identified
by PypeIt is with the `pypeit_chk_edges`_ script.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_chk_edges.rst

.. _pypeit_chk_edges:

pypeit_chk_edges
----------------

There are currently 2 options for viewing the slit edges on the image
used to construct them.  Each uses the `pypeit_chk_edges`_ script.

ginga
+++++

This is the default mode and it requires you to first launch
a `ginga` viewer in remote control mode::

    ginga --modules=RC

You should see `RC` in a tab on the upper right of the ginga viewer.

You may then run the script::

    pypeit_chk_edges Masters/MasterEdges_A_1_01.fits.gz


Here is an zoom-in screen shot for the :ref:`keck-lris-red` spectrograph.

.. image:: figures/slit_edges_ginga.png

A few things to note from this good-performing example:

 - The slits run nearly vertically across the image
 - The left/right edge of each slit identified is a green/red line
 - There were 13 slits identified (0 indexing)
 - The brightest `slit` is an alignment box and was discarded by the code

matplotlib
++++++++++

To avoid ginga, use the `--mpl` flag::

    pypeit_chk_edges Masters/MasterEdges_A_1_01.fits.gz --mpl

The color scheme is distinct and the labeling
now includes -1 or +1 for left/right edges.

SlitTrace
=========

TODO:
Describe that object here (or link to its own doc??)

.. _edges-trouble:

Edges Trouble Shooting
======================

See :doc:`slit_tracing` for a discussion of how to customize
and fuss with slit tracing.

Current EdgeTrace Data Model
============================

Internally, the image is held in EdgeTrace DataContainer.
The datamodel written to disk is:

TO APPEAR HERE

