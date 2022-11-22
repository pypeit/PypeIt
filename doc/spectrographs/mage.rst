.. highlight:: rest

*************
Magellan Mage
*************


Overview
========

This file summarizes several instrument specific
settings that are related to Magellan/Mage.


Short slits
===========

There are several issues related to the very short
slits of Magellan/Mage  (34 pixels or 10" unbinned).

Find Objects
------------

To have enough slit to 'properly' find objects,
we restrict the ``find_trim_edge`` parameter, i.e.:

.. code-block:: ini

    [scienceimage]
        find_trim_edge = 4, 4    # Slit is too short to trim 5,5 especially with 2x binning

For spatial binning, we recommend you to further reduce
this by the binning factor.

