.. include:: ../include/links.rst

.. _magellan_ldss3:

**************
Magellan LDSS3
**************

Overview
========

This file summarizes several instrument specific
items for the Magellan/LDSS3 spectrograph. 

The pipeline now only supports the VPH-ALL grism, 
and is only tested for cases when the science object is on chip 1 (c1) and 1x1 binning.

Prepare the .pypeit file
========================

To limit spurious object detections, we recommend restricting the object
finding in the ``.pypeit`` file, e.g.::

[reduce]
    [[findobj]]
        find_trim_edge = 100, 100  # ignore 100 pix at each slit edge

Sensitivity function
====================

The default sensitivity-function algorithm is ``IR``, which has been tested on
LDSS3 data and is robust.  ``UVIS`` is also an option: it is faster and gives
nearly identical results in our tests, but is less thoroughly validated for
this instrument setup.


