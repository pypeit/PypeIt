
.. include:: include/links.rst

.. _mosaic:

===============================================
Processing and Construction of Detector Mosaics
===============================================

If the geometry of each detector in a detector array is known *a priori*,
``PypeIt`` can construct a mosaic of the detector data for further processing,
instead of processing the detector data individually.  This is currently the
default approach for Gemini/GMOS only (soon we will add Keck/DEIMOS).

Coordinate Conventions
----------------------

Coordinate conventions are critical to understanding how ``PypeIt`` constructs
the detector mosaics.  These conventions are explained below.

The construction of the detector mosaic is done *after* the :ref:`image_proc`,
meaning that the images to be included in the mosaic obey the ``PypeIt``
convention of spectral pixels along the first axis and spatial pixels along the
second axis; i.e., the shape of each image is :math:`(N_{\rm spec}, N_{\rm
spat})`.  This is regardless of the 

*IN PROGRESS*



