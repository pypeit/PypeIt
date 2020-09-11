=================
Manual Extraction
=================

Overview
========

This document describes how to perform so-called Manual
Extraction in PypeIt.  This is generally for cases where the
object continuum is too faint to trigger the auto-magical
:doc:`object_finding` algorithm.

Process
=======

Here is the standard recipe:

1. Reduce the spectral image(s)
2. Examine the spec2d images with :ref:`pypeit_show_2dspec`
3. Record the spatial and spectral pixel where the trace should cross
4. Modify the PypeIt file as described below
5. :ref:`run-pypeit` again

Tracing
-------

The code will lay down a new trace and perform extraction
at each input location.  The trace used will be, in order
of decreasing preference:

1. The brightest object on the slit, offset to the input position
2. The standard star
3. The slit edges

Multi-Slit
----------

If you are running in multi-slit mode, you will add the spatial-spectral
pixel pair for each object to extract for each detector to the PypeIt file.

Here is an example::

    [reduce]
      [[extraction]]
        [[[manual]]]
           spat_spec = 212.3:1200,742.3:1200
           det = 1,2
           fwhm = 3.,3.

The above will lay down a new trace at spatial=212.3, spectral=1200
pixel on detector 1 and use a FWHM of 3 pixels.


Echelle
-------

For echelle, you only have to specify the object location in a single
order and the code will use its fractional position on all other orders.

Here is an example from the DevSuite (VLT_manual)::

    [reduce
      [[extraction]]
        [[[manual]]]
           spat_spec = 1181.8:3820.6
           det = 1
           fwhm = 3.

