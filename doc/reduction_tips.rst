
.. TODO: Is it useful to consolidate these here instead of putting them in the
   relevant docs?

==============
Reduction Tips
==============

Overview
========

This document describes commonly adjusted
:doc:`pypeit_par` related to
object finding, sky-subtraction, and extraction of the spectra.

It assumes (but does not require):

1. You have already generated a set of :doc:`out_spec2D` files.

2. You have inspected the spectra wish to improve aspects of them.


.. _object_finding_tips:

Object Finding
==============

Refer to :doc:`object_finding` for full details on the algorithm.
This process is guided by the :ref:`findobjpar`.

The most common to modify is ``snr_thresh``, which sets the
search for any source with *peak* flux in excess of ``snr_thresh``
times the RMS.  The default is 10 and you may wish to
reduce this parameter.   Add the following to the
:ref:`parameter_block`:

.. code-block:: ini

    [reduce]
        [[findobj]]
            snr_thresh = 3.

This will search for any source with peak flux 3-sigma above the
estimated RMS in the smashed slit profile.

No 1D Spectra Extracted
-----------------------

If you are missing 1D spectra, this means that PypeIt did not find any objects
in the corresponding frame.  Here are some common modifications that you can
make to your :ref:`parameter_block` to remedy this.  It may help to
run the reduction on a few of your exposures in ``-s`` mode to make sure that your
modifications are being implemented the way you intend.

- Ensure that the slits are correctly identified.  See :ref:`slit_tracing` for
  tips on how to adjust your slit edges (for example, modifying ``edge_thresh``
  or ``minimum_slit_length``) and add or remove a slit.

- Modify the significance threshold for object finding (see above).

- Modify the maximum number of objects that you expect to see in your frames.
  The default is 10, but, for example, if your exposure only contains one object
  of interest, you will want to change your :ref:`parameter_block`
  to look like this:

  .. code-block:: ini
	
    [reduce]
        [[findobj]]
            maxnumber = 1

- Modify the FWHM of your object of interest in pixels.  The default is 10, but
  you may want to increase or decrease this, depending on the number of pixels per
  spatial PSF:

  .. code-block:: ini

    [reduce]
        [[findobj]]
            find_fwhm = 8.


Extraction
==========

Emission Lines
--------------

It is common for bright emission lines to spatially extend
beyond the source continuum, especially for galaxies.  In
these cases, the code may reject the emission lines because
they present a different spatial profile from the majority
of the flux.

While this is a desired behavior for optimal extraction of
the continuum, it leads to incorrect and non-optimal fluxes
for the emission lines.

The current mitigation is to allow the code to reject the
pixels for profile estimation but then to include them in
extraction.  This may mean the incurrence of cosmic rays
in the extraction.

Here is the expert move.  Add the following to the
:ref:`parameter_block`:

.. code-block:: ini

    [reduce]
        [[extraction]]
            use_2dmodel_mask = False

And it is likely you will want to use the ``BOXCAR`` extractions
instead of ``OPTIMAL``.  But do a comparison.


