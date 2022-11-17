
.. include:: ../include/links.rst

.. _gnirs_howto:

================
Keck/NIRES HOWTO
================

Overview
========

This doc goes through a full run of PypeIt on one of the Keck/NIRES datasets
(specifically the ``NIRES`` dataset) in the `PypeIt Development Suite`_ (see
:ref:`dev-suite`).

PypeIt File
===========

To setup the pypeit file, first run :ref:`pypeit_setup`:

.. code-block:: bash

    pypeit_setup -r $HOME/Work/packages/PypeIt-development-suite/RAW_DATA/keck_nires/NIRES/ -s keck_nires -b -c A 

where ``-b`` indicates that the data uses background images and includes the
``calib``, ``comb_id``, ``bkg_id`` in the pypeit file.  This directory only has
one instrument configuration (NIRES only has one configuration anyway), so
setting ``-c A`` and ``-c all`` is equivalent.

The resulting pypeit file looks like this:

.. include:: ../include/keck_nires_A.pypeit.rst

For this example dataset, the details of the default pypeit file are not
correct; see :ref:`here <nires_config_report>` for an example where the frame
typing and dither pattern determinations are correct.

The issue with this database is that only two pointings in the ABBA pattern for
the standard star observations are available, the exposure time of the standard
observations is longer than the default, and the dither pattern of the science
frames are set to MANUAL.  All of this means you need to edit the pypeit file to
correct these errors.  The corrections needed are:

    - Set the frametypes for ``s190519_0059.fits`` and ``s190519_0060.fits`` to
      ``arc,tilt,standard``.

    - Set the ``bkg_id`` for the ``standard`` and ``science`` frames; see
      :ref:`a-b_differencing`.

The corrected version looks like this:

.. include:: ../include/keck_nires_A_corrected.pypeit.rst

Finally, note that by setting the two science frames to also be ``arc`` and
``tilt`` frames, we're indicating that the sky lines should be used to perform
the wavelength calibration.


Reliable image typing and sequence generation based on header cards is not yet implemented for GNIRS.
Hence, several modifications to the PypeIt file need to be made before executing :ref:`run-pypeit`.
Wavelength solutions can be robustly determined from OH sky lines, which due to 
flexure, is preferable to using the Ar lamp arcs. 
Typically the telluric star chosen is bright and will have high signal in a 
short amount of time resulting in weak sky lines relative to the telluric star 
itself. For this reason we cannot determine the tilt and and wavelength solutions from
the telluric standard (t_exp = 10s), thus, arc and tilt should be removed from the telluric
star's frame types.

Instead we will set its calibration ``calib`` ID to match that of the 
nearest (in time) science image, which instructs pypeit to use the tilt and 
wavelength solution from that sequence for the telluric standard.
**Note that even if you don't plan to telluric correct your data, it is advantageous
to always nevertheless include the telluric star in your PypeIt file.**
The reason is that PypeIt uses the telluric star as a crutch for object tracing, 
and if you have faint objects this will likely produce better results than if there 
is no telluric star in the PypeIt file, in which case PypeIt uses the slit/order boundaries as the
tracing crutch.

The ``calib`` ID should be set to the same number for the frames that will be
combined together to generate the calibration image in question. The best
strategy for choosing these frames and setting the associated ``calib`` IDs
depends on the exposure time, dither pattern, and the instrument in question.
The example PypeIt file block below is for an ABBA sequence with GNIRS. To
better understand how to set the ``calib`` IDs, review how to set the
``comb_id`` and ``bkg_id`` for an ABBA sequence, as described
:doc:`here<../A-B_differencing>`.

In the edited pypeit file below, we instruct PypeIt to average two A images
together as the science (e.g., frames ``*216*`` and ``*219*``) and subtract from
this the average of two B images (e.g., frames ``*217*`` and ``*218*``), which
will generate one set of spec2d and spec1d outputs. We do the same thing for the
B images, i.e. combine them and subtract from them the combine of the A images.

This entire ABBA sequence has the same OH lines. If we really want a distinct
wavelength and tilt solution to be generated from the average of the As and
another one to be generated from the average of the Bs, then we would set the
``calib`` IDs to be 0110, where the 0 and 1 are arbitrary numbers (i.e. it could
also be 4554). However, if the instrument is not expected to flex much in an
ABBA sequence, it is actually advantageous to combine the entire ABBA sequence
into one Master frame for wavelength calibration and tilt generation. The reason
for this is that by averaging four images, the flux from the science object gets
diluted. This is desirable for OH wavelength and tilt calibrations because the
object counts, particularly if the object is bright, is actually a contaminant.
In other words, we extract a 1d OH sky spectrum for the wavelength calibration
and trace the arc lines trajectories across the detector for the tilts.

Obviously a bright object can mess this up.  (For example in the optical you
would not turn the arc lamps on and take arcs while simultaneously observing a
star through the slit). In the IR often the sky is so bright relative to the
objects and contaminates so few spatial pixels that this is not much of a worry,
but it still good practice to average away the object flux by combining the
entire ABBA sequence into set of calibration frames. For this reason, we set the
``calib`` IDs to be, e.g., 0 for all the images in a given ABBA sequence.  The
``calib`` IDs in the next ABBA sequence is 1, etc.

The edited pypeit file (exactly the one to reduce our example Gemini/GNIRS data
set in the `PypeIt Development Suite`_) is:

.. include:: ../include/gemini_gnirs_A_corrected.pypeit.rst

Note that the telluric standard has its ``calib`` IDs set to all 0s, which
corresponds to the ``calib`` ID of the nearest science ABBA sequence in time.

Try pulling the data from the `PypeIt Development Suite`_ (see
:ref:`devsuite-raw-data`) and reducing it!

