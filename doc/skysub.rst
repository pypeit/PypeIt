
.. include:: include/links.rst

.. _skysub:

===============
Sky Subtraction
===============

Overview
========

This document describes how PypeIt performs sky subtraction.

See :ref:`skysubpar` for the complete
list of options related to sky subtraction.

.. _skysub-global:

Global
======

Phase I of sky subtraction is to perform a fit to the sky across
the entire slit.  By default, this is done twice:  once without
any knowledge of objects in the slit and then again after object
detection has taken place (these are masked).

Default masking of objects is relatively benign.  The FWHM
of each object is estimated and then those pixels above a
set threshold in the profile are masked.

One can enforce more aggressive masking by
setting ``mask_by_boxcar`` which will mask each object by the
``boxcar_radius`` set in :ref:`extractionpar`:

.. code-block:: ini

    [reduce]
       [[extraction]]
          boxcar_radius = 2.5  # arcsec
       [[skysub]]
          mask_by_boxcar = True


Local
=====

Assuming you perform :ref:`extraction-optimal` extraction,
the default is to refine the sky subtraction in tandem.

To turn this off (e.g. recommended for bright extended emission
lines on faint galaxy continua), set
``no_local_sky`` in :ref:`skysubpar`:

.. code-block:: ini

    [reduce]
       [[skysub]]
          no_local_sky = True

.. _skysub-regions:

Interactively defining the sky regions
======================================

PypeIt has an automatic algorithm (described above) to define
the sky regions, but this may not work in
your specific science case. There are several ways to define
the sky regions. The first option is to define the locations
on the slits where there is sky in your :ref:`pypeit_file`. The
command is a comma separated list of regions that represent
the locations on the slit (0 is the left edge, 100 is the
right edge):

.. code-block:: ini

    [reduce]
      [[skysub]]
           user_regions = :20,65:

where in the example above, the sky regions are defined as all
pixels in all slices that are in the leftmost 20 percent of the
slit (i.e. :20), and the rightmost 35 percent of the slit (65:).
You can specify as many regions as you like. For example, 45:55
would indicate that the innermost 10 percent of pixels contains
sky.

An alternative approach is to set the sky regions interactively.
This is the preferred approach if you want to set different sky
regions for every slit. Remember, you really should assign some
sky regions in every slit, otherwise the relative spectral
sensitivity correction will not work. To interactively define
the sky regions, you must first run through the reduction once,
and then use the following command:

.. code-block:: console

    pypeit_skysub_regions spec2d_file.fits

You will need to manually define the sky regions for each spec2d file.
You will see a GUI where you can click and drag regions on each
slit to define the sky regions. Hover the mouse over the window
and press the ``?`` key. This will print a list of options in the
terminal window, so that you know how to operate the GUI. A left
(right) mouse button click and drag will add (remove) pixels to
(from) the sky regions mask. Once you have defined some regions,
the red shaded regions represent the sky pixels. If you want to
set the sky regions for multiple slits, use the
"Assign sky regions to all slits"
bar on the right hand side of the GUI. The gray region represents
the slit, and the black regions represent outside the slit. You
need to click and drag only on the gray regions, or you can click
and drag from the gray to the black regions (i.e. you must click
and drag within this small window for it to work).

Alternatively, you can click the "Enter regions" button, which
will request input from the command line. You should now enter
the regions in the same format as above for the ``user_regions``.

If you're happy with the sky regions, press the
"Continue (and save changes)" button. If you do not wish to save
the sky regions, press the "Continue (and don't save changes)" button.
The menu bar at the top of the screen will prompt you if you
wish to save these sky regions (click on either YES or NO).
If you chose to save the regions file, the regions will be
saved in your ``Calibrations/`` folder, with a prefix ``SkyRegions``.
A given ``SkyRegions`` file is linked to a science frame
based on the name of the ``SkyRegions`` file.

Once you have defined all of the sky regions manually, you will need to explicitly
tell PypeIt to use the manually defined sky regions file by adding the following
lines to your :ref:`pypeit_file`:

.. code-block:: ini

    [reduce]
       [[skysub]]
          user_regions = user

and then re-run the reduction.

