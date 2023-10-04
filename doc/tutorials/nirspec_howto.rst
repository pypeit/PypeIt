
.. include:: ../include/links.rst

.. _nirspec_high_howto:

================
Keck/NIRSPEC High Resolution HOWTO
================

Overview
========

This tutorial is written to serve as both and guide and documentation for reducing 
Keck/NIRSPEC high resolution spectra using PypeIt. We will begin with a reduction of the sample
J Band dataset and include recommendations for reducing other sample datasets (Y, H, and K bands).

See :ref:`here <dev-suite>` to find the example dataset, please join our `PypeIt Users Slack <pypeit-users.slack.com>`__ (using
`this invitation link
<https://join.slack.com/t/pypeit-users/shared_invite/zt-1kc4rxhsj-vKU1JnUA~8PZE~tPlu~aTg>`__)
to ask for help, and/or `Submit an issue`_ to Github if you find a bug!


----

Directory Organization
=====

Before starting the reduction, we recommend creating a directory in which to store all of the PypeIt output files. PypeIt will 
create its own directory structure within that, but it is good to have a dedicated place to work. PypeIt will also reuse directory names
so creating a dedicated directory for a given reduction run is recommended. 


----

Setup
=====

PypeIt uses a special input file, the :ref:`pypeit_file`, to organize any user supplied keywords, the locations of the input data files, and 
the metadata corresponding to each file. This file is created and automatically populated by running :ref:`pypeit_setup`:

.. code-block:: bash

    pypeit_setup -r $HOME/Work/packages/PypeIt-development-suite/RAW_DATA/keck_nirspec_high/Jband/ -s keck_nirspec_high -b -c all

where the ``-r`` argument should be replaced by your local directory and the
``-b`` indicates that the data uses background images and should include the
``calib``, ``comb_id``, ``bkg_id`` in the pypeit file.  In the event that you are 
reducing data with multiple configurations, using ``-c all`` will create an output directory for each configuration. In this case,
since the dataset only has one configuration, using ``-c A`` would be equivalent. 

This will make a directory called ``keck_nirspec_high_A`` that holds a pypeit file
called ``keck_nirspec_high_A.pypeit`` that looks like this:

.. include:: ../include/keck_nires_A.pypeit.rst


At the moment, NIRSPEC does not keep track of the nod pattern that was used in observing and where in the nod pattern the exposure 
was taken. This means that setting AB pairs in the data for background subtraction must be done manually. 

The corrections we'll need to make are:
    - Set the frametypes for ``nspec***.fits`` to only ``science`` (i.e. remove both the ``arcs`` and ``tilts``. Since we already have 
    arclamp spectra for this dataset, we want to use those, rather than the sky lines. If you would like to use the sky lines, leave the 
    ``arc`` and ``tilts`` in and remove the lamp spectrum from the file list. We will discuss other changes you'll need to make to the
    parameter file below, in the Wavelength Calibration section.

    - Assign the appropriate A/B nod combinations by setting the ``bkg_id`` of each A exposure to the ``comb_id`` of its B counterpart 
    and vice versa. See :ref:`a-b_differencing` for other options for differencing and background subtraction.


The corrected version looks like this (pulled directly from the :ref:`dev-suite`):

.. include:: ../include/keck_nires_A_corrected.pypeit.rst



The importance of standard frames
---------------------------------

We recommend the observer always take standard frames during the night, ideally 
frequently enough to sample the airmass throughout the night, which will make the
standards more useful for telluric removal later on. Note that the main data-reduction
script (:ref:`run-pypeit`) does *not* perform the telluric correction. For that, users
will use the :ref:`pypeit_tellfit` script or perform the removal using their standards,
following the procedure discussed below. 

However, even if you don't intend to telluric-correct or flux-calibrate your data, it's useful to include
the standard star observations along with the reductions of your main science
target, particularly if the science target is faint.  If your object is faint,
tracing the object spectrum for extraction can be difficult using only the
signal from the source itself.  PypeIt will resort to a "tracing crutch" if the
source signal becomes too weak.  Without the bright standard star trace, the
tracing crutch used is the slit/order edge, which will not include the effects
of differential atmospheric refraction on the object centroid and therefore
yield a poorer spectral extraction.



Dither sequence
---------------

In this example dataset, the science object and the standard star are both only
observed at two offset positions. Since NIRSPEC does not identify which exposure
is taken in which position, we have to fill that information in manually, as discussed
above.

By setting ``comb_id=3`` and ``bkg_id=4`` for frame ``s190519_0067.fits``, we
are indicating that this frame should be treated as frame A and that frame
B should be used as its background.  We reverse the values of ``comb_id``
and ``bkg_id`` for frame ``s190519_0068.fits``, which indicates that this frame
should be treated as frame B and that frame A should be used as its
background.  I.e., the ``comb_id`` column effectively sets the numeric identity
of each frame, and the ``bkg_id`` column selects the numeric identity of the
frame that should be used as the background image.

When PypeIt reduces the frames in this example, it constructs two images, A-B
and B-A, such that the positive residuals in each image are the observed
source flux for that observation, which can be combined using :ref:`2D coadding
<coadd2d>`.

The use of the ``comb_id`` and ``bkg_id`` integers is very flexible, allowing
for many, more complicated dithering scenarios; see :ref:`a-b_differencing`.

----


Core Processing
===============

To perform the core processing of the NIRES data, use :ref:`run-pypeit`:

.. code-block:: bash

    run_pypeit keck_nires_A.pypeit

The code will run uninterrupted until the basic data-reduction procedures
(wavelength calibration, field flattening, object finding, sky subtraction, and
spectral extraction) are complete; see :doc:`../running`.  Processing of this
example dataset takes roughly 10 minutes (on a 2020 MacBook Pro with 16 GB of
RAM and a 2GHz i5 processor).

As the code processes your data, it will produce a number of files and QA plots
for you to inspect:

Order Edges
-----------

The code first uses the ``trace`` frames to find the order edges.  NIRES is a
fixed-format echelle, meaning that the trace results should always look the
same.  To show the results of the trace, run, e.g.:

.. code-block:: bash
    
    pypeit_chk_edges Calibrations/Edges_A_7_DET01.fits.gz

which will show the image and overlay the traces (green is the left edge;
magenta is the right edge); this should open a `ginga`_ window for you if one
is not open already.  Here is the result from this example dataset:

.. image:: ../figures/nires_trace.png
   :scale: 60%

An important check is to ensure that the code has correctly traced the bluest
(left-most) order.  PypeIt currently expects to find all 5 orders and will fault
if it does not.  
    
.. tip::

    If PypeIt faults because it did not find all 5 orders, try adjusting the
    ``edge_thresh`` parameter; see the :ref:`parameters` and specifically the
    :ref:`edgetracepar`.

Wavelength Calibration
----------------------

Next the code performs the wavelength calibration.  Via the :ref:`pypeit_file`,
we designated two sets of wavelength calibration frames, one for the standard
star and one for the science frame.  You should inspect the results for both.

First, it's important to understand PypeIt's :ref:`calib-naming` convention,
specifically the calibration group bit identities used in the output file names.
In this example, two :ref:`arc` files are produced:
``Calibrations/Arc_A_2_DET01.fits`` and ``Calibrations/Arc_A_4_DET01.fits``.
The ``2`` and ``4`` are the bits associated with the calibration group and link
back to which files are associated with each frame type.  The
:ref:`calibrations-calibfile` provides the direct association of input frame
with output calibration file.

Combined Arc Frame
++++++++++++++++++

You can view the combined arc frame used for, e.g., the standard star
observations with `ginga`_:

.. code-block:: bash

    ginga Calibrations/sArc_A_2_DET01.fits

1D Wavelength Solution
++++++++++++++++++++++

More importantly, you should check the result of the wavelength calibration
using the automatically generated QA file; see :ref:`qa-wave-fit`.  Below is the
wavelength-calibration QA plot for the reddest order (order=3).  The RMS of the
wavelength solution should be of order 0.1-0.2 pixels.  Such a plot is produced
for each order of each the combined arc image used for each calibration group.

.. figure:: ../figures/Arc_1dfit_A_2_DET01_S0003.png
   :width: 70%

   The wavelength-calibration QA plot for the reddest Keck/NIRES order
   (order=3), called ``Arc_1dfit_A_2_DET01_S0003.png``.  The left panel shows
   the arc spectrum extracted down the center of the order, with green text and
   lines marking lines used by the wavelength calibration.  Gray lines mark
   detected features that were *not* included in the wavelength solution.  The
   top-right panel shows the fit (red) to the observed trend in wavelength as a
   function of spectral pixel (blue crosses); gray circles are features that
   were rejected by the wavelength solution.  The bottom-right panel shows the
   fit residuals (i.e., data - model).

In addition, the script :ref:`pypeit-chk-wavecalib` provides a summary of the
wavelength calibration for all orders. We can run it with this simple call:

.. code-block:: bash

    pypeit_chk_wavecalib Calibrations/WaveCalib_A_2_DET01.fits

and it prints on screen the following (you may need to expand the width of your
terminal to see the full output):

