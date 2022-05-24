.. _object_finding:

==============
Object Finding
==============

This document describes how the code identifies objects within the slits/orders.

Overview
========

Object identification is a challenging process to code in general. The
challenges are that one must consider a large dynamic range between bright
telluric and spectrophotometric standards and faint continuum or emission-line
sources. When this is coupled with very short slits, the presence of both
positive and negative traces (when image differencing), and the detection of
objects on specific orders (for echelle spectra), many corner cases can arise.
Our general philosophy has been to try to create an algorithm that is fully
automated and detects the most probable sources (:math:`{\rm SNR} > 10\sigma`).
The user has the option to lower this threshold and/or request that specific
objects be manually extracted.

FindObj Parameters
==================

This reduction step is guided by the :ref:`pypeit_par:FindObjPar Keywords`.

Algorithms
==========

The algorithms described below attempts to identify the peak location of objects
on the slit/order as well as measure an initial FWHM of the object. This FWHM
then informs the width of the mask around this object trace in global
sky-subtraction, as well as the width of the region that ``PypeIt`` uses for
local sky-subtraction.

In a standard run, the approach to object finding adopted is:

    #. Perform a first global sky-subtraction without any masking of objects.
    #. Run the object finding code for the first time.
    #. Create a mask indicating object free sky pixels using the locations of
       the objects identified.
    #. Perform a second global sky-subtraction with the objects masked.
    #. Run the object finding code a second time.

These steps are performed by :func:`pypeit.find_objects.FindObjects.run`.

automatic
---------

The automated object finding algorithm is
:func:`~pypeit.core.findobj_skymask.objs_in_slit`. It performs the following
steps:

    #. Rectify the sky-subtracted frame.

    #. Create a mask by sigma clipping (median stat) down the spectral dimension
       to further reject CRs.  This could mask bright emission lines.

    #. Sum over the spectral direction but using the sigma clipping resulting in
       a vector representing counts as a function of spatial slit/order. Compute
       the formal S/N ratio of this collapsed image using the noise model.

    #. Smooth this 1d vector which represents S/N vs spatial slit/order position.

    #. Search for peaks in this S/N ratio vector above the threshold
       ``snr_thresh``. This is the quantity that is shown as the histogram in
       the object finding QA plots (link an example!). The default value is
       ``snr_thresh = 10``.

    #. Objects that are within ``find_trim_edge`` of the slit edges will be removed.

    #. Of the good objects (i.e. not close to the edge), the number of objects
       returned will be further restricted to a maximum number with the highest
       S/N ratio. Note that this maximum is set by ``maxnumber_sci`` parameter
       for science exposures, versus the ``maxnumber_std`` for standard
       exposures. For multislit the defaults are ``maxnumber_sci = 10`` and
       ``maxnumber_std = 5``. For echelle spectrographs with short slits they
       are ``maxnumber_sci = 2`` and ``maxnumber_std = 1``.

The most common parameter modifications we recommend are to adjust
``snr_thresh`` to enable the identification of fainter sources (at the risk of
false positives).  Reasonable results have been obtained with ``snr_thresh`` as
low as 5.0.

To make this modification, add the following to your
:doc:`pypeit_file`:

.. code-block:: ini

    [reduce]
        [[findobj]]
            snr_thresh = 5.0

The QA plot (see linked) clearly indicates the objects that are being excluded
due to proximity to the edges, and/or the maxnumber parameter. These can also be
adjusted, for example for an echelle spectrum with three objects on the slit
where one wants to work closer to slit edges:

.. code-block:: ini

    [reduce]
        [[findobj]]
            maxnumber_sci = 3
            find_trim_edge = 3,3

Finally, if your objects have spectral breaks, like with high-z quasars or
galaxies, then you can restrict the spectral region that the automated object
finding collapses out to search for objects via

.. code-block:: ini

    [reduce]
        [[findobj]]
            find_min_max = 1600, 2048

Which will only collapse out spectral pixels 1600-2048 when computing the 1d SNR
vs spatial position vector. The best way to choose these pixels is to run pypeit
without it set. Then run :ref:`pypeit_show_2dspec` to view the sky-subtracted
image and decide which pixels to use for object finding. Then re-run ``PypeIt``.

.. todo::

    WIP: Put in manual extraction stuff here in place of this stuff on
    interactive object finding.

Interactive object finding/tracing
----------------------------------

.. warning::

    The ``pypeit_find_objects`` script has been deprecated until it can be
    updated.  The following description is outdated!!

THE FOLLOWING IS UNDER DEVELOPMENT.

In some cases, the code may not find the object that you're after,
or may find several spurious objects. To add/remove/modify object
traces interactively, there is an interactive GUI utility.

pypeit_find_objects Science/spec2d.fits

and this will launch an interactive GUI that will allow you to perform
several simple operations on the object tracing. The
tool will produce a few lines of text that you can insert
into your .pypeit file, and this will allow for a
reproducible data reduction.

Using this tool, you will be able to delete spurious traces, add new object traces,
and manually set the FWHM of the object profile. To view a complete list of
the supported functions, press the '?' key on your keyboard when the
mouse is hovering over the panel displaying the 2D image. The detailed
information will be printed to the terminal (i.e. it is not displayed
on the GUI). Below we discuss some of the operations you can perform
with this GUI.

The main panel displays the 2D sky subtracted image of the data.
Darker shades correspond to higher flux. The green/blue lines display
the left/right slit edges. Dashed red lines indicate the object traces
currently stored. You can select an object trace by clicking
(left mouse button) near an object trace; the selected object trace
will be highlighted by a solid thick red line.

The bottom right panel displays the object profile (the profile is
only displayed when an object is selected). By clicking (left mouse
button) on this panel, you can set the FWHM of the object trace. The
FWHM is indicated by the vertical red lines in this panel.

The top (information) panel will provide information about the current
status of the object tracing, and will sometimes prompt the user for
a yes/no response (e.g. "Are you sure you want to delete this trace?").
You can select the answer by clicking on the yes/no button when they
appear.

Finally, there are two buttons on the right hand side of the GUI that
allow you to exit the tracing and print out a script for you to
include in your .pypeit file. **Please use these exit buttons instead of killing the window
from the menu bar**. The button labelled "Continue (and save changes)"
will exit the session and print to screen the relevent text needed
for inclusion in the .pypeit file. The button labelled
"Continue (don't save changes)" will exit the interactive session and
all of your interactive changes will be ignored.

Just below these exit buttons there are four radio buttons that allow
you to select a method to trace the object profiles. Below is a
description of each of these models:

+ *Object* - If there is a bright object on the same slit,
  the trace of the bright object will be used.
+ *Standard Star* - If a standard star has been acquired,
  the trace defined by the standard star will be used.
+ *Slit Edges* - The object trace will follow the same functional
  form as the function that defines the slit edge tracing.
+ *Manual* - Allows the user to manually define an arbitrary slit.

You can use the matplotlib tools to zoom in on the data frame (e.g.
using the rectangular selection tool). To toggle the panning and
zoom feature with the mouse button, press the 'p' key. To return
back to the original plotting extent, press the 'h' or the 'r' key.

To define a new object trace, select one of the first three methods
above, hover the mouse to the location you would like to lay down an
object trace, and press the 'a' key on the keyboard.

When using the "manual" object trace method, you need to define the
anchor points of the object trace. To define the anchor points, hover
the mouse to a location where you see data for the object and press
the 'm' key. This will add a point that helps to identify the object
trace. Add as many points as needed to accurately define the object
trace (a green curve displays the fitted object trace, while single
bullet points define the anchor points). To increase/decrease the
fitting order of the polynomial, press the '+/-' keys on the keyboard.
To delete an individual anchor point, hover near the anchor point
you wish to delete and press the 'n' key. Alternatively, if you want
to clear all anchor points and start again, press the 'c' key. Once
you are satisfied with the green curve defining your object trace,
press the 'a' key to add this to the object tracing.

The delete an object trace, select the object trace by clicking the
left mouse button near the object trace. Once selected, press the
'd' key. If you're sure you want to delete this trace, select "Yes"
from the information panel.

.. The following lines are commented out.
.. The script usage can be displayed by calling the script with the
.. ``-h`` option:

.. .. include:: help/pypeit_find_objects.rst


