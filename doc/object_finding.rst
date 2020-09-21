==============
Object Finding
==============

This document describes how the code identifies
objects within the slits/orders.

Overview
========

Object identification is a challenging process to
code, especially to allow for a large dynamic range
between bright continuum sources and faint emission
line sources.

Our general philosophy has been to err on the
bright side, i.e. detect more probable sources.

FindObj Parameters
==================

This reduction step is guided by the
:ref:`pypeit_par:FindObjPar Keywords`.

Algorithms
==========

Each of the algorithms described below attempt to
identify the peak location of objects in the slit
and then defines a left and right edge for each source.

In a standard run, these are performed twice.
Once before any sky subtraction and then again
after the global sky subtraction has been performed.

automatic
---------

The standard algorithm performs the following steps:

1. Rectify the sky-subtracted frame

2. Smooth this 2D image

3. Perform sigma clipping (median stat) down the wavelength dimension to further reject CRs.  This may eliminate bright emission lines.

4.  Smash the 2D image along the spectral dimension, to get a 1D array that represents the spatial profile of the exposure.

5. If **find_cont_fit** is True, fit a continuum to the profile and subtract it.

6.  Estimate the RMS scatter in the profile array.

7.  Find all objects with peak flux in excess of **sig_thresh** times the RMS.

8.  Eliminate any objects within a few percent of the slit edge. Parameterized by `trace object xedge`.

9.  Determine edges and background regions for each object.

10.  Restrict to **maxnumber** of objects, ordered by flux.

The most common modification that we recommend is to reduce
**sig_thresh** to enable the identification of fainter sources
(at the risk of false positives).  Reasonable results have been
obtained with **sig_thresh** as low as 3.0.

To make this modification, add the following to your
:doc:`pypeit_file`::

    [reduce]
      [[findobj]]
        sig_thresh = 3.



Interactive object finding/tracing
----------------------------------

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

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_find_objects.rst
