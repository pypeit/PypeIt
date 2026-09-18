.. include:: ../include/links.rst

.. _arc_tspec_howto:

=================
ARC-TSPEC-HOWTO
=================

Overview
========

This document provides a guide to using the TripleSpec spectrograph on the
ARC 3.5m with PypeIt.  It covers the specific configurations, calibration
procedures, and data reduction steps necessary.

TripleSpec is a fixed format spectrograph that records spectra in the
near infrared. It records five spectral orders for a slit of length
43 arcsec.

External (truss) lamps are available for flat fielding and wavelength 
calibration, but wavelength calibration is generally achieved using
sky lines in the near-IR. Near-IR flats are generally taken both
with truss lights on and off, so that thermal contributions from the
telescope can be subtracted.

Most near-IR observations are obtained with the object being "nodded"
along the slit to allow for pixel-to-pixel sky subtraction. PypeIt
has the ability to process difference images, but the user needs
to manually set up the difference image pairing.

Setup
=====

Since TripleSpec has a fixed configuration, there should only
be one setup file per night.

Run pypeit_setup
----------------

You create PypeIt setup files using the pypeit_setup command, giving
it the name of spectrograph and the directory with the raw data. You
probably want to run this from a separate directory in which you want
the reduced data to reside.

.. code-block:: bash

       pypeit_setup -s arc_tspec -r PATH_TO_RAW_DATA/ --gui

The --gui argument will make this open a graphical interface from which
you can inspect, edit, and/or save the different configurations. You can
also run pypeit_setup without this argument and let it create the 
configuration files by itself without showing you the details. 

pypeit_setup attempts to sort files into separate configurations, and
automatically identify which files are to be used for pixelflats,
illumflats, bias, arcs, tilts, standards, and science frames; see
:ref:`setup_doc` for additional details.

Since ARC TripleSpec headers do not contain information about the
truss lamps states (as of April 2026), pypeit_setup cannot determine
which images are flats, and which flats have the lamps on and which
have lamps off. As a result, you will need to manually identify these
in the .pypeit file, either through the GUI or after the file is saved.
For the illuminated flats, add  ``illumflat,pixelflat,trace`` to the frametype
column, and for the unilluminated flats, add the keyword ``lampoffflats``
to the frametype column.

You will also need to identify which frames to use for wavelength
calibration (frametype arc, tilt); you probably just want to use a single
science frame for this.

Finally, you will need to identify frame groups for difference imaging.
See the PypeIt documentation on :ref:`a-b_differencing`.

Note that the default setup for Triplespec turns off cosmic ray detection
since this seems to have issues over-detecting CRs. If you want to use
this, you will need to enable CR detection using mask_cr = True, and you 
will likely need to play with the LAcosmic CR parameters; see :ref:`image_proc`.

PypeIt has MANY parameters that can be used to adjust the reduction; see
:ref:`parameters` for details if you want/need to modify the default
behavior.

Non-standard binning
--------------------

The pipeline has not been tested for non-standard (non 1x1) 
binning. Please contact the developers if there are issues.


Main Run
========

Once the :doc:`../pypeit_file` is ready, the main call is
simply:

.. code-block:: bash

    run_pypeit arc_tspec_A.pypeit -o

The ``-o`` indicates that any existing output files should be overwritten.  As
there are none, it is superfluous but we recommend (almost) always using it.

The :doc:`../running` doc describes the process in some more detail.

Note that the PypeIt reduction process can take quite a while!

Inspecting Files
================

Calibrations
------------

The first set are :doc:`../calibrations/calibrations`.

Slit Edges
++++++++++

PypeIt will map the slit edges using the trace frames.

Wavelengths
+++++++++++

One should inspect the :doc:`../qa` for the wavelength
calibration.  These are PNGs in the QA/PNG/ folder.

Note:  there are multiple files generated for every order.
When the reduction is complete, you may prefer to scan
through them by opening the MF_[configuration].html file under ``QA/``;
note, however, that you may find other useful PNG files in the QA/PNG
directory that are not listed in the  MF_[configuration].html file.

The final wavelength solution is a two dimensional fit 
with pixel along the order one axis and the order itself
being the second. 

Remember, the default calibration is in vacuum wavelengths.


Spectra
-------

The code will generate 2D and 1D spectra outputs.  One per 
science frame, located in the ``Science/`` folder.

One can inspect the one dimensional spectra 
with :ref:`pypeit_show_1dspec`.

Using that tool, you can examine the individual orders. 
