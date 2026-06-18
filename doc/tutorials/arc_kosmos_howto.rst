.. include:: ../include/links.rst

.. _arc_kosmos_howto:

=================
ARC-KOSMOS-HOWTO
=================

Overview
========

This document provides a guide to using the KOSMOS spectrograph on the
ARC 3.5m with PypeIt.  It covers the specific configurations, calibration
procedures, and data reduction steps necessary.

KOSMOS has multiple configurations. It can be used in longslit mode
or with slit masks. In longslit mode, three different slit positions 
(low, center, high) can be used to adjust the wavelength coverage.
In addition, two grisms, for blue and red coverage, are available.

Both internal and external (truss) lamps are available for flat
fielding and wavelength calibration. The external lamps should provide
better calibration data, but require longer exposure times. The
current default KOSMOS setup is for the truss lamps, but this can
be overridden, see below.

There is non-negligible flexure in the instruments. Some users take
both internal and external wavelength calibrations at the beginning
and/or end of night, and then internal wavelength calibrations at the 
location of each object to do differential wavelength adjustments. 
PypeIt is not set up to do this, so current options are 1) to do
spectral flexure correction using sky lines (should be good in red
at least), or 2) use the internal lamps taken at the location
of the object for wavelength calibration.


Setup
=====

PypeIt will be run separately for each configurations used during
a night of data. 

Run pypeit_setup
----------------

You create PypeIt setup files using the pypeit_setup command, giving
it the name of spectrograph and the directory with the raw data. You
probably want to run this from a separate directory in which you want
the reduced data to reside.

.. code-block:: bash

    pypeit_setup -s arc_kosmos -r PATH_TO_RAW_DATA/ --gui

The --gui argument will make this open a graphical interface from which
you can inspect, edit, and/or save the different configurations. You can
also run pypeit_setup without this argument and let it create the 
configuration files by itself without showing you the details. 

pypeit_setup attempts to sort files into separate configurations, and
automatically identify which files are to be used for pixelflats,
illumflats, bias, arcs, tilts, standards, and science frames; see
:ref:`setup_doc` for additional details.

Different configurations are defined as those with different grism,
slit, and binning combinations. Note that, in the setup gui and the
configuration files, the grating is the grism, and the decker is the
slit.

If you want to use different calibration frames for different science
images (e.g., for internal lamps at the position of each object), 
you should create a separate .pypeit file for each group.

After running this, you should have one or more .pypeit configuration
files, e.g. arc_kosmos_A.pypeit. You should inspect these files to
make sure that the correct classifications have been made. You can
also remove or comment out individual frames.

Wavelength calibration frames (listed as arcs), **should only be
for the HeNeAr truss lamps**; they can be in separate or a single
exposure. The default line lists are for HeI, NeI, and ArI; if not
all of these lamps are used, add a section to the .pypeit file
to specify which lamps were used, e.g.:

.. code-block:: ini
  [calibrations]
   [[wavelengths]]
     lamps = NeI, ArI

Wavelength calibration works by cross-correlating a combined
arc frame with a reference frame. The reference frames have been
constructed from observations with He, Ne, and Ar lamps on, which
can lead to issues if the combined frame only has a subset, if
the missing lamp(s) have a lot of lines. There is a Ne only
reference frame available for the red configuration, to use that
you would add :

.. code-block:: ini

  [calibrations]
   [[wavelengths]]
     lamps = NeI
     reid_arxiv = 'arc_kosmos_red_ctr_ne.fits'

If you used the internal lamps in the blue, you can
use a reference frame made with Kr, Ne, and Ar:

.. code-block:: ini

  [calibrations]
   [[wavelengths]]
     lamps = KrI,NeI,ArI
     reid_arxiv = 'arc_kosmos_red_int_ctr.fits'

Since lines in the red are largely from Ne and Ar, the standard
red reference frame should work with either the internal or the
truss lamps.

Finally, if you are using slit masks and have slits far from the
center position, the wavelength identification may fail. In some cases,
we have found that this can be improved using:

.. code-block:: ini
  [calibrations]
   [[wavelengths]]
     nsnippet = 5

PypeIt has MANY other parameters that can be used to adjust the reduction; see
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

    run_pypeit arc_kosmos_A.pypeit -o

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
calibration.  These are PNGs in the QA/PNG/ directory.

Note:  there are multiple files generated for every slit.
When the reduction is complete, you may prefer to scan
through them by opening the MF_[configuration].html file under ``QA/``;
note, however, that you may find other useful PNG files in the QA/PNG
directory that are not listed in the  MF_[configuration].html file.

Note that it is possible to have apparently small residuals
but still have an incorrect wavelength solution if lines are
misidentified! This often leads to only having relatively few lines.
Output science spectra can be the final arbiter of problems, e.g. by
looking at sky lines.

Remember, the default calibration is in vacuum wavelengths.


Spectra
-------

The code will generate 2D and 1D spectra outputs.  One per 
science frame, located in the ``Science/`` folder.

One can inspect the one dimensional spectra 
with :ref:`pypeit_show_1dspec`.

