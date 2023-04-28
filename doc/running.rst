
.. include:: include/links.rst

====================================================
PypeIt's Core Data Reduction Executable and Workflow
====================================================

Overview
========

This document describes PypeIt's fundamental data-reduction executable and
summarizes its workflow.  Links to more detail on specific algorithms and advice
on adjustments that can/should be made to better reduce your data are included
throughout.

----

.. _run-pypeit:

run_pypeit
==========

The main script to run the PypeIt reduction is ``run_pypeit`` and should be
executed from the command-line (once you've activated the relevant PypeIt python
environment; see :ref:`environment`).  Here, we briefly describe how to use the
script and some of its options.

Before executing ``run_pypeit``, you must have

    #. Inspected and fussed with your setups; see :doc:`setup`.

    #. Set your current working directory to one of the setup sub-folders (e.g.,
       ``keck_deimos_A``).

    #. Edited the :doc:`pypeit_file` in that directory as recommended for a
       successful reduction.

    #. (optional) Removed any calibration files in the ``Calibrations/`` folder.
       This is particularly necessary if you're re-reducing data that was
       previously reduced by an older PypeIt version.  However, when in doubt,
       it's good practice to perform a fresh reduction by removing these files
       instead of attempting to re-use them.

.. _run-pypeit-usage:

usage
-----

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/run_pypeit.rst

Standard Call
-------------

A typical run of PypeIt is initiated with a command like::

    run_pypeit keck_lris_blue_multi_600_4000_d560.pypeit -o

The code launches, reads the :doc:`pypeit_file`, initiates a few internals,
and then proceeds to generate *a lot* of messages in the terminal window.
We fear that some of those comments are outdated or even misleading.
In short, only a PypeIt *developer* is likely to make too much sense of them.

Options
-------

There are a few standard options that you should consider.
These are listed in the :ref:`run-pypeit-usage` and we
describe them in a bit more detail here with guidance on
when to (or not to) use them.

-o
++

The ``-o`` or ``--overwrite`` command will over-write any existing
files and directories.  We recommend this be used the majority of the
time.  But if you know you only want to re-reduce a few science frames,
then remove them and run without ``-o``.

-m
++

This ``-m`` or ``--do_not_use_calibs`` flag tells PypeIt to **avoid** using any
existing calibration frames instead of loading from disk.

Using this can *greatly* slow down the code.

-s
++

This is the main debugging mode of PypeIt.  It will generate *a lot* of plots to
the screen.  It is probably too overwhelming for most users, i.e. best for
*developers*.  Instead, aim to inspect the output
:doc:`calibrations/calibrations` and :doc:`QA<qa>` plots.

-c
++

This will cause PypeIt to only reduce the calibrations.  In fact, if your
:ref:`pypeit_file` only has calibrations and you execute ``run_pypeit`` without
this option, an error will be raised that will suggest you run the code in this
mode.

----

.. _run-pypeit-workflow:

Workflow
========

A core design principle that allows PypeIt to reduce data from many different
:doc:`spectrographs/spectrographs` is that it isolates code that is
spectrograph-specific from code that is generally applicable to any slit-based
spectrograph.  In this section, we briefly describe the *general* workflow used
for *all* spectrographs.  For caveats and advice specific to a given
spectrograph, see :ref:`spec_details`.

The most basic workflow used by PypeIt is:

    #. *Initialization*: Parse the :ref:`pypeit_file` and consolidate with
       additional :doc:`dev/metadata` read from the file headers.  Metadata in
       the :ref:`pypeit_file` takes precedence over any metadata read from the
       file headers.  Setup the output directory structure.

    #. *Reduce Standard Star Frames*: If any are included, first reduce any
       standard star exposures.  As part of the ``run_pypeit`` workflow, the
       *only* use of the standard stars is to guide object tracing when the
       object flux is too faint to be traced directly.  Standard stars are used
       for :ref:`fluxing` via separate scripts.

    #. *Reduce Science Frames*: Perform the main reduction procedures required
       for all the science exposures.  Any standard stars associated with the
       science frames are used as a crutch for object tracing.

.. _run-pypeit-reduce-one:

Data Reduction Workflow
-----------------------

For each of the two "reduction" steps listed above, the only difference is the
crutch used for the object tracing:  The standard star frames always use the
slit-edge traces as the crutch, whereas the science frames will use the standard
star trace as the crutch, if one is available.  Whenever the code falls back to
a tracing crutch, beware of differential atmospheric refraction!

From the outer-most to the inner-most loop, the data reduction iterates through:

    #. the set of defined "calibration groups",

    #. the set of "frame combination" groups (or individual frames if no groups
       are defined),

    #. and the set of spectrograph :ref:`detectors<detectors>` or :ref:`mosaic`.

For each detector or detector mosaic, PypeIt:

    #. processes the relevant :doc:`calibration<calibrations/calibrations>`
       frames

    #. performs a first round of estimating a global sky model (see
       :ref:`skysub-global`) and :ref:`object_finding`

All frames undergo :ref:`image_proc` before more frame-specific processing.

For instruments with appropriate slit-mask metadata (e.g., Keck/DEIMOS), PypeIt
will then attempt to match the detected objects to the expected locations of the
objects from the slit-mask metadata, as well as determine which objects from the
slit-mask design were not detected and which detected objects are serendipitous.

Finally, PypeIt loops through all of the detectors or detector mosaics again to
perform the spectral :ref:`extraction`.

Workflow Alterations
--------------------

PypeIt defines three main algorithms for reducing data, depending on the type of
instrument used: ``MultiSlit`` (used for both long-slit and MOS observations),
``Echelle``, and ``IFU``.  The pipeline used for each spectrograph is provided
:doc:`here<spectrographs/spectrographs>`.

Key differences between the three algorithms include:

    - Object finding and extraction is performed independently for
      ``MultiSlit``, but ``Echelle`` reductions match object traces found on
      each order

    - Calibrations are performed in a slightly different order for ``IFU``
      observations

    - Sky-subtraction and flat-fielding are treated differently for ``IFU``
      observations

Beyond this, more fine-grained control over the PypeIt workflow and the
parameters used by many of the lower-level algorithms is enabled by
:ref:`parameters` that can be changed via the :ref:`pypeit_file`.






