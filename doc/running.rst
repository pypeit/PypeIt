
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

The docs below describe the standard, end-to-end reduction with 
PypeIt.  We point users to the :ref:`pypeit-run-to-calibstep` 
and :doc:`reduce_by_step` documentation for a discussion of how to 
run individual steps of PypeIt (calibrations and primary reductions).
We caution, however, that the code was primarily developed to perform 
end-to-end reductions and we anticipate it will perform best when 
executed in that fashion.

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

    #. (optional) Remove any calibration files in the ``Calibrations/`` folder.
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

.. _run-pypeit-ncpu:

Running on multiple CPUs
------------------------

The number of CPUs PypeIt may use is controlled by the ``[rdx] ncpu``
:ref:`parameter <parameters>` or, equivalently, the ``run_pypeit --ncpu``
command-line option (the command-line value takes precedence).  The default is
``ncpu = 1``, which runs the code serially.

What is parallelized
++++++++++++++++++++

With ``ncpu > 1``, PypeIt reduces **detectors (or detector mosaics)
concurrently**, one worker process per detector.  Three stages of the reduction
are distributed this way:

#. building (or loading) the calibrations,
#. processing the science images, including mosaicing, and
#. extracting the objects, including local sky subtraction.

Together these are roughly 70-77% of the run time of a Keck/DEIMOS reduction.

**Object finding and global sky subtraction remain serial.**  They cannot be
distributed as easily, because they fan out over the detectors, then hit a
cross-detector barrier -- the slitmask offset is computed by aggregating the
objects and slits found on *all* detectors -- and then fan out again.

Independently of the above, ``ncpu > 1`` also hands the QA PNG encoding to a
small pool of background threads, which recovers some of the QA overhead on
slitmask reductions that generate many hundreds of QA files.  That threading is
used only in the serial stages; inside a detector worker the QA figures are
written in line, so that four workers do not each spawn their own pool.

How many CPUs to ask for
++++++++++++++++++++++++

The number of workers is ``min(ncpu, ndetectors, ncores-1)``, so asking for more
CPUs than the data has detectors does nothing.  A single-detector longslit
reduction gains nothing from ``ncpu > 1`` and will be marginally slower, because
of the cost of starting the worker processes.

**In practice the limit is memory, not cores.**  Each worker holds its own
detector's images, so peak memory scales roughly linearly with the number of
detectors reduced at once.  For example, a Keck/DEIMOS reduction of four
mosaics peaks at about 36 GB at ``--ncpu 4``, against about 12 GB serially, and
a single MMT/Binospec long-slit detector needs 10-14 GB, so ``--ncpu 2`` on
Binospec needs 20-28 GB.  On a machine with 32 GB or less, ``--ncpu 2`` is the
safer choice for such data.  If a reduction is killed or the machine starts
swapping, lower ``ncpu`` first.

How the workers are started
+++++++++++++++++++++++++++

Detector parallelism works on Linux, macOS and Windows.  The workers are
started with the ``forkserver`` method of :mod:`multiprocessing` on Linux and
with ``spawn`` elsewhere.  Plain ``fork`` is deliberately not used: by the time
the detectors are distributed, the main process already has threads running,
and forking a multi-threaded process is unsafe on every platform.  Each worker
receives a copy of the inputs that are common to all detectors once, when it
starts, and the inputs specific to its detector with each job.

The start method can be forced by setting the environment variable
``PYPEIT_MP_START`` to any method your platform supports, e.g.:

.. code-block:: console

    PYPEIT_MP_START=spawn run_pypeit my.pypeit --ncpu 2

This is mainly useful for debugging, with one practical exception on macOS.
Under ``spawn`` each worker starts a fresh Python interpreter and imports
PypeIt, which costs a second or two every time the workers are started.  That
happens once for the calibrations and twice per exposure, so for reductions of
many exposures, or of a single detector, ``--ncpu 2`` can end up slower than a serial run.  On macOS,
``PYPEIT_MP_START=forkserver`` avoids that cost, because its workers are
forked from a server that has already imported PypeIt.  It is not the default
there because Apple's system libraries are not guaranteed to be safe to use
after a ``fork``; it has worked in our tests, but use it at your own
discretion.

Reproducibility
+++++++++++++++

Detector parallelism is designed so that the products of a reduction do **not**
depend on ``ncpu``.  How far that holds depends on the platform and the data:

- **On Linux**, where the BLAS/OpenMP thread count can be controlled (see
  below), running with ``--ncpu 4`` produces ``spec1d`` and ``spec2d`` files
  that are identical, bit for bit, to those from ``--ncpu 1`` -- provided the
  reduction is reproducible run-to-run in the first place.  This is verified
  on long-slit data.
- **Elsewhere**, the serial and parallel products are scientifically
  equivalent, but they are not guaranteed to be bit-identical.  In particular,
  on macOS with numpy linked against Apple's Accelerate library,
  ``threadpoolctl`` cannot see the BLAS, so the serial path is not pinned to one
  thread (the workers are, through environment variables they receive when
  they start), and long-slit reductions there have been seen to differ even
  between two identical serial runs.

.. warning::

    On **multi-slit** data, PypeIt is currently not reproducible run-to-run at
    all, independently of ``ncpu``: two identical serial reductions of the same
    data already differ, because edge centroiding returns a slightly different
    answer each time and the sigma-clipped fits downstream amplify it.  See
    `Issue 2208 <https://github.com/pypeit/PypeIt/issues/2208>`__.  Until that
    is fixed, products from two runs of such data will differ whether or not
    ``ncpu`` is changed, and the difference is **not** attributable to the
    parallelism.

Making the products independent of ``ncpu`` requires one deliberate trade-off.  PypeIt pins the BLAS/OpenMP
thread count to one while it reduces a detector, in *both* the serial and the
parallel paths.  It has to pin both: numerical libraries sum a reduction in a
different order when they split it across threads, and although the resulting
difference is at the level of the last bit, the sigma-clipped iterative fits in
the sky subtraction amplify it into a visibly different sky model.  Pinning only
the workers would therefore make the output depend on ``ncpu``.

The cost is that a default ``ncpu = 1`` reduction no longer uses multi-threaded
linear algebra, which makes it a few percent slower than it was in previous
versions.  In exchange, where the pinning takes effect, results are
reproducible across ``ncpu`` settings and across machines with different core
counts.

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


Step-by-step Processing and Monitoring
======================================

Beyond the end-to-end execution described above, PypeIt provides a script to
run the reduction one step at a time and a GUI to monitor and inspect a
reduction as it progresses:

.. toctree::
   :maxdepth: 1

   Step-by-step Processing <reduce_by_step>
   Dashboard <dashboard/dashboard>

See also
========

 - :ref:`pypeit-run-to-calibstep`: Process the calibration step-by-step
 - :ref:`step-by-step-reductions`: Process the science frames step-by-step
 - :ref:`dashboard`: Monitor and inspect a reduction in a GUI
 - :ref:`state`: The on-disk record of a reduction's progress



