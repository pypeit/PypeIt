
.. include:: include/links.rst

.. _setup_doc:

*****
Setup
*****

Overview
========

The :doc:`pypeit_file` is *the* critical component to any successful run of
PypeIt, and :ref:`pypeit_setup` provides an automated means of generating this
file based on a set of fits files.

Below, we describe how PypeIt defines unique instrument configurations
("setups") and sorts data in preparation for the data reduction.  Then we
describe our recommended method for generating the :ref:`pypeit_file`, which
includes multiple executions of :ref:`pypeit_setup`.  We also describe another
preparatory script, :ref:`pypeit_obslog`, which provides a simple listing of the
available data files; however, use of this script is optional in terms of setup
for reducing your data.

.. _setup-metadata:

Use of Metadata to Identify Instrument Configurations
=====================================================

PypeIt distinguishes between instrument configurations using metadata pulled
from the fits file headers, as defined specifically for each spectrograph. To
list *all* of the metadata pulled from each data file and how it maps to the
keywords in the data file headers, run, e.g.:

.. code-block:: bash

    pypeit_obslog keck_deimos -k

(see :ref:`pypeit_obslog` and
:func:`~pypeit.spectrographs.spectrograph.Spectrograph.meta_key_map`).  For Keck/DEIMOS, this
prints:

.. include:: include/deimos_meta_key_map.rst

where the left column provides the PypeIt-specific metadata
keyword and the right column is the associated instrument-specific
header card. Any metadata element with a header card set to ``None``
means that the metadata keyword or value is conditioned on multiple
header values. In this example, the header card that provides the
central wavelength given the grating angle (``dispangle``) for DEIMOS
is dependent on the grating used (``dispname``); see
:func:`~pypeit.spectrographs.keck_deimos.KeckDEIMOSSpectrograph.compound_meta`.

A subset of this metadata is used to determine the unique instrument
configurations used and associate all of the data files to the relevant
configuration.  Generally speaking, most instrument configurations are set by
the name (``dispname``) and central wavelength (``dispangle``) of the dispersing
element, the name of the decker or slit mask (``decker``), the name of the
beam-splitting dichroic (``dichroic``), and/or the on-chip binning
(``binning``).  The specific set of metadata keys used for each instrument is
set by the ``configuration_keys`` method; for DEIMOS, see
:func:`~pypeit.spectrographs.keck_deimos.KeckDEIMOSSpectrograph.configuration_keys`.

The following table provides the list of metadata used to define the instrument
configuration for *any* supported spectrograph (as of PypeIt version
``1.10.1dev``):

====================  ======= =============== ===================================================================
Metadata Key          Type    Example         Description
====================  ======= =============== ===================================================================
``amp``               str     SINGLE:B        Name of the amplifier used to read the detector
``arm``               str     VIS             Name of the spectrograph arm used to collect the data
``binning``           str     1,1             On-chip binning
``cenwave``           str     7500.0          Central wavelength of the disperser
``datasec``           str     [1:256,1:512]   The science region of the detector
``decker``            str     long_1.0        Name of the decker or slit mask
``decker_secondary``  str     long_1.0        Partial Slitmask/decker name; see :ref:`mosfire_config_report`
``detector``          str     CHIP1           Name of the detector
``dichroic``          str     560             Name of the dichroic
``dispangle``         float   7500.0          Central wavelength for the dispersing element at the observed angle
``dispname``          str     830G            Name of the dispersing element
``filter1``           str     J               Name of the order-sorting filter
``slitlength``        str     ..              Slit length; see :ref:`mosfire_config_report`
``slitwid``           str     ..              Slit width; see :ref:`mosfire_config_report`
====================  ======= =============== ===================================================================

*Every unique combination* of the relevant metadata found in any of
the fits files to be reduced represents a unique instrument
configuration for PypeIt to consider; these configurations can be
determined automatically by :ref:`pypeit_setup` and will be
identified by a capital letter, e.g., **A**. When executing
:ref:`run-pypeit`, however, the instrument configuration is set by
the :ref:`setup_block` in the :ref:`pypeit_file` and not redetermined
by the fits files. The latter allows the user flexibility to override
PypeIt's automated configuration settings. **Currently, each**
:ref:`pypeit_file` **should only provide data from one instrument
configuration.**

If you tend to observe with one instrument configuration and with a
simple set of calibrations, then the setup to run PypeIt should
be straightforward. However, if you use multiple configurations (e.g.
gratings, grating tilts, slitmasks), then you must pay more careful
attention to the setups.

Below, we describe how PypeIt automatically determines instrument
configurations for a set of files and constructs auto-generated
pypeit files.

----

.. _pypeit_obslog:

pypeit_obslog
=============

The ``pypeit_obslog`` script allows you to see a simple listing of the data
files in a given directory (or directories) and the metadata that PypeIt
will pull from their headers. As always, the script usage can be displayed by
calling the script with the ``-h`` option:

.. include:: help/pypeit_obslog.rst

For example, if you've been observing with Keck DEIMOS, you can go into the
directory with the raw data and execute:

.. code-block:: console

    pypeit_obslog keck_deimos

Or you can point to the directory that you want to list:

.. code-block:: console

    pypeit_obslog keck_deimos -r /path/to/raw/data

The listing is printed to ``stdout``, but you can also specify a file for the
output. You can also select the metadata by which to sort the output, change
the columns that are printed, etc.

By default, the columns provided in the output table are identical to the
columns in the :ref:`data_block` of a :ref:`pypeit_file`; to print columns with
all the metadata collected, add ``-c all`` to the command line.

----

.. _pypeit_setup:

pypeit_setup
============

``pypeit_setup`` is the script one executes to prepare for the data reduction by
automatically associating fits files to specific :ref:`frame_types` and
collecting groups of frames taken using unique instrument configurations. The
script usage can be displayed by calling the script with the ``-h`` option:

.. include:: help/pypeit_setup.rst

The following is a step-by-step procedure for preparing for the data
reduction.

0. Prepare
----------

Change your working directory to the one where you wish the reduced data
products to appear. Any directory will do; however, you likely *don't* want this
to be the same directory that holds the raw data.

1. First Execution
------------------

We recommend you first execute ``pypeit_setup`` like this::

    pypeit_setup -r path_to_your_raw_data/LB -s keck_lris_blue

where the two command-line options are *required* and provide:

  - ``-s``: Sets the spectrograph to be reduce. This is camera specific (e.g.,
    you need to reduce LRIS red data separately from LRIS blue data).

  - ``-r``: Sets the full path to the raw data and the root prefix of
    the FITS files

Specifically note that this first run does *not* include the ``-c``
argument.

This execution of ``pypeit_setup`` searches for all `*.fits` and
`*.fits.gz` files with the provided root directory. Generally, the
provided path should **not** contain a wild-card and it is best if
you provide the *full* path; however, you can search through multiple
directories as follows::

    pypeit_setup -r "/Users/xavier/Keck/LRIS/data/2016apr06/Raw/*/LB" -s keck_lris_blue

2. Inspect the outputs
----------------------

The call above creates two files in the ``setup_files/`` folder, created where
you executed ``pypeit_setup``:

  - ``{spectrograph}.sorted``: This shows the unique configurations and the list
    of frames associated with each configuration as determine by the automated
    procedures in PypeIt. Each unique configuration is given a capital letter
    identifier (e.g., A,B,C,D...).

  - ``{spectrograph}.obslog``: This provides the default log file, identical to
    the result of running :ref:`pypeit_obslog`.

Here is an example ``.sorted`` file for some example Keck DEIMOS data:

.. include:: include/keck_deimos.sorted.rst

This ``sorted`` file contains two configurations, ``A`` and ``B``.
The "setup block" (the information between the ``###...`` and
``#--...`` lines) describes the instrument configuration specific to
this setup as determined using the relevant metadata. Each "setup
block" is followed by a table listing the files and relevant metadata
for all files matched to that instrument configuration. The data
provided is specific to each instrument, as defined by, e.g.,
:func:`~pypeit.spectrographs.keck_deimos.KeckDEIMOSSpectrograph.pypeit_file_keys`.

The ``sorted`` file is only provided as a means of assessing the
automated setup identification and file sorting, and we encourage you
to briefly review the results. You may recognize that you are missing
calibrations or you may be surprised to see more configurations than
you were expecting. In the latter case, you can help us improve the
automated procedures in PypeIt by submitting an issue on GitHub
(`Submit an issue`_) and including the sorted file in the issue
description.

Importantly, you should use the ``sorted`` file to decide which
configuration (as selected by its letter) you wish to reduce. Also
note that PypeIt cannot interpret any edits you make to this
file; all user-level edits to the frame-typing, association of frames
with given configurations, etc., *must* be done via the :doc:`pypeit_file`.

3. Second execution: Write the pypeit file for one or more setups
-----------------------------------------------------------------

Provided you are happy with the ``sorted`` file, you should execute
:ref:`pypeit_setup` a second time with the `--cfg_split` (shortcut
`-c`) option. This will generate one or more sub-folders and populate
each with a :doc:`pypeit_file`. Some example uses of the ``-c``
option are:

    - ``-c A``: This will generate one folder+file for the chosen
      configuration
    - ``-c A,C``: This will generate one folder+file for each input
      configuration
    - ``-c all``: This will generate folders+files for all
      configurations

An example execution that only produces the :ref:`pypeit_file` for
the A configuration is::

    pypeit_setup -r path_to_your_raw_data/LB -s keck_lris_blue -c A

This example will generate a new folder named ``keck_lris_blue_A``
and within it will be a file named ``keck_lris_blue_A.pypeit``.

4. Edit the pypeit file
-----------------------

You will likely need to edit your :ref:`pypeit_file`.  This includes adding any
desired changes to the :ref:`parameter_block`, adding/removing files from the
:ref:`data_block`, correcting any erroneous frametypes (see :ref:`frame_types`),
specifying any desired association between science and calibration frames (see
:ref:`calibration-groups`), which science frames to combine (see
:ref:`2d_combine`), and which frames to use as backgrounds and/or part of an AB
offset sequence (see :ref:`a-b_differencing`).

.. warning:: 

    Any execution of :ref:`run-pypeit` will *crash* if your :doc:`pypeit_file`
    includes entries with ``None`` frame types defined.  You must either remove
    or edit those entries in the pypeit file by-hand after running
    :ref:`pypeit_setup`.

----

Options
-------

Consider your need for the following additional options:

-b option
+++++++++

You should consider using the ``-b`` option if you need to specify:

- the calibration frames that should be used with each science frame (the same
  can be achieved by dividing your dataset into multiple ``.pypeit`` files that
  all use the same calibrations),

- groups of science frames that should be combined, and/or

- groups of frames that should be treated as background frames in an on-off
  (e.g., ABBA) observation sequence.

This simply adds the relevant columns to the :ref:`pypeit_file` that you *must 
edit by hand*: ``calib``, ``comb_id``, and ``bkg_id``.  These columns are added
by default for most near-IR spectrographs.  See :ref:`calibrations` and
:doc:`A-B_differencing` for the syntax used for the data in these columns and
how PypeIt uses them; see also a worked example in the :ref:`gnirs_howto`.

-m option
+++++++++

If you wish to include a column where you can include
input for :doc:`manual`, use the ``-m`` option.

