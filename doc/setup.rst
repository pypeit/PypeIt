
.. include:: include/links.rst

.. _setup_doc:

*****
Setup
*****

Overview
========

The :doc:`pypeit_file` is *the* critical component to any successful
run of ``PypeIt``, and :ref:`pypeit_setup` provides an automated
means of generating this file based on a set of fits files. Below, we
describe how ``PypeIt`` defines unique instrument configurations
("setups") and sorts data in preparation for the data reduction.

``PypeIt`` distinguishes between instrument configurations using
metadata pulled from the fits file headers, as defined specifically
for each spectrograph. To list the metadata used to establish the
instrument configuration, e.g.:

.. code-block:: python

    from pypeit.spectrographs.util import load_spectrograph
    spec = load_spectrograph('keck_deimos')
    cfg_keys = spec.configuration_keys()
    for key in cfg_keys:
        print('{0:>10} {1:>10}'.format(key, str(spec.meta[key]['card']))) 

which will print::

      dispname   GRATENAM
        decker   SLMSKNAM
       binning       None
     dispangle       None
           amp    AMPMODE

where the left column provides the ``PypeIt``-specific metadata
keyword and the right column is the associated instrument-specific
header card. Any metadata element with a header card set to ``None``
means that the metadata keyword or value is conditioned on other
header values. In this example, the header card that provides the
central wavelength given the grating angle (``dispangle``) for DEIMOS
is dependent on the grating used (``dispname``); see
:func:`~pypeit.spectrographs.keck_deimos.compound_meta`. Generally
speaking, most instrument configurations are set by the name
(``dispname``) and central wavelength (``dispangle``) of the
dispersing element, the name of the decker or slit mask (``decker``),
the name of the beam-splitting dichroic (``dichroic``), and/or the
on-chip binning (``binning``). As of ``PypeIt`` version ``1.2.1dev``,
the following table provides the list of metadata used to define the
instrument configuration for *any* supported spectrograph.

=============== ======= =============== ===================================================================
Metadata Key    Type    Example         Description
=============== ======= =============== ===================================================================
``amp``         str     SINGLE:B        Name of the amplifier used to read the detector
``arm``         str     VIS             Name of the spectrograph arm used to collect the data
``binning``     str     1,1             On-chip binning
``datasec``     str     [1:256,1:512]   The science region of the detector
``decker``      str     long_1.0        Name of the decker or slit mask
``detector``    str     CHIP1           Name of the detector
``dichroic``    str     560             Name of the dichroic
``dispangle``   float   7500.0          Central wavelength for the dispersing element at the observed angle
``dispname``    str     830G            Name of the dispersing element
``filter1``     str     J               Name of the order-sorting filter
=============== ======= =============== ===================================================================

*Every unique combination* of the relevant metadata found in any of
the fits files to be reduced represents a unique instrument
configuration for ``PypeIt`` to consider; these configurations can be
determined automatically by :ref:`pypeit_setup` and will be
identified by a capital letter, e.g., **A**. When executing
:ref:`run-pypeit`, however, the instrument configuration is set by
the :ref:`setup_block` in the :ref:`pypeit_file` and not redetermined
by the fits files. The latter allows the user flexibility to override
``PypeIt``'s automated configuration settings. Currently, each
:ref:`pypeit_file` should only provide data from *one* instrument
configuration.

If you tend to observe with one instrument configuration and with a
simple set of calibrations, then the setup to run ``PypeIt`` should
be straightforward. However, if you use multiple configurations (e.g.
gratings, grating tilts, slitmasks), then you must pay more careful
attention to the setups.

Below, we describe how ``PypeIt`` automatically determines instrument
configurations for a set of files and constructs auto-generated
pypeit files.

.. _pypeit_setup:

pypeit_setup
============

``PypeIt`` includes the :ref:`pypeit_setup` script that one executes
to prepare for the data reduction by automatically associating fits
files to specific :ref:`frame_types` and collecting groups of frames
collected in a unique instrument configuration. The script usage can
be displayed by calling the script with the ``-h`` option:

.. include:: help/pypeit_setup.rst

The following is a step-by-step procedure for preparing for the data
reduction.

0. Prepare
----------

Change your working directory to the one where you wish the reduced
data products to appear. Any directory will do.

1. First Execution
------------------

We recommend you first execute ``pypeit_setup`` like this::

    pypeit_setup -r path_to_your_raw_data/LB -s keck_lris_blue

where the two command-line options are *required* and provide:

  - ``-s``: Sets the spectrograph to be reduce. This is camera
    specific.
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

The call above creates two files in a ``setup_files/`` folder:

  - ``{spectrograph}_{date}.pypeit``: This is a dummy file that can
    be ignored.
  - ``{spectrograph}_{date}.sorted``: This shows the unique
    configurations and the list of frames associated with each
    configuration as determine by the automated procedures in
    ``PypeIt``. Each unique configuration is given a capital letter
    identifier (e.g., A,B,C,D...).

Here is an example ``.sorted`` file for some example Keck DEIMOS data:

.. include:: include/keck_deimos.sorted.rst

This ``sorted`` file contains two configurations, ``A`` and ``B``.
The "setup block" (the information between the ``###...`` and
``#--...`` lines) describes the instrument configuration specific to
this setup as determined using the relevant metadata. Each "setup
block" is followed by a table listing the files and relevant metadata
for all files matched to that instrument configuration. The data
provided is specific to each instrument, as defined by, e.g.,
:func:`~pypeit.spectrographs.keck_deimos.pypeit_file_keys`.

The ``sorted`` file is only provided as a means of assessing the
automated setup identification and file sorting, and we encourage you
to briefly review the results. You may recognize that you are missing
calibrations or you may be surprised to see more configurations than
you were expecting. In the latter case, you can help us improve the
automated procedures in ``PypeIt`` by submitting an issue on GitHub
(`Submit an issue`_) and including the sorted file in the issue
description.

Importantly, you should use the ``sorted`` file to decide which
configuration (as selected by its letter) you wish to reduce. Also
note that ``PypeIt`` cannot interpret any edits you make to this
file; all user-level edits to the frame-typing, association of frames
with given configurations, etc., must be done via the :doc:`pypeit_file`.

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

-b option
+++++++++

If you wish to specify pairs (or groups) of files to use for
background subtraction (e.g. A-B), then include the `-b` option. This
simply adds the relevant columns to the :ref:`pypeit_file` that you
will need to edit by hand. The two columns added, ``comb_id`` and
``bkg_id``, are added by default for most near-IR spectrographs. See
:doc:`A-B_differencing` for the syntax used for the data in these
columns and how ``PypeIt`` uses them.

