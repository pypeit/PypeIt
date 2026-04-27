.. include:: ../include/links.rst

.. _deimos_config_report:

Automated sorting of DEIMOS frames by instrument configuration
==============================================================

Version History
---------------

=========   ================   =========== ===========
*Version*   *Author*           *Date*      ``PypeIt``
=========   ================   =========== ===========
1.0         Kyle Westfall      13 Oct 2020 1.1.2dev
1.1         Debora Pelliccia   12 Jul 2021 1.9.2.dev
=========   ================   =========== ===========

----

Basics
------

Sorting frames by the configuration of the instrument, ensuring that
this configuration is the same for all coupled sets of calibration and
science frames, is performed by the :ref:`pypeit_setup` script; see
:ref:`setup_doc`. :ref:`pypeit_setup` uses automated procedures to sort
the frames and write a :ref:`pypeit_file` for each unique
configuration (or for some down-selected set).

DEIMOS configuration identification
-----------------------------------

The DEIMOS instrument configuration is determined by a unique
combination of the following keywords:

===============     ====================================================================
``fitstbl`` key     Header Key
===============     ====================================================================
``dispname``        ``GRATENAM``
``dispangle``       ``G3TLTWAV`` or ``G4TLTWAV``, depending on the value of ``GRATEPOS``
``decker``          ``SLMSKNAM``
``binning``         ``BINNING``
``amp``             ``AMPMODE``
``filter1``         ``DWFILNAM``
===============     ====================================================================

as determined by
:func:`pypeit.metadata.PypeItMetaData.unique_configurations`. The
``AMPMODE`` value is included, even though ``PypeIt`` (currently)
restricts itself to only attempting to reduce frames read by the B
and A amplifiers; see
:func:`~pypeit.spectrographs.keck_deimos.KeckDEIMOSSpectrograph.valid_configuration_values`.
Additionally, ``PypeIt`` requires all frames to have ``MOSMODE ==
'Spectral'``. Frames that do not match these header keyword
restrictions will not be included in the automatically generated
:ref:`pypeit_file` created by :ref:`pypeit_setup`.

For DEIMOS, the unique configurations are determined by collating the
relevant metadata from the headers of all frames found by a run of
:ref:`pypeit_setup`, *except* those that are designated as bias or
dark frames. The reason is that bias and darks can have header data
(e.g., ``dispangle``) that do not match the instrument configuration
that an observer intended for their use; e.g., the frames were taken
before the instrument was fully configured for the night's
observations. To match these frames to a specific configuration,
``PypeIt`` uses the ``DATE-OBS`` header keyword to match the frames
to the configurations with frames taken on the same date.

..
    warning [DP: this is not valid anymore with the new changes. The bias/darks frames are assigned to every configuration taken on the same date]
    The fact that the bias and dark frames are matched by date to a
    *single* configuration leads to a complication if the DEIMOS
    configuration is changed during the night. I.e., any biases/darks
    taken on the same date will only be associated with the first
    configuration. It will be up to the observer to edit the
    auto-generated pypeit files constructed for each configuration by
    :ref:`pypeit_setup` to ensure that the correct biases and darks
    are included for each configuration. It is possible to include
    the same biases/darks in multiple configurations in the case
    that, e.g., biases/darks were taken in the evening but not in the
    morning after the instrument configuration change.

DEIMOS calibration groups
-------------------------

``PypeIt`` uses the concept of a "calibration group" to define a
complete set of calibration frames (e.g., arcs, flats, biases) and
the science frame to which these calibration frames should be
applied. By default, :ref:`pypeit_setup` uses the configuration
identifier (e.g., ``A``) to assign frames to a single calibration
group. No automated procedure exists to do anything except this.
However, a user can edit the :ref:`pypeit_file` to, within a given
configuration, assign specific calibration frames to specific science
frames using the data in the ``calib`` column. For example, if the
observer takes both evening and morning calibration frames, the
``PypeIt`` user can use the :ref:`pypeit_file` to associate
calibrations taken at the end of the night with the morning
calibrations and vice versa. For DEIMOS specifically, columns with
the ``DATE-OBS`` and ``UTC`` header data for each frame are provided
to assist with this.

Testing
-------

Requirement PD-3 states: "As a user, I want the pipeline to
automatically and correctly associate calibrations with science
data."

``PypeIt`` meets this requirement in the majority of use cases. One
exception to this is given in the warning above with respect to how
bias and dark frames are assigned to "calibration groups".

The test used to demonstrate PD-1 is satisfied
(:ref:`deimos_frames_report`) is also relevant here in that each
directory in the ``PypeIt`` dev suite with DEIMOS data correctly
identifies the frame types and associates them with a single
configuration, all written to a single pypeit file.

To test that ``PypeIt`` can successfully identify multiple
configurations among a set of files, we have added the
``test_setup_keck_deimos_multiconfig`` and
``test_setup_keck_deimos_multiconfig_clean`` tests to
``${PYPEIT_DEV}/unit_tests/test_setups.py``.

To run these tests:

.. code-block:: bash

    cd ${PYPEIT_DEV}/unit_tests
    pytest test_setup.py::test_setup_keck_deimos_multiconfig -W ignore
    pytest test_setup.py::test_setup_keck_deimos_multiconfig_clean -W ignore

The tests require that you have downloaded the ``PypeIt``
:ref:`dev-suite` and defined the ``PYPEIT_DEV`` environmental
variable that points to the relevant directory.

Both tests collect the names of all files in the following two
directories::

    ${PYPEIT_DEV}/RAW_DATA/keck_deimos/830G_L_8100
    ${PYPEIT_DEV}/RAW_DATA/keck_deimos/830G_L_8400

test_setup_keck_deimos_multiconfig
----------------------------------

The algorithm for this test is as follows:

    1. Use :class:`~pypeit.pypeitsetup.PypeItSetup` to automatically
       identify the configurations for these files.

    2. Check that the code found two configurations and wrote the
       pypeit files for each.

    3. For each configuration:

        a. Read the pypeit file

        b. Check that the name for the setup is correct ('A' or 'B')

        c. Check that the calibration group is the same for all frames ('0' or '1')

test_setup_keck_deimos_multiconfig_clean
----------------------------------------

The algorithm for this test is as follows:

    1. Similar to the previous test, use
       :class:`~pypeit.pypeitsetup.PypeItSetup` to automatically
       identify the configurations for these files.

    2. Check that the code found two configurations.

    3. Check that the only bias frame in the list was correctly
       identified.

    4. Check that the full table containing both configurations
       contains the correct number of files (25).

    5. Check that "cleaning" the configurations of frames that cannot
       be reduced by ``PypeIt`` (those with ``MOSMODE != 'Spectral'``
       or ``AMPMODE != SINGLE:B`` or ``AMPMODE != SINGLE:A``), using
       :func:`~pypeit.metadata.PypeItMetaData.clean_configurations`
       does not remove any file because all of the dev-suite files
       are valid.

    6. After artificially changing the metadata for two files so that
       they *do not* satisfy the necessary metadata restrictions,
       check that
       :func:`~pypeit.metadata.PypeItMetaData.clean_configurations`
       removes both of these files.

Because these tests are now included in the ``PypeIt``
:ref:`unit-tests`, these configuration checks are performed by the
developers for every new version of the code.
