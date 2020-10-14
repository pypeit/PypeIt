.. include:: ../include/links.rst

.. _deimos_config_report:

Automated sorting of DEIMOS frames by instrument configuration
==============================================================

Version History
---------------

=========   =============   =========== ===========
*Version*   *Author*        *Date*      ``PypeIt``
=========   =============   =========== ===========
1.0         Kyle Westfall   13 Oct 2020 1.1.2dev
=========   =============   =========== ===========

----

Basics
------

Sorting frames by the configuration of the instrument, ensuring that
this configuration is the same for all coupled sets of calibration and
science frames, is performed by the :ref:`pypeit_setup` script; see
:doc:`setup`. :ref:`pypeit_setup` uses automated procedures to sort
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
===============     ====================================================================

as determined by
:func:`pypeit.metadata.PypeItMetaData.unique_configurations`. The
``AMPMODE`` value is included, even though ``PypeIt`` (currently)
restricts itself to only attempting to reduce frames read by the B
amplifier; see
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
to the first configuration with frames taken on the same date.

.. warning::

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


as demonstrated by the test at
``pypeit/tests/test_frametype.py``.  To run the test:

.. code-block:: bash

    cd pypeit/tests
    pytest test_frametype.py::test_deimos -W ignore

The test requires that you have downloaded the ``PypeIt``
:ref:`dev-suite` and defined the ``PYPEIT_DEV`` environmental
variable that points to the relevant directory. The algorithm of the
test is as follows:

    1. Find all the directories in the :ref:`dev-suite` with Keck
       DEIMOS data.

    2. For each directory (i.e., instrument setup):

        a. Make sure there is a "by-hand" version of the pypeit file
           for this setup where a human (one of the pypeit
           developers) has ensured the frame types are correct.

        b. Effectively run :ref:`pypeit_setup` on each of the
           instrument setups to construct a new pypeit file with the
           automatically generated frame types.
           
        c. Read both the by-hand and automatically generated frame
           types from these two pypeit files and check that they are
           identical. This check is *only* performed for the
           calibration frames, not any ``science`` or ``standard``
           frames.

Because this test is now included in the ``PypeIt``
:ref:`unit-tests`, this frame-typing check is performed by the
developers for every new version of the code.
