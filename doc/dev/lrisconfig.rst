.. include:: ../include/links.rst

.. _lris_config_report:

Automated sorting of LRIS frames by instrument configuration
==============================================================

Version History
---------------


=========   ================   =========== ===========
*Version*   *Author*           *Date*      ``PypeIt``
=========   ================   =========== ===========
1.0         Debora Pelliccia    6 Sep 2023 1.13.1.dev
=========   ================   =========== ===========

----

Basics
------

To prepare for the data reduction, PypeIt, first, automatically associates fits
files to specific :ref:`frame_types` (see :ref:`lris_frames_report`) and, then,
collects groups of frames in unique instrument configurations (see below). This is performed
by the :ref:`pypeit_setup` script, which sorts the frames and writes a
:ref:`pypeit_file` for each unique configuration. See :ref:`setup_doc`.


LRIS configuration identification
---------------------------------

The LRIS instrument configurations are determined, in the same
way for ``keck_lris_red``, ``keck_lris_red_orig``, ``keck_lris_red_mark4``,
``keck_lris_blue``, and ``keck_lris_blue_orig`` (unless otherwise noted),
by the function :func:`pypeit.metadata.PypeItMetaData.unique_configurations`,
which finds unique combinations of the following keywords:

===============     ====================================================================
``fitstbl`` key     Header Key
===============     ====================================================================
``dispname``        ``GRANAME`` (LRIS RED) or ``GRISNAME`` (LRIS BLUE)
``dichroic``        ``DICHNAME``
``decker``          ``SLITNAME``
``binning``         ``BINNING``
``amp``             ``NUMAMPS`` (``TAPLINES`` for ``keck_lris_red_mark4``)
``dispangle``       ``GRANGLE`` (only LRIS RED)
``cenwave``         ``WAVELEN`` (only LRIS RED)
===============     ====================================================================

The unique configurations are determined by collating the relevant metadata from the headers
of all frames found by a run of :ref:`pypeit_setup`, *except* those that are designated as
bias frames. The reason is that bias frames can have header data (e.g., ``dispangle``)
that do not match the instrument configuration that an observer intended for their use;
e.g., the frames were taken before the instrument was fully configured for the night's
observations. Therefore, PypeIt uses the ``dateobs``, ``binning``, ``amp`` keys to match
the bias frames to the configurations with frames taken on the same date, with
the same binning and on the same amplifier.

After that, :func:`pypeit.metadata.PypeItMetaData.set_configurations` associates each frame
to the relevant unique configuration ("setup"), by assigning a setup identifier
(e.g., A,B,C,D...) to every frames for which the values of the above keywords match the
values of the specific unique configuration.

LRIS calibration groups
-----------------------

PypeIt uses the concept of a "calibration group" to define a complete set of
calibration frames (e.g., arcs, flats) and the science frames to which these calibration
frames should be applied.

By default, :ref:`pypeit_setup` uses the setup identifier to assign frames to a single
calibration group. Frames that are in the same calibration group will have the same PypeIt
keyword ``calib``. No automated procedure exists to do anything except this.
However, the user can edit the :ref:`pypeit_file` to, within a given configuration, assign
specific calibration frames to specific science frames using the data in the ``calib`` column
of the :ref:`data_block`.

Testing
-------

Requirement PLL-17 states: "As a user, I expect the pipeline to automatically
and correctly associate calibrations with science frames."

PypeIt meets this requirement in the majority of use cases, as shown by the tests below.

Note that the tests described in :ref:`lris_frames_report` are also relevant here
since they show that PypeIt correctly identifies LRIS data and associates
them with a single configuration, all written to a single pypeit file.

To test that PypeIt can successfully identify multiple
configurations among a set of files, we have added five tests
``${PYPEIT_DEV}/unit_tests/test_setups.py``. They are:

- ``test_setup_keck_lris_blue_multiconfig()``,
- ``test_setup_keck_lris_blue_orig_multiconfig()``
- ``test_setup_keck_lris_red_multiconfig()``
- ``test_setup_keck_lris_red_orig_multiconfig()``
- ``test_setup_keck_lris_red_mark4_multiconfig()``

Here is an example of how to run the tests:

.. code-block:: bash

    cd ${PYPEIT_DEV}/unit_tests
    pytest test_setup.py::test_setup_keck_lris_blue_multiconfig -W ignore

The tests require that you have downloaded the PypeIt
:ref:`dev-suite` and defined the ``PYPEIT_DEV`` environmental
variable that points to the relevant directory.

The algorithm for all these tests is the same and is as follows:

    1. Collect the names of all files in selected LRIS directories
       (separately for ``keck_lris_blue``, ``keck_lris_blue_orig``, ``keck_lris_red``,
       ``keck_lris_red_orig``, ``keck_lris_red_mark4``).

    2. Use :class:`~pypeit.pypeitsetup.PypeItSetup` to automatically
       identify the configurations for these files.

    3. Check that the code found two configurations and wrote the
       pypeit files for each.

    4. For each configuration:

        a. Read the pypeit file

        b. Check that the name for the setup is correct ('A' or 'B')

        c. Check that the calibration group is the same for all frames ('0' or '1')


Because these tests are now included in the PypeIt
:ref:`unit-tests`, these configuration checks are performed by the
developers for every new version of the code.
