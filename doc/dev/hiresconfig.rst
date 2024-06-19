.. include:: ../include/links.rst

.. _hires_config:

Automated sorting of HIRES frames by instrument configuration
=============================================================

Version History
---------------


=========   ================   =========== ===========
*Version*   *Author*           *Date*      ``PypeIt``
=========   ================   =========== ===========
1.0         Debora Pelliccia   17 Jun 2024 1.15.1.dev
=========   ================   =========== ===========

----

Basics
------

To prepare for the data reduction, PypeIt, first, automatically associates fits
files to specific :ref:`frame_types` (see :ref:`hires_frames`) and, then,
collects groups of frames in unique instrument configurations (see below). This is performed
by the :ref:`pypeit_setup` script, which sorts the frames and writes a
:ref:`pypeit_file` for each unique configuration. See :ref:`setup_doc`.


HIRES configuration identification
---------------------------------

The HIRES instrument configurations are determined by the function
:func:`pypeit.metadata.PypeItMetaData.unique_configurations`,
which finds unique combinations of the following keywords:

===============     ============
``fitstbl`` key     Header Key
===============     ============
``dispname``        ``XDISPERS``
``decker``          ``DECKNAME``
``binning``         ``BINNING``
``filter1``         ``FIL1NAME``
``echangle``        ``ECHANGL``
``xdangle``         ``XDANGL``
===============     ============

The unique configurations are determined by collating the relevant metadata from the headers
of all frames found by a run of :ref:`pypeit_setup`, *except* those that are designated as
bias and slitless_pixflat frames. Bias and slitless_pixflat frames can have header data (e.g., ``filter1``)
that do not match the instrument configuration that an observer intended for their use.
Therefore, PypeIt uses the ``dispname`` and ``binning`` keys to match the bias and
slitless_pixflat frames to the configurations with frames taken with the same cross-disperser
and same binning.

After that, :func:`pypeit.metadata.PypeItMetaData.set_configurations` associates each frame
to the relevant unique configuration ("setup"), by assigning a setup identifier
(e.g., A,B,C,D...) to every frames for which the values of the above keywords match the
values of the specific unique configuration.

HIRES calibration groups
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

To test that PypeIt can successfully identify multiple
configurations among a set of files, we have added the
``test_setup_keck_hires_multiconfig()`` test to
``${PYPEIT_DEV}/unit_tests/test_setups.py``.

Here is an example of how to run the test:

.. code-block:: bash

    cd ${PYPEIT_DEV}/unit_tests
    pytest test_setup.py::test_setup_keck_hires_multiconfig -W ignore

The tests require that you have downloaded the PypeIt
:ref:`dev-suite` and defined the ``PYPEIT_DEV`` environmental
variable that points to the relevant directory.

The algorithm for this test is as follows:

    1. Collect the names of all files in selected HIRES directories.

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
