.. include:: ../include/links.rst

.. _radec_object_report:

RA, Dec and object name assignment to 1D extracted spectra
==========================================================

Version History
---------------

=========   ================   =========== ===========
*Version*   *Author*           *Date*      ``PypeIt``
=========   ================   =========== ===========
1.0         Debora Pelliccia   25 Jan 2021  1.3.1dev
1.1         Debora Pelliccia   02 Apr 2021  1.3.4dev
1.2         Debora Pelliccia   28 Jul 2021  1.4.3dev
1.3         Debora Pelliccia   21 Oct 2021  1.6.1dev
1.4         Debora Pelliccia    6 Sep 2023  1.13.1dev
=========   ================   =========== ===========

----

Basics
------

The procedure used to assign RA, Dec and object name to each 1D extracted spectrum
is currently available for these :ref:`slitmask_info_instruments` only and is performed right after
the object finding procedure described in :ref:`object_finding`.


Procedure
---------

RA, Dec and object name assignment is primarily performed by
:func:`pypeit.slittrace.assign_addobjs_alldets`.

First, ``Pypeit``, using :func:`pypeit.slittrace.get_maskdef_objpos_offset_alldets`,
computes for every detectors that the user wants to reduce the offset of the observed
slitmask from the position expected by the design file (see :ref:`slitmask_ids_report`
for a description on how the slitmask design matching is performed), and stores it in
the :class:`~pypeit.slittrace.SlitTraceSet` datamodel (see :ref:`slits` for a
description of the provided information and for a way to visualize them). There are four
different options that the user can choose to compute the offset: using objects with high
significance detections, using only one bright object in a selected slit, using alignment stars,
or using the dither offset recorded in the science frame header (available only for MOSFIRE).
Or, alternatively, the user can input a known offset in pixels (see `Application`_).

Then, :func:`pypeit.slittrace.average_maskdef_offset` determines an average slitmask offset over all
the detectors (this does not happen for MOSFIRE observations, which have only one detector), which is
used by :func:`pypeit.slittrace.assign_addobjs_alldets` to assign RA, Dec and object name to detected
objects and to force extract objects that were not detected (for the latter see
:ref:`add_missing_obj_report`).

The function :func:`~pypeit.slittrace.assign_addobjs_alldets`, for each detector, goes through
all the slits and checks if the measured distance of the detected objects from the left edge
of the slit (corrected for the offset computed in the previous step) is within a certain tolerance
(see `Application`_ for details on how to control the value of this parameter) of the distance
expected from the slitmask design (differences between the expected and the measured slit length
are taken into account). Correcting for the slitmask offset allows ``PypeIt`` to deal also with
dithered observations.
If the measured distance is within the set tolerance, the ``RA``, ``DEC`` and ``MASKDEF_OBJNAME``
of the detected object are updated with the coordinates and name of the targeted object.
If the measured distance is not within the tolerance, the detected object is considered a
serendipitous object.
Using the coordinates of the slit center available from the slitmask design and the distance in pixels
(converted then in arcsec) between the traced object and the center of the slit, the coordinates
of the serendipitous object are estimates and recorded in  ``RA`` and ``DEC`` while
``MASKDEF_OBJNAME`` is set to "SERENDIP". For both cases, the ``MASKDEF_EXTRACT`` attribute
is set to **False** (see :ref:`add_missing_obj_report`).


Application
-----------

To perform the RA, Dec and object name assignment to extracted spectra, the parameters
described in the *Application* section of :ref:`slitmask_ids_report` must be set.
Moreover, the **assign_obj** flag in :ref:`slitmaskpar` must be **True**.
This is the default for DEIMOS (except when the *LongMirr* or the *LVM* mask is used) and
MOSFIRE (except when the *LONGSLIT* mask is used).
Seven other parameters control this procedure. Six are for the slitmask offset determination
and one is for the RA, Dec and object name assignment. They are the following.

- **snr_thrshd**: objects detected above this S/N threshold are used to
  compute the slitmask offset. This is the default behaviour for DEIMOS unless **slitmask_offset**,
  **bright_maskdef_id** or **use_alignbox** is set. Default value is **snr_thrshd=50**.

- **bright_maskdef_id**: ``maskdef_id`` (corresponding to ``dSlitId`` and ``Slit_Number`` in the DEIMOS/LRIS
  and MOSFIRE slitmask design, respectively) of a slit containing a bright object that will be used
  to compute the slitmask offset. This parameter is optional (default value is **bright_maskdef_id=None**)
  and is ignored if **slitmask_offset** is provided or **use_dither_offset = True**. However,
  this parameter is highly recommended for MOSFIRE observations if a bright object is present in
  the slitmask, since it allows to trace the small drifts of the objects position that have been
  typically seen in MOSFIRE data.

- **use_alignbox**: flag to use stars in alignment boxes to compute the slitmask offset. This is
  available only for DEIMOS observations and it is set as the default for this instrument.
  If **use_alignbox = True** PypeIt will NOT compute the offset using **snr_thrshd** or
  **bright_maskdef_id**.

- **use_dither_offset**: flag to use the dither offset recorded in the header of science frames as the
  value of the slitmask offset. This is currently only available for MOSFIRE reduction and
  it is set as the default for this instrument. If **use_dither_offset = True** PypeIt will NOT compute the
  offset using `snr_thrshd` or `bright_maskdef_id`. However, it is ignored if ``slitmask_offset`` is provided.

- **slitmask_offset**: user-provided slitmask offset (pixels) from the position expected
  by the slitmask design. This is optional (default value is **slitmask_offset=None**),
  and if set PypeIt will NOT compute the offset, i.e., the above parameters will be ignored.

- **obj_toler**: sets the tolerance in arcsec for the matching process between the measured
  coordinates of the extracted spectrum and the expected coordinates of the targeted object.
  The default value is **obj_toler = 1**.


Access
------

``PypeIt`` users can access the objects information in several ways.

- Ra, Dec, object name and ``maskdef_id`` are visible in the .txt file with a list of all extracted spectra,
  generated at the end the ``PypeIt`` reduction.
- Ra, Dec, object name are also visible by running `pypeit_show_1d --list` (see :ref:`pypeit_show_1dspec`)
- Object names are visible in `ginga` when running `pypeit_show_2d` (see :ref:`pypeit_show_2dspec`)


Testing
-------
- PYPEIT-DEIMOS (PD)

    Requirement PD-8 states: "As a user, I expect products that are associated with my mask
    definition (object names, object positions)."

    Requirement PD-35 states: "As a user, I expect to be able to reduce data taken with a nodding pattern."

- PYPEIT-MOSFIRE (PM)

    Requirement PM-9 states: "As a user, I expect products that are associated with my mask
    definition (object names, object positions)."

    Requirement PM-12 states: "Use slitmask information to determine object position in the slit."

``PypeIt`` meets this requirement as demonstrated by the tests at ``pypeit/tests/test_slitmask.py``.
To run the test:

.. code-block:: bash

    cd pypeit/tests
    pytest test_slitmask.py::test_assign_maskinfo_add_missing -W ignore
    pytest test_slitmask.py::test_dith_obs -W ignore

The tests require that you have downloaded the ``PypeIt`` :ref:`dev-suite` (including the folder compressed in 
"Cooked_pypeit_dev_vX.XX.X.tar.gz", where vX.XX.X indicates the version of the file) and defined
the ``PYPEIT_DEV`` environmental variable that points to the relevant directory. These tests are
run using only one instrument setup and for DEIMOS only one detector.

First test ``pytest test_slitmask::test_assign_maskinfo -W ignore``.
The algorithm of the test is repeated twice (once for a DEIMOS dataset and once for a MOSFIRE dataset)
and is as follows:

    1. Load the information relative to the specific instrument (DEIMOS, MOSFIRE).

    2. Load the :ref:`instr_par` parameters and select the detector.

    3. Build a trace image using three flat-field images from a specific dataset in the :ref:`dev-suite`.

    4. Update the instrument configuration parameters to include configurations specific for the
       used instrument setup. Among others, this step sets the **assign_obj** flag in
       :ref:`slitmaskpar` to **True**.

    5. Run the slit tracing procedure using :class:`~pypeit.edgetrace.EdgeTraceSet`, during which
       the slitmask ID assignment is performed (see :ref:`slitmask_ids_report`), and the ``maskdef_id``
       and object information associated to each slit are recorded in the
       :class:`~pypeit.slittrace.SlitTraceSet` datamodel.

    6. A file containing previously extracted 1D spectra is loaded from the **Cooked** folder in the :ref:`dev-suite` 
       and the objects information re-initialized, i.e.,  ``RA``, ``DEC`` and ``MASKDEF_OBJNAME`` are set to ``None``
       for each spectrum.

    7. :func:`pypeit.slittrace.average_maskdef_offset` is run to determine an average slitmask offset over
       all the detectors, which is used by :func:`pypeit.slittrace.assign_addobjs_alldets` to assign
       the object RA, Dec and object name to each extracted spectrum.

    8. Read  ``RA``, ``DEC`` and ``MASKDEF_OBJNAME`` for a selected slit and check if those correspond to
       the expected values. The expected values are taken from the :class:`~pypeit.slittrace.SlitTraceSet`
       datamodel. See :ref:`slits` for a description of the provided information and for a way
       to visualize them.


Second test ``pytest test_slitmask::test_dith_obs -W ignore`` is performed only on DEIMOS data.
The algorithm of this test is identical to the first test, but using a different DEIMOS dataset obtained with
dithered observations.

Because these tests are now included in the ``PypeIt`` :ref:`unit-tests`, this check is performed by the
developers for every new version of the code.