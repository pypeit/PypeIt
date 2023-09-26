.. include:: ../include/links.rst

.. _add_missing_obj_report:

Object position on the slit from slitmask design
=======================================================

Version History
---------------

=========   =================== =========== ===========
*Version*   *Author*            *Date*      ``PypeIt``
=========   =================== =========== ===========
1.0         Debora Pelliccia    02 Apr 2021  1.3.4dev
1.1         Debora Pelliccia    28 Jul 2021  1.4.3dev
1.3         Debora Pelliccia    21 Oct 2021  1.6.1dev
1.4         J. Xavier Prochaska 25 Jan 2022  1.7.1dev
1.5         Debora Pelliccia     6 Sep 2023  1.13.1dev
=========   =================== =========== ===========

----

Basics
------

The procedure to determine the location on the slit of undetected objects
using the slitmask design information is currently available for these :ref:`slitmask_info_instruments` only
and it is performed right after the object finding (see :ref:`object_finding`)
and the RA, Dec and object name assignment procedures (see :ref:`radec_object_report`) have been completed.


Procedure
---------

The determination of the position on the slit of non detected objects is primarily performed by
:func:`pypeit.slittrace.assign_addobjs_alldets`. This function, first, assigns RA, Dec and object name
to the detected objects (see :ref:`radec_object_report`).

Then, ``PypeIt``, for each detector, goes through all the slits and for each slit checks if the
target object was detected (this is done by checking if ``MASKDEF_OBJNAME`` corresponds to the object
name of the target). If the answer is yes, it goes to the next slit. If the answer is no, a new
:class:`pypeit.specobj.SpecObj` is added to the :class:`pypeit.specobjs.SpecObjs` class. The expected
position on the slit (corrected for **slitmask_offset**, see :ref:`radec_object_report` for how
to compute this offset) is recorded in the :class:`~pypeit.specobjs.SpecObjs`'s attribute
``SPAT_PIXPOS``. Correcting for **slitmask_offset** allows ``PypeIt`` to deal also with dithered
observations. Other relevant attributes are also updated, i.e., ``TRACE_SPAT``, ``SPAT_FRACPOS``,
``OBJID``, ``FWHM``, ``RA``, ``DEC``, ``MASKDEF_OBJNAME``, ``MASKDEF_ID`` (see spec1D
:ref:`spec1D-datamodel` for a description of these parameters).
The attribute ``MASKDEF_EXTRACT`` is set to **True** to flag the spectra that have been extracted
from undetected objects.


Application
-----------

To perform the determination of the location on the slit of undetected objects, the parameters
described in the *Application* section of :ref:`slitmask_ids_report` and
:ref:`radec_object_report` must be set.
Moreover, **extract_missing_objs** flag in :ref:`slitmaskpar` must be **True**.
This is the default for DEIMOS (except when the *LongMirr* or the *LVM* mask is used) and
MOSFIRE (except when the *LONGSLIT* or the *long2pos* mask is used).

See :ref:`slitmaskpar` for more details.

For LRIS one needs to add these explicitly to the :ref:`pypeit_file`, e.g.::

    [reduce]
    [[slitmask]]
        use_alignbox = False
        assign_obj = True
        extract_missing_objs = True

For faint alignment stars, one may wish to set `use_alignbox = True`.

Access
------

``PypeIt`` users can access the objects information resulted from this procedure in several ways.

- A flag that identifies the undetected objects for which the extraction was "forced" is visible in
  the .txt file with a list of all extracted spectra, generated at the end the ``PypeIt`` reduction.
- The same flag is also visible when running `pypeit_show_1d --list` (see :ref:`pypeit_show_1dspec`)
- The forced extraction are shown in a different color than the detected objects (yellow vs. orange)
  in `ginga` when running `pypeit_show_2d` (see :ref:`pypeit_show_2dspec`)
- the **slitmask_offset** value is reported when running ``pypeit_parse_slits Science/spec2d_XXX.fits``.


Testing
-------

- PYPEIT-DEIMOS (PD)

    Requirement PD-11 states: "Use slitmask information to determine object position in the slit."

    Requirement PD-12 states: "As a user, I expect the pipeline to extract 1d spectra using standard extraction
    method such as optimal extraction, even for sources that don’t have continuum using information contained
    in the mask definition"

- PYPEIT-MOSFIRE (PM)

    Requirement PM-9 states: "As a user, I expect products that are associated with my mask
    definition (object names, object positions)."

    Requirement PM-12 states: "Use slitmask information to determine object position in the slit."

    Requirement PM-14 states: "As a user, I expect the pipeline to extract 1d spectra using standard extraction
    method such as optimal extraction, even for sources that don’t have continuum using information contained
    in the mask definition"

    Requirement PM-16 states: "Use slitmask information to determine the location of the object in the slit."

``PypeIt`` meets these requirements as demonstrated by the test at ``pypeit/tests/test_slitmask.py``.
To run the test:

.. code-block:: bash

    cd pypeit/tests
    pytest test_slitmask.py::test_assign_maskinfo_add_missing -W ignore

The test requires that you have downloaded the ``PypeIt`` :ref:`dev-suite` (including the folder compressed in
"Cooked_pypeit_dev_vX.XX.X.tar.gz", where vX.XX.X indicates the version of the file) and defined
the ``PYPEIT_DEV`` environmental variable that points to the relevant directory. This test is
run using only one instrument setup and for DEIMOS only one detector.

The algorithm of the test is repeated twice (once for a DEIMOS dataset and once for a MOSFIRE dataset)
and is as follows:

    1. Load the information relative to the specific instrument (DEIMOS, MOSFIRE).

    2. Load the :ref:`instr_par` parameters and select the detector.

    3. Build a trace image using three flat-field images from a specific dataset in the :ref:`dev-suite`.

    4. Update the instrument configuration parameters to include configurations specific for the
       used instrument setup. Among others, this step sets the **extract_missing_objs** flag in
       :ref:`slitmaskpar` to **True**.

    5. Run the slit tracing procedure using :class:`~pypeit.edgetrace.EdgeTraceSet`, during which
       the slitmask ID assignment is performed (see :ref:`slitmask_ids_report`), and the ``maskdef_id``
       and object information associated to each slit are recorded in the :class:`~pypeit.slittrace.SlitTraceSet`
       datamodel.

    6. A file containing previously extracted 1D spectra is loaded from the **Cooked** folder in the :ref:`dev-suite` 
       and the objects information re-initialized, i.e.,  ``RA``, ``DEC``, ``MASKDEF_OBJNAME`` and ``MASKDEF_EXTRACT``
       are set to ``None`` for the detected objects, while the spectra of undetected object are removed.

    7. :func:`pypeit.slittrace.average_maskdef_offset` is run to determine an average slitmask offset over
       all the detectors, which is used by :func:`pypeit.slittrace.assign_addobjs_alldets` to assign
       the object RA, Dec and object name to each extracted spectrum and to add a new
       :class:`~pypeit.specobj.SpecObj` to the :class:`~pypeit.specobjs.SpecObjs` class for each
       undetected object, updating all the relevant attributes (see `Procedure`_).

    9. Read ``SPAT_PIXPOS`` for two undetected objects and check if those correspond to
       the expected position of the object on the slit. The expected positions are verified by visual inspection
       of the 2D spectrum.


Because these tests are now included in the ``PypeIt`` :ref:`unit-tests`, this check is performed by the
developers for every new version of the code.