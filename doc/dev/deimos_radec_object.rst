.. include:: ../include/links.rst

.. _deimos_radec_object_report:

RA, Dec and object name assignment to DEIMOS extracted spectra
==============================================================

Version History
---------------

=========   ================   =========== ===========
*Version*   *Author*           *Date*      ``PypeIt``
=========   ================   =========== ===========
1.0         Debora Pelliccia   25 Jan 2021  1.3.1dev
1.1         Debora Pelliccia   02 Apr 2021  1.3.4dev
=========   ================   =========== ===========

----

Basics
------

The procedure used to assign RA, Dec and object name to each DEIMOS 1D extracted spectrum
is performed right after the object finding procedure described in :ref:`object_finding`.


Procedure
---------

RA, Dec and object name assignment is primarily performed by
:func:`pypeit.slittrace.SlitTraceSet.assign_maskinfo`. This function, first, reads in the slitmask
design information stored in the :class:`~pypeit.slittrace.SlitTraceSet` datamodel (see
:ref:`deimos_slitmask_ids_report` for a description on how the slitmask design matching is performed,
and :ref:`master_slits` for a description of the provided information and for a way to visualize them).

Then, :func:`~pypeit.slittrace.SlitTraceSet.assign_maskinfo` goes through all the slits in a
selected detector and for each detected object checks if the measured distance of the object from the
left edge of the slit is within a certain tolerance (see `Application`_ for details on how to control the value of
this parameter) of the distance expected from the slitmask design (differences between the expected and
the measured slit length are taken into account).
If the answer is yes, the ``RA``, ``DEC`` and ``MASKDEF_OBJNAME`` of the detected object are updated
with the coordinates and name of the targeted object. If the answer is no, the detected object is
considered a serendipitous object. Using the coordinates of the slit center available from the slitmask
design and the distance in pixels (converted then in arcsec) between the traced object and the center
of the slit, the coordinates of the serendipitous object are estimates and recorded in  ``RA``
and ``DEC`` while ``MASKDEF_OBJNAME`` is set to "SERENDIP". For both cases, the ``MASKDEF_EXTRACT`` attribute
is set to **False** (see :ref:`deimos_add_missing_obj_report`).


Application
-----------

To perform the RA, Dec and object name assignment to DEIMOS extracted spectra, the parameters described in
the *Application* section of :ref:`deimos_slitmask_ids_report` must be set. Moreover, the **assign_obj** flag in
:ref:`pypeit_par:SlitMaskPar Keywords` must be **True**.  This is the default for DEIMOS,
except when the *LongMirr* or the *LVM* mask is used. One other keyword controls this procedure and it is **obj_toler**.
This keyword sets the tolerance in arcsec for the matching process between
the measured coordinates of the extracted spectrum and the expected coordinates of the targeted object.
The default value is **obj_toler = 5**.


Access
------

``PypeIt`` users can access the DEIMOS objects information in several ways.

- Ra, Dec, object name and ``maskdef_id`` are visible in the .txt file with a list of all extracted spectra,
  generated at the end the ``PypeIt`` reduction.
- Ra, Dec, object name are also visible by running `pypeit_show_1d --list` (see :ref:`out_spec1D:pypeit_show_1dspec`)
- Object names are visible in `ginga` when running `pypeit_show_2d` (see :ref:`out_spec2D:pypeit_show_2dspec`)


Testing
-------

Requirement PD-8 states: "As a user, I expect products that are associated with my mask
definition (object names, object positions)."

``PypeIt`` meets this requirement as demonstrated by the test at ``pypeit/tests/test_slitmask.py``.
To run the test:

.. code-block:: bash

    cd pypeit/tests
    pytest test_slitmask::test_assign_maskinfo -W ignore

The tests require that you have downloaded the ``PypeIt`` :ref:`dev-suite` (including the folder compressed in 
"Cooked_pypeit_dev_vX.XX.X.tar.gz", where vX.XX.X indicates the version of the file) and defined
the ``PYPEIT_DEV`` environmental variable that points to the relevant directory. These tests are
run using only one instrument setup and only one DEIMOS detector.

The algorithm of the test is as follows:

    1. Load the information relative to the specific instrument (DEIMOS).

    2. Load the :ref:`pypeit_par:Instrument-Specific Default Configuration` parameters and select the detector.

    3. Build a trace image using three flat-field images from a specific DEIMOS dataset in the :ref:`dev-suite`.

    4. Update the DEIMOS configuration parameters to include configurations specific for the
       used instrument setup. Among others, this step sets the **assign_obj** flag in
       :ref:`pypeit_par:SlitMaskPar Keywords` to **True**.

    5. Run the slit tracing procedure using :class:`~pypeit.edgetrace.EdgeTraceSet`, during which
       the slitmask ID assignment is performed (see :ref:`deimos_slitmask_ids_report`), and the ``maskdef_id``
       and object information associated to each slit are recorded in the
       :class:`~pypeit.slittrace.SlitTraceSet` datamodel.

    6. A file containing previously extracted 1D spectra is loaded from the **Cooked** folder in the :ref:`dev-suite` 
       and the objects information re-initialized, i.e.,  ``RA``, ``DEC`` and ``MASKDEF_OBJNAME`` are set to ``None``
       for each spectrum.

    7. :func:`~pypeit.slittrace.SlitTraceSet.assign_maskinfo` is then run and the object RA, Dec and object
       name are assigned to each extracted spectrum.

    8. Read  ``RA``, ``DEC`` and ``MASKDEF_OBJNAME`` for a selected slit and check if those correspond to
       the expected values. The expected values are taken from the :class:`~pypeit.slittrace.SlitTraceSet`
       datamodel. See :ref:`master_slits` for a description of the provided information and for a way
       to visualize them.


Because these tests are now included in the ``PypeIt`` :ref:`unit-tests`, this check is performed by the
developers for every new version of the code.