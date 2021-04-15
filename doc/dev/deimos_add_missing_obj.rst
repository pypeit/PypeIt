.. include:: ../include/links.rst

.. _deimos_add_missing_obj_report:

DEIMOS object position on the slit from slitmask design
=======================================================

Version History
---------------

=========   ================   =========== ===========
*Version*   *Author*           *Date*      ``PypeIt``
=========   ================   =========== ===========
1.0         Debora Pelliccia   02 Apr 2021  1.3.4dev
=========   ================   =========== ===========

----

Basics
------

The procedure to determine the location on the slit of DEIMOS undetected objects
using the slitmask design information is performed right after the object finding
(see :ref:`object_finding`) and the RA, Dec and object name assignment procedures
(see :ref:`deimos_radec_object_report`) have been completed.


Procedure
---------

The determination of the position on the slit of DEIMOS non detected objects is primarily performed by
:func:`pypeit.slittrace.SlitTraceSet.mask_add_missing_obj`. This function relies on the output
from :func:`pypeit.slittrace.SlitTraceSet.assign_maskinfo` which assign RA, Dec and object name
to the detected objects and provides an array of all the expected position corrected for the differences
between the expected and measured slit length (see :ref:`deimos_radec_object_report`).

``PypeIt`` goes through all the slits in the selected detector and for each slit checks if the
target object was detected (this is done by checking if ``MASKDEF_OBJNAME`` corresponds to the object name
of the target). If the answer is yes, it goes to the next slit. If the answer is no, a new
:class:`pypeit.specobj.SpecObj` is added to the :class:`pypeit.specobjs.SpecObjs` class. The expected
position is provided as an output of :func:`~pypeit.slittrace.SlitTraceSet.assign_maskinfo` and is recorded
in the :class:`~pypeit.specobjs.SpecObjs`'s attribute ``SPAT_PIXPOS``. The user can provide an additional
offset between expected and measured position if needed (see `Application`_ for details on how to control
the value of this parameter). Other relevant attributes are also updated, i.e., ``TRACE_SPAT``, ``SPAT_FRACPOS``,
``OBJID``, ``FWHM``, ``RA``, ``DEC``, ``MASKDEF_OBJNAME``, ``MASKDEF_ID`` (see spec1D
:ref:`out_spec1D:Current Data Model` for a description of these parameters). The attribute ``MASKDEF_EXTRACT``
is set to **True** to flag the spectra that have been extracted from undetected objects.


Application
-----------

To perform the determination of the location on the slit of undetected objects, the parameters described in
the *Application* section of :ref:`deimos_slitmask_ids_report` and :ref:`deimos_radec_object_report` must be set.
Moreover, **extract_missing_objs** flag in :ref:`pypeit_par:SlitMaskPar Keywords` must be **True**.  This is the
default for DEIMOS, except when the *LongMirr* or the *LVM* mask is used. One other keyword control this procedure.
It is **slitmask_offset**, which sets a user provided offset in pixels between the measured and expected
position of the slitmask. The default is zero.

See :ref:`pypeit_par:SlitMaskPar Keywords` for more details.

Access
------

``PypeIt`` users can access the DEIMOS objects information resulted from this procedure in several ways.

- A flag that identifies the undetected objects for which the extraction was "forced" is visible in
  the .txt file with a list of all extracted spectra, generated at the end the ``PypeIt`` reduction.
- The same flag is also visible when running `pypeit_show_1d --list` (see :ref:`out_spec1D:pypeit_show_1dspec`)
- The forced extraction are shown in a different color than the detected objects (yellow vs. orange)
  in `ginga` when running `pypeit_show_2d` (see :ref:`out_spec2D:pypeit_show_2dspec`)
- the **slitmask_offset** value is reported when running ``pypeit_chk_2dslits Science/spec2d_XXX.fits``.


Testing
-------

Requirement PD-11 states: "Use slitmask information to determine object position in the slit."

Requirement PD-12 states: "As a user, I expect the pipeline to extract 1d spectra using standard extraction
method such as optimal extraction, even for sources that don’t have continuum using information contained
in the mask definition"

``PypeIt`` meets these requirements as demonstrated by the test at ``pypeit/tests/test_slitmask.py``.
To run the test:

.. code-block:: bash

    cd pypeit/tests
    pytest test_slitmask::test_add_missing_obj -W ignore

The tests require that you have downloaded the ``PypeIt`` :ref:`dev-suite` (including the folder compressed in 
"Cooked_pypeit_dev_vX.XX.X.tar.gz", where vX.XX.X indicates the version of the file) and defined
the ``PYPEIT_DEV`` environmental variable that points to the relevant directory. These tests are
run using only one instrument setup and only one DEIMOS detector.

The algorithm of the test is as follows:

    1. Load the information relative to the specific instrument (DEIMOS).

    2. Load the :ref:`pypeit_par:Instrument-Specific Default Configuration` parameters and select the detector.

    3. Build a trace image using three flat-field images from a specific DEIMOS dataset in the :ref:`dev-suite`.

    4. Update the DEIMOS configuration parameters to include configurations specific for the
       used instrument setup. Among others, this step sets the **extract_missing_objs** flag in
       :ref:`pypeit_par:SlitMaskPar Keywords` to **True**.

    5. Run the slit tracing procedure using :class:`~pypeit.edgetrace.EdgeTraceSet`, during which
       the slitmask ID assignment is performed (see :ref:`deimos_slitmask_ids_report`), and the ``maskdef_id``
       and object information associated to each slit are recorded in the :class:`~pypeit.slittrace.SlitTraceSet`
       datamodel.

    6. A file containing previously extracted 1D spectra is loaded from the **Cooked** folder in the :ref:`dev-suite` 
       and the objects information re-initialized, i.e.,  ``RA``, ``DEC``, ``MASKDEF_OBJNAME`` and ``MASKDEF_EXTRACT``
       are set to ``None`` for the detected objects, while the spectra of undetected object are removed.

    7. :func:`~pypeit.slittrace.SlitTraceSet.assign_maskinfo` is then run and the object RA, Dec and object
       name are assigned to each detected objects.

    8. :func:`~pypeit.slittrace.SlitTraceSet.mask_add_missing_obj` is also run and adds a new
       :class:`~pypeit.specobj.SpecObj` to the :class:`~pypeit.specobjs.SpecObjs` class for each undetected
       undetected object, updating all the relevant attributes (see `Procedure`_).

    9. Read ``SPAT_PIXPOS`` for two undetected objects and check if those correspond to
       the expected position of the object on the slit. The expected positions are verified by visual inspection
       of the 2D spectrum.


Because these tests are now included in the ``PypeIt`` :ref:`unit-tests`, this check is performed by the
developers for every new version of the code.