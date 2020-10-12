.. include:: ../include/links.rst

.. _deimos_slitmask_ids_report:

DEIMOS slitmask ID assignment
=============================

Version History
---------------

=========   ================   =========== ===========
*Version*   *Author*           *Date*      ``PypeIt``
=========   ================   =========== ===========
1.0         Debora Pelliccia   9 Oct 2020  1.1.2dev
=========   ================   =========== ===========

----

Basics
------

The procedure used to assign DEIMOS slitmask ID to each slit 
is part of the more general slit tracing procedure described in :doc:`slit_tracing`.

Slitmask ID assignment
----------------------

DEIMOS slitmask ID assignment is primarily performed by
:func:`pypeit.edgetrace.EdgeTraceSet.maskdesign_matching`. This function matches the slit 
edges traced by ``PypeIt`` to the slit edges predicted by the DEIMOS optical model. These 
predictions are generated, using the slitmask design information stored in the DEIMOS data, 
by the function :func:`pypeit.spectrographs.keck_deimos.KeckDEIMOSSpectrograph.mask_to_pixel_coordinates`,
which relies on the existence of pre-generated detectors maps pre- and post-grating,
called `amap` and `bmap`, respectively.

The function :func:`~pypeit.edgetrace.EdgeTraceSet.maskdesign_matching` assigns to each slit 
a `maskdef_id`, which corresponds to `dSlitId` in the DEIMOS slitmask design, and records it
in the :class:`~pypeit.slittrace.SlitTraceSet` datamodel.


Application
-----------

To perform the slitmask ID assignment, the **use_maskdesign** flag in :ref:`pypeit_par:EdgeTracePar Keywords`
must be *True*.  This is the default for DEIMOS, except when *LongMirr* mask is used.

Three other :ref:`pypeit_par:EdgeTracePar Keywords` control the slitmask ID assignment;
these are: **maskdesign_maxsep**, **maskdesign_sigrej**, **maskdesign_step**.


Testing
-------

Requirement PD-9 states: "Ingest slitmask information into a data structure."

``PypeIt`` meets this requirement as demonstrated by the test at
``pypeit/tests/test_maskdesign_matching.py``.  To run the test:

.. code-block:: bash

    cd pypeit/tests
    pytest test_maskdesign_matching.py::test_maskdef_id -W ignore

The test requires that you have downloaded the ``PypeIt`` :ref:`dev-suite` and defined 
the ``PYPEIT_DEV`` environmental variable that points to the relevant directory. This test is 
run using only one instrument setup and only one DEIMOS detector.

The algorithm of the test is as follows:

    1. Load the information relative to the specific instrument (DEIMOS).

    2. Load the :ref:`pypeit_par:Instrument-Specific Default Configuration` parameters and select the detector.

    3. Build a trace image using three flat-field images from a specific DEIMOS dataset in the :ref:`dev-suite`.

    4. Update the DEIMOS configuration parameters to include configurations specific for the
       used instrument setup. Among others, this step sets the **use_maskdesign** flag in 
       :ref:`pypeit_par:EdgeTracePar Keywords` to *True*.

    5. Run the slit tracing procedure using :class:`~pypeit.edgetrace.EdgeTraceSet`, during which
       the slitmask ID assignment is performed, and record the `maskdef_id` associated to each slit
       in the :class:`~pypeit.slittrace.SlitTraceSet` datamodel.

    6. Read the `maskdef_id` for the first and last slits of the selected detector and check if
       those correspond to the expected values. The expected values are determined by previously 
       reducing the same dataset with the already existing *DEEP2 IDL-based pipeline*, which is 
       optimized to match the slits found on the DEIMOS detectors to the slitmask information.

Because this test is now included in the ``PypeIt`` :ref:`unit-tests`, this check is performed by the
developers for every new version of the code.
