.. include:: ../include/links.rst

.. _nires_frames_report:

Automated typing of NIRES frames
================================

Version History
---------------

=========   =============   =========== ===========
*Version*   *Author*        *Date*      ``PypeIt``
=========   =============   =========== ===========
1.0         Kyle Westfall    9 Aug 2022 1.10.1dev
=========   =============   =========== ===========

----

Basics
------

The general procedure used to assign frames a given type is described
here: :ref:`frame_types`.

NIRES frame typing
------------------

The primary typing of NIRES frames is performed by
:func:`pypeit.spectrographs.keck_nires.KeckNIRESSpectrograph.check_frame_type`.
This function checks the values of various header keywords against a
set of criteria used to classify the frame type. The headers cards
required for the frame-typing and their associated keyword in the
:class:`~pypeit.metadata.PypeItMetaData` object are:

===============     ============
``fitstbl`` key     Header Key
===============     ============
``idname``          ``OBSTYPE``
===============     ============

The criteria used to select each frame type are as follows:

=============   ==========================================
Frame           ``OBSTYPE``                               
=============   ==========================================
``science``     ``'Object'``, ``'object'``                             
``dark``        ``'dark'``                                
``pixelflat``   ``'domeflat'``                             
``trace``       ``'domeflat'``                             
``arc``         ``'Object'``, ``'object'``                             
``tilt``        ``'Object'``, ``'object'``                             
=============   ==========================================

.. note::

    - Frames with type ``bias``, ``illumflat``, and ``pinhole`` are *not* typed.

    - By default, the exposure time (``ELAPTIME``) is not used to distinguish
      frame types; however, this can be changed using the ``exprng`` parameter
      in the :ref:`pypeit_file`; see also :ref:`frame_types`.

Testing
-------

Requirement PN-14 states: "As a user, I expect the pipeline to automatically
recognize calibration files."

``PypeIt`` meets this requirement as demonstrated by the test at
``PypeIt-development-suite/unit_tests/test_frametype.py``.  To run the test:

.. code-block:: bash

    cd PypeIt-development-suite/unit_tests
    pytest test_frametype.py::test_nires -W ignore

The test requires that you have downloaded the ``PypeIt``
:ref:`dev-suite` and defined the ``PYPEIT_DEV`` environmental
variable that points to the relevant directory. The algorithm of the
test is as follows:

    #. Find all the directories in the :ref:`dev-suite` with Keck
       NIRES data.

    #. For each directory (i.e., instrument setup; currently there is only one):

        #. Make sure there is a "by-hand" version of the pypeit file for this
           setup where a human (one of the pypeit developers) has ensured the
           frame types are correct.

        #. Effectively run :ref:`pypeit_setup` on each of the instrument setups
           to construct a new pypeit file with the automatically generated frame
           types.
           
        #. Read both the by-hand and automatically generated frame types from
           these two pypeit files and check that they are identical. This check
           is *only* performed for the calibration frames, not any ``science``
           or ``standard`` frames.

Because this test is now included in the ``PypeIt``
:ref:`unit-tests`, this frame-typing check is performed by the
developers for every new version of the code.

