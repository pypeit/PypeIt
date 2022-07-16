.. include:: ../include/links.rst

.. _mosfire_frames_report:

Automated typing of Keck/MOSFIRE frames
=======================================

Version History
---------------

=========   ================   =========== ===========
*Version*   *Author*           *Date*      ``PypeIt``
=========   ================   =========== ===========
1.0         Debora Pelliccia   24 Oct 2021 1.6.1.dev
1.1         Debora Pelliccia   12 Jul 2021 1.9.2.dev
=========   ================   =========== ===========

----

Basics
------

The general procedure used to assign frames a given type is described
here: :ref:`frame_types`.

MOSFIRE frame typing
--------------------

The primary typing of MOSFIRE frames is performed by
:func:`pypeit.spectrographs.keck_mosfire.KeckMOSFIRESpectrograph.check_frame_type`.
This function checks the values of various header keywords against a
set of criteria used to classify the frame type. The header cards
required for the frame-typing and their associated keyword (when available)
in the :class:`~pypeit.metadata.PypeItMetaData` object are:

    ================   =================
       Header Key       ``fitstbl`` key
    ================   =================
      ``TRUITIME``        ``exptime``
      ``FLATSPEC``       ``lampstat01``
      ``PWSTATA7``          No key
      ``PWSTATA8``          No key
       ``OBJECT``         ``object``
       ``FILTER``         ``filter1``
         No key           ``idname``
      ``MASKNAME``        ``decker``
    ================   =================

``idname`` is defined using a combination of header keys, including ``PWSTATA7`` and ``PWSTATA8`` (see below).


The criteria used to select each frame type are as follows:

================   =====================   ============   ============   ============   ============   =============   ======================   ===================================
Frame              ``idname``              ``TRUITIME``   ``FLATSPEC``   ``PWSTATA7``   ``PWSTATA8``   ``FILTER``      ``OBJECT``               ``MASKNAME``
================   =====================   ============   ============   ============   ============   =============   ======================   ===================================
``science``        ``'object'``              ``>20s``        ``0``          ``0``          ``0``       ``!= 'Dark'``   ``not include 'Flat'``   Not used
``standard``       ``'object'``              ``<20s``        ``0``          ``0``          ``0``       ``!= 'Dark'``   ``not include 'Flat'``   Not used
``dark``           ``'dark'``                Not used        ``0``         Not used       Not used      ``'Dark'``     Not used                 Not used
``pixelflat``      ``'flatlamp'``            Not used        ``1``          ``0``          ``0``       ``!= 'Dark'``   Not used                 Not used
``trace``          ``'flatlamp'``            Not used        ``1``          ``0``          ``0``       ``!= 'Dark'``   Not used                 Not used
``illumflat``      ``'flatlamp'``            Not used        ``1``          ``0``          ``0``       ``!= 'Dark'``   Not used                 Not used
``lampoffflats``   ``'flatlampoff'``         Not used        ``0``          ``0``          ``0``       ``!= 'Dark'``   ``include 'Flat'``       Not used
``arc``            ``'arclamp'``             Not used        ``0``          ``1``          ``1``       ``!= 'Dark'``   Not used                 Not used
``arc``            ``'object'``              Not used        ``0``          ``0``          ``0``       ``!= 'Dark'``   ``not include 'Flat'``   ``not include 'long2pos_specphot'``
``tilt``           ``'arclamp'``             Not used        ``0``          ``1``          ``1``       ``!= 'Dark'``   Not used                 Not used
``tilt``           ``'object'``              Not used        ``0``          ``0``          ``0``       ``!= 'Dark'``   ``not include 'Flat'``   ``not include 'long2pos_specphot'``
================   =====================   ============   ============   ============   ============   =============   ======================   ===================================

Note that, by default, the exposure time (``TRUITIME``) is only used
to distinguish between ``science`` and ``standard`` frames; the criteria
for ``TRUITIME`` can be changed using the ``exprng``
parameter in the :ref:`pypeit_file`; see also :ref:`frame_types`.

Importantly, note that a MOSFIRE frame is never given a ``pinhole``
type. Also note that the criteria used to select ``arc`` and ``tilt``
frames are identical. The same is true for ``pixelflat``, ``trace``,
and ``illumflat`` frames.


Testing
-------
- PYPEIT-MOSFIRE
    Requirement PM-4 states: "As a user, I want the pipeline to
    automatically classify my calibrations."

``PypeIt`` meets this requirement as demonstrated by the test at
``${PYPEIT_DEV}/unit_tests/test_frametype.py``.  To run the test:

.. code-block:: bash

    cd ${PYPEIT_DEV}/unit_tests
    pytest test_frametype.py::test_mosfire -W ignore

The test requires that you have downloaded the ``PypeIt``
:ref:`dev-suite` and defined the ``PYPEIT_DEV`` environmental
variable that points to the relevant directory. The algorithm of the
test is as follows:

    1. Find all the directories in the :ref:`dev-suite` where a selected
       sample of Keck MOSFIRE data is.

    2. For each directory (i.e., instrument setup):

        a. Make sure there is a "by-hand" version of the pypeit file
           for this setup where a human (one of the pypeit
           developers) has ensured the frame types are correct.

        b. Effectively run :ref:`pypeit_setup` on each of the
           instrument setups to construct a new pypeit file with the
           automatically generated frame types.

        c. Read both the by-hand and automatically generated frame
           types from these two pypeit files and check that they are
           identical.

Because this test is now included in the ``PypeIt``
:ref:`unit-tests`, this frame-typing check is performed by the
developers for every new version of the code.
