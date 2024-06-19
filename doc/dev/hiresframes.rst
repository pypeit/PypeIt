.. include:: ../include/links.rst

.. _hires_frames:

Automated typing of HIRES frames
================================

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

The general procedure used to assign frames a given type is described
here: :ref:`frame_types`.

HIRES frame typing
-----------------

The primary typing of HIRES frames is performed by
:func:`pypeit.spectrographs.keck_hires.KECKHIRESSpectrograph.check_frame_type`.
This function checks the values of various header keywords against a
set of criteria used to classify the frame type.
The header cards required for the frame-typing and their associated keyword in the
:class:`~pypeit.metadata.PypeItMetaData` object are:

===============     ============
``fitstbl`` key     Header Key
===============     ============
``exptime``         ``ELAPTIME``
``hatch``           ``HATOPEN``
``lampstat01``      See below
No key              ``XCOVOPEN``
No key              ``AUTOSHUT``
===============     ============

``lampstat01`` is defined using a combination of header keywords, which include
``LAMPCAT1``, ``LAMPCAT2``, ``LAMPQTZ2``, ``LAMPNAME``. If ``LAMPCAT1 = True`` or
``LAMPCAT2 = True``, ``lampstat01`` will be equal to ``'ThAr1'`` or ``'ThAr2'``, respectively.
If ``LAMPQTZ2 = True`` or ``LAMPNAME = 'quartz1'``, ``lampstat01`` will be equal to ``'on'``.


The criteria used to select each frame type are as follows:

====================   ============   ============   ============   ======================================   ======================================================
Frame                  ``hatch``      ``AUTOSHUT``   ``XCOVOPEN``   ``lampstat01``                           ``exptime``
====================   ============   ============   ============   ======================================   ======================================================
``science``            ``True``       ``True``       ``True``       ``'off'``                                ``>601s``
``standard``           ``'open'``     ``True``       ``True``       ``'off'``                                ``>1s`` & ``<600s``
``bias``               ``False``      ``False``      ``True``       ``'off'``                                ``<0.001s``
``dark``               ``False``      ``True``       ``True``       ``'off'``                                Not used
``slitless_pixflat``   ``False``      ``True``       ``False``      ``'off'``                                ``<60s``
``pixelflat``          ``False``      ``True``       ``True``       ``'on'``                                 ``<60s``
``trace``              ``False``      ``True``       ``True``       ``'on'``                                 ``<60s``
``illumflat``          ``False``      ``True``       ``True``       ``'on'``                                 ``<60s``
``arc``                ``False``      ``True``       ``True``       ``'ThAr1'`` or ``'ThAr2'``               Not used
``tilt``               ``False``      ``True``       ``True``       ``'ThAr1'`` or ``'ThAr2'``               Not used
====================   ============   ============   ============   ======================================   ======================================================

Note that PypeIt employs commonly used value of ``exptime`` to distinguish frame type;
however, if needed, the user can specify a different value by
using the ``exprng`` parameter in the :ref:`pypeit_file`; see also :ref:`frame_types`.

The ``science`` and ``standard`` frames have identical selection criteria, except for the
``exptime`` value. In order to better distinguish between the two types, the ``RA`` and ``DEC`` header
keywords are also used to assign the ``standard`` type to frames with ``RA`` and ``DEC`` values that are
within 10 arcmin of one of the standard stars available in PypeIt (see :ref:`standards`).

The criteria used to select ``arc`` and ``tilt`` frames are identical; the same is true for
``pixelflat``, ``trace``, and ``illumflat`` frames. Note that if both ``pixelflat`` and
``slitless_pixflat`` frames are identified, the ``pixelflat`` assignment will be removed
so that the ``slitless_pixflat`` frames will be used for the flat fielding.

Finally, note that a HIRES frame is never given a ``pinhole`` type.


Testing
-------

To test that PypeIt can successfully identify HIRES framt types
among a set of files, we have added the
``test_hires()`` test to ``${PYPEIT_DEV}/unit_tests/test_frametype.py``.

Here is an example of how to run the test:

.. code-block:: bash

    cd ${PYPEIT_DEV}/unit_tests
    pytest test_frametype.py::test_hires -W ignore

The tests requires that you have downloaded the PypeIt
:ref:`dev-suite` and defined the ``PYPEIT_DEV`` environmental
variable that points to the relevant directory. The algorithm for
all these tests is the same and is as follows:

    1. Find the directories in the :ref:`dev-suite` with Keck
       HIRES data.

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


