.. include:: ../include/links.rst

.. _deimos_frames_report:

Automated typing of DEIMOS frames
=================================

Version History
---------------

=========   =============   =========== ===========
*Version*   *Author*        *Date*      ``PypeIt``
=========   =============   =========== ===========
1.0         Kyle Westfall   11 Sep 2020 1.1.1
1.1         Kyle Westfall   18 Sep 2020 1.1.2dev
1.2         Kyle Westfall   28 Sep 2020 1.1.2dev
1.3         Kyle Westfall   14 Oct 2020 1.2.1dev
1.4         Kyle Westfall   28 Jan 2021 1.3.1dev
=========   =============   =========== ===========

----

Basics
------

The general procedure used to assign frames a given type is described
here: :ref:`frame_types`.

DEIMOS frame typing
-------------------

The primary typing of DEIMOS frames is performed by
:func:`pypeit.spectrographs.keck_deimos.KeckDEIMOSSpectrograph.check_frame_type`.
This function checks the values of various header keywords against a
set of criteria used to classify the frame type. The headers cards
required for the frame-typing and their associated keyword in the
:class:`~pypeit.metadata.PypeItMetaData` object are:

===============     ============
``fitstbl`` key     Header Key
===============     ============
``mode``            ``MOSMODE``
``exptime``         ``ELAPTIME``
``lampstat01``      ``LAMPS``
``hatch``           ``HATCHPOS``
``idname``          ``OBSTYPE``
===============     ============

The criteria used to select each frame type are as follows:

=============   ==========================================  ===========   ============    ============
Frame           ``OBSTYPE``                                 ``LAMPS``     ``HATCHPOS``    ``MOSMODE``
=============   ==========================================  ===========   ============    ============
``science``     ``'Object'``                                ``'Off'``     ``'open'``      ``'Spectral'``
``bias``        ``'Bias'``                                  ``'Off'``     ``'closed'``    ``'Spectral'``
``dark``        ``'Dark'``                                  ``'Off'``     ``'closed'``    ``'Spectral'``
``pixelflat``   ``'IntFlat'``                               ``!='Off'``   ``'closed'``    ``'Spectral'``
``pixelflat``   ``'DmFlat'``, ``'SkyFlat'``                 ``!='Off'``   ``'open'``      ``'Spectral'``
``trace``       ``'IntFlat'``                               ``!='Off'``   ``'closed'``    ``'Spectral'``
``trace``       ``'DmFlat'``, ``'SkyFlat'``                 ``!='Off'``   ``'open'``      ``'Spectral'``
``illumflat``   ``'IntFlat'``                               ``!='Off'``   ``'closed'``    ``'Spectral'``
``illumflat``   ``'DmFlat'``, ``'SkyFlat'``                 ``!='Off'``   ``'open'``      ``'Spectral'``
``arc``         ``'Line'``                                  ``!='Off'``   ``'closed'``    ``'Spectral'``
``tilt``        ``'Line'``                                  ``!='Off'``   ``'closed'``    ``'Spectral'``
=============   ==========================================  ===========   ============    ============

Note that, by default, the exposure time (``ELAPTIME``) is not used to
distinguish frame types; however, this can be changed using the ``exprng``
parameter in the :ref:`pypeit_file`; see also :ref:`frame_types`.

Importantly, note that a DEIMOS frame is never given a ``pinhole``
type. Also note that the criteria used to select ``pixelflat``,
``trace``, and ``illumflat`` are identical; the same is true for
``arc`` and ``tilt`` frames. No distinction is currently made between
``IntFlat``, ``DmFlat``, or ``SkyFlat`` and how they're used in the
code.

Also note that ``PypeIt`` will ignore frames observed under
conditions that do not meet the restricted values set by
:func:`~pypeit.spectrographs.keck_deimos.KeckDEIMOSSpectrograph.valid_configuration_values`.
This currently requires all frames to have ``MOSMODE == 'Spectral'``
and ``AMPMODE == 'SINGLE:B'`` or ``AMPMODE == 'SINGLE:A'``. Frames that do not meet these criteria
will not be included in the automatically generated
:ref:`pypeit_file` created by :ref:`pypeit_setup`.

Testing
-------

Requirement PD-1 states: "As a user, I want the pipeline to
automatically classify my calibrations."

``PypeIt`` meets this requirement as demonstrated by the test at
``${PYPEIT_DEV}/unit_tests/test_frametype.py``.  To run the test:

.. code-block:: bash

    cd ${PYPEIT_DEV}/unit_tests
    pytest test_frametype.py::test_deimos -W ignore

The test requires that you have downloaded the ``PypeIt``
:ref:`dev-suite` and defined the ``PYPEIT_DEV`` environmental
variable that points to the relevant directory. The algorithm of the
test is as follows:

    1. Find all the directories in the :ref:`dev-suite` with Keck
       DEIMOS data.

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
