.. include:: ../include/links.rst

.. _lris_frames_report:

Automated typing of LRIS (BLUE and RED) frames
==============================================

Version History
---------------


=========   ================   =========== ===========
*Version*   *Author*           *Date*      ``PypeIt``
=========   ================   =========== ===========
1.0         Debora Pelliccia    6 Sep 2023 1.13.1.dev
=========   ================   =========== ===========

----

Basics
------

The general procedure used to assign frames a given type is described
here: :ref:`frame_types`.

LRIS frame typing
-----------------

The primary typing of LRIS frames is performed by
:func:`pypeit.spectrographs.keck_lris.KeckLRISSpectrograph.check_frame_type`.
This function checks the values of various header keywords against a
set of criteria used to classify the frame type. The same criteria are used for
``keck_lris_red``, ``keck_lris_red_orig``, ``keck_lris_red_mark4``, ``keck_lris_blue``,
and ``keck_lris_blue_orig``, unless otherwise noted.
The header cards required for the frame-typing and their associated keyword in the
:class:`~pypeit.metadata.PypeItMetaData` object are:

===============     ======================================================
``fitstbl`` key     Header Key
===============     ======================================================
``exptime``         ``ELAPTIME`` (``TELAPSE`` for ``keck_lris_red_mark4``)
``hatch``           ``TRAPDOOR``
``lampstat01``      See below
===============     ======================================================

``lampstat01`` is defined using a combination of header keywords, which include
``LAMPS``, ``MERCURY``, ``NEON``, ``ARGON``, ``CADMIUM``, ``ZINC``, ``HALOGEN``,
``KRYPTON``, ``XENON``, ``FEARGON``, ``DEUTERI``, ``FLAMP1``, ``FLAMP2``, ``FLIMAGIN``,
``FLSPECTR``. Since LRIS header keywords have changed over time, the exact combination
of keywords used to define ``lampstat01`` varies depending on the available header keywords.

The criteria used to select each frame type are as follows:

=============   ============   ======================================   ===========================================
Frame           ``hatch``      ``lampstat01``                           ``exptime``
=============   ============   ======================================   ===========================================
``science``     ``'open'``     ``'off'``                                ``>61s``
``standard``    ``'open'``     ``'off'``                                ``>1s`` & ``<61s``
``bias``        ``'closed'``   ``'off'``                                ``<0.001s``
``pixelflat``   ``'closed'``   ``'Halogen' or '2H'``                    ``<60s``(LRIS RED) or ``<300s`` (LRIS BLUE)
``pixelflat``   ``'open'``     ``'on'``                                 ``<60s``(LRIS RED) or ``<300s`` (LRIS BLUE)
``trace``       ``'closed'``   ``'Halogen'`` or ``'2H'``                ``<60s``(LRIS RED) or ``<300s`` (LRIS BLUE)
``trace``       ``'open'``     ``'on'``                                 ``<60s``(LRIS RED) or ``<300s`` (LRIS BLUE)
``illumflat``   ``'closed'``   ``'Halogen'`` or ``'2H'``                ``<60s``(LRIS RED) or ``<300s`` (LRIS BLUE)
``illumflat``   ``'open'``     ``'on'``                                 ``<60s``(LRIS RED) or ``<300s`` (LRIS BLUE)
``arc``         ``'closed'``   ``!= 'Halogen', '2H', 'on', 'off'``      Not used
``tilt``        ``'closed'``   ``!= 'Halogen', '2H', 'on', 'off'``      Not used
=============   ============   ======================================   ===========================================

Note that PypeIt employs commonly used value of ``exptime`` to distinguish frame type;
however, if needed, the user can specify a different value by
using the ``exprng`` parameter in the :ref:`pypeit_file`; see also :ref:`frame_types`.

The criteria used to select ``arc`` and ``tilt`` frames are identical; the same is true for
``pixelflat``, ``trace``, and ``illumflat`` frames. However, it's important to note that
PypeIt is able to correctly assign the ``pixelflat``, ``trace``, and ``illumflat`` types
to the internal and dome flat frames, but the twilight flats will generally have the
``science`` or ``standard`` type. Therefore, the user should manually change their frame type
in the :ref:`pypeit_file`.

Finally, note that a LRIS frame is never given a ``pinhole`` or ``dark`` type.


Testing
-------

Requirement PLL-16 states: "As a user, I expect the pipeline to automatically classify my data."

``PypeIt`` meets this requirement as demonstrated by the tests at
``${PYPEIT_DEV}/unit_tests/test_frametype.py``. There is one test
per spectrograph:

- ``test_lris_blue()``
- ``test_lris_blue_orig()``
- ``test_lris_red()``
- ``test_lris_red_orig()``
- ``test_lris_red_mark4()``

Here is an example of how to run the tests:

.. code-block:: bash

    cd ${PYPEIT_DEV}/unit_tests
    pytest test_frametype.py::test_lris_blue -W ignore

The tests requires that you have downloaded the PypeIt
:ref:`dev-suite` and defined the ``PYPEIT_DEV`` environmental
variable that points to the relevant directory. The algorithm for
all these tests is the same and is as follows:

    1. Find the directories in the :ref:`dev-suite` with Keck
       LRIS data (separately for ``keck_lris_blue``,
       ``keck_lris_blue_orig``, ``keck_lris_red``,
       ``keck_lris_red_orig``, ``keck_lris_red_mark4``).

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


