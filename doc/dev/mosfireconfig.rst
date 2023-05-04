.. include:: ../include/links.rst

.. _mosfire_config_report:

Automated sorting of Keck/MOSFIRE frames by instrument configuration
====================================================================

Version History
---------------

=========   ================   =========== ===========
*Version*   *Author*           *Date*      ``PypeIt``
=========   ================   =========== ===========
1.0         Debora Pelliccia   12 Jul 2021 1.9.2.dev
=========   ================   =========== ===========

----

Basics
------

To prepare for the data reduction, ``PypeIt`` automatically associates fits
files to specific :ref:`frame_types` (see :ref:`mosfire_frames_report`) and
collects groups of frames in unique instrument configurations (see below). This is performed
by the :ref:`pypeit_setup` script, which sorts the frames and write a
:ref:`pypeit_file` for each unique configuration. See :ref:`setup_doc`.


MOSFIRE configuration identification
------------------------------------

The MOSFIRE instrument configurations are determined by :func:`pypeit.metadata.PypeItMetaData.unique_configurations`,
which finds unique combinations of the following keywords:

====================  ===========
``fitstbl`` key       Header Key
====================  ===========
``dispname``          ``OBSMODE``
``decker_secondary``  No Key
``filter1``           ``FILTER``
``slitlength``        No Key
``slitwid``           No Key
====================  ===========

``decker_secondary``, ``slitlength``, and ``slitwid`` are defined as follow:

    - ``decker_secondary`` is by default equal to the header key ``MASKNAME``. However, when "LONGSLIT" is included
      in ``MASKNAME``, ``decker_secondary`` is equal to "LONGSLIT".
    - ``slitlength`` is the length of the slit, expressed in number of CSUs (Configurable Slit Units). This value is
      included in the header key ``MASKNAME`` for *LONGSLIT* mask. Therefore, ``slitlength`` is
      available only for *LONGSLIT* mask, and is ``None`` for *multi-object* and *long2pos* masks.
    - ``slitwid`` is the width of the slit, expressed in arcsec. This value is included in the header
      key ``MASKNAME`` for *LONGSLIT* mask. Therefore, ``slitwid`` is available only for
      *LONGSLIT* mask, and is ``None`` for *multi-object* and *long2pos* masks.


The unique configurations are determined by collating the relevant metadata from the headers
of all frames found by a run of :ref:`pypeit_setup`. Each unique configuration is given a capital letter
identifier (e.g., A,B,C,D...).

After that, :func:`pypeit.metadata.PypeItMetaData.set_configurations` associates each frame to the relevant
unique configuration ("setup"), by assigning a setup identifier (e.g., A,B,C,D...) to every frames for which
the values of the above keywords match the values of the specific unique configuration.

This process is slightly modified for *LONGSLIT* and *long2pos* masks. It is common observation practise with
MOSFIRE to use *LONGSLIT* masks that are 46 CSUs long for the calibration of science frames obtained with
shorter *LONGSLIT* masks of the same width.
For this reason, when *LONGSLIT* masks are used, we don't require that the value of ``slitlength`` for the
calibration frames (arcs, flats) matches the value of the specific unique configuration. This allows to
associate the same calibration frames (only the ones with *LONGSLIT* masks that are 46 CSUs long) to different
unique configurations, i.e., where the science/standard frames are obtained with shorter *LONGSLIT*.

Similarly, when *long2pos* masks are reduced, calibration frames taken with masks that have ``MASKNAME``
equal to "long2pos" are used to calibrated science/standard frames that have ``MASKNAME`` equal to "long2pos_specphot".
In this case, the default behaviour would identify the calibration frames and the science/standard frames
as part of different unique configurations. Instead, when *long2pos* masks are used, we modify in place the requirement
for the calibration frames. Specifically, only when we are assigning the setup identifier to the calibration frames,
we temporarily change the required value of ``decker_secondary`` in the specific unique configuration from
"long2pos_specphot" to "long2pos". This allows to identify "long2pos" calibration frames and "long2pos_specphot"
science/standard frames as part of the same unique configuration.


MOSFIRE calibration groups
--------------------------

``PypeIt`` uses the concept of a "calibration group" to define a complete set of
calibration frames (e.g., arcs, flats) and the science frames to which these calibration
frames should be applied.

By default, :ref:`pypeit_setup` uses the setup identifier (e.g., A,B,C,D...) to assign frames to a single calibration
group. Frames that are in the same calibration group will have the same ``PypeIt`` keyword ``calib``.
No automated procedure exists to do anything except this. However, the user can edit the :ref:`pypeit_file` to,
within a given configuration, assign specific calibration frames to specific science frames using the data in
the ``calib`` column of the :ref:`data_block`.

.. _mosfire_combid_bkgid:

MOSFIRE combination and background groups
-----------------------------------------

``PypeIt`` is able to reduce data taken with a nodding pattern, by grouping the science frames into combination
and background groups. Science frames that should be combined together are assigned the same combination ID
(``comb_id``), while a background ID (``bkg_id``) identifies frames that are used as background images.
Frames with the same value of ``bkg_id`` will be combined together. The values of ``comb_id`` and ``bkg_id``
are provided in the :ref:`pypeit_file` as two columns in the :ref:`data_block`, so that users
can modify them according to their preferred reduction. See more detail in :ref:`a-b_differencing`.

For MOSFIRE, ``PypeIt`` attempts to automatically assign ``comb_id`` and ``bkg_id`` to the science frames, by
using the information on the nodding pattern available in the files headers. Specifically, the keywords used are:

===============     ============
``fitstbl`` key     Header Key
===============     ============
``dithoff``         ``YOFFSET``
``dithpat``         ``PATTERN``
``dithpos``         ``FRAMEID``
===============     ============

which are also provided in the :ref:`data_block`.

If the observations were taken with a "Slit Nod"/"Mask Nod" ``dithpat`` or using the *long2pos* slitmask,
the :ref:`data_block` will look like::

                  filename | frametype | ... |  dithpat | dithpos | dithoff | frameno | calib | comb_id | bkg_id
    MF.20141126.17372.fits |   science | ... | Mask Nod |       A |    10.0 |     182 |     4 |      21 |     22
    MF.20141126.17526.fits |   science | ... | Mask Nod |       B |   -10.0 |     183 |     4 |      22 |     21
    MF.20141126.17686.fits |   science | ... | Mask Nod |       A |    10.0 |     184 |     4 |      23 |     24
    MF.20141126.17842.fits |   science | ... | Mask Nod |       B |   -10.0 |     185 |     4 |      24 |     23


where all the science frames have different ``comb_id`` (i.e., no frames will be combined), while the ``bkg_id``
for the frame at the "A" ``dithpos`` is equal to the ``comb_id`` of the frame at the "B" ``dithpos`` and vice versa.
This combination of ``comb_id`` and ``bkg_id`` will create four reduced frames::

    MF.20141126.17372.fits - MF.20141126.17526.fits (A-B)
    MF.20141126.17526.fits - MF.20141126.17372.fits (B-A)
    MF.20141126.17686.fits - MF.20141126.17842.fits (A-B)
    MF.20141126.17842.fits - MF.20141126.17686.fits (B-A)

If the observations were taken with an "ABAB" or "ABBA" ``dithpat``, the frames in the same dither
sequence will be combined. Here is an example for "ABBA"::

                  filename | frametype | ... |  dithpat | dithpos | dithoff | frameno | calib | comb_id | bkg_id
    MF.20130903.35554.fits |  standard | ... |     ABBA |       A |     6.0 |     276 |     0 |       1 |      2
    MF.20130903.35576.fits |  standard | ... |     ABBA |       B |    -6.0 |     277 |     0 |       2 |      1
    MF.20130903.35593.fits |  standard | ... |     ABBA |       B |    -6.0 |     278 |     0 |       2 |      1
    MF.20130903.35620.fits |  standard | ... |     ABBA |       A |     6.0 |     279 |     0 |       1 |      2
    MF.20130903.35679.fits |  standard | ... |     ABBA |       A |     6.0 |     280 |     0 |       5 |      6
    MF.20130903.35697.fits |  standard | ... |     ABBA |       B |    -6.0 |     281 |     0 |       6 |      5
    MF.20130903.35710.fits |  standard | ... |     ABBA |       B |    -6.0 |     282 |     0 |       6 |      5
    MF.20130903.35726.fits |  standard | ... |     ABBA |       A |     6.0 |     283 |     0 |       5 |      6

This combination of ``comb_id`` and ``bkg_id`` will create four reduced frames::

    MF.20130903.35554.fits+MF.20130903.35620.fits - MF.20130903.35576.fits+MF.20130903.35593.fits (AA-BB)
    MF.20130903.35576.fits+MF.20130903.35593.fits - MF.20130903.35554.fits+MF.20130903.35620.fits (BB-AA)
    MF.20130903.35679.fits+MF.20130903.35726.fits - MF.20130903.35697.fits+MF.20130903.35710.fits (AA-BB)
    MF.20130903.35697.fits+MF.20130903.35710.fits - MF.20130903.35679.fits+MF.20130903.35726.fits (BB-AA)

Lastly, if observations were taken with the *long2pos_specphot* slitmask, where two frames are taken at the
"A" and "B" ``dithpos`` and one frame is taken at the center (``dithoff = 0``) of a wider slit (see
https://www2.keck.hawaii.edu/inst/mosfire/long2pos.html ), one of the two "A" and "B" frames are used as
background image for the frame taken at the center of the wider slit. Here is an example::

                  filename | frametype | ... |  dithpat | dithpos | dithoff | frameno | calib | comb_id | bkg_id
    MF.20181217.55882.fits |  standard | ... | long2pos |       A |     0.0 |     286 |     7 |      27 |     28
    MF.20181217.55901.fits |  standard | ... | long2pos |       B |    -7.0 |     287 |     7 |      28 |     29
    MF.20181217.55920.fits |  standard | ... | long2pos |       A |     7.0 |     288 |     7 |      29 |     28

This combination of ``comb_id`` and ``bkg_id`` will create three reduced frames::

    MF.20181217.55882.fits - MF.20181217.55901.fits (center pos -B)
    MF.20181217.55901.fits - MF.20181217.55920.fits (A-B)
    MF.20181217.55920.fits - MF.20181217.55901.fits (B-A)




Testing
-------

- Requirement PM-5 states: "As a user, I want the pipeline to automatically correctly associate
  calibrations with science data."

- Requirement PM-34 states: "As a user, I expect the calibration association system to use
  calibrations from full long slits (46 bars long) for short “long slits” of the same width."

- Requirement PM-7 states: "Use longslit calibrations for longslit calibrations even if the
  length of the slit is different."

- Requirement PM-8 states: "Use narrow slit long2pos calibration for long2pos reduction."

- Requirement PM-13 states: "As a user, I expect my dithering pattern to automatically
  turn into the proper reduction sequence. If I used ABBA, I expect this to be taken into account"

``PypeIt`` meets these requirements in the majority of use cases.

The test used to demonstrate that PM-4 is satisfied (:ref:`mosfire_frames_report`)
is also relevant here since it shows that ``PypeIt`` correctly
identifies MOSFIRE data frame types and associates them with a single
configuration, all written to a single pypeit file.

To test that ``PypeIt`` can successfully identify multiple configurations among a set of files,
and can assign ``comb_id`` and ``bkg_id`` to science frames following the information on the
dither pattern, we have added the ``test_setup_keck_mosfire_multiconfig`` test to
``${PYPEIT_DEV}/unit_tests/test_setups.py``.

To run this test:

.. code-block:: bash

    cd ${PYPEIT_DEV}/unit_tests
    pytest test_setups.py::test_setup_keck_mosfire_multiconfig -W ignore

The test requires that you have downloaded the ``PypeIt``
:ref:`dev-suite` and defined the ``PYPEIT_DEV`` environmental
variable that points to the relevant directory.

The algorithm for this test is as follows:

    1. Collect the names of all files in the following two directories::

        ${PYPEIT_DEV}/RAW_DATA/keck_mosfire/K_long
        ${PYPEIT_DEV}/RAW_DATA/keck_mosfire/long2pos1_H

    2. Use :class:`~pypeit.pypeitsetup.PypeItSetup` to automatically
       identify the configurations for these files.

    3. Check that the code found four configurations and wrote the
       pypeit files for each.

    4. For each configuration:

        a. Read the pypeit file

        b. Check that the name for the setup is correct ('A', 'B', 'C', or 'D')

        c. Check that the calibration group is the same for all frames ('0', '1', '2', '3')

        d. Check that ``comb_id`` and ``bkg_id`` for the science frames are what expected. The
           dither sequences used here are: "ABBA", "long2pos", "ABA'B'", "Mask Nod".


Because this test is now included in the ``PypeIt`` :ref:`unit-tests`, these configuration checks
are performed by the developers for every new version of the code.
