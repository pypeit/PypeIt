.. include:: ../include/links.rst

.. _nires_config_report:

Automated sorting of Keck/NIRES frames by instrument configuration
====================================================================

Version History
---------------

=========   ================   =========== ===========
*Version*   *Author*           *Date*      ``PypeIt``
=========   ================   =========== ===========
1.0         Debora Pelliccia   11 Nov 2022 1.11.1dev
1.1         Debora Pelliccia   24 Mar 2023 1.12.2dev
=========   ================   =========== ===========

----

Basics
------

To prepare for the data reduction, ``PypeIt`` automatically associates fits
files to specific :ref:`frame_types` (see :ref:`nires_frames_report`) and
collects groups of frames in unique instrument configurations. This is performed
by the :ref:`pypeit_setup` script, which sorts the frames and write a
:ref:`pypeit_file` for each unique configuration. See :ref:`setup_doc`.


NIRES configuration identification
------------------------------------

The instrument configurations are determined by :func:`pypeit.metadata.PypeItMetaData.unique_configurations`,
using a combination of header keys. However, since NIRES has a fixed configuration, the only header key used
to identify the desired configuration is ``INSTR``, which corresponds to the PypeIt ``fitstbl`` key ``dispname``.
``INSTR`` can be equal to ``spec`` when NIRES is used in spectroscopy mode, or ``imag`` if it is used in imaging
mode. Therefore, for our purpose of spectroscopic reduction, only one configuration is available.

NIRES calibration groups
--------------------------

PypeIt uses the concept of a "calibration group" to define a complete set of
calibration frames (e.g., arcs, flats) and the science frames to which these calibration
frames should be applied.


By default, :ref:`pypeit_setup` uses the setup identifier (e.g., A,B,C,D...) to assign frames to a single calibration
group. Since NIRES has only one configuration, i.e., ony one setup identifier, all the frames would have
the same PypeIt keyword ``calib``. However, since it is likely, during an observing night, to observe different
targets, PypeIt automatically assigns different ``calib`` values to the NIRES science frames of different targets.
This is done because PypeIt currently only uses OH sky lines in the science frames to perform the wavelength
calibration (instead of using the actual arc frames), therefore having different ``calib`` values for different
targets will allow each target to have its own wavelength calibration.
The standard frames, on the contrary, by default will not be used as arc/tilt frames since they are generally
taken with short exposures (i.e., the sky line may not be bright enough), therefore PypeIt will assign to them
the ``calib`` value of the closest science frames, i.e., they will share the wavelength calibration with the
science frame.
Finally, since usually only one set of flat observations are taken for the different targets, PypeIt
automatically sets ``calib = all`` for the flat frames, so that it can use them for the calibration
of all the different targets. See :ref:`calibration-groups`.
The user can always edit the :ref:`pypeit_file` to assign specific calibration frames to specific science
frames using the data in the ``calib`` column of the :ref:`data_block`.


.. _nires_combid_bkgid:

NIRES combination and background groups
-----------------------------------------

PypeIt is able to reduce data taken with a nodding pattern, by grouping the science frames into combination
and background groups. Science frames that should be combined together are assigned the same combination ID
(``comb_id``), while a background ID (``bkg_id``) identifies frames that are used as background images.
Frames with the same value of ``bkg_id`` will be combined together. The values of ``comb_id`` and ``bkg_id``
are provided in the :ref:`pypeit_file` as two columns in the :ref:`data_block`, so that users
can modify them according to their preferred reduction. See more detail in :ref:`a-b_differencing`.

For NIRES, ``PypeIt`` attempts to automatically assign ``comb_id`` and ``bkg_id`` to the science frames, by
using the information on the nodding pattern available in the files headers. Specifically, the keywords used are:

===============     ============
``fitstbl`` key     Header Key
===============     ============
``dithoff``         ``YOFFSET``
``dithpat``            no Key
``dithpos``         ``DPATIPOS``
===============     ============

The header key that provides information on the dither positions is ``DPATIPOS``, which defines those as positions
1, 2, 3, or 4 of a given dither pattern (``DPATNAME``). PypeIt, using the information from both
``DPATNAME`` and ``DPATIPOS``, expresses the dither positions (``dithpat``) as A, B, or C frames.
``dithoff``, ``dithpat``, and ``dithpos`` are provided in the :ref:`data_block`.

The dither patterns parsed by PypeIt are: "ABAB", "ABBA", "ABBAprime", "ABpat", "ABC" see examples below.
``comb_id`` and ``bkg_id`` will not be assigned if:

- ``dithoff`` is zero for every frames of a dither sequence;
- ``dithpat`` is NONE or MANUAL, or is none of the above patterns.

In these cases, the user should manually input the ``comb_id`` and ``bkg_id`` values.


If the observations were taken with a "ABpat" ``dithpat`` the :ref:`data_block` will look like:

.. code-block:: console

                filename |        frametype | ... | dithpat | dithpos | dithoff | frameno | calib | comb_id | bkg_id
    s181127_0076.fits.gz | arc,science,tilt | ... |   ABpat |       A |     2.5 |      76 |     1 |      57 |     58
    s181127_0077.fits.gz | arc,science,tilt | ... |   ABpat |       B |    -2.5 |      77 |     1 |      58 |     57


where the science frames have different ``comb_id`` (i.e., no frames will be combined), while the ``bkg_id``
for the A frame is equal to the ``comb_id`` of the B frame and vice versa.
This combination of ``comb_id`` and ``bkg_id`` will create two reduced frames:

.. code-block:: ini

    s181127_0076.fits.gz - s181127_0077.fits.gz (A-B)
    s181127_0077.fits.gz - s181127_0076.fits.gz (B-A)

If the observations were taken with an "ABAB", "ABBA", or "ABBAprime" ``dithpat``, the frames in the same dither
sequence will be combined. Here is an example for "ABBA":

.. code-block:: console

                filename |        frametype | ... | dithpat | dithpos | dithoff | frameno | calib | comb_id | bkg_id
    s181127_0020.fits.gz | arc,science,tilt | ... |    ABBA |       A |     2.5 |      20 |     1 |       5 |      6
    s181127_0021.fits.gz | arc,science,tilt | ... |    ABBA |       B |    -2.5 |      21 |     1 |       6 |      5
    s181127_0022.fits.gz | arc,science,tilt | ... |    ABBA |       B |    -2.5 |      22 |     1 |       6 |      5
    s181127_0023.fits.gz | arc,science,tilt | ... |    ABBA |       A |     2.5 |      23 |     1 |       5 |      6

This combination of ``comb_id`` and ``bkg_id`` will create two reduced frames:

.. code-block:: ini

    s181127_0020.fits.gz+s181127_0023.fits.gz - s181127_0021.fits.gz+s181127_0022.fits.gz (AA-BB)
    s181127_0021.fits.gz+s181127_0022.fits.gz - s181127_0020.fits.gz+s181127_0023.fits.gz (BB-AA)

Lastly, if observations were taken with an "ABC" ``dithpat``, where the A frame is taken at the center of the slit
(``dithoff = 0``) while the B and C frames are taken at a +/- offset, the B frames will be used as
background image for the frame taken at the center. Here is an example:

.. code-block:: console

                  filename |        frametype | ... | dithpat | dithpos | dithoff | frameno | calib | comb_id | bkg_id
    NR.20181126.38930.fits | arc,science,tilt | ... |     ABC |       A |     0.0 |      31 |     1 |       1 |      2
    NR.20181126.39604.fits | arc,science,tilt | ... |     ABC |       B |     5.0 |      32 |     1 |       2 |      3
    NR.20181126.40277.fits | arc,science,tilt | ... |     ABC |       C |    -5.0 |      33 |     1 |       3 |      2


This combination of ``comb_id`` and ``bkg_id`` will create three reduced frames:

.. code-block:: ini

    NR.20181126.38930.fits - NR.20181126.39604.fits (A-B)
    NR.20181126.39604.fits - NR.20181126.40277.fits (B-C)
    NR.20181126.40277.fits - NR.20181126.39604.fits (C-B)


Testing
-------

- Requirement PN-15 states: "As a user, I expect the pipeline to recognize dither positions from the header."

- Requirement PM-16 states: "As a user, I expect the pipeline to associate a pair of observations for sky
  subtraction and also allow for a manual selection and association. ABBA should associate A-B and B-A."


``PypeIt`` meets these requirements in the majority of use cases.

The test used to demonstrate that PN-14 is satisfied (:ref:`nires_frames_report`)
is also relevant here since it shows that ``PypeIt`` correctly
identifies NIRES data frame types and associates them with a single
configuration, all written to a single pypeit file.

To test that ``PypeIt`` can successfully assign ``comb_id`` and ``bkg_id``
to science frames following the information on the
dither pattern, we have added the ``test_setup_keck_nires_comb`` test to
``${PYPEIT_DEV}/unit_tests/test_setups.py``.

To run this test:

.. code-block:: bash

    cd ${PYPEIT_DEV}/unit_tests
    pytest test_setups.py::test_setup_keck_nires_comb -W ignore

The test requires that you have downloaded the ``PypeIt``
:ref:`dev-suite` and defined the ``PYPEIT_DEV`` environmental
variable that points to the relevant directory.

The algorithm for this test is run on three datasets,
'ABpat_wstandard', 'ABC_nostandard', 'ABBA_nostandard', and is as follows:

1. Collect the names of all files from the following directory:

.. code-block:: ini

    ${PYPEIT_DEV}/RAW_DATA/keck_nires/{dataset}

2. Use :class:`~pypeit.pypeitsetup.PypeItSetup` to automatically
   identify the configurations for these files.

3. Check that the code found one setup.

4. Read in a pre-generated .pypeit file with the correct calibration, combination, and background ids.

5. Check that the ``calib``, ``comb_id``, and ``bkg_id`` values for science frames in the automatically
identified files are the same as the ones in the pre-generated .pypeit file.


The dither sequences used here are: "ABBA", "ABC", "ABpat".
Because this test is now included in the ``PypeIt`` :ref:`unit-tests`, these configuration checks
are performed by the developers for every new version of the code.
