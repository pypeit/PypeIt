.. include:: ../include/links.rst

.. _lris:

=========
Keck LRIS
=========


Overview
========

This file summarizes several instrument specific
settings that are related to the Keck/LRIS (RED and BLUE) spectrograph.

Common Items to LRISb and LRISr
===============================

.. _lris_FITS_format:

FITS format
+++++++++++

Before May 2009, both LRISb and LRISr observations were stored in standard simple
`FITS <https://fits.gsfc.nasa.gov/fits_home.html>`__ format consisting of a
single primary HDU without any extensions. Subsequently, the data were stored in
multi-extension FITS files, including four extensions, one for each amplifier.
To handle both formats, PypeIt defined 2 sets of spectrographs for each LRIS
camera, one for the pre-May 2009 data (``keck_lris_blue_orig`` and ``keck_lris_red_orig``)
and one for the post-May 2009 data (``keck_lris_blue`` and ``keck_lris_red``).
Note that this change in FITS format coincided with the installation of the LBNL
detectors (2kx4k) in LRISr (see :ref:`keck-lris-red`).

.. _lris_slitmask:

Slit-masks
++++++++++

PypeIt can now incorporate slitmask information in the reduction
routine for ``keck_lris_red``, ``keck_lris_red_mark4``, and ``keck_lris_blue``
similar to its DEIMOS capabilities (see :ref:`deimos-mask-matching`).
That is, if the trace calibration files with mask information are used, PypeIt is
capable of using said information to match the traced slit edges to those predicted
by the slitmask design, assign to each extracted spectrum the corresponding RA, Dec and
object name information, force the extraction of undetected objects at the location
expected from the slitmask design, identify serendipitous sources and, subsequently,
collate by RA/Dec. See `Additional Reading`_ for more information.

Unfortunately, LRIS raw frames do not come ready with slitmask
data and thus this information needs to be inserted by the user before
processing with PypeIt if they are desirous of incorporating
the above mentioned features into their reduction. To do so, the user
must first obtain the mask design files and then process them with
the software `TILSOTUA <https://github.com/jsulli27/tilsotua>`__.
Here are the steps to follow:

#. Obtain the mask design files, which include:

    #. Output files produced by `autoslit <https://www2.keck.hawaii.edu/inst/lris/autoslit_WMKO.html>`__.

        - One file with the ``".file3"`` extension containing milling information.

        - One file with the ``".file1"`` extension containing the object catalog matched to the slitmask slits.

    #. The ASCII object list file fed as input to
       `autoslit <https://www2.keck.hawaii.edu/inst/lris/autoslit_WMKO.html>`__ to generate the files above.

   .. note::

         ``".file3"`` is mandatory while the other two files can be optionally excluded to
         debug `TILSOTUA <https://github.com/jsulli27/tilsotua>`__.

#. Process the design files with `TILSOTUA <https://github.com/jsulli27/tilsotua>`__ :
   The design files contain the milling blueprint (the ``BluSlits`` table).
   When using the ``".file3"`` design file, TILSOTUA creates a FITS file with the ``BluSlits`` table
   following the UCO/Lick template. If the ``".file1"`` file and the object list are provided,
   the FITS mask design file will also includes the ``DesiSlits``, ``ObjectCat`` and ``SlitObjMap``
   binary tables, otherwise they will be empty. These tables include information on the
   slitmask design, the object catalog and the mapping between the two, similar to the
   binary tables in DEIMOS raw frames. TILSOTUA populates these tables using its ``xytowcs``
   function (in `LRIS_Mask_Coords_to_WCS.py
   <https://github.com/jsulli27/tilsotua/blob/master/tilsotua/LRIS_Mask_Coords_to_WCS.py>`__).
   This function can be run by providing two parameters:

        - the ``data_input_name``, which is either the FITS or ``".file3"`` mask
          design file (be sure the name includes the extension);

        - the ``output_file``, which is the name of the output file that
          TILSOTUA will generate. Do not include any extension such as ``.fits``.

   If only the ``".file3"`` file is provided, the calling sequence is:

   .. code-block:: python

        from tilsotua import xytowcs

        xytowcs(data_input_name="yourmaskname.file3",output_file="yourmaskname_output")

   Although the other parameters are optional for `xytowcs` as a standalone code, users interested in applying the slitmask information to their PypeIt reduction **must provide the `obj_file` and `file1` files to ensure that object names are assigned to the extracted spectra**.

   If the ``".file1"`` file and the object list are provided, the calling sequence is:

    .. code-block:: python

        from tilsotua import xytowcs

        xytowcs(data_input_name="yourmaskname.file3",output_file="yourmaskname_output",
                obj_file="yourtargets.obj", file1="yourmaskname.file1")

   It is assumed that the entries in ``file1`` and ``obj_file`` have unique ``Name`` values, i.e., make
   sure you have a unique identifier for each object. Without this, it is not possible to correctly
   reconcile the two tables.

#. Add the TILSOTUA-generated slitmask design information to your raw trace FITS files:
   The user must first verify that TILSOTUA has indeed processed the files correctly. This implies:

   - TILSOTUA has correctly identified the alignment stars (see the QA plot it generates).

   - TILSOTUA has estimated the values of the ``TopDist`` and ``BotDist`` columns in the ``SlitObjMap`` table correctly.

   Once satisfied with the processed FITS file from TILSOTUA, the user can append the binary tables
   populated by TILSOTUA to the LRIS trace FITS files as additional HDUs, e.g.:

   .. code-block:: python

        from astropy.io import fits

        tracehdus = fits.open("trace_rawframe.fits")
        autoslithdus = fits.open("yourmaskname_output.fits")

        for hdu in autoslithdus[1:]:
            tracehdus.append(hdu)
        tracehdus.writeto("trace_with_maskinfo.fits")

If processed correctly, PypeIt should now fully utilize its arsenal of slitmask processing tools
to reduce and coadd spectra with the WCS information incorporated.

Flexure
+++++++

There is substantial flexure in the LRIS instrument and
the default settings attempts to characterize both the spectral
and spatial effects.

See the :doc:`../calibrations/flexure` notes if you wish
to turn either of these off.

.. _lrisb:

LRIS BLUE
=========

This section provides information on both ``keck_lris_blue_orig`` and ``keck_lris_blue``.
When not specified, the information applies to both.

Default Settings
++++++++++++++++

See :ref:`instr_par-keck_lris_blue` and :ref:`instr_par-keck_lris_blue_orig` for
a listing of modifications to the default settings.
*You do not have to add these changes to your PypeIt reduction file!*  This is just a listing of
how the parameters used for LRISb differ from the defaults listed in the preceding tables on that page.
Moreover, additional modifications may have been made for specific setups, e.g, for different grisms,
or different slitmasks, etc. You can see a list of all the used parameters in the ``keck_lris_blue_XXX.par``
file generated by PypeIt at the beginning of the reduction.

Calibrations
++++++++++++

Wavelength calibration
^^^^^^^^^^^^^^^^^^^^^^

Arcs
----

We recommend that you turn on *most* of the standard
arc lamps, including those slated for the red side.

These are::

    Ne,Ar,Cd,Kr,Xe,Zn,Hg

The archived solutions expect most of these lamps. FeAr lamp can also be used,
since a FeAr line list is available for reduction process, although it appears
to be less useful to obtain good wavelength solutions.

Wavelength Solution
-------------------

As default, the wavelength calibration is performed using the :ref:`wvcalib-fulltemplate`
algorithm. The templates are created in the same way as done for Keck/DEIMOS
(see :ref:`deimos_wavecalib`) and kept in the ``data/arc_lines/reid_arxiv`` directory.
There are four templates, one per each LRISb grism:

===========  =======================================================
   GRISM                         template
===========  =======================================================
 300/5000     keck_lris_blue_B300_5000_d680_ArCdHgKrNeXeZnFeAr.fits
 400/3400       keck_lris_blue_B400_3400_d560_ArCdHgNeZnFeAr.fits
 600/4000       keck_lris_blue_B600_4000_d560_ArCdHgKrNeXeZn.fits
 1200/3400      keck_lris_blue_B1200_3400_d560_ArCdHgNeZn.fits
===========  =======================================================

PypeIt will automatically choose the right template according to the specific dataset.
These templates work for both ``keck_lris_blue_orig`` and ``keck_lris_blue``.

Flat Fielding
^^^^^^^^^^^^^

Pixel Flat
----------

It is recommend to correct for pixel-to-pixel variations using a slitless
flat.  If you did not take such calibration frames or cannot process them,
you may wish to use an archival.
`This link <https://drive.google.com/drive/folders/1YmDgCgXrsRbkuH_Pc_MLShWVdSrMkoFP?usp=sharing>`__
has the existing ones staged by the PypeIt team.

And then set the following in your :ref:`pypeit_file`:

.. code-block:: ini

    [calibrations]
        [[flatfield]]
            frame = path_to_the_file/PYPEIT_LRISb_pixflat_B600_2x2_17sep2009.fits.gz

.. warning::

    Internal flats may be too bright and need to be tested.

Trace Flat
----------

We strongly recommend on-the-sky trace flats through full instrument
setup.  Aim for 1000 counts per pixel above bias.
These are best achieved by taking twilight flats within 15 minutes
of sunset/sunrise.

.. warning::

    Internal/dome flats are likely to be too faint in the very blue.

.. TODO: Can somebody comment on this. This warning says that the internal flats
   are too faint in the very blue for tracing, while the warning above says that they may be
   too bright for pixelflat.  Is it true?

.. _400-3400-grism:

400/3400 grism
**************

If you are using this grism, you are likely aware there are
strong ghosts.  We have found these complicate edge tracing
for dome flats (sky flats appear ok).  Therefore, you may
need to increase the ``edge_thresh`` parameter to 
40 for successful performance, i.e.:

.. code-block:: ini

    [calibrations]
        [[slitedges]]
            edge_thresh = 40

.. _keck-lris-red:

LRIS RED
========

This section provides information on ``keck_lris_red_orig``, ``keck_lris_red``,
and ``keck_lris_red_mark4``. When not specified, the information applies to all.

.. _lrisr-detectors:

Detectors
+++++++++

There have been 3 (or is it 4?!) generations of detectors in the LRISr camera.
The first detector change (from Tektronik 2kx2k to LBNL 2kx4k) happened in May 2009,
concurrently with the change in the FITS file format (see :ref:`lris_FITS_format`).
To reduce the data taken before May 2009, the user should use the
``keck_lris_red_orig`` spectrograph. For data taken with the LBNL detectors,
after May 2009, the user should use the ``keck_lris_red`` spectrograph.
The newest detector was installed around May 2022 and is referred to as
the Mark4 detector. To reduce the data taken with this detector,
the user should use the ``keck_lris_red_mark4`` spectrograph.


Default Settings
++++++++++++++++

See :ref:`instr_par-keck_lris_red`, :ref:`instr_par-keck_lris_red_orig`,
and :ref:`instr_par-keck_lris_red_mark4`  for a listing of modifications to the default settings.
*You do not have to add these changes to your PypeIt reduction file!* This is just a listing of how
the parameters used for LRISr differ from the defaults listed in the preceding tables on that page.
Moreover, additional modifications may have been made for specific setups, e.g, for different gratings,
or different slitmasks, etc. You can see a list of all the used parameters in the ``keck_lris_red_XXX.par``
file generated by PypeIt at the beginning of the reduction.

Calibrations
++++++++++++

Wavelength calibration
^^^^^^^^^^^^^^^^^^^^^^

Arcs
----

We recommend that you turn on *most* of the standard
arc lamps::

    Ne,Ar,Cd,Kr,Xe,Zn,Hg

The archived solutions expect most (or all) of these lamps.

.. _lrisr_wavesol:

Wavelength Solution
-------------------
As default, the wavelength calibration is performed using the :ref:`wvcalib-fulltemplate`
algorithm. The templates are created in the same way as done for Keck/DEIMOS
(see :ref:`deimos_wavecalib`) and kept in the ``data/arc_lines/reid_arxiv`` directory.
When the first detector change happened (from Tektronik 2kx2k to LBNL 2kx4k) in May 2009
(see :ref:`lrisr-detectors`), the pixel size changed from 24 to 15 microns. Because of this,
different templates are used for ``keck_lris_red_orig`` and ``keck_lris_red``. Luckily,
the pixel size did not change when the Mark4 detector was installed, so the same templates
are used for ``keck_lris_red`` and ``keck_lris_red_mark4``.
These are the templates available, one per each LRISr grating:

============  ================================================  =============================================================
  GRATING            template for ``keck_lris_red_orig``          template for ``keck_lris_red`` and ``keck_lris_red_mark4``
============  ================================================  =============================================================
  150/7500        keck_lris_red_orig_R150_7500_ArHgNe.fits                keck_lris_red_R150_7500_ArCdHgNeZn.fits
  300/5000      keck_lris_red_orig_R300_5000_ArCdHgNeZn.fits            keck_lris_red_R300_5000_ArCdHgKrNeXeZn.fits
  400/8500      keck_lris_red_orig_R400_8500_ArCdHgNeZn.fits            keck_lris_red_R400_8500_ArCdHgKrNeXeZn.fits
  600/5000      keck_lris_red_orig_R600_5000_ArCdHgNeZn.fits            keck_lris_red_R600_5000_ArCdHgKrNeXeZn.fits
  600/7500      keck_lris_red_orig_R600_7500_ArCdHgNeZn.fits            keck_lris_red_R600_7500_ArCdHgKrNeXeZn.fits
  600/10000     keck_lris_red_orig_R600_10000_ArCdHgNeZn.fits           keck_lris_red_R600_10000_ArCdHgKrNeXeZn.fits
  831/8200      keck_lris_red_orig_R831_8200_ArCdHgNeZn.fits            keck_lris_red_R831_8200_ArCdHgKrNeXeZn.fits
  900/5500      keck_lris_red_orig_R900_5500_ArCdHgNeZn.fits              keck_lris_red_R900_5500_ArCdHgNeZn.fits
 1200/7500      keck_lris_red_orig_R1200_7500_ArCdHgNeZn.fits           keck_lris_red_R1200_7500_ArCdHgKrNeXeZn.fits
 1200/9000                                                                     keck_lris_red_R1200_9000.fits
============  ================================================  =============================================================

PypeIt will automatically choose the right template according to the specific dataset.
Note that the 1200/9000 grating was first released in 2013, so it was not available for
``keck_lris_red_orig``.

Known issues
============

LRISb Slit Edges
++++++++++++++++

When observing in long-slit mode, PypeIt might set the slit incorrectly
for detector 2.  This may occur if the counts from the flat field
are too low (e.g., using internal flats rather than twilight
flats with a higher signal in the blue).
Therefore, if you use internal flats, be careful to inspect the
slits defined by PypeIt as described in :doc:`../calibrations/edges`.

If the defined slit(s) does not cover the portion of
the illuminated detector where your source falls, you
can manually define the slit position as described
in :ref:`slit-tracing-missing-slit`.

Here is an example for the PypeIt file:

.. code-block:: ini

    [calibrations]
        [[slitedges]]
            add_slits = 2:788:10:650
            sync_predict = nearest

This will force a slit onto the detector for reduction.

Multi-slit
++++++++++

The code may identify a 'ghost' slit in empty detector real
estate if your mask does not fill most of the field.  Be prepared
to ignore it.


Additional Reading
==================

Here are additional docs related to Keck/LRIS.  Note many of them are related
to the development of PypeIt for use with LRIS data:

.. TODO: Generally useful information in these dev docs should be moved into
.. user-level doc pages, even if that means repeating information.

.. toctree::
   :maxdepth: 1

   ../dev/lrisframes
   ../dev/lrisconfig
   ../dev/slitmask_ids
   ../dev/radec_object
   ../dev/add_missing_obj
