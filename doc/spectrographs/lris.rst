=========
Keck LRIS
=========


Overview
========

This file summarizes several instrument specific
settings that are related to the Keck/LRIS spectrograph.

Common Items
============

Flexure
+++++++

There is substantial flexure in the LRIS instrument and
the default settings attempts to characterize both the spectral
and spatial effects.

See the :doc:`../calibrations/flexure` notes if you wish
to turn either of these off.

.. _lrisb:

keck_lris_blue
==============

LRISb Default Settings
++++++++++++++++++++++

See :ref:`instr_par-keck_lris_blue` for
a listing of modifications to the default settings.  *You do not have to add these changes to
your PypeIt reduction file!*  This is just a listing of how the parameters used
for Keck/LRIS differ from the defaults listed in the preceding tables on
that page.

Taking Calibrations for LRISb
+++++++++++++++++++++++++++++

Arcs
----

We recommend that you turn on *all* of the standard
arc lamps,  including those slated for the red side.

These are::

    Ne,Ar,Cd,Kr,Xe,Zn,Hg

The archived solutions expect all of these lamps.

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

.. _400-3400-grism:


400/3400 grism
++++++++++++++

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

keck_lris_red
=============

Detectors
+++++++++

There have been 3 (or is it 4?!) generations of detectors
in the LRISr camera.  In PypeIt parlance, the original is named ``keck_lris_red_orig``,
the LBNL detectors (2kx4k) are ``keck_lris_red``, and the newest
Mark4 detector is ``keck_lris_red_mark4``.   

For the latter (Mark4), the wavelengths have been incorporated for the 
R400 grating only so far but the arxiv solutions from the LBNL detector
may work ok.  Check the outputs!

LRISr Default Settings
++++++++++++++++++++++

See :ref:`instr_par-keck_lris_red` for
a listing of modifications to the default settings.  *You do not have to add these changes to
your PypeIt reduction file!*  This is just a listing of how the parameters used
for Keck/LRIS differ from the defaults listed in the preceding tables on
that page.

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

Slit-masks
++++++++++

PypeIt can now incorporate slitmask information in the reduction
routine for LRIS similar to its DEIMOS capabilities.  That is, if the trace
calibrations files with mask information are fed to PypeIt, it is 
capable of using said information to determine object coordinates, 
identify targeted and serendipitous source and subsequently, collate by
ra/dec. 

Unfortunately, LRIS raw frames do not come ready with slitmask
data and thus this information needs to be inserted by the user before
processing with PypeIt if they are desirous of incorporating
the abovementioned features into their reduction. 
Here are the steps to do so:

#. Obtain the mask design files. The design files include:

    #. The AUTOSLIT-generated mask design files.

        #. One file with the ".file3" extension containing milling information.

        #. One file with the ".file1" extension containing the object catalog
           corresponding to the mask slits.

    #. The ASCII object list file fed as input to AUTOSLIT to generate the files
       above.

   ".file3" is mandatory while the other two files can be optionally excluded to
   debug `TILSOTUA <https://github.com/jsulli27/tilsotua>`__.
    
#. Process the design files with `TILSOTUA
   <https://github.com/jsulli27/tilsotua>`__ : The design files contain the
   milling blueprint (the `BluSlits` table).  When using the ".file3" design
   files, TILSOTUA creates FITS files based on the UCO/Lick template. The FITS
   mask design files have empty `DesiSlits`, `ObjectCat` and `SlitObjMap` binary
   tables. DEIMOS users may be familiar with these tables from their raw frames.
   TILSOTUA populates these tables using its ``xytowcs`` function (in
   ``LRIS_Mask_Coords_to_WCS.py``). One provides the code with two parameters:
   ``input_file_name`` is either the FITS or ".file3" mask design file (be sure
   the name includes the extension), and ``output_file_base`` is the prefix for
   the the four files that get created by the code. The calling sequence is:

   .. code-block:: python
    
        xytowcs(input_file_name="yourmaskname.file3",output_file_base="yourmaskname_output.fits")

#. The `ObjectCat` and `SlitObjMap` are only populated if ".file1" and the object list are provided.
    e.g.

    .. code-block:: python

        xytowcs(input_file_name="yourmaskname.file3",output_file_base="yourmaskname_output.fits",
                obj_file="yourtargets.obj", file1="yourmaskname.file1")

    It is assumed that the entries in `file1` and `obj_file` have unique `Name` values. i.e. Make
    sure you have a unique identifier for each object. Without this, it is not possible to correctly
    reconcile the two tables.
    
#. Append TILSOTUA's output to your raw trace files: Once the user is satisfied
   with the processed FITS file from TILSOTUA, append the binary tables to the
   trace FITS files. The user must first verify that TILSOTUA has indeed
   processed the files correctly. This implies:

    #. TILSOTUA has correctly identified the alignment stars (see the QA plot it generates).

    #. TILSOTUA has estimated the `TopDist` and `BotDist` in the `SlitObjMap` table correctly.

One may append the binary tables from the outputs as additional `HDUs` in the LRIS trace files. e.g.

    .. code-block:: python

        from astropy.io import fits

        tracehdus = fits.open("trace_rawframe.fits")
        autoslithdus = fits.open("yourmaskname_output.fits")

        for hdu in autoslithdus[1:]:
            tracehdus.append(hdu)
        tracehdus.writeto("trace_with_maskinfo.fits")
        
If processed correctly, PypeIt should now fully utilize its 
arsenal of slitmask processing tools to reduce and coadd spectra 
with the WCS information incorporated.

.. TODO: be specific about what you mean by "append the binary tables"

.. TODO: Does the above mean that LRIS should be included in lists of
   instruments that use mask design information.  Most relevant places claim we
   can only do this for DEIMOS and MOSFIRE.

