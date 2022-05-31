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
the default settings attemps to characterize both the spectral
and spatial effects.

See the :doc:`flexure` notes if you wish
to turn either of these off.

.. _LRISb:

keck_lris_blue
==============

LRISb Default Settings
++++++++++++++++++++++

See :ref:`pypeit_par:Instrument-Specific Default Configuration` for
a listing of modifications to the default settings.

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
`This link <https://drive.google.com/drive/folders/1YmDgCgXrsRbkuH_Pc_MLShWVdSrMkoFP?usp=sharing>`_
has the existing ones staged by the PypeIt team.

And then set the following in your :doc:`pypeit_file`::

    [calibrations]
      [[flatfield]]
           frame = path_to_the_file/PYPEIT_LRISb_pixflat_B600_2x2_17sep2009.fits.gz

WARNING: Internal flats may be too bright and need to be tested.

Trace Flat
----------

We strongly recommend on-the-sky trace flats through full instrument
setup.  Aim for 1000 counts per pixel above bias.
These are best achieved by taking twilight flats within 15 minutes
of sunset/sunrise.

WARNING: Internal/dome flats are likely to be too faint in the
very blue.

.. _400-3400-grism:


400/3400 grism
++++++++++++++

If you are using this grism, you are likely aware there are
strong ghosts.  We have found these complicate edge tracing
for dome flats (sky flats appear ok).  Therefore, you may
need to increase the `edge_thresh` parameter to 
40 for successful performance, i.e.::

    [calibrations]
        [[slitedges]]
            edge_thresh = 40

.. _keck-lris-red:

keck_lris_red
=============

Detectors
+++++++++

There have been 3 (or is it 4?!) generations of detectors
in the LRISr camera.  The original is named `keck_lris_red_orig`,
the LBNL detectors (2kx4k) are `keck_lris_red` and the newest
Mark4 detector is `keck_lris_red_mark4`.   

For the latter (Mark4), the wavelengths have been incorporated for the 
R400 grating only so far but the arxiv solutions from the LBNL detector
may work ok.  Check the outputs!

LRISr Default Settings
++++++++++++++++++++++

Here are the deviations from the default settings
for LRISr::

    settings trace slits sigdetect 50.0   # Good for relatively bright dome flats
    settings trace slits pca params [3,2,1,0]

Known issues
============

LRISb Slit Edges
++++++++++++++++

When observing in long-slit mode, PypeIt might set the slit incorrectly
for detector 2.  This may occur if the counts from the flat field
are too low (e.g., using internal flats rather than twilight
flats with a higher signal in the blue).
Therefore, if you use internal flats, be careful to inspect the
slits defined by PypeIt as described in :doc:`master_edges`.

If the defined slit(s) does not cover the portion of
the illuminated detector where your source falls, you
can manually define the slit position as described
in :ref:`slit_tracing:Missing A Slit`.


Here is an example for the PypeIt file::

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
routine for LRIS similar to its DEIMOS capabilities. i.e. If the trace
calibrations files with mask information are fed to PypeIt, it is 
capable of using said information to determine object coordinates, 
identify targeted and serendipitous source and subsequently, collate by
ra/dec. 

Unfortunately, LRIS raw frames do not come ready with slitmask
data and thus this information needs to be inserted by the user before
processing with PypeIt if they are desirous of incorporating
the abovementioned features into their reduction. 
Here are the steps to do so:

#. Obtain the mask design files. The design files can be in one of two forms:

    #. The AUTOSLIT-generated mask design files (those with the ".file3" extension). 
    #. FITS files of the mask designs from UCO/Lick: As of 2022 Jan 27th, when the AUTOSLIT mask design files (ascii files that end with ".file3" by default) are fed to the slitmask database, a FITS file is generated with the milling blueprint. The FITS files have HDUs akin to DEIMOS raw files (sans the raw image of course). Please contact Steve Allen of UCO/Lick (UCSC) to procure these files.
    
#. Process the design files with `TILSOTUA <https://github.com/jsulli27/tilsotua>`_ : The design files contain the milling blueprint (the `BluSlits` table).  When using the ".file3" design files, TILSOTUA creates FITS files based on the UCO/Lick template. The FITS mask design files have empty `DesiSlits`, `ObjectCat` and `SlitObjMap` binary tables. DEIMOS users may be familiar with these tables from their raw frames. TILSOTUA populates these tables using its `xytowcs` function (in `LRIS_Mask_Coords_to_WCS.py`). One provides the code with two parameters: `input_file_name` is either the FITS or ".file3" mask design file (be sure the name includes the extension), and `output_file_base` is the prefix for the the four files that get created by the code. The calling sequence is:: 
    
    xytowcs(input_file_name,output_file_base)

#. Append TILSOTUA's output to your raw trace files: Once the user is satisfied with the processed FITS file from TILSOTUA, append the binary tables to the trace FITS files. The user must first verify that TILSOTUA has indeed processed the files correctly. This implies:

    #. TILSOTUA has correctly identified the alignment stars (see the QA plot it generates).
    #. TILSOTUA has estimated the `TopDist` and `BotDist` in the `SlitObjMap` table correctly.

If processed correctly, PypeIt should now fully utilize its 
arsenal of slitmask processing tools to reduce and coadd spectra 
with the WCS information incorporated. 
