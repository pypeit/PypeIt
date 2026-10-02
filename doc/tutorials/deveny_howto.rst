.. include:: ../include/links.rst

.. _deveny_howto:

=================
LDT-DeVeny HOWTO
=================

Overview
========

This doc goes through a full run of PypeIt on a night of single-slit
observations with LDT/DeVeny.
The following was performed on a Macbook Pro with 16 GB RAM and took ~25min
for the whole set.

.. _deveny_setup:

Setting Up and Running a PypeIt Reduction
=========================================

This section outlines the highlights of how to use PypeIt with LDT/DeVeny data.
It is a condensed and paraphrased version of the :ref:`cookbook`, which is 
routinely updated and should be referenced for complete and detailed
instructions.

.. note::

    Before you get too far, it is important to understand that PypeIt
    reorients all 2D image data (from any spectrograph) so that the spectral
    axis is vertical with increasing wavelength corresponding to increasing
    pixel number. In the case of DeVeny data, this amounts to a 90\
    :math:`^\circ` CW rotation of the images with respect to the original
    files. Don't Panic!

.. tip::

   At the :ref:`bottom of this page <deveny_workflow>` there is a "cheat sheet"
   of common DeVeny PypeIt workflows.


Planning Your Observations for Reduction with PypeIt
----------------------------------------------------

Because PypeIt is an "end-to-end" data reduction pipeline with minimal
opportunity to interact with the reduction in progress, pre-telescope
planning is required to obtain the proper calibration frames. While most
observing programs will already collect all of the frames necessary for
smooth operation of the pipeline, several items bear pointing out:

-  Bias frames are used to remove fixed-pattern noise in the data and
   generate the default bad pixel mask for reductions.

-  Dome flats are used for the dual purposes of removing pixel-to-pixel
   variations in sensitivity and tracing the edges of the slit.  (Slit edges 
   can vary from grating to grating, and will be more apparent following
   the future installation of the decker). Dome flats (or, optionally,
   sky flats) also may be used for correcting for the variable illumination
   function along the slit (generally < 1% variation), **but this feature
   must be explicitly turned on during reduction**.

-  Wavelength calibration is the piece most likely to cause headaches
   for any spectroscopy program. The user needs to decide which
   combination of lamps will provide suitable calibration for their program.
   PypeIt performs an unclipped mean combine of specified arc frames into a
   single ``Arc`` calibration file. As of ``v1.13.0``, the default DeVeny
   parameters allow for combination of frames taken with individual lamps
   (*e.g.*, separate Hg- and Ar-only frames), as well as multilamp frames
   **but you must modify parameters following** :ref:`deveny_wavecalib`.

   The selected slit width also plays into how well PypeIt matches the
   calibration spectrum with the corresponding line lists. While it is
   sometimes possible to attempt calibration on arc frames taken with a
   wide slit opening (>2"), for best results use arc spectra taken with an
   optimal slit width (*i.e.*, projected slit width on the detector of 2.5 -
   3.0 pixels) to ensure matching by the automated algorithms.

   Additionally, because of the spectral-direction flexure of the DeVeny
   camera, **do not attempt** to combine comparison arc frames from **different
   telescope positions**. The shift in line positions between positions will
   create a hot-mess calibration frame and the wavelength calibration will
   fail. PypeIt's flexure-correction algorithm (see :ref:`deveny_flexure`)
   uses night-sky lines to adjust the wavelength calibration for individual
   science frames, so the use of *in situ* arcs may not be necessary. It
   is, however, possible to correlate individual science frames with
   individual arc images, and this is discussed under :ref:`deveny_groups`
   -- if you go this route, we suggest you also compare it with the
   single-pointing arcs and PypeIt's flexure correction and let LDT staff know
   how well they do.

   .. As of ``v1.14.1``, complete wavelength templates using the Hg, Cd, and Ar
   .. lamps are available for the 150g/mm (DV1), 300g/mm (DV2, DV3),
   .. 500g/mm (DV5), 600g/mm (DV6, DV7), and 1200g/mm (DV9) gratings.
   .. PypeIt will automatically match calibration spectra from these
   .. gratings against the appropriate template using its ``full_template``
   .. method. For the other two gratings (DV4 and DV8), PypeIt attempts to
   .. identify the lines in the spectrum *ex nihilo* using its ``holy-grail``
   .. method.

-  All 9 in-service gratings have been tested with PypeIt and appropriate
   grating-specific parameters have been included in the ``v1.15.0`` release.
   If you have issues with the pipeline crashing or incorrect reduction of your
   data, please contact LDT staff for troubleshooting.

-  *To ensure your calibrations will work with PypeIt*, test the pipeline on a
   preexisting data set whose calibration frames were taken in the same way you
   expect to take them. If this testing is done ahead of time, it will save
   much frustration later. It is also possible to run the pipeline on-the-fly
   on your observing night to ensure you have collected a workable calibration
   set.

.. _`deveny_organize`:

Organize the Data to be Reduced
-------------------------------

Download a single night's data from the site computers to your reduction
machine, as described in the DeVeny User Manual. The easiest method is using
secure copy (``scp``), but feel free to use whatever method you prefer\ [2]_.

Be sure your data directory includes calibration frames taken using the same
grating, tilt, and rear filter (order blocking) settings as your science data.
Focus frames may be present or deleted -- PypeIt will ignore them. You should have:

-  Bias frames (to remove fixed-pattern noise from the data)

-  Dome Flat frames (to remove pixel-to-pixel sensitivity variations and trace
   the slit edges)

-  Comparison Arc frames (for wavelength calibration)

-  Science frames (the whole point)

Optionally, depending on the science requirements of your program, you
may also include:

-  Spectrophotometric Standard Star frames (for flux calibration)

-  Sky Flat frames (to correct for variations in illumination along the
   slit, seldom required but may be applicable to certain science programs)

This raw data directory is the root of the directory tree PypeIt uses
for organizing the processing files and processed data (see
:ref:`deveny_filestructure`). Make sure it is on a local drive (rather than
network storage) for speed and efficiency, since PypeIt reads and writes files
frequently during the reduction process.


Setup
-----

.. tip::

   All PypeIt command-line scripts (*e.g.*, ``pypeit_setup``) will display
   script usage by calling the script with the ``-h`` option.

Running PypeIt on a set of data is controlled by a :ref:`pypeit_file` that
details what the software should do to each file along the way to producing
reduced and calibrated data. PypeIt determines unique instrument
configurations, and sorts data in preparation for the data reduction.

The package provides two setup methods (:ref:`automated <pypeit_setup>` and
GUI) -- both available through the ``pypeit_setup`` script -- that read through
the FITS headers in the raw directory to generate the Reduction File (and
directory tree) based on what it finds.  The Setup GUI (call using
``pypeit_setup -G``) provides the ability to interactively produce a PypeIt
Reduction File, which is very helpful for first-time users.

For DeVeny data, instrument configurations are defined by unique combinations
of the grating (FITS keyword ``GRATING``), grating tilt angle (``GRANGLE``),
rear order-blocking filter (``FILTREAR``), and CCD binning (``CCDSUM``). PypeIt
maps the various DeVeny FITS keywords onto a set of internal metadata keys for
processing. The relevant PypeIt metadata keys for DeVeny configurations (which
you will see in your reduction files) are:

::

       Metadata Key   FITS Header
       ------------   -----------
           dispname       GRATING
            cenwave       GRANGLE
            filter1      FILTREAR
            binning        CCDSUM

The PypeIt metadata key ``cenwave`` is the computed central wavelength of the
spectrum in Angstroms, derived from the grating and tilt angle,  rounded to the
nearest 5\ :math:`\mathring{A}`.

#. **Run** ``pypeit_setup``

   The first run will produce the setup files that should be inspected to
   ensure the code has properly divvied up the FITS files into the proper
   configuration(s). For most DeVeny programs (a single grating tilt and rear
   filter used with the installed grating), should find a single instrument
   configuration. Run the script:

   ::

      $ pypeit_setup -s ldt_deveny

   where the required command-line option ``-s`` sets the spectrograph
   configuration parameters.


#. **Inspect the Outputs**

   The ``ldt_deveny.obslog`` file should somewhat resemble your own
   time-ordered observing log for this set of data, with the relevant FITS
   keywords mapped to their PypeIt metadata keys. This is a good time to ensure
   that all the files you expect to see are in fact present.

   Any collimator focus frames (which you should have identified with
   FITS header keyword ``IMAGETYP = FOCUS``) will have a ``frametype`` listed
   as ``None`` in this file and are commented out. If there are
   non-focus frames with ``frametype None`` listed, this indicates the FITS
   keyword was not correctly set. You should note the affected frames so that
   you can later edit the relevant PypeIt Reduction File(s)
   (:ref:`deveny_edit`) with the correct frame type.

   The ``ldt_deveny.sorted`` file is divided into sections enumerating the
   unique instrument configurations and the list of frames associated
   therewith. Each unique configuration is given a capital letter identifier
   (A, B, C, D...).  Below are example headers from a file for LDT/DeVeny data taken with
   two different order-blocking filters on the same night:

   ::

           ##########################################################
           Setup A
                dispname: DV1 (150/5000)
                 cenwave: 7220.0
                 filter1: Clear
                 binning: 1,1

   ::

           ##########################################################
           Setup B
                dispname: DV1 (150/5000)
                 cenwave: 7220.0
                 filter1: OG570
                 binning: 1,1

   PypeIt does not use this file to guide reductions, but it is provided **as a
   means for the user to assess the automated setup, identification, and file
   sorting**. If, at the start of your observing session, you did not select the
   grating or rear filter in the LOUI before taking exposures, those frames
   will have ``UNKNOWN`` listed in the associated header field. In this case,
   you should go back and edit the FITS headers with the proper values and
   rerun step #1 above.

   The ``ldt_deveny.calib`` file enumerates all of the PypeIt frame types found,
   the calibration files associated therewith, and the raw data frames combined
   to produce them. This version of the calibration association file is
   informational only, but it may be helpful for thinking about grouping frames
   into separate calibration groups, if necessary (see :ref:`deveny_groups`).

#. **Run** ``pypeit_setup`` **again**

   Provided you are happy with the ``ldt_deveny.sorted`` file, you are ready to
   write the ``.pypeit`` file(s) for one or more setups. Executing the
   ``pypeit_setup`` script a second time with the ``-c`` option will create one
   or more sub-folders and populate each with a :ref:`pypeit_file`. See
   :ref:`the setup documentation <setup_doc>` for details on the various options
   available for use with this script.

   An example execution that only produces setup files for the ``A``
   configuration is:

   ::

      $ pypeit_setup -s ldt_deveny -c A

   This will generate a subfolder ``ldt_deveny_A`` containing two files: the
   base PypeIt Reduction File ``ldt_deveny_A.pypeit``, and its calibration
   association file ``ldt_deveny_A.calib``.


.. _`deveny_edit`:

Edit Your PypeIt File
---------------------

The :ref:`pypeit_file` dictates how the pipeline is executed on your raw data
files. While you just generated the file automatically (above), it can (and
should) be edited by the user to ensure the reduction proceeds as expected.

Each unique instrument configuration will have its own PypeIt Reduction
File. In the case of DeVeny, this means different rear filters, grating
tilt angles, binning schemes, or even different gratings used on different
nights.  See :ref:`the relevant documentation <pypeit_file>` for descriptions
of the file format and common edits a user may wish to make.  

Here is the ``.pypeit`` file for this example:

.. code-block::

    # Auto-generated PypeIt input file using PypeIt version: 1.18.2.dev613+gf221c3831
    # UTC 2026-02-05T17:44:04.030+00:00

    # User-defined execution parameters
    [rdx]
        spectrograph = ldt_deveny

    # Setup
    setup read
    Setup A:
    binning: 1,2
    cenwave: 4900.0
    dispname: DV6 (600/4900)
    filter1: CLEAR
    setup end

    # Data block 
    data read
    path .
                filename |       frametype |                 ra |                dec |                     target |       dispname | binning |                mjd | airmass | exptime | filter1 | dispangle | slitwid | lampstat01 | calib
    # 20221102.0011.fits |            None | 283.12604166666665 |  34.68172222222222 | Collimator Focus  (HgCdAr) | DV6 (600/4900) |     1,2 | 59885.007247800924 |     1.0 |     5.0 |   CLEAR |     27.04 |     1.2 | Cd, Ar, Hg |     0
    # 20221102.0012.fits |            None | 283.12604166666665 |  34.68172222222222 | Collimator Focus  (HgCdAr) | DV6 (600/4900) |     1,2 |  59885.00744490741 |     1.0 |     5.0 |   CLEAR |     27.04 |     1.2 | Cd, Ar, Hg |     0
    # 20221102.0013.fits |            None | 283.12604166666665 |  34.68172222222222 | Collimator Focus  (HgCdAr) | DV6 (600/4900) |     1,2 |  59885.00764201389 |     1.0 |     5.0 |   CLEAR |     27.04 |     1.2 | Cd, Ar, Hg |     0
    # 20221102.0014.fits |            None | 283.12604166666665 |  34.68172222222222 | Collimator Focus  (HgCdAr) | DV6 (600/4900) |     1,2 |  59885.00783900463 |     1.0 |     5.0 |   CLEAR |     27.04 |     1.2 | Cd, Ar, Hg |     0
    # 20221102.0015.fits |            None | 283.12604166666665 |  34.68172222222222 | Collimator Focus  (HgCdAr) | DV6 (600/4900) |     1,2 |  59885.00803611111 |     1.0 |     5.0 |   CLEAR |     27.04 |     1.2 | Cd, Ar, Hg |     0
    # 20221102.0016.fits |            None | 283.12604166666665 |  34.68172222222222 | Collimator Focus  (HgCdAr) | DV6 (600/4900) |     1,2 | 59885.008233333334 |     1.0 |     5.0 |   CLEAR |     27.04 |     1.2 | Cd, Ar, Hg |     0
    # 20221102.0017.fits |            None | 283.12604166666665 |  34.68172222222222 | Collimator Focus  (HgCdAr) | DV6 (600/4900) |     1,2 |  59885.00843043981 |     1.0 |     5.0 |   CLEAR |     27.04 |     1.2 | Cd, Ar, Hg |     0
    # 20221102.0018.fits |            None | 283.12604166666665 |  34.68172222222222 | Collimator Focus  (HgCdAr) | DV6 (600/4900) |     1,2 |  59885.00862789352 |     1.0 |     5.0 |   CLEAR |     27.04 |     1.2 | Cd, Ar, Hg |     0
    # 20221102.0019.fits |            None | 283.12604166666665 |  34.68172222222222 | Collimator Focus  (HgCdAr) | DV6 (600/4900) |     1,2 |  59885.00882488426 |     1.0 |     5.0 |   CLEAR |     27.04 |     1.2 | Cd, Ar, Hg |     0
    20221102.0020.fits |        arc,tilt | 283.12604166666665 |  34.68172222222222 |                HgCdAr Arcs | DV6 (600/4900) |     1,2 | 59885.010723611114 |     1.0 |     5.0 |   CLEAR |     27.04 |     1.2 | Cd, Ar, Hg |     0
    20221102.0021.fits |        arc,tilt | 283.12604166666665 |  34.68172222222222 |                HgCdAr Arcs | DV6 (600/4900) |     1,2 | 59885.010830439816 |     1.0 |     5.0 |   CLEAR |     27.04 |     1.2 | Cd, Ar, Hg |     0
    20221102.0022.fits |        arc,tilt | 283.12604166666665 |  34.68172222222222 |                HgCdAr Arcs | DV6 (600/4900) |     1,2 |  59885.01093726852 |     1.0 |     5.0 |   CLEAR |     27.04 |     1.2 | Cd, Ar, Hg |     0
    20221102.0023.fits |        arc,tilt | 283.12604166666665 |  34.68172222222222 |                HgCdAr Arcs | DV6 (600/4900) |     1,2 | 59885.011044212966 |     1.0 |     5.0 |   CLEAR |     27.04 |     1.2 | Cd, Ar, Hg |     0
    20221102.0024.fits |        arc,tilt | 283.12604166666665 |  34.68172222222222 |                HgCdAr Arcs | DV6 (600/4900) |     1,2 |  59885.01115104167 |     1.0 |     5.0 |   CLEAR |     27.04 |     1.2 | Cd, Ar, Hg |     0
    20221102.0025.fits |        arc,tilt | 283.12604166666665 |  34.68172222222222 |                HgCdAr Arcs | DV6 (600/4900) |     1,2 |  59885.01125787037 |     1.0 |     5.0 |   CLEAR |     27.04 |     1.2 | Cd, Ar, Hg |     0
    20221102.0026.fits |        arc,tilt | 283.12604166666665 |  34.68172222222222 |                HgCdAr Arcs | DV6 (600/4900) |     1,2 |  59885.01136469907 |     1.0 |     5.0 |   CLEAR |     27.04 |     1.2 | Cd, Ar, Hg |     0
    20221102.0027.fits |        arc,tilt | 283.12604166666665 |  34.68172222222222 |                HgCdAr Arcs | DV6 (600/4900) |     1,2 |  59885.01147152778 |     1.0 |     5.0 |   CLEAR |     27.04 |     1.2 | Cd, Ar, Hg |     0
    20221102.0028.fits |        arc,tilt | 283.12604166666665 |  34.68172222222222 |                HgCdAr Arcs | DV6 (600/4900) |     1,2 |  59885.01157835648 |     1.0 |     5.0 |   CLEAR |     27.04 |     1.2 | Cd, Ar, Hg |     0
    20221102.0029.fits |        arc,tilt | 283.12604166666665 |  34.68172222222222 |                HgCdAr Arcs | DV6 (600/4900) |     1,2 | 59885.011685185185 |     1.0 |     5.0 |   CLEAR |     27.04 |     1.2 | Cd, Ar, Hg |     0
    20221102.0001.fits |            bias | 283.12604166666665 |  34.68172222222222 |                       Bias | DV6 (600/4900) |     1,2 | 59885.002777546295 |     1.0 |     0.0 |   CLEAR |      16.5 |  55.853 |        off |     0
    20221102.0002.fits |            bias | 283.12604166666665 |  34.68172222222222 |                       Bias | DV6 (600/4900) |     1,2 |  59885.00282650463 |     1.0 |     0.0 |   CLEAR |      16.5 |  55.853 |        off |     0
    20221102.0003.fits |            bias | 283.12604166666665 |  34.68172222222222 |                       Bias | DV6 (600/4900) |     1,2 |  59885.00287546296 |     1.0 |     0.0 |   CLEAR |      16.5 |  55.853 |        off |     0
    20221102.0004.fits |            bias | 283.12604166666665 |  34.68172222222222 |                       Bias | DV6 (600/4900) |     1,2 |  59885.00292615741 |     1.0 |     0.0 |   CLEAR |      16.5 |  55.853 |        off |     0
    20221102.0005.fits |            bias | 283.12604166666665 |  34.68172222222222 |                       Bias | DV6 (600/4900) |     1,2 |  59885.00297511574 |     1.0 |     0.0 |   CLEAR |      16.5 |  55.853 |        off |     0
    20221102.0006.fits |            bias | 283.12604166666665 |  34.68172222222222 |                       Bias | DV6 (600/4900) |     1,2 |  59885.00302418981 |     1.0 |     0.0 |   CLEAR |      16.5 |  55.853 |        off |     0
    20221102.0007.fits |            bias | 283.12604166666665 |  34.68172222222222 |                       Bias | DV6 (600/4900) |     1,2 | 59885.003073148146 |     1.0 |     0.0 |   CLEAR |      16.5 |  55.853 |        off |     0
    20221102.0008.fits |            bias | 283.12604166666665 |  34.68172222222222 |                       Bias | DV6 (600/4900) |     1,2 | 59885.003122222224 |     1.0 |     0.0 |   CLEAR |      16.5 |  55.853 |        off |     0
    20221102.0009.fits |            bias | 283.12604166666665 |  34.68172222222222 |                       Bias | DV6 (600/4900) |     1,2 |  59885.00317118056 |     1.0 |     0.0 |   CLEAR |      16.5 |     1.2 |        off |     0
    20221102.0010.fits |            bias | 283.12604166666665 |  34.68172222222222 |                       Bias | DV6 (600/4900) |     1,2 |  59885.00322013889 |     1.0 |     0.0 |   CLEAR |      16.5 |     1.2 |        off |     0
    20221102.0072.fits |            bias |  76.38124999999998 |  52.83255555555556 |                       Bias | DV6 (600/4900) |     1,2 |   59885.2970068287 |     1.2 |     0.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0073.fits |            bias |  76.38124999999998 |  52.83255555555556 |                       Bias | DV6 (600/4900) |     1,2 | 59885.297055902774 |     1.2 |     0.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0074.fits |            bias |  76.38124999999998 |  52.83255555555556 |                       Bias | DV6 (600/4900) |     1,2 | 59885.297104861114 |     1.2 |     0.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0075.fits |            bias |  76.38124999999998 |  52.83255555555556 |                       Bias | DV6 (600/4900) |     1,2 |  59885.29715381945 |     1.2 |     0.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0076.fits |            bias |  76.38124999999998 |  52.83255555555556 |                       Bias | DV6 (600/4900) |     1,2 |  59885.29720277778 |     1.2 |     0.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0077.fits |            bias |  76.38124999999998 |  52.83255555555556 |                       Bias | DV6 (600/4900) |     1,2 |  59885.29725173611 |    1.19 |     0.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0078.fits |            bias |  76.38124999999998 |  52.83255555555556 |                       Bias | DV6 (600/4900) |     1,2 |  59885.29730081018 |    1.19 |     0.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0079.fits |            bias |  76.38124999999998 |  52.83255555555556 |                       Bias | DV6 (600/4900) |     1,2 | 59885.297349768516 |    1.19 |     0.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0080.fits |            bias |  76.38124999999998 |  52.83255555555556 |                       Bias | DV6 (600/4900) |     1,2 |  59885.29739872685 |    1.19 |     0.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0081.fits |            bias |  76.38124999999998 |  52.83255555555556 |                       Bias | DV6 (600/4900) |     1,2 |  59885.29744768519 |    1.19 |     0.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0030.fits | pixelflat,trace | 283.12604166666665 |  34.68172222222222 |                 Dome Flats | DV6 (600/4900) |     1,2 | 59885.014693171295 |    1.47 |    20.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0031.fits | pixelflat,trace | 283.12604166666665 |  34.68172222222222 |                 Dome Flats | DV6 (600/4900) |     1,2 |  59885.01497372685 |    1.47 |    20.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0032.fits | pixelflat,trace | 283.12604166666665 |  34.68172222222222 |                 Dome Flats | DV6 (600/4900) |     1,2 |  59885.01525416667 |    1.47 |    20.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0033.fits | pixelflat,trace | 283.12604166666665 |  34.68172222222222 |                 Dome Flats | DV6 (600/4900) |     1,2 | 59885.015534722224 |    1.47 |    20.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0034.fits | pixelflat,trace | 283.12604166666665 |  34.68172222222222 |                 Dome Flats | DV6 (600/4900) |     1,2 | 59885.015815162034 |    1.47 |    20.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0035.fits | pixelflat,trace | 283.12604166666665 |  34.68172222222222 |                 Dome Flats | DV6 (600/4900) |     1,2 |  59885.01609560185 |    1.47 |    20.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0036.fits | pixelflat,trace | 283.12604166666665 |  34.68172222222222 |                 Dome Flats | DV6 (600/4900) |     1,2 |  59885.01637604167 |    1.47 |    20.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0037.fits | pixelflat,trace | 283.12604166666665 |  34.68172222222222 |                 Dome Flats | DV6 (600/4900) |     1,2 |  59885.01665648148 |    1.47 |    20.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0038.fits | pixelflat,trace | 283.12604166666665 |  34.68172222222222 |                 Dome Flats | DV6 (600/4900) |     1,2 | 59885.016936921296 |    1.47 |    20.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0039.fits | pixelflat,trace | 283.12604166666665 |  34.68172222222222 |                 Dome Flats | DV6 (600/4900) |     1,2 |  59885.01721747685 |    1.47 |    20.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0040.fits | pixelflat,trace | 283.12604166666665 |  34.68172222222222 |                 Dome Flats | DV6 (600/4900) |     1,2 |  59885.01749791667 |    1.47 |    20.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0041.fits | pixelflat,trace | 283.12604166666665 |  34.68172222222222 |                 Dome Flats | DV6 (600/4900) |     1,2 |  59885.01777835648 |    1.47 |    20.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0042.fits | pixelflat,trace | 283.12604166666665 |  34.68172222222222 |                 Dome Flats | DV6 (600/4900) |     1,2 |   59885.0180587963 |    1.47 |    20.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0043.fits | pixelflat,trace | 283.12604166666665 |  34.68172222222222 |                 Dome Flats | DV6 (600/4900) |     1,2 |  59885.01833935185 |    1.47 |    20.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0044.fits | pixelflat,trace | 283.12604166666665 |  34.68172222222222 |                 Dome Flats | DV6 (600/4900) |     1,2 |  59885.01861979167 |    1.47 |    20.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0045.fits | pixelflat,trace | 283.12604166666665 |  34.68172222222222 |                 Dome Flats | DV6 (600/4900) |     1,2 |  59885.01890023148 |    1.47 |    20.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0046.fits | pixelflat,trace | 283.12604166666665 |  34.68172222222222 |                 Dome Flats | DV6 (600/4900) |     1,2 |   59885.0191806713 |    1.47 |    20.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0047.fits | pixelflat,trace | 283.12604166666665 |  34.68172222222222 |                 Dome Flats | DV6 (600/4900) |     1,2 | 59885.019461111115 |    1.47 |    20.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0048.fits | pixelflat,trace | 283.12604166666665 |  34.68172222222222 |                 Dome Flats | DV6 (600/4900) |     1,2 |  59885.01974166666 |    1.47 |    20.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0049.fits | pixelflat,trace | 283.12604166666665 |  34.68172222222222 |                 Dome Flats | DV6 (600/4900) |     1,2 |  59885.02002210648 |    1.47 |    20.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0050.fits |         science |          327.79725 | 28.861611111111113 |                 BD+28 4211 | DV6 (600/4900) |     1,2 |  59885.06604016203 |    1.03 |   600.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0051.fits |         science |  307.3576666666666 |  40.23319444444444 |                  B03_159_A | DV6 (600/4900) |     1,2 |  59885.07566712963 |    1.01 |   600.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0052.fits |         science |         307.366125 |  40.21288888888889 |                  B03_159_C | DV6 (600/4900) |     1,2 |  59885.08816435185 |    1.02 |   600.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0053.fits |         science |  307.2952916666666 |  40.16397222222222 |                  B03_159_E | DV6 (600/4900) |     1,2 |  59885.09651574074 |    1.03 |   600.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0054.fits |         science |  307.3106666666667 |  40.16686111111111 |                  B03_159_I | DV6 (600/4900) |     1,2 |  59885.10378668981 |    1.04 |   600.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0055.fits |         science |          307.32025 |  40.21702777777778 |                  B03_159_O | DV6 (600/4900) |     1,2 |  59885.11452175926 |    1.06 |   600.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0056.fits |         science |  307.3915416666666 |  40.21613888888889 |                  B03_159_P | DV6 (600/4900) |     1,2 |  59885.12204780093 |    1.07 |   600.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0057.fits |         science | 307.38879166666663 |           40.17725 |                  B03_159_Q | DV6 (600/4900) |     1,2 | 59885.179968055556 |    1.26 |   600.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0058.fits |         science |  307.3976666666666 | 40.169555555555554 |                  B03_159_R | DV6 (600/4900) |     1,2 | 59885.187330671295 |    1.29 |   600.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0059.fits |         science |         307.302875 |  40.20963888888889 |                  B03_159_T | DV6 (600/4900) |     1,2 | 59885.194894444445 |    1.33 |   600.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0060.fits |         science |  307.2955416666666 |  40.15897222222222 |                  B03_159_U | DV6 (600/4900) |     1,2 |  59885.20247199074 |    1.38 |   600.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0061.fits |         science |          36.354375 |  62.00716666666667 |                   B03_54_G | DV6 (600/4900) |     1,2 |  59885.21155717593 |     1.2 |   600.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0062.fits |         science | 36.325041666666664 | 62.135527777777774 |                   B03_54_H | DV6 (600/4900) |     1,2 |  59885.21962025463 |    1.18 |   600.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0063.fits |         science | 36.420958333333324 |  62.07527777777778 |                   B03_54_L | DV6 (600/4900) |     1,2 | 59885.226822453704 |    1.17 |   600.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0064.fits |         science | 36.323499999999996 |  62.09730555555556 |                   B03_54_M | DV6 (600/4900) |     1,2 | 59885.234117708336 |    1.16 |   600.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0065.fits |         science | 36.451333333333324 |  62.04697222222222 |                   B03_54_P | DV6 (600/4900) |     1,2 |  59885.24142002315 |    1.16 |   600.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0066.fits |         science |  36.37491666666666 |  62.11938888888889 |                   B03_54_Q | DV6 (600/4900) |     1,2 | 59885.248722337965 |    1.15 |   600.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0067.fits |         science |  349.0227083333333 |  60.04413888888889 |                   B03_43_A | DV6 (600/4900) |     1,2 |  59885.25696423611 |    1.19 |   600.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0068.fits |         science | 348.91866666666664 |  60.05127777777778 |                   B03_43_B | DV6 (600/4900) |     1,2 |  59885.26450902778 |    1.21 |   600.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0069.fits |         science |  349.1099166666666 |  60.02152777777778 |                   B03_43_C | DV6 (600/4900) |     1,2 |  59885.27215428241 |    1.23 |   600.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0070.fits |         science |  348.9552083333333 | 60.027055555555556 |                   B03_43_F | DV6 (600/4900) |     1,2 | 59885.279801157405 |    1.25 |   600.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    20221102.0071.fits |         science |  76.38133333333332 |  52.83241666666667 |                   G191-B2B | DV6 (600/4900) |     1,2 |  59885.28974270833 |    1.22 |   600.0 |   CLEAR |     27.04 |     1.2 |        off |     0
    data end




Specific LDT/DeVeny considerations:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#. The DeVeny-specific modifications to default PypeIt reduction parameters are
   already included in the
   :class:`~pypeit.spectrographs.ldt_deveny.LDTDeVenySpectrograph` class and
   loaded using the ``spectrograph = ldt_deveny`` line at the top of the
   :ref:`parameter_block` -- it is not necessary to reproduce all those
   parameters in the :ref:`parameter_block` of your file. What do go here are
   changes **away** from the *DeVeny default configuration* you wish to use for
   reducing a particular data set. For instance, to specify that PypeIt should
   use the ``illumflat`` files to correct for illumination variations along the
   slit and that it should only find and extract the one brightest object in
   each science frame, you would add the following to your parameter block:

   .. code-block:: ini

      [baseprocess]
         use_illumflat = True
      [reduce]
         [[findobj]]
            maxnumber_sci = 1

   A discussion of typical parameter changes that may apply to DeVeny data
   is given at :ref:`deveny_parmods`, and an exhaustive discussion of all
   parameters may be found at :ref:`parameters`.

#. Here is yet another reminder to **not include bad calibration frames** in
   the reduction (frames that you do not want to use, frames with incorrectly
   identified types, or frames that could not be automatically classified and
   have a ``None`` type). Check them now and remove or comment them out if
   they are bad.

   You may need to add configuration-independent files from one setup to the
   Data Block of another, but PypeIt is getting better at including
   setup-independent files in all configurations.  In any event, it is
   important to double-check that all files needed for the reduction are
   present in your ``.pypeit`` file. If needed, simply copy the needed lines
   from one file to the other so that both setups have access to, *e.g.*,
   the bias frames. The ordering of table rows in the :ref:`pypeit_file` does
   not matter, so don't worry about adding lines in the "proper" location.

#. Check the ``frametype`` of all files. For DeVeny reductions, you need at
   least one file with each of the following Frame Types (see
   :ref:`deveny_organize`):

   -  ``bias``: Bias frames (removing fixed-pattern noise)

   -  ``pixelflat,trace``: Flat fielding (removing pixel-to-pixel sensitivity variations)
      and edge tracing

   -  ``arc,tilt``: Two-dimensional wavelength calibration (colorizing the
      black-and-white spectrum)

   -  ``science``: Science exposure (answering the grand questions of the universe)

   Remove / comment out all images with a ``frametype`` of ``None``, or correct
   the value.  PypeIt will NOT run if any of the uncommented frames have
   ``None`` under ``frametype``.

   Additionally, frames may have type ``illumflat`` if you are doing
   illumination corrections along the slit.  While other spectrographs support
   ``standard`` star frames, at the moment, DeVeny does not need anything of
   this type, and spectrophotometric standards should be marked as ``science``
   frames.

   .. tip::
      A given image can have multiple frame types (*e.g.*, ``arc,tilt``).
      Simply enter the types as a comma-separated list without spaces.


#. Check ``target`` names for all files for both accuracy and for illegal
   characters.  The ``target`` name is used as part of the reduced data
   filename -- accurate names help identify objects later.

   .. important::

      Because PypeIt uses the ``target`` name (pulled from the ``OBJNAME`` FITS
      keyword, entered by the observer in the DeVeny LOUI) as part of the
      reduced data filename, this column **must** include only legal characters
      for your filesystem. In general, forward slash (``/``) is always
      disallowed (sorry, comet and interstellar object observers), but other
      characters may be a concern on your particular filesystem. Additionally,
      parentheses or other characters in ``target`` names may cause issues if
      such characters are not escaped in shell environments.  Editing the name
      in the PypeIt Reduction File (and not in the actual FITS file itself) is
      sufficient for the limitations mentioned here.


#. Adjust the ``calib`` groupings for calibration associations. See
   :ref:`calibration-groups` for an exhaustive discussion or
   :ref:`deveny_groups` for a more tailored outline.  For LDT/DeVeny,
   care must be exercised in grouping arc frames for wavelength calibration.
   Given the large shifts along the spectral axis of the DeVeny CCD caused by
   flexure (:math:`\sim\pm` 10 pixels), some observers prefer to take *in situ* arcs at the
   location of each object rather than rely upon PypeIt's flexure correction
   based on night sky lines (see :ref:`deveny_flexure` for a discussion of
   flexure corrections). The ensemble of *in situ* arcs should definitely
   **not** be grouped together for PypeIt wavelength calibration, as the
   flexure-induced shifts between frames can produce an unusable mess of
   multiple, shifted lines. The safe move here is to assign each set of frames
   at a given pointing (zenith, object A, standard C, etc.) to a unique
   calibration group.


Run the Reduction
-----------------

PypeIt is designed (and currently only able) to do end-to-end reductions,
resulting in a fully processed 2D spectral image and extracted 1D spectra (if
any objects were found) from each science frames. Once you have completed the
setup steps above, you are just about ready to run the pipeline.

The script to run a reduction is :ref:`run-pypeit`.  See that documentation
page for all relevant script options and workflows.

.. caution::

   When you upgrade PypeIt versions, changes to the underlying data models
   (which are largely not backwards compatible) may cause errors if you try to
   use calibration files processed with an earlier version. The safe move is to
   completely reprocess all data currently being used when PypeIt is upgraded,
   including deleting and recreating all processed calibrations.  Your
   currently installed version of PypeIt may be checked using
   ``pypeit_version``, and the version used to create any output file is listed
   in the FITS header with the keyword ``VERSPYP``.


.. _deveny_outputs:

Primary Output Files and Post-Processing Scripts
================================================

.. _deveny_calibrations:

Examine the Calibration Files
-----------------------------

As PypeIt begins churning through your reduction, it will create and write to
disk calibration frames in the ``Calibrations/`` subfolder of the, *e.g.*,
``ldt_deveny_A/`` directory (see :ref:`deveny_filestructure`). Additional
Quality Assurance files will be written to the ``QA/`` subfolder for some
types of Calibration frames. It is important to take the time to inspect these
calibration outputs as they are generated.

.. tip::
   
   To process just the calibrations without trying to process the ``science``
   data, use the command::

      run_pypeit -c ldt_deveny_<setup>.pypeit

The naming convention for :ref:`calibrations` frames is a bit cumbersome, but
follows a regular pattern.  Here is a brief listing of the Calibration frames
produced (in the order in which they are created):

-  :ref:`bias` -- Processed combined bias frame used to remove
   fixed-pattern noise from all other images.

-  :ref:`edges` -- Collection of images and FITS binary tables
   describing the slit traces. While this file is primarily of interest
   for multislit or echelle spectrographs (DeVeny has but one slit and
   no cross-disperser, after all), it is instructive to quickly peek at
   this file to ensure the code correctly identified the slit (and not
   some artifact at the edge of the CCD):

   ::

      $ pypeit_chk_edges Calibrations/Edges_A_0_DET01.fits.gz

   This command will launch a :ref:`GUI viewer <pypeit_chk_edges>` to display
   the combined trace image along with a sobel-filtered version used to
   identify illumination discontinuities in the spatial direction (see figure
   below). For DeVeny data, it should identify a single, long slit with a
   (spatial ID) approximately half the spatial extent of the CCD image
   (mid-200s for spatially unbinned data). The exact number will vary from
   grating to grating due to differing small roll angles about the dispersion
   axis when the gratings were installed in their cells.

.. _deveny_edges_DV2:
.. figure:: ../figures/deveny_edges_DV2.png
   :alt: DV2 edges example
   :scale: 50
   :class: with-shadow
         
   Example of output from the ``pypeit_chk_edges`` script for data taken with
   the DV2 grating. The green and magenta lines in the center panels mark
   the left and right edges of the detected slit, respectively. The CCD is
   about 2.9' wide, so at least one edge of the 2.5' slit should be
   visible.

-  :ref:`slits` -- This file contains the distilled PypeIt-internal information
   on the traced slit edges, derived from the ``Edges`` file and organized in
   FITS binary tables. The best way to assess these data is in the
   ``pypeit_chk_edges`` GUI. Once again, there should only be one slit for
   DeVeny data.

-  :ref:`arc` -- Processed combined arc spectral image, where the frames are
   combined using an unclipped mean combine algorithm. Closely examine this
   image in a tool like ``ds9`` to ensure it will be suitable for generating a
   wavelength solution. If not, try editing the calibration group information
   in the PypeIt Reduction File to include only a subset of the arc frames
   taken at the same telescope position and rerunning ``run_pypeit``.

-  :ref:`tiltimg` -- Image used to trace the tilting of spectral lines across
   the slit traces to produce an accurate 2D wavelength solution for the
   detector. For the case of DeVeny (single slit trace on the sole detector),
   this is identical to the ``Arc`` image.

-  :ref:`wave_calib` -- Contains the 1D wavelength solution for this setup.
   Inspect the wavelength solution using the ``pypeit_chk_wavecalib`` script.
   Below is an example output from data taken with the DV2 grating
   (:math:`\theta_{\rm grangle} = 22.54^\circ`, :math:`\lambda_c = 5195\mathring{A}`):

   ::

      $ pypeit_chk_wavecalib Calibrations/WaveCalib_A_0_DET01.fits

         N. SpatID minWave Wave_cen maxWave dWave Nlin     IDs_Wave_range    IDs_Wave_cov(%) mesured_fwhm  RMS
        --- ------ ------- -------- ------- ----- ---- --------------------- --------------- ------------ -----
          0    276  2924.1   5151.2  7385.8 2.173   19  3132.752 -  7274.940            92.8          4.8 0.141

   The central wavelength and wavelength range should be close to what you set
   using values from the LOUI and `obstools
   <https://lowellobservatory.github.io/LDTObserverTools/deveny_grangle.html>`__
   package.  The dispersion (``dWave``) should be close to the value listed in
   the DeVeny Users Manual for the selected grating. Note that the ``SpatID``
   listed here should match that from ``pypeit_chk_edges``.

-  :ref:`tilts` -- Contains the 2D mapping of the slit to lines of constant
   wavelength. The quality of this step is shown in the images of the
   ``QA/PNGs`` directory (examples below), and should rarely need much scrutiny
   for DeVeny data if you have strong arc lines and a good wavelength solution.

.. grid:: 3

   .. grid-item::
      :columns: 4

      .. image:: ../figures/deveny_tilts_2d.png
         :alt: 2D arc tilts
         :class: with-shadow

   .. grid-item::
      :columns: 7

      .. image:: ../figures/deveny_tilts_spat.png
         :alt: Spatial tilt residuals
         :class: with-shadow

   .. grid-item::
      :columns: 12

      .. image:: ../figures/deveny_tilts_spec.png
         :alt: Spectral tilt residuals
         :class: with-shadow

   .. grid-item::
      :columns: 12

      Example PypeIt QA plots for the ``Tilts`` file associated with the
      example DV2 data set.

-  :ref:`flat` -- Processed combined dome flat fields for removing
   pixel-to-pixel sensitivity variations. PypeIt fits a basis spline
   (``bspline``) to the spectral direction to remove the structure in the flat
   lamp spectra, and should yield a normalized image with all values close to
   unity. Examine the normalized flat field frame using the
   ``pypeit_chk_flats`` utility.  The GUI also shows the 2D wavelength solution
   derived from when you mouse over the various images. This is a good guide
   for determining whether artifacts seen in the flats are caused by low
   signal at extreme wavelengths.


Examine the Science Spectra
---------------------------

As PypeIt runs, it will begin generating 2D and 1D spectra outputs in the
``Science/`` folder for each science frame in the PypeIt Reduction File. Feel
free to examine the files as they are created, even while the code continues to
process the other raw frames.

Examine the 2D Spectral Images
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   During the data-reduction process, PypeIt will create a reduced 2D spectral
   image product for each science frame prior to the extraction of 1D spectra.
   These products are stored in multi-extension FITS files with names like:

   ::

      spec2d_20221102.0069-B03_43_C_DeVeny_20221102T063154.130.fits

   The complete description of these files is given at :ref:`spec-2d-output`,
   including how to use the viewing tool ``pypeit_show_2dspec``.  An example
   GUI view of an object observed with the DV6 grating is shown below.

   .. figure:: ../figures/deveny_spec2d.png
      :alt: DV6 2D reduced spectrum example
      :width: 100%
      :class: with-shadow

      Example of the PypeIt reduced 2D spectrum for an object observed with DV6
      displayed with the ``pypeit_show_2dspec`` script. *Top Left:* the
      calibrated science image, *top right:* sky-subtracted and masked image
      along the slit bounds (green and magenta lines), *bottom left:* the
      sky-subtracted image divided by the pixel-by-pixel uncertainty to yield a
      residual map including the object, *bottom right:* the same residual map
      but with the object subtracted.  Note that three objects have been
      identified and extracted (orange traces and labels).

   PypeIt names each extracted object by its spatial position on the reduced image
   [``SPAT``], slit position on the reduced image [``SLIT``] and the detector
   number [``DET``]. For instance, the three objects shown above have the
   labels ``SPAT0033-SLIT0126-DET01``, ``SPAT0128-SLIT0126-DET01``, and
   ``SPAT0231-SLIT0126-DET01``. The single-slit nature of DeVeny means that
   multiple objects extracted from a given image will have names differing only
   in the ``SPAT`` code.

Examine the Extracted 1D Spectra
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   If one or more objects have been automatically or manually identified in the
   reduced 2D spectral image, 1D data products will be produced. These 1D
   products are the primary outputs of PypeIt, and consist of a series of
   1-dimensional arrays: vacuum wavelength, extracted flux (using one or more
   methods), and associated error arrays for each identified object. These
   arrays are packaged into multi-extension FITS files, and are accompanied by
   ``.txt`` files with extraction information (*read*: table of contents) for
   each 1D spectrum.

   The 1D spectral FITS files have names like:

   ::

      spec1d_20221102.0069-B03_43_C_DeVeny_20221102T063154.130.fits

   The complete description of these files is given at :ref:`spec-1d-output`,
   including how to use the viewing tool ``pypeit_show_1dspec``.  An example
   GUI view of the object described above is shown below.

   .. figure:: ../figures/deveny_spec1d.png
      :alt: DV6 1D reduced spectrum example
      :width: 100%
      :class: with-shadow

      Example of the PypeIt reduced and extracted 1D spectrum for the brightest
      object shown in the 2D spectrum above. The red dotted line indicates the
      1-:math:`\sigma` uncertainty in the flux values.

   The accompanying ``.txt`` file contains information about the extracted
   object(s), including FWHM of the optimal extraction in arcseconds (this
   should be similar to the seeing on the observing night, convolved with
   jitter in the star position along the slit), the SNR of the extracted
   spectrum (useful in identifying spurious objects), and the RMS in pixels
   of the wavelength solution (for DeVeny should be the same for every object):

   ::

      | slit |                    name | spat_pixpos | spat_fracpos | box_width | opt_fwhm |    s2n | wv_rms |
      |  120 | SPAT0033-SLIT0120-DET01 |        32.9 |        0.117 |      3.80 |    2.188 |   9.13 |  0.065 |
      |  120 | SPAT0128-SLIT0120-DET01 |       128.2 |        0.535 |      3.80 |    2.121 | 146.12 |  0.065 |
      |  120 | SPAT0231-SLIT0120-DET01 |       231.2 |        0.984 |      3.80 |    1.641 |   5.17 |  0.065 |

   By default, ``pypeit_show_1dspec`` loads the first (lowest ``SPAT`` code)
   object extracted from the 2D spectrum. Examination of the spectral image with
   ``pypeit_show_2dspec`` or printing the ``.txt`` file will help you identify
   which extracted object(s) corresponds to your desired target(s). If there
   are spurious low-signal objects identified, you may re-run the reduction
   with adjusted object-finding parameters (see :ref:`deveny_objfind`). A
   particular extracted object may be loaded by using the ``--obj`` option to
   :ref:`pypeit_show_1dspec`.

   By default, PypeIt performs both a boxcar (top-hat) extraction around the
   trace and a Horne optimal extraction\ [3]_ using the fitted spatial profile.
   The boxcar-extracted spectrum may be displayed using the ``--extract BOX``
   option to ``pypeit_show_1dspec``, otherwise the optimal extraction is
   displayed (if available).

.. _deveny_missing1d:

Missing 1D Spectra
^^^^^^^^^^^^^^^^^^

   Sometimes PypeIt will not extract all (or any) of the objects you expect to
   be in a given frame. This can look like either:

   -  some, but not all, of the expected objects were found and extracted
      orange traces on the images of ``pypeit_show_2dspec``) and the
      ``spec1d`` file has fewer entries than expected, or

   -  no objects were found and no ``spec1d`` file was created.

   In either of these cases, the steps for attempting to extract such
   missing objects are the same:

   #. You may modify the object finding parameters in your PypeIt Reduction
      File (see :ref:`deveny_objfind`), remove this ``spec2d_*.fits`` file, and
      rerun ``run_pypeit`` **without the** ``-o`` **option**. This will have
      the effect of processing only the one frame, and should run fairly
      quickly. If the missing objects are found, you're done.

   #. If the objects are still not extracted with repeated parameter
      modification, you can attempt to :ref:`manually identify and extract
      <manual>` the object.  For the example 2D spectrum described above, to 
      manually extract the faint object between the left two identified
      objects, the ``manual`` column for this frame would read
      ``1:77.5:1000:3.1``, where the FWHM (*in pixels*) is the value of
      extracted objects listed in the ``spec1d`` text file divided by the
      spatial plate scale of the spectral image (the DeVeny plate scale of
      0.34"/pixel times the spatial binning).
      
      The resulting 2D spectral image with the manual trace and 1D spectrum
      of the manually extracted object are shown below.  In this case, even
      though the object is detectable to the human eye, it does not contain
      enough signal to produce a useable spectrum (SNR ~ 1).

   .. grid:: 2

      .. grid-item::
         :columns: 6

         .. image:: ../figures/deveny_spec2d_manual.png
            :alt: Manual extraction in a spec2d file
            :class: with-shadow

      .. grid-item::
         :columns: 6

         .. image:: ../figures/deveny_spec1d_manual.png
            :alt: Manual extraction in a spec1d file
            :class: with-shadow

      .. grid-item::
         :columns: 12

         Example of PypeIt manual extraction. The left panel is the 2D spectrum
         with the manually object identified in blue, and the right panel is its
         extracted spectrum.


Post-Processing the Files
-------------------------

While the main PypeIt run ends with ``spec1d`` files, this is not the end of
the processing available with the package. There are several
:ref:`post-processing steps <further_proc_scripts>` that may be considered,
depending on the needs of your particular science program:

-  :ref:`Coadding 2D spectral images <coadd2d>` of the same target to increase S/N in the
   extracted spectra.

-  :ref:`Flux calibration <fluxing>` of extracted 1D spectra.

-  :ref:`Coadding / collating flux-calibrated 1D spectra <coadd1d>` of the same object
   into separate files.

-  :ref:`Telluric correction <telluric_correction>` for NIR spectra (only relevant for the very red
   end of DeVeny's range).

.. _deveny_coadd2d:

Coadding 2D Spectral Images
^^^^^^^^^^^^^^^^^^^^^^^^^^^

PypeIt has the ability to :ref:`coadd 2D spectral images <coadd2d>` of the same
object to increase signal-to-noise prior to object finding and extraction.
While it is possible to simply combine (without weighting) individual exposures
by using the ``comb_id`` column in the PypeIt Reduction File, 2D coadding
accounts for spectral and/or spatial shifts in the spectrum on the CCD.  The
former is important given the spectral flexure seen in DeVeny's camera, and the
latter can help with jitter in the position of the object along the slit due to
manual guiding or imperfect replacement of the object on the slit between
observations.  Coadding aligns the frames spectrally and spatially before
running the object finding and extraction routines.

Coadding is done after the main PypeIt run (as it requires the wavelength
calibration and slit definitions produced during the reduction) and is executed
with the :ref:`pypeit-coadd-2dspec` script. Because the input file format for
this script can be a bit cumbersome, there is a :ref:`setup script
<pypeit_setup_coadd2d>` available that ingests the ``.pypeit`` file or reads
FITS headers in a directory as a starting point. 

In a case of astronomical meta observation, LDT/DeVeny took some spectra of the
JWST spacecraft during its operational mission at the Earth-Sun L2 point.  Four
300-second spectra were taken, and manual guiding was undertaken to keep the
object on the slit as a consequence of the quality of the ephemeris.  

To illustrate the difference between a straight combination and coadding, these
spectra were subject to both procedures and the results shown below.

Single-Frame Spectra of JWST
++++++++++++++++++++++++++++

The extracted 1D spectra from the individual frames are shown below.

.. figure:: ../figures/deveny_1d_jwst_spec.png
   :alt: Spectra of JWST spacecraft
   :width: 100%
   :class: with-shadow

   Four 300-second 1D spectra of the JWST spacecraft.  Note the variation in
   object position along the slit moves by ~5 pixels from the first (top-left)
   to last (bottom-right) frames.  These spectra were displayed in Ginga using
   the ``--ginga`` option to ``pypeit_show_1dspec``.

The contents of the associated ``.txt`` files are (listed together for
clarity):

::

   | slit |                    name | obj_id | spat_pixpos | spat_fracpos | box_width | opt_fwhm |  s2n | wv_rms |
   |  241 | SPAT0250-SLIT0241-DET01 |    250 |       250.1 |        0.519 |      3.80 |    2.144 | 6.33 |  0.088 |
   |  241 | SPAT0249-SLIT0241-DET01 |    249 |       248.7 |        0.516 |      3.80 |    1.887 | 7.92 |  0.088 |
   |  241 | SPAT0247-SLIT0241-DET01 |    247 |       246.6 |        0.511 |      3.80 |    1.960 | 4.65 |  0.088 |
   |  241 | SPAT0245-SLIT0241-DET01 |    245 |       245.2 |        0.508 |      3.80 |    2.055 | 6.77 |  0.088 |

We see that the individual spectra range from an integrated S/N of 4.7 to 6.8,
and all have OPT FWHM between 1.9" and 2.1"


Combining the Frames Directly in the Main PypeIt Run
++++++++++++++++++++++++++++++++++++++++++++++++++++

The simplest way to combine these frames to increase signal-to-noise is to
perform a straight combination during the main PypeIt run.  To do this, you
would need to call ``pypeit_setup`` with the ``-b`` flag to include the
"background pair" columns at the far right of the PypeIt Reduction File.
By default, all calibration frames are given the value ``-1``, and science
frames are numbered sequentially.  To combine frames directly, assign the same
``comb_id`` value to all frames.  In the example here, (portions of) the Data
block our PypeIt Reduction File would look like:

.. code-block::

             filename |       frametype |       ra |    dec | target | airmass | exptime | slitwid | lampstat01 | calib | comb_id | bkg_id
   20250909.0066.fits |         science | 330.8203 | 0.0115 |   JWST |    1.24 |   300.0 |     1.0 |        off |     0 |       7 |     -1
   20250909.0067.fits |         science | 330.8220 | 0.0171 |   JWST |    1.23 |   300.0 |     1.0 |        off |     0 |       7 |     -1
   20250909.0068.fits |         science | 330.8239 | 0.0214 |   JWST |    1.23 |   300.0 |     1.0 |        off |     0 |       7 |     -1
   20250909.0069.fits |         science | 330.8255 | 0.0254 |   JWST |    1.22 |   300.0 |     1.0 |        off |     0 |       7 |     -1

The resulting spectrum is named for the first file in the group, and includes
in the header information about the frames that went into the combination.

The resulting spec1d ``.txt`` file for this combination is:

::

   | slit |                    name | obj_id | spat_pixpos | spat_fracpos | box_width | opt_fwhm |  s2n | wv_rms |
   |  241 | SPAT0248-SLIT0241-DET01 |    248 |       247.6 |        0.514 |      3.80 |    2.330 | 8.90 |  0.088 |

The total S/N has improved to 8.9, but note that the OPT FWHM of the extracted
has increased to 2.3" as a result of JWST moving along the slit.


Coadding the 2D Spectra
+++++++++++++++++++++++

To prepare for coadding the processed 2D spectra of the individual frames, we
can use the :ref:`pypeit_setup_coadd2d` script.  Since we are interested in
coadding just the frames of JWST spacecraft spectra, we can use the ``--obj``
option, in addition to specifying the location of the science spectra.  If
we run the script in the ``ldt_deveny_A`` directory (which contains the
PypeIt Reduction File), this looks like:

::

   $ pypeit_setup_coadd2d -d Science/ --obj JWST

The processing input file ``ldt_deveny_JWST.coadd2d`` is created in the working
directory, and has the contents:

.. code-block:: ini

   # Auto-generated Coadd2D input file using PypeIt version: 1.18.0
   # UTC 2025-09-26T21:27:12.718+00:00

   # User-defined execution parameters
   [rdx]
      spectrograph = ldt_deveny
      redux_path = .
      scidir = Science
      qadir = QA
   [calibrations]
      calib_dir = Calibrations
      [[wavelengths]]
         refframe = observed
   [coadd2d]
      offsets = auto
      weights = auto
      spat_toler = 5
      spec_samp_fact = 1.0
      spat_samp_fact = 1.0
   [flexure]
      spec_method = skip
   [reduce]
      [[findobj]]
         skip_skysub = True

   # Data block 
   spec2d read
   path ./Science
                                                    filename
   spec2d_20250909.0066-JWST_DeVeny_20250909T052930.310.fits
   spec2d_20250909.0067-JWST_DeVeny_20250909T053627.110.fits
   spec2d_20250909.0068-JWST_DeVeny_20250909T054151.040.fits
   spec2d_20250909.0069-JWST_DeVeny_20250909T054703.360.fits
   spec2d end

The setup script adds many of the parameter knobs you might need to turn to
make your 2D coadd successful.  See the :ref:`Coadd2D parameters <coadd2dpar>`
for a detailed listing of what these and others can do for your data.

Once you have created the input file, run the coadd:

::

   $ pypeit_coadd_2dspec ldt_deveny_JWST.coadd2d

The resulting spec1d ``.txt`` file for the coadd is:

::

   | slit |                    name | obj_id | spat_pixpos | spat_fracpos | box_width | opt_fwhm |   s2n |
   |  248 | SPAT0254-SLIT0248-DET01 |    254 |       253.6 |        0.512 |      3.80 |    1.988 | 11.82 |

The total S/N is now 11.8, and the OPT FWHM is 2.0" -- right in line with the
individual frame extractions.

A visual comparison of the straight-combined spectrum (left) and 2D coadded
spectrum (right) is shown below.

.. figure:: ../figures/deveny_1d_jwst_comb_coadd.png
   :alt: Comparison of 2D combination methods
   :width: 100%
   :class: with-shadow

   The straight combined (using ``comb_id`` -- left) and coadded (using
   ``pypeit_coadd_2dspec`` -- right) versions of the 4 spectra of the JWST
   spacecraft.  The ``comb_id`` version has an integrated S/N of 8.9 and an
   optimally extracted profile FWHM of 2.3", whereas the ``pypeit_coadd_2spec``
   version has an integrated S/N of 11.8 and an optimally extracted profile
   FWHM of 2.0".


Which to Use for DeVeny Data?
+++++++++++++++++++++++++++++

It depends.  If you have autoguiding set up for a series of spectra and expect
the object to remain at the same location on the slit for all exposures, then
the straight combination is fine.  If you have a wandering object (like this
example), or have variable S/N on the individual frames (*e.g.*, from clouds),
then coadding might be the better path forward.


.. _deveny_flux:

Flux Calibrating 1D Spectra
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The main PypeIt run returns extracted 1D spectra, measured in instrumental
units (namely, *electrons*). For some science programs, this is sufficient,
and further processing is unnecessary prior to analysis (skip ahead to
:ref:`deveny_specutils`). Other programs either benefit from or require
correcting for the relative spectral sensitivity of the instrument and
converting the instrumental intensity into physical flux units before the
spectra can be analyzed. PypeIt provides routines for :ref:`creating a
sensitivity function <pypeit_sensfunc>` for your data set from observations
of `spectrophotometric standard stars
<https://www.eso.org/sci/observing/tools/standards/spectra/stanlis.html>`__,
and applying that to the remainder of the science data.

If you plan to flux calibrate your spectra, it is imperative to include one
or more `spectrophotometric standard stars 
<https://www.eso.org/sci/observing/tools/standards/spectra/stanlis.html>`__
in your observing program. Exactly which stars and when to observe them depend
on the specific requirements of your science program.  Please see the
:ref:`fluxing` documentation for a description of how to perform this step.

.. important::

   When performing the flux calibration, ``spec1d`` files are modified in
   place, adding the additional components of the data model (*e.g.*,
   ``OPT_FLAM``, ``BOX_FLAM``, etc.) as FITS extensions.  Running
   ``run_pypeit`` with the ``-o`` overwrite flag will cause flux calibration
   information for a given object to be lost requiring the re-running of
   ``pypeit_flux_calibrate``.

-  The first step is to build a sensitivity function using ``pypeit_sensfunc``
   from your observed spectrophotometric standard star that translates the
   count rate (in :math:`{\rm e}^- / {\rm s}`) on the detector as a function of
   wavelength into a flux density (in units of
   :math:`10^{-17} {\rm erg} / {\rm s} / {\rm cm}^2 / \mathring{A}`). Due to factors such
   as grating blaze and the transmission function of the optics in the
   telescope and spectrograph, this sensitivity function will not be uniform
   and requires careful fitting.

   The script will produce an output sensitivity function file in the working
   directory -- you may name the output file anything you like, but it is
   generally helpful to use  something identifiable to the setup and/or date of
   the observation. The figure below shows the throughput plot for the
   spectrophotometric standard star G191-B2B taken on 2022-11-02UT (the same
   night as the other DV6 data shown above).

   .. figure:: ../figures/deveny_sensfunc_throughput.png
      :alt: DV6 sensitivity function throughput
      :width: 50%
      :class: with-shadow

      Example of PypeIt sensitivity function throughput. This observation was
      taken of G191-B2B with DV6 with a 1.2" slit on the night of 2022-11-02UT.

-  Once you are satisfied with with the sensitivity function, the next step is
   to use ``pypeit_flux_setup`` to create a ``.flux`` input file that drives
   the actual flux calibration process. As with the Pypeit Reduction File, you
   will need to edit the ``ldt_deveny.flux`` file to ensure the flux
   calibration proceeds as expected.  See :ref:`pypeit_flux_setup` for a
   description of necessary edits.  The most common for DeVeny users will be to
   specify the sensitivity function file(s) to be used and specify the ``UVIS``
   algorithm be used (for observations blueward of ~9000\ :math:`\mathring{A}`):

   .. code-block:: ini

      [fluxcalib]
         extinct_correct = True  # Set to True if your SENSFUNC derived with the UVIS algorithm

-  After all of the setup work above, the actual flux calibration execution
   is quite straightforward with a call to ``pypeit_flux_calibrate``. All of
   the file information and parameter adjustments are in the
   ``ldt_deveny.flux`` file, and this script requires no additional
   information. Examples of flux-calibrated spectra for the two objects
   described above (automatically identified and manually identified) are shown
   below.

   .. grid:: 2

      .. grid-item::
         :columns: 6

         .. image:: ../figures/deveny_spec1d_fluxed_A.png
            :alt: Automatically identified object, flux calibrated
            :class: with-shadow

      .. grid-item::
         :columns: 6

         .. image:: ../figures/deveny_spec1d_fluxed_B.png
            :alt: Manually identified object, flux calibrated
            :class: with-shadow

      .. grid-item::
         :columns: 12

         Example of flux-calibrated spectra for the objects shown above. As with
         the uncalibrated spectra, the red dashed line indicates the 1-:math:`\sigma`
         uncertainty in the data.


Coadding / Collating Flux-Calibrated 1D Spectra
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

PypeIt has the ability to :ref:`coadd flux-calibrated 1D spectra <coadd1d>` of
the same object. This may be because you have exposures of the same object from
different nights or the object was placed in different locations along the slit
in different frames, either of which precludes coadding the processed 2D
spectral images. In this case, you may use the ``pypeit_coadd_1dspec`` script
for coadding these individual flux-calibrated extracted spectra. This step is
less common for DeVeny users; read the :ref:`coadd1d` documentation if you wish
to perform this action.


Performing a Telluric Correction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For observations done at the extreme red end of the DeVeny's range
(:math:`\gtrsim 9000 \mathring{A}`), you will want to perform a telluric correction to
minimize the effects of atmospheric emission on your data. If you need to
perform this step, please read through the :ref:`telluric_correction`
documentation, and let LDT staff know the use case and how well it worked.


.. _deveny_specutils:

Loading PypeIt 1D Spectra into ``specutils`` for Analysis
---------------------------------------------------------

PypeIt is a package for reducing spectroscopic data from raw frames collected
at the telescope to 1D spectra, ready for analysis. To do the actual analysis
in service of your particular science program, you will need to employ other
tools. One possibility is the AstroPy-coordinated package ``specutils``\ [4]_.

As of ``v1.12.2``, PypeIt includes a loader for importing pipeline outputs into
``specutils``, and can import either the ``spec1d`` (all objects in a frame) or
the ``OneSpec`` (output of ``pypeit_collate_1d``) PypeIt 1D spectral format. 
These loaders automatically recognize PypeIt data from the FITS headers and
properly parse the data into class instance(s).

See the :ref:`spec1D-specutils` for details of implementation.  What you do
with the loaded object(s) will be defined by the requirements of your science
program and is beyond the scope of this documentation.




.. [2]
   `<https://tools.ietf.org/html/rfc1149>`__

.. [3]
   `Horne, K. 1986, PASP, 98,
   609 <https://ui.adsabs.harvard.edu/abs/1986PASP...98..609H/abstract>`__

.. [4]
   `<https://specutils.readthedocs.io/en/stable/index.html>`__ In
   contrast to PypeIt itself, the use of ``specutils`` *does* require knowledge
   of Python for use. This is but one possible analysis tool, and the reader is
   encouraged to seek out the best tool for their particular work.
