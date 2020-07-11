**************
PypeIt scripts
**************

PypeIt is packaged with several scripts that should have
been installed directly into your path (e.g. ~/anaconda/bin).

Pipeline Scripts
++++++++++++++++


pypeit_chk_for_calibs
=====================

This script, which is similar to :ref:`pypeit-setup`, examines a set
of files for an input spectrograph and scans for the standard calibrations.
It raises warnings when these are not found.

Here is the usage::

    usage: pypeit_chk_for_calibs [-h] [-s SPECTROGRAPH] [-e EXTENSION] root

    Script to check for calibrations [v1]

    positional arguments:
      root                  File path+root, e.g. /data/Kast/b

    optional arguments:
      -h, --help            show this help message and exit
      -s SPECTROGRAPH, --spectrograph SPECTROGRAPH
                            A valid spectrograph identifier: keck_deimos,
                            keck_lris_blue, keck_lris_red, keck_nires,
                            keck_nirspec_low, keck_mosfire, keck_hires_red,
                            keck_kcwi, shane_kast_blue, shane_kast_red,
                            shane_kast_red_ret, tng_dolores, wht_isis_blue,
                            wht_isis_red, vlt_xshooter_uvb, vlt_xshooter_vis,
                            vlt_xshooter_nir, vlt_fors2, gemini_gnirs,
                            gemini_flamingos1, gemini_flamingos2,
                            gemini_gmos_south_ham, gemini_gmos_north_e2v,
                            gemini_gmos_north_ham, magellan_fire,
                            magellan_fire_long, magellan_mage, lbt_mods1r,
                            lbt_mods1b, lbt_mods2r, lbt_mods2b, lbt_luci1,
                            lbt_luci2, mmt_binospec, mdm_osmos_mdm4k, not_alfosc
      -e EXTENSION, --extension EXTENSION
                            File extension; compression indicators (e.g. .gz) not
                            required.



And a typical call::

    pypeit_chk_calibs /PypeIt-development-suite/RAW_DATA/not_alfosc/grism4/ALDc2 -s not_alfosc

After a running stream of detailed notes, it prints a table of results
to the screen::

    setups pass     scifiles
    ------ -------- ---------------
         A False ALDc200205.fits
      None True


.. _pypeit-setup:

pypeit_setup
============

This setups files for data reduction.  See :doc:`setup` for details

run_pypeit
==========

This is the main executable for PypeIt.  See :doc:`running` for details.

pypeit_view_fits
================

This is a wrapper to the Ginga image viewer.  It is a bit of a kludge
in that it writes a dummy tmp.fits file to the harddrive and sends
that into Ginga.  The dummy file is deleted afterwards.::

    unix> pyp_view_fits -h
    usage: pyp_view_fits [-h] [--list] [--raw_lris] [--exten EXTEN] file

    positional arguments:
      file           FITS file

    optional arguments:
      -h, --help     show this help message and exit
      --list         List the extensions only? (default: False)
      --raw_lris
      --exten EXTEN  FITS extension (default: None)



Data Processing Scripts
+++++++++++++++++++++++

pypeit_coadd_1dspec
===================

See :doc:`coadd1d` for further details.

Calibration Scripts
+++++++++++++++++++

pypeit_arcid_plot
=================

Generate a PDF plot from a MasterFrame_WaveCalib.json file.
This may be useful to ID lines in other data.::

    unix> pypeit_arcid_plot -h
    usage: pypeit_arcid_plot [-h] wave_soln title outfile

    positional arguments:
      wave_soln   MasterWaveSoln file [JSON]
      title       Title for the plot
      outfile     Output PDF file

    optional arguments:
      -h, --help  show this help message and exit

pypeit_lowrdx_pixflat
=====================

Convert a LowRedux pixel flat into a PypeIt ready file::

    unix> pypeit_lowrdx_pixflat -h
    usage: pypeit_lowrdx_pixflat [-h] lowrdx_file new_file

    positional arguments:
      lowrdx_file  LowRedux Pixel Flat FITS file
      new_file     PypeIt FITS file

    optional arguments:
      -h, --help   show this help message and exit


pypeit_chk_edges
================

Inspect the slit/order edges identified by PypeIt in a RC Ginga
window::

    wolverine> pypeit_chk_edges -h
    usage: pypeit_chk_edges [-h] [--chname CHNAME] [--dumb_ids] root

    Display MasterTrace image in a previously launched RC Ginga viewer

    positional arguments:
      root             PypeIt Master Trace file root [e.g.
                       MasterTrace_A_01_aa.fits]

    optional arguments:
      -h, --help       show this help message and exit
      --chname CHNAME  Channel name for image in Ginga (default: MTrace)
      --dumb_ids       Slit ID just by order? (default: False)

pypeit_chk_flats
================

Inspect the flat field images produced by PypeIt in a RC Ginga
window.  This includes the stacked 'raw' image, the pixel flat,
the illumination flat, and the flat model::

    wolverine> pypeit_chk_flats -h
    usage: pypeit_chk_flats [-h] master_file

    Display MasterFlat images in a previously launched RC Ginga viewer

    positional arguments:
      master_file  PypeIt MasterFlat file [e.g. MasterFlat_A_1_01.fits]

    optional arguments:
      -h, --help   show this help message and exit


pypeit_flux_setup
=================

This sets up files for fluxing, coadding and telluric corrections.
Note the pypeit files generated by this scripts need your changes:

    - Give sensfunc file name in the fluxing pypeit file
    - Give sensfunc file name in the coadding pypeit file
    - The coadding pypeit file includes all objects extracted from
      your main reduction, so you need to pick up the one you are
      interested in and remove all others in the coadding pypeit file
      (between coadd1d read and coadd1d end)

See :doc:`fluxing`, :doc:`coadd1d`, and :doc:`telluric` for details::

    unix> pypeit_flux_setup -h
    usage: pypeit_flux_setup sci_path [-h] [--objmodel]

    positional arguments:
      sci_path           the path for your Science folder

    optional arguments:
      -h, --help     show this help message and exit
      --objmodel     set objmodel for telluric fitting (default: qso)
