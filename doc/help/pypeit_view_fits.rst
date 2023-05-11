.. code-block:: console

    $ pypeit_view_fits -h
    usage: pypeit_view_fits [-h] [--list] [--proc] [--bkg_file BKG_FILE]
                            [--exten EXTEN] [--det [DET ...]] [--chname CHNAME]
                            [--embed]
                            spectrograph file
    
    View FITS files with ginga
    
    positional arguments:
      spectrograph         A valid spectrograph identifier: bok_bc,
                           gemini_flamingos1, gemini_flamingos2,
                           gemini_gmos_north_e2v, gemini_gmos_north_ham,
                           gemini_gmos_north_ham_ns, gemini_gmos_south_ham,
                           gemini_gnirs, gtc_maat, gtc_osiris, gtc_osiris_plus,
                           jwst_nircam, jwst_nirspec, keck_deimos, keck_hires,
                           keck_kcwi, keck_lris_blue, keck_lris_blue_orig,
                           keck_lris_red, keck_lris_red_mark4, keck_lris_red_orig,
                           keck_mosfire, keck_nires, keck_nirspec_low, lbt_luci1,
                           lbt_luci2, lbt_mods1b, lbt_mods1r, lbt_mods2b,
                           lbt_mods2r, ldt_deveny, magellan_fire,
                           magellan_fire_long, magellan_mage, mdm_osmos_mdm4k,
                           mmt_binospec, mmt_bluechannel, mmt_mmirs, not_alfosc,
                           not_alfosc_vert, ntt_efosc2, p200_dbsp_blue,
                           p200_dbsp_red, p200_tspec, shane_kast_blue,
                           shane_kast_red, shane_kast_red_ret, soar_goodman_blue,
                           soar_goodman_red, tng_dolores, vlt_fors2, vlt_sinfoni,
                           vlt_xshooter_nir, vlt_xshooter_uvb, vlt_xshooter_vis,
                           wht_isis_blue, wht_isis_red
      file                 FITS file
    
    optional arguments:
      -h, --help           show this help message and exit
      --list               List the extensions only? (default: False)
      --proc               Process the image (i.e. orient, overscan subtract,
                           multiply by gain) using pypeit.images.buildimage. Note
                           det=mosaic will not work with this option (default:
                           False)
      --bkg_file BKG_FILE  FITS file to be subtracted from the image in file.--proc
                           must be set in order for this option to work. (default:
                           None)
      --exten EXTEN        Show a FITS extension in the raw file. Note --proc and
                           --mosaic will not work with this option. (default: None)
      --det [DET ...]      Detector(s) to show. If more than one, the list of
                           detectors, i.e. --det 4 8 to show detectors 4 and 8. This
                           combination must be one of the allowed mosaics hard-coded
                           for the selected spectrograph. Using "mosaic" for
                           gemini_gmos, keck_deimos, or keck_lris will show the
                           mosaic of all detectors. (default: 1)
      --chname CHNAME      Name of Ginga tab (default: Image)
      --embed              Upon completion embed in ipython shell (default: False)
    