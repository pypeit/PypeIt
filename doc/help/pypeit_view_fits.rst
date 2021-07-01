.. code-block:: console

    $ pypeit_view_fits -h
    usage: pypeit_view_fits [-h] [--list] [--proc] [--exten EXTEN] [--det DET]
                            [--chname CHNAME]
                            spectrograph file
    
    View FITS files with ginga
    
    positional arguments:
      spectrograph     A valid spectrograph identifier: bok_bc, gemini_flamingos1,
                       gemini_flamingos2, gemini_gmos_north_e2v,
                       gemini_gmos_north_ham, gemini_gmos_north_ham_ns,
                       gemini_gmos_south_ham, gemini_gnirs, gtc_osiris, keck_deimos,
                       keck_hires_red, keck_kcwi, keck_lris_blue,
                       keck_lris_blue_orig, keck_lris_red, keck_lris_red_orig,
                       keck_mosfire, keck_nires, keck_nirspec_low, lbt_luci1,
                       lbt_luci2, lbt_mods1b, lbt_mods1r, lbt_mods2b, lbt_mods2r,
                       ldt_deveny, magellan_fire, magellan_fire_long, magellan_mage,
                       mdm_osmos_mdm4k, mmt_binospec, mmt_bluechannel, mmt_mmirs,
                       not_alfosc, ntt_efosc2, p200_dbsp_blue, p200_dbsp_red,
                       p200_tspec, shane_kast_blue, shane_kast_red,
                       shane_kast_red_ret, soar_goodman_red, tng_dolores, vlt_fors2,
                       vlt_sinfoni, vlt_xshooter_nir, vlt_xshooter_uvb,
                       vlt_xshooter_vis, wht_isis_blue, wht_isis_red
      file             FITS file
    
    optional arguments:
      -h, --help       show this help message and exit
      --list           List the extensions only? (default: False)
      --proc           Process the image (i.e. orient, overscan subtract, multiply
                       by gain) using pypeit.images.buildimage. Note det=mosaic will
                       not work with this option (default: False)
      --exten EXTEN    Show a FITS extension in the raw file. Note --proc and
                       --mosaic will not work with this option. (default: None)
      --det DET        Detector number. To mosaic keck_deimos or keck_lris images,
                       set equal to mosaic. (default: 1)
      --chname CHNAME  Name of Ginga tab (default: Image)
    