.. code-block:: console

    $ run_pypeit -h
    usage: run_pypeit [-h] [-v VERBOSITY] [-r REDUX_PATH] [-m] [-s] [-o] [-c]
                      pypeit_file
    
    ##  [1;37;42mPypeIt : The Python Spectroscopic Data Reduction Pipeline v1.15.1.dev151+g015ee23ce.d20240423[0m
    ##  
    ##  Available spectrographs include:
    ##   bok_bc, gemini_flamingos1, gemini_flamingos2, gemini_gmos_north_e2v,
    ##   gemini_gmos_north_ham, gemini_gmos_north_ham_ns,
    ##   gemini_gmos_south_ham, gemini_gnirs_echelle, gemini_gnirs_ifu,
    ##   gtc_maat, gtc_osiris, gtc_osiris_plus, jwst_nircam, jwst_nirspec,
    ##   keck_deimos, keck_esi, keck_hires, keck_kcrm, keck_kcwi,
    ##   keck_lris_blue, keck_lris_blue_orig, keck_lris_red,
    ##   keck_lris_red_mark4, keck_lris_red_orig, keck_mosfire, keck_nires,
    ##   keck_nirspec_high, keck_nirspec_high_old, keck_nirspec_low, lbt_luci1,
    ##   lbt_luci2, lbt_mods1b, lbt_mods1r, lbt_mods2b, lbt_mods2r, ldt_deveny,
    ##   magellan_fire, magellan_fire_long, magellan_mage, mdm_modspec,
    ##   mdm_osmos_mdm4k, mdm_osmos_r4k, mmt_binospec, mmt_bluechannel,
    ##   mmt_mmirs, not_alfosc, not_alfosc_vert, ntt_efosc2, p200_dbsp_blue,
    ##   p200_dbsp_red, p200_tspec, shane_kast_blue, shane_kast_red,
    ##   shane_kast_red_ret, soar_goodman_blue, soar_goodman_red, tng_dolores,
    ##   vlt_fors2, vlt_sinfoni, vlt_xshooter_nir, vlt_xshooter_uvb,
    ##   vlt_xshooter_vis, wht_isis_blue, wht_isis_red
    
    positional arguments:
      pypeit_file           PypeIt reduction file (must have .pypeit extension)
    
    options:
      -h, --help            show this help message and exit
      -v VERBOSITY, --verbosity VERBOSITY
                            Verbosity level between 0 [none] and 2 [all]
      -r REDUX_PATH, --redux_path REDUX_PATH
                            Path to directory for the reduction. Only advised for
                            testing
      -m, --do_not_reuse_calibs
                            Do not load previously generated calibrations, even ones
                            made during the run.
      -s, --show            Show reduction steps via plots (which will block further
                            execution until clicked on) and outputs to ginga.
                            Requires remote control ginga session via "ginga
                            --modules=RC,SlitWavelength &"
      -o, --overwrite       Overwrite any existing files/directories
      -c, --calib_only      Only run on calibrations
    