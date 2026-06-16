.. code-block:: console

    $ run_pypeit -h
    [1;34musage: [0m[1;35mrun_pypeit[0m [[32m-h[0m] [[32m-v [33mVERBOSITY[0m] [[36m--log_file [33mLOG_FILE[0m]
                      [[36m--log_level [33mLOG_LEVEL[0m] [[32m-r [33mREDUX_PATH[0m] [[32m-m[0m] [[32m-s[0m] [[32m-o[0m] [[32m-c[0m]
                      [32mpypeit_file[0m
    
    PypeIt: The Python Spectroscopic Data Reduction Pipeline
    Version 2.0.2.dev420+gc9fb1e8f6.d20260611
    
    Available spectrographs include:
        aat_uhrf, apf_levy, bok_bc, gemini_flamingos1, gemini_flamingos2,
        gemini_gmos_north_e2v, gemini_gmos_north_ham,
        gemini_gmos_north_ham_ns, gemini_gmos_south_ham, gemini_gnirs_echelle,
        gemini_gnirs_ifu, gtc_maat, gtc_osiris, gtc_osiris_plus, jwst_nircam,
        jwst_nirspec, keck_deimos, keck_esi, keck_hires, keck_kcrm, keck_kcwi,
        keck_lris_blue, keck_lris_blue_orig, keck_lris_red,
        keck_lris_red_mark4, keck_lris_red_orig, keck_mosfire, keck_nires,
        keck_nirspec_high, keck_nirspec_high_old, keck_nirspec_low, lbt_luci1,
        lbt_luci2, lbt_mods1b, lbt_mods1b_proc, lbt_mods1r, lbt_mods1r_proc,
        lbt_mods2b, lbt_mods2b_proc, lbt_mods2r, lbt_mods2r_proc, ldt_deveny,
        magellan_fire, magellan_fire_long, magellan_mage, mdm_modspec,
        mdm_osmos_mdm4k, mdm_osmos_r4k, mmt_binospec, mmt_bluechannel,
        mmt_mmirs, not_alfosc, not_alfosc_vert, ntt_efosc2, p200_dbsp_blue,
        p200_dbsp_red, p200_ngps_i, p200_ngps_r, p200_tspec, shane_kast_blue,
        shane_kast_red, shane_kast_red_ret, soar_goodman_blue,
        soar_goodman_red, subaru_focas, tng_dolores, vlt_fors2, vlt_sinfoni,
        vlt_uves_blue, vlt_uves_red, vlt_xshooter_nir, vlt_xshooter_uvb,
        vlt_xshooter_vis, wht_isis_blue, wht_isis_red
    
    [1;34mpositional arguments:[0m
      [1;32mpypeit_file[0m           PypeIt reduction file (must have .pypeit extension)
    
    [1;34moptions:[0m
      [1;32m-h[0m, [1;36m--help[0m            show this help message and exit
      [1;32m-v[0m, [1;36m--verbosity[0m [1;33mVERBOSITY[0m
                            Verbosity level, which must be 0, 1, or 2. Level 0
                            includes warning and error messages, level 1 adds
                            informational messages, and level 2 adds debugging
                            messages and the calling sequence.
      [1;36m--log_file[0m [1;33mLOG_FILE[0m   Name for the log file. If set to "default", a default
                            name is used. If None, a log file is not produced.
      [1;36m--log_level[0m [1;33mLOG_LEVEL[0m
                            Verbosity level for the log file. If a log file is
                            produce and this is None, the file log will match the
                            console stream log.
      [1;32m-r[0m, [1;36m--redux_path[0m [1;33mREDUX_PATH[0m
                            Path to directory for the reduction. Only advised for
                            testing
      [1;32m-m[0m, [1;36m--do_not_reuse_calibs[0m
                            Do not load previously generated calibrations, even ones
                            made during the run.
      [1;32m-s[0m, [1;36m--show[0m            Show reduction steps via plots (which will block further
                            execution until clicked on) and outputs to ginga.
                            Requires remote control ginga session via "ginga
                            --modules=RC,SlitWavelength &"
      [1;32m-o[0m, [1;36m--overwrite[0m       Overwrite any existing files/directories
      [1;32m-c[0m, [1;36m--calib_only[0m      Only run on calibrations
    