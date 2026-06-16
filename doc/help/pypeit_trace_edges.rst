.. code-block:: console

    $ pypeit_trace_edges -h
    [1;34musage: [0m[1;35mpypeit_trace_edges[0m [[32m-h[0m] [[32m-v [33mVERBOSITY[0m] [[36m--log_file [33mLOG_FILE[0m]
                              [[36m--log_level [33mLOG_LEVEL[0m] ([32m-f [33mPYPEIT_FILE[0m |
                              [32m-t [33mTRACE_FILE[0m) [[32m-g [33mGROUP[0m] [[32m-d [33m[DETECTOR ...][0m]
                              [[32m-s [33mSPECTROGRAPH[0m] [[32m-b [33mBINNING[0m] [[32m-p [33mREDUX_PATH[0m]
                              [[32m-c [33mCALIB_DIR[0m] [[32m-o[0m] [[36m--debug [33mDEBUG[0m]
    
    Trace slit edges
    
    [1;34moptions:[0m
      [1;32m-h[0m, [1;36m--help[0m            show this help message and exit
      [1;32m-v[0m, [1;36m--verbosity[0m [1;33mVERBOSITY[0m
                            Verbosity level, which must be 0, 1, or 2. Level 0
                            includes warning and error messages, level 1 adds
                            informational messages, and level 2 adds debugging
                            messages and the calling sequence. (default: 2)
      [1;36m--log_file[0m [1;33mLOG_FILE[0m   Name for the log file. If set to "default", a default
                            name is used. If None, a log file is not produced.
                            (default: default)
      [1;36m--log_level[0m [1;33mLOG_LEVEL[0m
                            Verbosity level for the log file. If a log file is
                            produce and this is None, the file log will match the
                            console stream log. (default: None)
      [1;32m-f[0m, [1;36m--pypeit_file[0m [1;33mPYPEIT_FILE[0m
                            PypeIt reduction file (default: None)
      [1;32m-t[0m, [1;36m--trace_file[0m [1;33mTRACE_FILE[0m
                            Image to trace (default: None)
      [1;32m-g[0m, [1;36m--group[0m [1;33mGROUP[0m     If providing a pypeit file, use the trace images for
                            this calibration group. If None, use the first
                            calibration group. (default: None)
      [1;32m-d[0m, [1;36m--detector[0m [1;33m[DETECTOR ...][0m
                            Detector(s) to process. If more than one, the list of
                            detectors must be one of the allowed mosaics hard-coded
                            for the selected spectrograph. Using "mosaic" for
                            gemini_gmos, keck_deimos, or keck_lris will use the
                            default mosaic. (default: None)
      [1;32m-s[0m, [1;36m--spectrograph[0m [1;33mSPECTROGRAPH[0m
                            A valid spectrograph identifier, which is only used if
                            providing files directly: aat_uhrf, apf_levy, bok_bc,
                            gemini_flamingos1, gemini_flamingos2,
                            gemini_gmos_north_e2v, gemini_gmos_north_ham,
                            gemini_gmos_north_ham_ns, gemini_gmos_south_ham,
                            gemini_gnirs_echelle, gemini_gnirs_ifu, gtc_maat,
                            gtc_osiris, gtc_osiris_plus, jwst_nircam, jwst_nirspec,
                            keck_deimos, keck_esi, keck_hires, keck_kcrm, keck_kcwi,
                            keck_lris_blue, keck_lris_blue_orig, keck_lris_red,
                            keck_lris_red_mark4, keck_lris_red_orig, keck_mosfire,
                            keck_nires, keck_nirspec_high, keck_nirspec_high_old,
                            keck_nirspec_low, lbt_luci1, lbt_luci2, lbt_mods1b,
                            lbt_mods1b_proc, lbt_mods1r, lbt_mods1r_proc,
                            lbt_mods2b, lbt_mods2b_proc, lbt_mods2r,
                            lbt_mods2r_proc, ldt_deveny, magellan_fire,
                            magellan_fire_long, magellan_mage, mdm_modspec,
                            mdm_osmos_mdm4k, mdm_osmos_r4k, mmt_binospec,
                            mmt_bluechannel, mmt_mmirs, not_alfosc, not_alfosc_vert,
                            ntt_efosc2, p200_dbsp_blue, p200_dbsp_red, p200_ngps_i,
                            p200_ngps_r, p200_tspec, shane_kast_blue,
                            shane_kast_red, shane_kast_red_ret, soar_goodman_blue,
                            soar_goodman_red, subaru_focas, tng_dolores, vlt_fors2,
                            vlt_sinfoni, vlt_uves_blue, vlt_uves_red,
                            vlt_xshooter_nir, vlt_xshooter_uvb, vlt_xshooter_vis,
                            wht_isis_blue, wht_isis_red (default: None)
      [1;32m-b[0m, [1;36m--binning[0m [1;33mBINNING[0m
                            Image binning in spectral and spatial directions. Only
                            used if providing files directly; default is 1,1.
                            (default: None)
      [1;32m-p[0m, [1;36m--redux_path[0m [1;33mREDUX_PATH[0m
                            Path to top-level output directory. (default: current
                            working directory)
      [1;32m-c[0m, [1;36m--calib_dir[0m [1;33mCALIB_DIR[0m
                            Name for directory in output path for calibration
                            file(s) relative to the top-level directory. (default:
                            Calibrations)
      [1;32m-o[0m, [1;36m--overwrite[0m       Overwrite any existing files/directories (default:
                            False)
      [1;36m--debug[0m [1;33mDEBUG[0m         Debug level. (1) Show the result of each stage of the
                            tracing algorithm (previously the --show option). (2)
                            Also show summary plots related to the PCA decomposition
                            and the slit and order matching. (3) Also show the
                            individual polynomial fits to the detected edges.
                            (default: 0)
    