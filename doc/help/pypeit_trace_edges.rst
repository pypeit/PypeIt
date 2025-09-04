.. code-block:: console

    $ pypeit_trace_edges -h
    usage: pypeit_trace_edges [-h] (-f PYPEIT_FILE | -t TRACE_FILE) [-g GROUP]
                              [-d [DETECTOR ...]] [-s SPECTROGRAPH] [-b BINNING]
                              [-p REDUX_PATH] [-c CALIB_DIR] [-o] [--debug DEBUG]
                              [--show] [-v VERBOSITY]
    
    Trace slit edges
    
    options:
      -h, --help            show this help message and exit
      -f, --pypeit_file PYPEIT_FILE
                            PypeIt reduction file (default: None)
      -t, --trace_file TRACE_FILE
                            Image to trace (default: None)
      -g, --group GROUP     If providing a pypeit file, use the trace images for
                            this calibration group. If None, use the first
                            calibration group. (default: None)
      -d, --detector [DETECTOR ...]
                            Detector(s) to process. If more than one, the list of
                            detectors must be one of the allowed mosaics hard-coded
                            for the selected spectrograph. Using "mosaic" for
                            gemini_gmos, keck_deimos, or keck_lris will use the
                            default mosaic. (default: None)
      -s, --spectrograph SPECTROGRAPH
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
                            vlt_sinfoni, vlt_xshooter_nir, vlt_xshooter_uvb,
                            vlt_xshooter_vis, wht_isis_blue, wht_isis_red (default:
                            None)
      -b, --binning BINNING
                            Image binning in spectral and spatial directions. Only
                            used if providing files directly; default is 1,1.
                            (default: None)
      -p, --redux_path REDUX_PATH
                            Path to top-level output directory. (default: current
                            working directory)
      -c, --calib_dir CALIB_DIR
                            Name for directory in output path for calibration
                            file(s) relative to the top-level directory. (default:
                            Calibrations)
      -o, --overwrite       Overwrite any existing files/directories (default:
                            False)
      --debug DEBUG         Debug level. (1) Show the result of each stage of the
                            tracing algorithm (previously the --show option). (2)
                            Also show summary plots related to the PCA decomposition
                            and the slit and order matching. (3) Also show the
                            individual polynomial fits to the detected edges.
                            (default: 0)
      --show                DEPRECATED! If set, the code will assume you mean to set
                            --debug 1. (default: False)
      -v, --verbosity VERBOSITY
                            Verbosity level between 0 [none] and 2 [all]. Default:
                            1. Level 2 writes a log with filename
                            trace_edges_YYYYMMDD-HHMM.log (default: 1)
    