.. code-block:: console

    $ pypeit_trace_edges -h
    usage: pypeit_trace_edges [-h] (-f PYPEIT_FILE | -t TRACE_FILE) [-g GROUP]
                              [-d DETECTOR] [-s SPECTROGRAPH] [-b BINNING]
                              [-p REDUX_PATH] [-m MASTER_DIR] [-o] [--debug]
                              [--show]
    
    optional arguments:
      -h, --help            show this help message and exit
      -f PYPEIT_FILE, --pypeit_file PYPEIT_FILE
                            PypeIt reduction file (default: None)
      -t TRACE_FILE, --trace_file TRACE_FILE
                            Image to trace (default: None)
      -g GROUP, --group GROUP
                            If providing a pypeit file, use the trace images for
                            this calibration group. If None, use the first
                            calibration group. (default: None)
      -d DETECTOR, --detector DETECTOR
                            Only analyze the specified detector; otherwise analyze
                            all or detectors selected by the pypeit file, if
                            provided. (default: None)
      -s SPECTROGRAPH, --spectrograph SPECTROGRAPH
                            A valid spectrograph identifier, which is only used if
                            providing files directly: gemini_flamingos1,
                            gemini_flamingos2, gemini_gmos_north_e2v,
                            gemini_gmos_north_ham, gemini_gmos_north_ham_ns,
                            gemini_gmos_south_ham, gemini_gnirs, keck_deimos,
                            keck_hires_red, keck_kcwi, keck_lris_blue,
                            keck_lris_blue_orig, keck_lris_red,
                            keck_lris_red_orig, keck_mosfire, keck_nires,
                            keck_nirspec_low, lbt_luci1, lbt_luci2, lbt_mods1b,
                            lbt_mods1r, lbt_mods2b, lbt_mods2r, magellan_fire,
                            magellan_fire_long, magellan_mage, mdm_osmos_mdm4k,
                            mmt_binospec, mmt_bluechannel, mmt_mmirs, not_alfosc,
                            p200_dbsp_blue, p200_dbsp_red, p200_tspec,
                            shane_kast_blue, shane_kast_red, shane_kast_red_ret,
                            tng_dolores, vlt_fors2, vlt_xshooter_nir,
                            vlt_xshooter_uvb, vlt_xshooter_vis, wht_isis_blue,
                            wht_isis_red (default: None)
      -b BINNING, --binning BINNING
                            Image binning in spectral and spatial directions. Only
                            used if providing files directly; default is 1,1.
                            (default: None)
      -p REDUX_PATH, --redux_path REDUX_PATH
                            Path to top-level output directory. Default is the
                            current working directory. (default: None)
      -m MASTER_DIR, --master_dir MASTER_DIR
                            Name for directory in output path for Master file(s)
                            relative to the top-level directory. (default:
                            Masters)
      -o, --overwrite       Overwrite any existing files/directories (default:
                            False)
      --debug               Run in debug mode. (default: False)
      --show                Show the stages of trace refinements (only for the new
                            code). (default: False)
    