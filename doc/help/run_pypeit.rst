.. code-block:: console

    $ run_pypeit -h
    usage: run_pypeit [-h] [-v VERBOSITY] [-t] [-r REDUX_PATH] [-m] [-s] [-o]
                      [-d DETECTOR] [-c]
                      pypeit_file
    
    ##  [1;37;42mPypeIt : The Python Spectroscopic Data Reduction Pipeline v1.3.1dev[0m
    ##  
    ##  Available spectrographs include:
    ##   gemini_flamingos1, gemini_flamingos2, gemini_gmos_north_e2v,
    ##   gemini_gmos_north_ham, gemini_gmos_north_ham_ns,
    ##   gemini_gmos_south_ham, gemini_gnirs, keck_deimos, keck_hires_red,
    ##   keck_kcwi, keck_lris_blue, keck_lris_blue_orig, keck_lris_red,
    ##   keck_lris_red_orig, keck_mosfire, keck_nires, keck_nirspec_low,
    ##   lbt_luci1, lbt_luci2, lbt_mods1b, lbt_mods1r, lbt_mods2b, lbt_mods2r,
    ##   magellan_fire, magellan_fire_long, magellan_mage, mdm_osmos_mdm4k,
    ##   mmt_binospec, mmt_bluechannel, mmt_mmirs, not_alfosc, p200_dbsp_blue,
    ##   p200_dbsp_red, p200_tspec, shane_kast_blue, shane_kast_red,
    ##   shane_kast_red_ret, tng_dolores, vlt_fors2, vlt_xshooter_nir,
    ##   vlt_xshooter_uvb, vlt_xshooter_vis, wht_isis_blue, wht_isis_red
    
    positional arguments:
      pypeit_file           PypeIt reduction file (must have .pypeit extension)
    
    optional arguments:
      -h, --help            show this help message and exit
      -v VERBOSITY, --verbosity VERBOSITY
                            Verbosity level between 0 [none] and 2 [all]
      -t, --hdrframetype    Use file headers and the instument-specific keywords
                            to determinethe type of each frame
      -r REDUX_PATH, --redux_path REDUX_PATH
                            Path to directory for the reduction. Only advised for
                            testing
      -m, --do_not_reuse_masters
                            Do not load previously generated MasterFrames, even
                            ones made during the run.
      -s, --show            Show reduction steps via plots (which will block
                            further execution until clicked on) and outputs to
                            ginga. Requires remote control ginga session via
                            "ginga --modules=RC &"
      -o, --overwrite       Overwrite any existing files/directories
      -d DETECTOR, --detector DETECTOR
                            Detector to limit reductions on. If the output files
                            exist and -o is used, the outputs for the input
                            detector will be replaced.
      -c, --calib_only      Only run on calibrations
    