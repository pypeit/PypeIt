.. code-block:: console

    $ pypeit_setup -h
    [1;34musage: [0m[1;35mpypeit_setup[0m [[32m-h[0m] [[32m-v [33mVERBOSITY[0m] [[36m--log_file [33mLOG_FILE[0m]
                        [[36m--log_level [33mLOG_LEVEL[0m] [[32m-s [33mSPECTROGRAPH[0m]
                        [[32m-r [33mROOT [ROOT ...][0m] [[32m-e [33mEXTENSION[0m] [[32m-d [33mOUTPUT_PATH[0m] [[32m-o[0m]
                        [[32m-c [33mCFG_SPLIT[0m] [[32m-b[0m] [[32m-f[0m] [[32m-m[0m] [[32m-k[0m] [[32m-p [33mPARAM_BLOCK_FILE[0m]
                        [[32m-G[0m]
    
    Parse data files to construct a pypeit file in preparation for reduction using
    'run_pypeit'
    
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
      [1;32m-s[0m, [1;36m--spectrograph[0m [1;33mSPECTROGRAPH[0m
                            A valid spectrograph identifier: aat_uhrf, apf_levy,
                            bok_bc, gemini_flamingos1, gemini_flamingos2,
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
      [1;32m-r[0m, [1;36m--root[0m [1;33mROOT [ROOT ...][0m
                            Root to search for data files. You can provide the top-
                            level directory (e.g., /data/Kast) or the search string
                            up through the wildcard (.e.g, /data/Kast/b). Use the
                            --extension option to set the types of files to search
                            for. (default: current working directory)
      [1;32m-e[0m, [1;36m--extension[0m [1;33mEXTENSION[0m
                            File extension to use. Must include the period (e.g.,
                            ".fits") and it must be one of the allowed extensions
                            for this spectrograph. If None, root directory will be
                            searched for all files with any of the allowed
                            extensions. (default: None)
      [1;32m-d[0m, [1;36m--output_path[0m [1;33mOUTPUT_PATH[0m
                            Path to top-level output directory. (default: current
                            working directory)
      [1;32m-o[0m, [1;36m--overwrite[0m       Overwrite any existing files/directories (default:
                            False)
      [1;32m-c[0m, [1;36m--cfg_split[0m [1;33mCFG_SPLIT[0m
                            Generate the PypeIt files and folders by input
                            configuration. To write all unique configurations
                            identifed, use 'all', otherwise provide the list of
                            configuration letters; e.g., 'A,B' or 'B,D,E' or 'E'.
                            (default: None)
      [1;32m-b[0m, [1;36m--background[0m      Include the background-pair columns for the user to edit
                            (default: False)
      [1;32m-f[0m, [1;36m--flexure[0m         Include the manual spatial shift (flexure) column for
                            the user to edit (default: False)
      [1;32m-m[0m, [1;36m--manual_extraction[0m
                            Include the manual extraction column for the user to
                            edit (default: False)
      [1;32m-k[0m, [1;36m--keep_bad_frames[0m
                            Keep all frames, even if they are identified as having
                            bad/unrecognized configurations that cannot be reduced
                            by pypeit. (This is the opposite of the --bad_frames
                            option in pypeit_obslog; i.e., you have to tell
                            pypeit_setup to keep these frames, whereas you have to
                            tell pypeit_obslog to remove them. (default: False)
      [1;32m-p[0m, [1;36m--param_block_file[0m [1;33mPARAM_BLOCK_FILE[0m
                            File containing the additional PypeIt user parameters to
                            be added to the parameter block of the generated
                            reduction file (default: None)
      [1;32m-G[0m, [1;36m--gui[0m             Run setup in a GUI (default: False)
    