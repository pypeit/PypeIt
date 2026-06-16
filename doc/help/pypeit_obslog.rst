.. code-block:: console

    $ pypeit_obslog -h
    [1;34musage: [0m[1;35mpypeit_obslog[0m [[32m-h[0m] [[32m-v [33mVERBOSITY[0m] [[36m--log_file [33mLOG_FILE[0m]
                         [[36m--log_level [33mLOG_LEVEL[0m] [[32m-r [33mROOT[0m] [[32m-k[0m] [[32m-c [33mCOLUMNS[0m] [[32m-b[0m]
                         [[32m-t [33mBAD_TYPES[0m] [[32m-g[0m] [[32m-i[0m] [[32m-s [33mSORT[0m] [[32m-e [33mEXTENSION[0m]
                         [[32m-d [33mOUTPUT_PATH[0m] [[32m-o[0m] [[32m-f [33mFILE[0m]
                         [32mspec[0m
    
    Construct an observing log for a set of files from the provided spectrograph
    using PypeItMetaData.
    
    [1;34mpositional arguments:[0m
      [1;32mspec[0m                  A valid spectrograph identifier: aat_uhrf, apf_levy,
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
                            wht_isis_blue, wht_isis_red
    
    [1;34moptions:[0m
      [1;32m-h[0m, [1;36m--help[0m            show this help message and exit
      [1;32m-v[0m, [1;36m--verbosity[0m [1;33mVERBOSITY[0m
                            Verbosity level, which must be 0, 1, or 2. Level 0
                            includes warning and error messages, level 1 adds
                            informational messages, and level 2 adds debugging
                            messages and the calling sequence. (default: 2)
      [1;36m--log_file[0m [1;33mLOG_FILE[0m   Name for the log file. If set to "default", a default
                            name is used. If None, a log file is not produced.
                            (default: None)
      [1;36m--log_level[0m [1;33mLOG_LEVEL[0m
                            Verbosity level for the log file. If a log file is
                            produce and this is None, the file log will match the
                            console stream log. (default: None)
      [1;32m-r[0m, [1;36m--root[0m [1;33mROOT[0m       Root to search for data files. You can provide the top-
                            level directory (e.g., /data/Kast) or the search string
                            up through the wildcard (.e.g, /data/Kast/b). Use the
                            --extension option to set the types of files to search
                            for. (default: current working directory)
      [1;32m-k[0m, [1;36m--keys[0m            Do not produce the log; simply list the pypeit-specific
                            metadata keys available for this spectrograph and their
                            associated header cards. Metadata keys with header cards
                            that are None have no simple mapping between keyword and
                            header card. (default: False)
      [1;32m-c[0m, [1;36m--columns[0m [1;33mCOLUMNS[0m
                            A comma-separated list of columns to include in the
                            output table. Each column must be a valid pypeit
                            metadata keyword specific to this spectrograph (run
                            pypeit_obslog with the -k argument to see the valid
                            list). Additional valid keywords are directory,
                            filename, frametype, framebit, setup, calib, and
                            calibbit. If 'all', all columns collected for the pypeit
                            metadata table are included. If 'pypeit', the columns
                            are the same as those included in the pypeit file.
                            (default: pypeit)
      [1;32m-b[0m, [1;36m--bad_frames[0m      Clean the output of bad frames that cannot be reduced by
                            pypeit. (default: False)
      [1;32m-t[0m, [1;36m--bad_types[0m [1;33mBAD_TYPES[0m
                            Dictates how frames that could not be given a valid type
                            should be treated. Options are: "keep" to include them
                            in the output, "rm" to remove them from the output,
                            "only" to only include the frames with unknown types in
                            the output (i.e, the frames with determined types are
                            excluded). (default: keep)
      [1;32m-g[0m, [1;36m--groupings[0m       Use this option to only determine the frame type. By
                            default, the script groups frames into expected
                            configuration and calibration groups, and it adds the
                            default combination groups. (default: True)
      [1;32m-i[0m, [1;36m--interact[0m        Once the metadata table is created, start an embedded
                            IPython session that you can use to interact with the
                            table (an Astropy.Table called fitstbl) directly.
                            (default: False)
      [1;32m-s[0m, [1;36m--sort[0m [1;33mSORT[0m       Metadata keyword (pypeit-specific) to use to sort the
                            output table. (default: mjd)
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
      [1;32m-f[0m, [1;36m--file[0m [1;33mFILE[0m       Name for the ascii output file. Any leading directory
                            path is stripped; use -d to set the output directory. If
                            None, the table is just printed to stdout. If set to
                            'default', the file is set to [spectrograph].obslog.
                            Note the file will *not* be written if you also include
                            the -i option to embed and interact with the table (you
                            can write the table using the astropy.table.Table.write
                            method in the embedded IPython session). The table is
                            always written in ascii format using
                            format=ascii.fixed_with for the call to
                            Astropy.table.Table.write . (default: None)
    