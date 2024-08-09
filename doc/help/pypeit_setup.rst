.. code-block:: console

    $ pypeit_setup -h
    usage: pypeit_setup [-h] [-s SPECTROGRAPH] [-r ROOT [ROOT ...]] [-e EXTENSION]
                        [-d OUTPUT_PATH] [-o] [-c CFG_SPLIT] [-b] [-f] [-m]
                        [-v VERBOSITY] [-k] [-G]
    
    Parse data files to construct a pypeit file in preparation for reduction using
    'run_pypeit'
    
    options:
      -h, --help            show this help message and exit
      -s SPECTROGRAPH, --spectrograph SPECTROGRAPH
                            A valid spectrograph identifier: bok_bc,
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
                            lbt_mods1r, lbt_mods2b, lbt_mods2r, ldt_deveny,
                            magellan_fire, magellan_fire_long, magellan_mage,
                            mdm_modspec, mdm_osmos_mdm4k, mdm_osmos_r4k,
                            mmt_binospec, mmt_bluechannel, mmt_mmirs, not_alfosc,
                            not_alfosc_vert, ntt_efosc2, p200_dbsp_blue,
                            p200_dbsp_red, p200_tspec, shane_kast_blue,
                            shane_kast_red, shane_kast_red_ret, soar_goodman_blue,
                            soar_goodman_red, tng_dolores, vlt_fors2, vlt_sinfoni,
                            vlt_xshooter_nir, vlt_xshooter_uvb, vlt_xshooter_vis,
                            wht_isis_blue, wht_isis_red (default: None)
      -r ROOT [ROOT ...], --root ROOT [ROOT ...]
                            Root to search for data files. You can provide the top-
                            level directory (e.g., /data/Kast) or the search string
                            up through the wildcard (.e.g, /data/Kast/b). Use the
                            --extension option to set the types of files to search
                            for. (default: current working directory)
      -e EXTENSION, --extension EXTENSION
                            File extension; compression indicators (e.g. .gz) not
                            required. (default: .fits)
      -d OUTPUT_PATH, --output_path OUTPUT_PATH
                            Path to top-level output directory. (default: current
                            working directory)
      -o, --overwrite       Overwrite any existing files/directories (default:
                            False)
      -c CFG_SPLIT, --cfg_split CFG_SPLIT
                            Generate the PypeIt files and folders by input
                            configuration. To write all unique configurations
                            identifed, use 'all', otherwise provide the list of
                            configuration letters; e.g., 'A,B' or 'B,D,E' or 'E'.
                            (default: None)
      -b, --background      Include the background-pair columns for the user to edit
                            (default: False)
      -f, --flexure         Include the manual spatial shift (flexure) column for
                            the user to edit (default: False)
      -m, --manual_extraction
                            Include the manual extraction column for the user to
                            edit (default: False)
      -v VERBOSITY, --verbosity VERBOSITY
                            Verbosity level between 0 [none] and 2 [all]. Default:
                            1. Level 2 writes a log with filename setup_YYYYMMDD-
                            HHMM.log (default: 1)
      -k, --keep_bad_frames
                            Keep all frames, even if they are identified as having
                            bad/unrecognized configurations that cannot be reduced
                            by pypeit. (This is the opposite of the --bad_frames
                            option in pypeit_obslog; i.e., you have to tell
                            pypeit_setup to keep these frames, whereas you have to
                            tell pypeit_obslog to remove them. (default: False)
      -G, --gui             Run setup in a GUI (default: False)
    