.. code-block:: console

    $ pypeit_obslog -h
    usage: pypeit_obslog [-h] [-r ROOT] [-k] [-c COLUMNS] [-b] [-t BAD_TYPES] [-g]
                         [-i] [-s SORT] [-e EXTENSION] [-d OUTPUT_PATH] [-o]
                         [-f FILE] [-G]
                         spec
    
    Construct an observing log for a set of files from the provided spectrograph
    using PypeItMetaData.
    
    positional arguments:
      spec                  A valid spectrograph identifier: bok_bc,
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
                            wht_isis_blue, wht_isis_red
    
    options:
      -h, --help            show this help message and exit
      -r ROOT, --root ROOT  Root to search for data files. You can provide the top-
                            level directory (e.g., /data/Kast) or the search string
                            up through the wildcard (.e.g, /data/Kast/b). Use the
                            --extension option to set the types of files to search
                            for. (default: current working directory)
      -k, --keys            Do not produce the log; simply list the pypeit-specific
                            metadata keys available for this spectrograph and their
                            associated header cards. Metadata keys with header cards
                            that are None have no simple mapping between keyword and
                            header card. (default: False)
      -c COLUMNS, --columns COLUMNS
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
      -b, --bad_frames      Clean the output of bad frames that cannot be reduced by
                            pypeit. (default: False)
      -t BAD_TYPES, --bad_types BAD_TYPES
                            Dictates how frames that could not be given a valid type
                            should be treated. Options are: "keep" to include them
                            in the output, "rm" to remove them from the output,
                            "only" to only include the frames with unknown types in
                            the output (i.e, the frames with determined types are
                            excluded). (default: keep)
      -g, --groupings       Use this option to only determine the frame type. By
                            default, the script groups frames into expected
                            configuration and calibration groups, and it adds the
                            default combination groups. (default: True)
      -i, --interact        Once the metadata table is created, start an embedded
                            IPython session that you can use to interact with the
                            table (an Astropy.Table called fitstbl) directly.
                            (default: False)
      -s SORT, --sort SORT  Metadata keyword (pypeit-specific) to use to sort the
                            output table. (default: mjd)
      -e EXTENSION, --extension EXTENSION
                            File extension; compression indicators (e.g. .gz) not
                            required. (default: .fits)
      -d OUTPUT_PATH, --output_path OUTPUT_PATH
                            Path to top-level output directory. (default: current
                            working directory)
      -o, --overwrite       Overwrite any existing files/directories (default:
                            False)
      -f FILE, --file FILE  Name for the ascii output file. Any leading directory
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
      -G, --gui             View the obs log in a GUI (default: False)
    