.. code-block:: console

    $ pypeit_cache_github_data -h
    usage: pypeit_cache_github_data [-h] [--exclude EXCLUDE [EXCLUDE ...] |
                                    --include INCLUDE [INCLUDE ...]]
                                    [--spec_dependent_only] [--force_update]
                                    spectrograph [spectrograph ...]
    
    Script to download/cache PypeIt github data
    
    positional arguments:
      spectrograph          A valid spectrograph identifier: aat_uhrf, bok_bc,
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
      --exclude EXCLUDE [EXCLUDE ...]
                            A subset of the directories to *exclude* from the list
                            of files to download. Options are: tests, reid_arxiv,
                            nist, standards, skisim, sensfunc. This option is
                            mutually exclusive with --include. (default: ['tests'])
      --include INCLUDE [INCLUDE ...]
                            The directories to *include* in the list of files to
                            download. Use "--include all" to include all
                            directories. Options are: all, tests, reid_arxiv, nist,
                            standards, skisim, sensfunc. This option is mutually
                            exclusive with --exclude. (default: None)
      --spec_dependent_only
                            Only include files that are specific to the provided
                            list of spectrographs. By default, the script also
                            includes any files in the selected directories that are
                            *not* specific to a given spectrograph (e.g.,
                            atmospheric extinction curves). (default: False)
      --force_update        Force re-download of existing files (default: False)
    