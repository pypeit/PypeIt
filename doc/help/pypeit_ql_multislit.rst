.. code-block:: console

    $ pypeit_ql_multislit -h
    usage: pypeit_ql_multislit [-h] [--spec_samp_fact SPEC_SAMP_FACT]
                               [--spat_samp_fact SPAT_SAMP_FACT] [--bkg_redux]
                               [--flux] [--mask_cr] [--writefits] [--no_gui]
                               [--box_radius BOX_RADIUS] [--offset OFFSET]
                               [--redux_path REDUX_PATH] [--master_dir MASTER_DIR]
                               [--embed] [--show] [--det [DET ...]]
                               spectrograph full_rawpath files [files ...]
    
    Script to produce quick-look multislit PypeIt reductions
    
    positional arguments:
      spectrograph          A valid spectrograph identifier: bok_bc,
                            gemini_flamingos1, gemini_flamingos2,
                            gemini_gmos_north_e2v, gemini_gmos_north_ham,
                            gemini_gmos_north_ham_ns, gemini_gmos_south_ham,
                            gemini_gnirs, gtc_maat, gtc_osiris, gtc_osiris_plus,
                            jwst_nircam, jwst_nirspec, keck_deimos, keck_hires,
                            keck_kcwi, keck_lris_blue, keck_lris_blue_orig,
                            keck_lris_red, keck_lris_red_mark4, keck_lris_red_orig,
                            keck_mosfire, keck_nires, keck_nirspec_low, lbt_luci1,
                            lbt_luci2, lbt_mods1b, lbt_mods1r, lbt_mods2b,
                            lbt_mods2r, ldt_deveny, magellan_fire,
                            magellan_fire_long, magellan_mage, mdm_osmos_mdm4k,
                            mmt_binospec, mmt_bluechannel, mmt_mmirs, not_alfosc,
                            not_alfosc_vert, ntt_efosc2, p200_dbsp_blue,
                            p200_dbsp_red, p200_tspec, shane_kast_blue,
                            shane_kast_red, shane_kast_red_ret, soar_goodman_blue,
                            soar_goodman_red, tng_dolores, vlt_fors2, vlt_sinfoni,
                            vlt_xshooter_nir, vlt_xshooter_uvb, vlt_xshooter_vis,
                            wht_isis_blue, wht_isis_red
      full_rawpath          Full path to the raw files
      files                 list of frames i.e. img1.fits img2.fits
    
    options:
      -h, --help            show this help message and exit
      --spec_samp_fact SPEC_SAMP_FACT
                            Make the wavelength grid finer (spec_samp_fact < 1.0) or
                            coarser (spec_samp_fact > 1.0) by this sampling factor,
                            i.e. units of spec_samp_fact are pixels. (default: 1.0)
      --spat_samp_fact SPAT_SAMP_FACT
                            Make the spatial grid finer (spat_samp_fact < 1.0) or
                            coarser (spat_samp_fact > 1.0) by this sampling factor,
                            i.e. units of spat_samp_fact are pixels. (default: 1.0)
      --bkg_redux           If set the script will perform difference imaging
                            quicklook. Namely it will identify sequences of AB pairs
                            based on the dither pattern and perform difference
                            imaging sky subtraction and fit for residuals (default:
                            False)
      --flux                This option will multiply in sensitivity function to
                            obtain a flux calibrated 2d spectrum (default: False)
      --mask_cr             This option turns on cosmic ray rejection. This improves
                            the reduction but doubles runtime. (default: False)
      --writefits           Write the ouputs to a fits file (default: False)
      --no_gui              Do not display the results in a GUI (default: False)
      --box_radius BOX_RADIUS
                            Set the radius for the boxcar extraction (default: None)
      --offset OFFSET       Override the automatic offsets determined from the
                            headers. Offset is in pixels. This option is useful if a
                            standard dither pattern was not executed. The offset
                            convention is such that a negative offset will move the
                            (negative) B image to the left. (default: None)
      --redux_path REDUX_PATH
                            Location where reduction outputs should be stored.
                            (default: current working directory)
      --master_dir MASTER_DIR
                            Location of PypeIt Master files used for the reduction.
                            (default: None)
      --embed               Upon completion embed in ipython shell (default: False)
      --show                Show the reduction steps. Equivalent to the -s option
                            when running pypeit. (default: False)
      --det [DET ...]       Detector(s) to show. If more than one, the list of
                            detectors must be one of the allowed mosaics hard-coded
                            for the selected spectrograph. (default: 1)
    