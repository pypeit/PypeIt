.. code-block:: console

    $ pypeit_ql -h
    [1;34musage: [0m[1;35mpypeit_ql[0m [[32m-h[0m] [[32m-v [33mVERBOSITY[0m] [[36m--log_file [33mLOG_FILE[0m]
                     [[36m--log_level [33mLOG_LEVEL[0m] [[36m--raw_files [33mRAW_FILES [RAW_FILES ...][0m]
                     [[36m--raw_path [33mRAW_PATH[0m] [[36m--sci_files [33mSCI_FILES [SCI_FILES ...][0m]
                     [[36m--redux_path [33mREDUX_PATH[0m] [[36m--parent_calib_dir [33mPARENT_CALIB_DIR[0m]
                     [[36m--setup_calib_dir [33mSETUP_CALIB_DIR[0m] [[36m--clear_science[0m]
                     [[36m--calibs_only[0m] [[36m--overwrite_calibs[0m] [[36m--det [33mDET [DET ...][0m]
                     [[36m--slitspatnum [33mSLITSPATNUM[0m] [[36m--maskID [33mMASKID[0m]
                     [[36m--boxcar_radius [33mBOXCAR_RADIUS[0m] [[36m--snr_thresh [33mSNR_THRESH[0m]
                     [[36m--ignore_std[0m] [[36m--skip_display[0m] [[36m--removetrace[0m] [[36m--coadd2d[0m]
                     [[36m--spec_samp_fact [33mSPEC_SAMP_FACT[0m]
                     [[36m--spat_samp_fact [33mSPAT_SAMP_FACT[0m] [[36m--offsets [33mOFFSETS[0m]
                     [[36m--weights [33mWEIGHTS[0m] [[36m--only_slits [33mONLY_SLITS [ONLY_SLITS ...][0m]
                     [[36m--try_old[0m]
                     [32mspectrograph[0m
    
    Script to produce quick-look PypeIt reductions
    
    [1;34mpositional arguments:[0m
      [1;32mspectrograph[0m          A valid spectrograph identifier: aat_uhrf, apf_levy,
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
      [1;36m--raw_files[0m [1;33mRAW_FILES [RAW_FILES ...][0m
                            Either a PypeIt-formatted input file with the list of
                            raw images to process and the relevant path, or a space-
                            separated list of the filenames (e.g., "img1.fits
                            img2.fits"). For the latter entry mode, the path
                            containing the files is set using --raw_path. (default:
                            None)
      [1;36m--raw_path[0m [1;33mRAW_PATH[0m   Directory with the raw files to process. Ignored if a
                            PypeIt-formatted file is provided using the --rawfiles
                            option. (default: current working directory)
      [1;36m--sci_files[0m [1;33mSCI_FILES [SCI_FILES ...][0m
                            A space-separated list of raw file names that are
                            science exposures. These files must *also* be in the
                            list of raw files. Use of this option overrides the
                            automated PypeIt frame typing. Should only be used of
                            automatic frame typing fails or is undesirable.
                            (default: None)
      [1;36m--redux_path[0m [1;33mREDUX_PATH[0m
                            Path for the QL reduction outputs. (default: current
                            working directory)
      [1;36m--parent_calib_dir[0m [1;33mPARENT_CALIB_DIR[0m
                            Directory with/for calibrations for *all* instrument
                            configurations/setups. If provided, the data for your
                            instrument configuration will be placed or pulled from a
                            relevant sub-directory. If None, the redux_path is used.
                            (default: None)
      [1;36m--setup_calib_dir[0m [1;33mSETUP_CALIB_DIR[0m
                            Directory with/for calibrations specific to your
                            instrument configuration/setup. Use of this option
                            circumvents the automated naming system for the
                            configuration/setup sub-directories. If None, the code
                            will try to find relevant calibrations in the
                            parent_calib_dir. If no calibrations exist in that
                            directory that match the instrument setup/configuration
                            of the provided data, the code will construct new
                            calibrations (assuming relevant raw files are provided).
                            (default: None)
      [1;36m--clear_science[0m       Remove the existing output science directories to force
                            a fresh reduction. If False, any existing directory
                            structure will remain, and any alterations to existing
                            science files will follow the normal behavior of
                            run_pypeit. (default: False)
      [1;36m--calibs_only[0m         Reduce only the calibrations? (default: False)
      [1;36m--overwrite_calibs[0m    Re-process and overwrite any existing calibration files.
                            (default: False)
      [1;36m--det[0m [1;33mDET [DET ...][0m   A space-separated set of detectors or detector mosaics
                            to reduce. By default, *all* detectors or default
                            mosaics for this instrument will be reduced. Detectors
                            in a mosaic must be a mosaic "allowed" by PypeIt and
                            should be provided as comma-separated integers (with no
                            spaces). For example, to separately reduce detectors 1
                            and 5 for Keck/DEIMOS, you would use --det 1 5; to
                            reduce mosaics made up of detectors 1,5 and 3,7, you
                            would use --det 1,5 3,7 (default: None)
      [1;36m--slitspatnum[0m [1;33mSLITSPATNUM[0m
                            Reduce the slit(s) as specified by the slitspatnum
                            value(s) (default: None)
      [1;36m--maskID[0m [1;33mMASKID[0m       Reduce the slit(s) as specified by the maskID value(s)
                            (default: None)
      [1;36m--boxcar_radius[0m [1;33mBOXCAR_RADIUS[0m
                            Set the radius for the boxcar extraction in arcseconds
                            (default: None)
      [1;36m--snr_thresh[0m [1;33mSNR_THRESH[0m
                            Change the default S/N threshold used during source
                            detection (default: None)
      [1;36m--ignore_std[0m          If standard star observations are automatically
                            detected, ignore those frames. Otherwise, they are
                            included with the reduction of the science frames.
                            (default: False)
      [1;36m--skip_display[0m        Run the quicklook without displaying any results. The
                            default skip_display=False will show the results.
                            (default: False)
      [1;36m--removetrace[0m         When the image is shown, do not overplot traces in the
                            skysub, sky_resid, and resid channels (default: False)
      [1;36m--coadd2d[0m             Perform default 2D coadding. (default: False)
      [1;36m--spec_samp_fact[0m [1;33mSPEC_SAMP_FACT[0m
                            If coadding, adjust the wavelength grid sampling by this
                            factor. For a finer grid, set value to <1.0; for coarser
                            sampling, set value to >1.0). (default: 1.0)
      [1;36m--spat_samp_fact[0m [1;33mSPAT_SAMP_FACT[0m
                            If coadding, adjust the spatial grid sampling by this
                            factor. For a finer grid, set value to <1.0; for coarser
                            sampling, set value to >1.0). (default: 1.0)
      [1;36m--offsets[0m [1;33mOFFSETS[0m     If coadding, spatial offsets to apply to each image; see
                            the [coadd2d][offsets] parameter. Options are restricted
                            here to either maskdef_offsets or auto. If not
                            specified, the (spectrograph-specific) default is used.
                            (default: None)
      [1;36m--weights[0m [1;33mWEIGHTS[0m     If coadding, weights used to coadd images; see the
                            [coadd2d][weights] parameter. Options are restricted
                            here to either uniform or auto. If not specified, the
                            (spectrograph-specific) default is used. (default: None)
      [1;36m--only_slits[0m [1;33mONLY_SLITS [ONLY_SLITS ...][0m
                            If coadding, only coadd this space-separated set of
                            slits. If not provided, all slits are coadded. (default:
                            None)
      [1;36m--try_old[0m             Attempt to load old datamodel versions. A crash may
                            ensue.. (default: False)
    