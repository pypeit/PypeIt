.. code-block:: console

    $ pypeit_ql -h
    usage: pypeit_ql [-h] [--raw_files RAW_FILES [RAW_FILES ...]]
                     [--raw_path RAW_PATH] [--sci_files SCI_FILES [SCI_FILES ...]]
                     [--redux_path REDUX_PATH] [--parent_calib_dir PARENT_CALIB_DIR]
                     [--setup_calib_dir SETUP_CALIB_DIR] [--clear_science]
                     [--calibs_only] [--overwrite_calibs] [--det DET [DET ...]]
                     [--slitspatnum SLITSPATNUM] [--maskID MASKID]
                     [--boxcar_radius BOXCAR_RADIUS] [--snr_thresh SNR_THRESH]
                     [--ignore_std] [--skip_display] [--coadd2d]
                     [--only_slits ONLY_SLITS [ONLY_SLITS ...]] [--offsets OFFSETS]
                     [--weights WEIGHTS] [--spec_samp_fact SPEC_SAMP_FACT]
                     [--spat_samp_fact SPAT_SAMP_FACT] [--try_old]
                     spectrograph
    
    Script to produce quick-look PypeIt reductions
    
    positional arguments:
      spectrograph          A valid spectrograph identifier: bok_bc,
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
      --raw_files RAW_FILES [RAW_FILES ...]
                            Either a PypeIt-formatted input file with the list of
                            raw images to process and the relevant path, or a space-
                            separated list of the filenames (e.g., "img1.fits
                            img2.fits"). For the latter entry mode, the path
                            containing the files is set using --raw_path. (default:
                            None)
      --raw_path RAW_PATH   Directory with the raw files to process. Ignored if a
                            PypeIt-formatted file is provided using the --rawfiles
                            option. (default: current working directory)
      --sci_files SCI_FILES [SCI_FILES ...]
                            A space-separated list of raw file names that are
                            science exposures. These files must *also* be in the
                            list of raw files. Use of this option overrides the
                            automated PypeIt frame typing. Should only be used of
                            automatic frame typing fails or is undesirable.
                            (default: None)
      --redux_path REDUX_PATH
                            Path for the QL reduction outputs. (default: current
                            working directory)
      --parent_calib_dir PARENT_CALIB_DIR
                            Directory with/for calibrations for *all* instrument
                            configurations/setups. If provided, the data for your
                            instrument configuration will be placed or pulled from a
                            relevant sub-directory. If None, the redux_path is used.
                            (default: None)
      --setup_calib_dir SETUP_CALIB_DIR
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
      --clear_science       Remove the existing output science directories to force
                            a fresh reduction. If False, any existing directory
                            structure will remain, and any alterations to existing
                            science files will follow the normal behavior of
                            run_pypeit. (default: False)
      --calibs_only         Reduce only the calibrations? (default: False)
      --overwrite_calibs    Re-process and overwrite any existing calibration files.
                            (default: False)
      --det DET [DET ...]   A space-separated set of detectors or detector mosaics
                            to reduce. By default, *all* detectors or default
                            mosaics for this instrument will be reduced. Detectors
                            in a mosaic must be a mosaic "allowed" by PypeIt and
                            should be provided as comma-separated integers (with no
                            spaces). For example, to separately reduce detectors 1
                            and 5 for Keck/DEIMOS, you would use --det 1 5; to
                            reduce mosaics made up of detectors 1,5 and 3,7, you
                            would use --det 1,5 3,7 (default: None)
      --slitspatnum SLITSPATNUM
                            Reduce the slit(s) as specified by the slitspatnum
                            value(s) (default: None)
      --maskID MASKID       Reduce the slit(s) as specified by the maskID value(s)
                            (default: None)
      --boxcar_radius BOXCAR_RADIUS
                            Set the radius for the boxcar extraction in arcseconds
                            (default: None)
      --snr_thresh SNR_THRESH
                            Change the default S/N threshold used during source
                            detection (default: None)
      --ignore_std          If standard star observations are automatically
                            detected, ignore those frames. Otherwise, they are
                            included with the reduction of the science frames.
                            (default: False)
      --skip_display        Run the quicklook without displaying any results.
                            (default: True)
      --coadd2d             Perform default 2D coadding. (default: False)
      --only_slits ONLY_SLITS [ONLY_SLITS ...]
                            If coadding, only coadd this space-separated set of
                            slits. If not provided, all slits are coadded. (default:
                            None)
      --offsets OFFSETS     If coadding, spatial offsets to apply to each image; see
                            the [coadd2d][offsets] parameter. Options are restricted
                            here to either maskdef_offsets or auto. If not
                            specified, the (spectrograph-specific) default is used.
                            (default: None)
      --weights WEIGHTS     If coadding, weights used to coadd images; see the
                            [coadd2d][weights] parameter. Options are restricted
                            here to either uniform or auto. If not specified, the
                            (spectrograph-specific) default is used. (default: None)
      --spec_samp_fact SPEC_SAMP_FACT
                            If coadding, adjust the wavelength grid sampling by this
                            factor. For a finer grid, set value to <1.0; for coarser
                            sampling, set value to >1.0). (default: 1.0)
      --spat_samp_fact SPAT_SAMP_FACT
                            If coadding, adjust the spatial grid sampling by this
                            factor. For a finer grid, set value to <1.0; for coarser
                            sampling, set value to >1.0). (default: 1.0)
      --try_old             Attempt to load old datamodel versions. A crash may
                            ensue.. (default: False)
    