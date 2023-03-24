.. code-block:: console

    $ pypeit_ql -h
    usage: pypeit_ql [-h] [--rawfile_list RAWFILE_LIST]
                     [--full_rawpath FULL_RAWPATH] [--raw_extension RAW_EXTENSION]
                     [--rawfiles RAWFILES [RAWFILES ...]]
                     [--sci_files SCI_FILES [SCI_FILES ...]]
                     [--redux_path REDUX_PATH] [--calib_dir CALIB_DIR]
                     [--masters_dir MASTERS_DIR] [--calibs_only] [--clobber_calibs]
                     [--slitspatnum SLITSPATNUM] [--maskID MASKID]
                     [--boxcar_radius BOXCAR_RADIUS] [--det DET] [--no_stack]
                     spectrograph
    
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
    
    options:
      -h, --help            show this help message and exit
      --rawfile_list RAWFILE_LIST
                            File providing raw files to reduce including their
                            path(s) (default: None)
      --full_rawpath FULL_RAWPATH
                            Full path to the raw files. Used with --rawfiles or
                            --raw_extension (default: None)
      --raw_extension RAW_EXTENSION
                            Extension for raw files in full_rawpath. Only use if
                            --rawfile_list and --rawfiles are not provided (default:
                            .fits)
      --rawfiles RAWFILES [RAWFILES ...]
                            space separated list of raw frames e.g. img1.fits
                            img2.fits. These must exist within --full_rawpath
                            (default: None)
      --sci_files SCI_FILES [SCI_FILES ...]
                            space separated list of raw frames to be specified as
                            science exposures (over-rides PypeIt frame typing)
                            (default: None)
      --redux_path REDUX_PATH
                            Full path to where QL reduction should be run. (default:
                            current working directory)
      --calib_dir CALIB_DIR
                            Location folders of calibration reductions (default:
                            None)
      --masters_dir MASTERS_DIR
                            Location of PypeIt Master files used for the reduction.
                            (default: None)
      --calibs_only         Reduce only the calibrations? (default: False)
      --clobber_calibs      Clobber existing calibration files? (default: False)
      --slitspatnum SLITSPATNUM
                            Reduce the slit(s) as specified by the slitspatnum
                            value(s) (default: None)
      --maskID MASKID       Reduce the slit(s) as specified by the maskID value(s)
                            (default: None)
      --boxcar_radius BOXCAR_RADIUS
                            Set the radius for the boxcar extraction in arcseconds
                            (default: None)
      --det DET             Detector to reduce. Same format as detnum (default:
                            None)
      --no_stack            Do *not* stack multiple science frames (default: True)
    