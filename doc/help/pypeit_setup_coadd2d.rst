.. code-block:: console

    $ pypeit_setup_coadd2d -h
    [1;34musage: [0m[1;35mpypeit_setup_coadd2d[0m [[32m-h[0m] [[32m-v [33mVERBOSITY[0m] [[36m--log_file [33mLOG_FILE[0m]
                                [[36m--log_level [33mLOG_LEVEL[0m] ([32m-f [33mPYPEIT_FILE[0m |
                                [32m-d [33mSCIENCE_DIR [SCIENCE_DIR ...][0m) [[36m--keep_par[0m]
                                [[36m--obj [33mOBJ [OBJ ...][0m] [[36m--det [33mDET [DET ...][0m]
                                [[36m--only_slits [33mONLY_SLITS [ONLY_SLITS ...][0m]
                                [[36m--exclude_slits [33mEXCLUDE_SLITS [EXCLUDE_SLITS ...][0m]
                                [[36m--spat_toler [33mSPAT_TOLER[0m] [[36m--offsets [33mOFFSETS[0m]
                                [[36m--weights [33mWEIGHTS[0m]
                                [[36m--spec_samp_fact [33mSPEC_SAMP_FACT[0m]
                                [[36m--spat_samp_fact [33mSPAT_SAMP_FACT[0m]
    
    Prepare a configuration file for performing 2D coadds
    
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
      [1;32m-f[0m, [1;36m--pypeit_file[0m [1;33mPYPEIT_FILE[0m
                            PypeIt reduction file (default: None)
      [1;32m-d[0m, [1;36m--science_dir[0m [1;33mSCIENCE_DIR [SCIENCE_DIR ...][0m
                            One or more directories with spec2d files to stack (use
                            wildcard to specify multiple directories). (default:
                            None)
      [1;36m--keep_par[0m            Propagate all parameters from the pypeit file to the
                            coadd2d file(s). If not set, only the required
                            parameters and their default values are included in the
                            output file(s). (default: True)
      [1;36m--obj[0m [1;33mOBJ [OBJ ...][0m   Limit the coadd2d files created to observations of the
                            specified target. If not provided, a coadd2D file is
                            written for each target found in the science directory.
                            The target names are included in the PypeIt spec2d file
                            names.For example, the target for spec2d file "spec2d_cN
                            20170331S0216-pisco_GNIRS_20170331T085412.181.fits" is
                            "pisco". (default: None)
      [1;36m--det[0m [1;33mDET [DET ...][0m   A space-separated set of detectors or detector mosaics
                            to coadd. By default, *all* detectors or default mosaics
                            for this instrument will be coadded. Detectors in a
                            mosaic must be a mosaic "allowed" by PypeIt and should
                            be provided as comma-separated integers (with no
                            spaces). For example, to separately coadd detectors 1
                            and 5 for Keck/DEIMOS, you would use --det 1 5; to coadd
                            mosaics made up of detectors 1,5 and 3,7, you would use
                            --det 1,5 3,7 (default: None)
      [1;36m--only_slits[0m [1;33mONLY_SLITS [ONLY_SLITS ...][0m
                            A space-separated set of slits to coadd. Example syntax
                            for argument is DET01:175,DET02:205 or MSC02:2234. If
                            not provided, all slits are coadded. If both --det and
                            --only_slits are provided, --det will be ignored. This
                            and --exclude_slits are mutually exclusive. If both are
                            provided, --only_slits takes precedence. (default: None)
      [1;36m--exclude_slits[0m [1;33mEXCLUDE_SLITS [EXCLUDE_SLITS ...][0m
                            A space-separated set of slits to exclude in the
                            coaddition. This and --only_slits are mutually
                            exclusive. If both are provided, --only_slits takes
                            precedence. (default: None)
      [1;36m--spat_toler[0m [1;33mSPAT_TOLER[0m
                            Desired tolerance in spatial pixel used to identify
                            slits in different exposures. If not provided, the
                            default value for the specific instrument/configuration
                            is used. (default: None)
      [1;36m--offsets[0m [1;33mOFFSETS[0m     Spatial offsets to apply to each image; see the
                            [coadd2d][offsets] parameter. Options are restricted
                            here to either maskdef_offsets or auto. If not
                            specified, the (spectrograph-specific) default is used.
                            Other options exist but must be entered by directly
                            editing the coadd2d file. (default: None)
      [1;36m--weights[0m [1;33mWEIGHTS[0m     Weights used to coadd images; see the [coadd2d][weights]
                            parameter. Options are restricted here to either uniform
                            or auto. If not specified, the (spectrograph-specific)
                            default is used. Other options exist but must be entered
                            by directly editing the coadd2d file. (default: None)
      [1;36m--spec_samp_fact[0m [1;33mSPEC_SAMP_FACT[0m
                            Make the wavelength grid finer (spec_samp_fact < 1.0) or
                            coarser (spec_samp_fact > 1.0) by this sampling factor,
                            i.e. units of spec_samp_fact are pixels. (default: 1.0)
      [1;36m--spat_samp_fact[0m [1;33mSPAT_SAMP_FACT[0m
                            Make the spatial grid finer (spat_samp_fact < 1.0) or
                            coarser (spat_samp_fact > 1.0) by this sampling factor,
                            i.e. units of spat_samp_fact are pixels. (default: 1.0)
    