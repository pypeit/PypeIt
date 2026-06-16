.. code-block:: console

    $ pypeit_collate_1d -h
    [1;34musage: [0m[1;35mpypeit_collate_1d[0m [[32m-h[0m] [[32m-v [33mVERBOSITY[0m] [[36m--log_file [33mLOG_FILE[0m]
                             [[36m--log_level [33mLOG_LEVEL[0m]
                             [[36m--spec1d_files [33m[SPEC1D_FILES ...][0m]
                             [[36m--par_outfile [33mPAR_OUTFILE[0m] [[36m--outdir [33mOUTDIR[0m]
                             [[36m--spec1d_outdir [33mSPEC1D_OUTDIR[0m] [[36m--tolerance [33mTOLERANCE[0m]
                             [[36m--match_using [33mMATCH_USING[0m] [[36m--dry_run[0m] [[36m--ignore_flux[0m]
                             [[36m--flux[0m]
                             [[36m--exclude_slit_trace_bm [33mEXCLUDE_SLIT_TRACE_BM[0m]
                             [[36m--exclude_serendip[0m] [[36m--wv_rms_thresh [33mWV_RMS_THRESH[0m]
                             [[36m--refframe [33m{observed,heliocentric,barycentric}[0m]
                             [[36m--chk_version[0m]
                             [32m[input_file][0m
    
    Flux/Coadd multiple 1d spectra from multiple nights and prepare a directory for
    the KOA.
    
    [1;34mpositional arguments:[0m
      [1;32minput_file[0m            (Optional) File for guiding the collate process.
                            Parameters in this file are overidden by the command
                            line. The file must have the following format:
                             
                            [collate1d]
                              tolerance             <tolerance>
                              outdir                <directory to place output files>
                              spec1d_outdir         <directory to place modified spec1ds, if any>
                              exclude_slit_trace_bm <slit types to exclude>
                              exclude_serendip      If set serendipitous objects are skipped.
                              match_using           Whether to match using "pixel" or
                                                    "ra/dec"
                              dry_run               If set the matches are displayed
                                                    without any processing
                              flux                  Flux calibrate using archived sensfuncs.
                              ignore_flux           Ignore any flux calibration information in
                                                    spec1d files.
                              wv_rms_thresh         If set, any objects with a wavelength rms > than the input
                                                    value are skipped, else all wavelength rms values are accepted.
                              refframe              Perform reference frame correction prior to coadding.
                                                    Options are ['observed', 'heliocentric', 'barycentric']. Defaults to None.
                             
                            spec1d read
                            <path to spec1d files, wildcards allowed>
                            ...
                            end
    
    [1;34moptions:[0m
      [1;32m-h[0m, [1;36m--help[0m            show this help message and exit
      [1;32m-v[0m, [1;36m--verbosity[0m [1;33mVERBOSITY[0m
                            Verbosity level, which must be 0, 1, or 2. Level 0
                            includes warning and error messages, level 1 adds
                            informational messages, and level 2 adds debugging
                            messages and the calling sequence.
      [1;36m--log_file[0m [1;33mLOG_FILE[0m   Name for the log file. If set to "default", a default
                            name is used. If None, a log file is not produced.
      [1;36m--log_level[0m [1;33mLOG_LEVEL[0m
                            Verbosity level for the log file. If a log file is
                            produce and this is None, the file log will match the
                            console stream log.
      [1;36m--spec1d_files[0m [1;33m[SPEC1D_FILES ...][0m
                            One or more spec1d files to flux/coadd/archive. Can
                            contain wildcards
      [1;36m--par_outfile[0m [1;33mPAR_OUTFILE[0m
                            Output to save the parameters
      [1;36m--outdir[0m [1;33mOUTDIR[0m       The path where all coadded output files and report files
                            will be placed. Defaults to the current directory.
      [1;36m--spec1d_outdir[0m [1;33mSPEC1D_OUTDIR[0m
                            The path where all modified spec1d files are placed.
                            These are only created if flux calibration or refframe
                            correction are asked for. Defaults to overwriting
                            existing spec1ds.
      [1;36m--tolerance[0m [1;33mTOLERANCE[0m
                            The tolerance used when comparing the coordinates of
                            objects. If two objects are within this distance from
                            each other, they are considered the same object. If
                            match_using is 'ra/dec' (the default) this is an angular
                            distance. The defaults units are arcseconds but other
                            units supported by astropy.coordinates.Angle can be used
                            (`e.g.`, '0.003d' or '0h1m30s'). If match_using is
                            'pixel' this is a float.
      [1;36m--match_using[0m [1;33mMATCH_USING[0m
                            Determines how 1D spectra are matched as being the same
                            object. Must be either 'pixel' or 'ra/dec'.
      [1;36m--dry_run[0m             If set, the script will display the matching File and
                            Object Ids but will not flux, coadd or archive.
      [1;36m--ignore_flux[0m         If set, the script will only coadd non-fluxed spectra
                            even if flux data is present. Otherwise fluxed spectra
                            are coadded if all spec1ds have been fluxed calibrated.
      [1;36m--flux[0m                If set, the script will flux calibrate using archived
                            sensfuncs before coadding.
      [1;36m--exclude_slit_trace_bm[0m [1;33mEXCLUDE_SLIT_TRACE_BM[0m
                            A list of slit trace bitmask bits that should be
                            excluded. Comma separated.
      [1;36m--exclude_serendip[0m    Whether to exclude SERENDIP objects from collating.
      [1;36m--wv_rms_thresh[0m [1;33mWV_RMS_THRESH[0m
                            If set, any objects with a wavelength RMS > this value
                            are skipped, else all wavelength RMS values are
                            accepted.
      [1;36m--refframe[0m [1;33m{observed,heliocentric,barycentric}[0m
                            Perform reference frame correction prior to coadding.
                            Options are: observed, heliocentric, barycentric
      [1;36m--chk_version[0m         If True enforce strict PypeIt version checking to ensure
                            that all files were created with the current version of
                            PypeIt. If set to False, the code will attempt to read
                            out-of-date files and keep going. Beware (!!) that this
                            can lead to unforeseen bugs that either cause the code
                            to crash or lead to erroneous results. I.e., you really
                            need to know what you are doing if you set this to
                            False!
    