.. code-block:: console

    $ pypeit_sensfunc -h
    [1;34musage: [0m[1;35mpypeit_sensfunc[0m [[32m-h[0m] [[32m-v [33mVERBOSITY[0m] [[36m--log_file [33mLOG_FILE[0m]
                           [[36m--log_level [33mLOG_LEVEL[0m] [[36m--extr [33m{OPT,BOX}[0m]
                           [[36m--algorithm [33m{UVIS,IR}[0m] [[36m--multi [33mMULTI[0m] [[32m-o [33mOUTFILE[0m]
                           [[32m-s [33mSENS_FILE[0m] [[32m-f[0m] [[36m--debug[0m] [[36m--par_outfile [33mPAR_OUTFILE[0m]
                           [32mspec1dfiles [spec1dfiles ...][0m
    
    Compute a sensitivity function
    
    [1;34mpositional arguments:[0m
      [1;32mspec1dfiles[0m           file(s) of the reduced standard star spectrum. These can
                            be either spec1d*.fits files or the output of
                            `pypeit_coadd_1dspec` (except for cross-dispersed
                            echelle data). Multiple files can be provided, but they
                            are helpful onlyif they cover different wavelength
                            ranges, since thisscript will splice (not combine) them
                            together.
    
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
      [1;36m--extr[0m [1;33m{OPT,BOX}[0m      Override the default extraction method used for
                            computing the sensitivity function.  Note that it is not
                            possible to set --extr and simultaneously use a .sens
                            file with the --sens_file option. If you are using a
                            .sens file, set the algorithm there via:
                             
                                [sensfunc]
                                     extr = BOX
                             
                            The extraction options are: OPT or BOX
      [1;36m--algorithm[0m [1;33m{UVIS,IR}[0m
                            Override the default algorithm for computing the
                            sensitivity function.  Note that it is not possible to
                            set --algorithm and simultaneously use a .sens file with
                            the --sens_file option. If you are using a .sens file,
                            set the algorithm there via:
                             
                                [sensfunc]
                                     algorithm = IR
                             
                            The algorithm options are:
                             
                            UVIS = Should be used for data with lambda < 7000A.  No
                            detailed model of telluric absorption but corrects for
                            atmospheric extinction.
                             
                            IR = Should be used for data with lambbda > 7000A.
                            Performs joint fit for sensitivity function and telluric
                            absorption using HITRAN models.
                             
      [1;36m--multi[0m [1;33mMULTI[0m         List of detector numbers to splice together for
                            instruments with multiple detectors arranged in the
                            spectral direction, e.g. --multi = '3,7'.  Note that it
                            is not possible to set --multi and simultaneously use a
                            .sens file with the --sens_file option.  If you are
                            using a .sens file, set the multi_spec_det param there
                            via:
                             
                                [sensfunc]
                                    multi_spec_det = 3,7
                             
      [1;32m-o[0m, [1;36m--outfile[0m [1;33mOUTFILE[0m
                            Output file for sensitivity function. If not specified,
                            the sensitivity function will be written out to a
                            standard filename in the current working directory, i.e.
                            if the standard spec1d file is named
                            spec1d_b24-Feige66_KASTb_foo.fits the sensfunc will be
                            written to sens_b24-Feige66_KASTb_foo.fits. A QA file
                            will also be written as
                            sens_spec1d_b24-Feige66_KASTb_foo_QA.pdf and a file
                            showing throughput plots to
                            sens_spec1d_b24-Feige66_KASTb_foo_throughput.pdf. The
                            same extensions for QA and throughput will be used if
                            outfile is provided but with .fits trimmed off if it is
                            in the filename.
      [1;32m-s[0m, [1;36m--sens_file[0m [1;33mSENS_FILE[0m
                            Configuration file with sensitivity function parameters
      [1;32m-f[0m, [1;36m--use_flat[0m        Use the extracted spectrum of the flatfield calibration
                            to estimate the blaze function when generating the
                            sensitivity function. This is helpful to account for
                            small scale undulations in the sensitivity function. The
                            spec1dfile must contain the extracted flatfield response
                            in order to use this option. This spectrum is extracted
                            by default, unless you did not compute a pixelflat
                            frame. Note that it is not possible to set --use_flat
                            and simultaneously use a .sens file with the --sens_file
                            option. If you are using a .sens file, set the use_flat
                            flag with the argument:
                             
                                [sensfunc]
                                     use_flat = True
      [1;36m--debug[0m               show debug plots?
      [1;36m--par_outfile[0m [1;33mPAR_OUTFILE[0m
                            Name of output file to save the parameters used by the
                            fit
    