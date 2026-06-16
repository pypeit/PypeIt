.. code-block:: console

    $ pypeit_flux_calib -h
    [1;34musage: [0m[1;35mpypeit_flux_calib[0m [[32m-h[0m] [[32m-v [33mVERBOSITY[0m] [[36m--log_file [33mLOG_FILE[0m]
                             [[36m--log_level [33mLOG_LEVEL[0m] [[36m--par_outfile[0m] [[36m--try_old[0m]
                             [32mflux_file[0m
    
    Flux calibrate 1D spectra produced by PypeIt
    
    [1;34mpositional arguments:[0m
      [1;32mflux_file[0m             File to guide fluxing process.  This file must have the
                            following format:
                             
                            flux read
                                 filename | sensfile
                              spec1dfile1 | sensfile1
                              spec1dfile2 | 
                                 ...    
                            flux end
                             
                            OR
                             
                            flux read
                                 filename | sensfile
                              spec1dfile1 | sensfile1
                              spec1dfile2 | sensfile2
                              spec1dfile3 | sensfile3
                                 ...    
                            flux end
                             
                            OR
                             
                            [fluxcalib]
                              use_archived_sens = True
                            flux read
                                 filename
                              spec1dfile1
                              spec1dfile2
                              spec1dfile3
                                 ...    
                            flux end
                             
                            That is, you must specify either a sensfile for all
                            spec1dfiles on the first line, specify one sensfile for
                            each spec1dfile, or specify no sensfiles and use an
                            archived one.
                            Archived sensfiles are available for the following
                            spectrographs: keck_deimos
                             
    
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
      [1;36m--par_outfile[0m         Output to save the parameters
      [1;36m--try_old[0m             Attempt to load old datamodel versions. A crash may
                            ensue..
    