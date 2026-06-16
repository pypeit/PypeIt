.. code-block:: console

    $ pypeit_multislit_flexure -h
    [1;34musage: [0m[1;35mpypeit_multislit_flexure[0m [[32m-h[0m] [[32m-v [33mVERBOSITY[0m] [[36m--log_file [33mLOG_FILE[0m]
                                    [[36m--log_level [33mLOG_LEVEL[0m] [[36m--clobber[0m] [[36m--debug[0m]
                                    [32mflex_file[0m [32moutroot[0m
    
    Calculate and apply flexure corrections for 1D spectra produced by PypeIt.
    
    [1;34mpositional arguments:[0m
      [1;32mflex_file[0m             File to guide flexure corrections for this multi-slit
                            mode.  This file must have the following format:
                             
                            flexure read
                              filename
                              spec1dfile1
                              spec1dfile2
                                 ...    
                            flexure end
                             
                             
      [1;32moutroot[0m               Output fileroot for the flexure fits saved as FITS.
    
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
      [1;36m--clobber[0m             Clobber output files
      [1;36m--debug[0m               show debug plots?
    