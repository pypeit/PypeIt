.. code-block:: console

    $ pypeit_coadd_2dspec -h
    [1;34musage: [0m[1;35mpypeit_coadd_2dspec[0m [[32m-h[0m] [[32m-v [33mVERBOSITY[0m] [[36m--log_file [33mLOG_FILE[0m]
                               [[36m--log_level [33mLOG_LEVEL[0m] [[36m--show[0m] [[36m--debug_offsets[0m]
                               [[36m--peaks[0m] [[36m--basename [33mBASENAME[0m] [[36m--debug[0m]
                               [32mcoadd2d_file[0m
    
    Coadd 2D spectra produced by PypeIt
    
    [1;34mpositional arguments:[0m
      [1;32mcoadd2d_file[0m          File to guide 2d coadds
    
    [1;34moptions:[0m
      [1;32m-h[0m, [1;36m--help[0m            show this help message and exit
      [1;32m-v[0m, [1;36m--verbosity[0m [1;33mVERBOSITY[0m
                            Verbosity level, which must be 0, 1, or 2. Level 0
                            includes warning and error messages, level 1 adds
                            informational messages, and level 2 adds debugging
                            messages and the calling sequence. (default: 2)
      [1;36m--log_file[0m [1;33mLOG_FILE[0m   Name for the log file. If set to "default", a default
                            name is used. If None, a log file is not produced.
                            (default: default)
      [1;36m--log_level[0m [1;33mLOG_LEVEL[0m
                            Verbosity level for the log file. If a log file is
                            produce and this is None, the file log will match the
                            console stream log. (default: None)
      [1;36m--show[0m                Show the reduction steps. Equivalent to the -s option
                            when running pypeit. (default: False)
      [1;36m--debug_offsets[0m       Show QA plots useful for debugging automatic offset
                            determination (default: False)
      [1;36m--peaks[0m               Show the peaks found by the object finding algorithm.
                            (default: False)
      [1;36m--basename[0m [1;33mBASENAME[0m   Basename of files to save the parameters, spec1d, and
                            spec2d (default: None)
      [1;36m--debug[0m               show debug plots? (default: False)
    