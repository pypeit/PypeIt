.. code-block:: console

    $ pypeit_chk_edges -h
    [1;34musage: [0m[1;35mpypeit_chk_edges[0m [[32m-h[0m] [[32m-v [33mVERBOSITY[0m] [[36m--log_file [33mLOG_FILE[0m]
                            [[36m--log_level [33mLOG_LEVEL[0m] [[36m--slits_file [33mSLITS_FILE[0m]
                            [[36m--mpl[0m] [[36m--try_old[0m]
                            [32mtrace_file[0m
    
    Display trace image and edge traces
    
    [1;34mpositional arguments:[0m
      [1;32mtrace_file[0m            PypeIt Edges file [e.g. Edges_A_0_DET01.fits.gz]
    
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
      [1;36m--slits_file[0m [1;33mSLITS_FILE[0m
                            PypeIt Slits file [e.g. Slits_A_1_01.fits]. If this file
                            does not exist or is not provided, PypeIt will attempt
                            to read the default file name (in the Calibrations
                            directory). Ignored if plotting using a matplotlib
                            window instead of ginga. (default: None)
      [1;36m--mpl[0m                 Use a matplotlib window instead of ginga to show the
                            trace (default: False)
      [1;36m--try_old[0m             Attempt to load old datamodel versions. A crash may
                            ensue.. (default: False)
    