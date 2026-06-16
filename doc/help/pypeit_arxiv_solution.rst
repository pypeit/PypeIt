.. code-block:: console

    $ pypeit_arxiv_solution -h
    [1;34musage: [0m[1;35mpypeit_arxiv_solution[0m [[32m-h[0m] [[32m-v [33mVERBOSITY[0m] [[36m--log_file [33mLOG_FILE[0m]
                                 [[36m--log_level [33mLOG_LEVEL[0m] [[32m-s [33mSLIT[0m] [[36m--try_old[0m]
                                 [32mfile[0m [32mbinning[0m
    
    Read in a WaveCalib solution and convert it into the format required for the
    PypeIt full template archive
    
    [1;34mpositional arguments:[0m
      [1;32mfile[0m                  WaveCalib file
      [1;32mbinning[0m               Spectral binning
    
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
      [1;32m-s[0m, [1;36m--slit[0m [1;33mSLIT[0m       Slit number to use (default: 0)
      [1;36m--try_old[0m             Attempt to load old datamodel versions. A crash may
                            ensue.. (default: False)
    