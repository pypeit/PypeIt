.. code-block:: console

    $ pypeit_compile_wvarxiv -h
    [1;34musage: [0m[1;35mpypeit_compile_wvarxiv[0m [[32m-h[0m] [[32m-v [33mVERBOSITY[0m] [[36m--log_file [33mLOG_FILE[0m]
                                  [[36m--log_level [33mLOG_LEVEL[0m] [[36m--append[0m]
                                  [32mwvarxiv_folder[0m [32minstrument[0m [32mgrating[0m
    
    Read in a set of wxarxiv solutions from Identify and compile them into a single
    fits file to be used with the reidentify method.
    
    [1;34mpositional arguments:[0m
      [1;32mwvarxiv_folder[0m        Location of the WVarxiv files
      [1;32minstrument[0m            Name of instrument. e.g. keck_lris_blue, keck_deimos,
                            gemini_gmos_south_ham
      [1;32mgrating[0m               Instrument grating name. E.g. b600, r400, 600_10000.
    
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
      [1;36m--append[0m              Append to an existing file for this instrument.
                            (default: False)
    