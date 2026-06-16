.. code-block:: console

    $ pypeit_install_telluric -h
    [1;34musage: [0m[1;35mpypeit_install_telluric[0m [[32m-h[0m] [[32m-v [33mVERBOSITY[0m] [[36m--log_file [33mLOG_FILE[0m]
                                   [[36m--log_level [33mLOG_LEVEL[0m] [[36m--force_update[0m]
                                   [[36m--local_file[0m]
                                   [32mfiles [files ...][0m
    
    Script to download/install PypeIt telluric files
    
    [1;34mpositional arguments:[0m
      [1;32mfiles[0m                 Exact paths to TelFits files to be downloaded from the
                            Cloud and installed in the PypeIt cache
    
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
      [1;36m--force_update[0m        Force download of latest version of the telluric grid
                            (default: False)
      [1;36m--local_file[0m          This is a local file to be installed in the cache
                            (default: False)
    