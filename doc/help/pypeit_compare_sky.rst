.. code-block:: console

    $ pypeit_compare_sky -h
    [1;34musage: [0m[1;35mpypeit_compare_sky[0m [[32m-h[0m] [[32m-v [33mVERBOSITY[0m] [[36m--log_file [33mLOG_FILE[0m]
                              [[36m--log_level [33mLOG_LEVEL[0m] [[36m--exten [33mEXTEN[0m] [[36m--optimal[0m]
                              [[36m--scale_user [33mSCALE_USER[0m]
                              [32mfile[0m [32mskyfile[0m
    
    Compare the extracted sky spectrum against an archived sky model maintained by
    PypeIt.
    
    [1;34mpositional arguments:[0m
      [1;32mfile[0m                  spec1d Spectral file
      [1;32mskyfile[0m               Archived PypeIt sky file (e.g. paranal_sky.fits)
    
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
      [1;36m--exten[0m [1;33mEXTEN[0m         FITS extension (default: None)
      [1;36m--optimal[0m             Show Optimal? Default is boxcar (default: False)
      [1;36m--scale_user[0m [1;33mSCALE_USER[0m
                            Scale user spectrum by a factor (default: 1.0)
    