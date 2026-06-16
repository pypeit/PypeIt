.. code-block:: console

    $ pypeit_coadd_datacube -h
    [1;34musage: [0m[1;35mpypeit_coadd_datacube[0m [[32m-h[0m] [[32m-v [33mVERBOSITY[0m] [[36m--log_file [33mLOG_FILE[0m]
                                 [[36m--log_level [33mLOG_LEVEL[0m] [[36m--det [33mDET[0m] [[32m-o[0m]
                                 [32mfile[0m
    
    Read in an array of spec2D files and convert them into a datacube
    
    [1;34mpositional arguments:[0m
      [1;32mfile[0m                  filename.coadd3d file
    
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
      [1;36m--det[0m [1;33mDET[0m             Detector (default: 1)
      [1;32m-o[0m, [1;36m--overwrite[0m       Overwrite any existing files/directories (default:
                            False)
    