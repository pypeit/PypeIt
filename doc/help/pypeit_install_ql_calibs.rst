.. code-block:: console

    $ pypeit_install_ql_calibs -h
    [1;34musage: [0m[1;35mpypeit_install_ql_calibs[0m [[32m-h[0m] [[32m-v [33mVERBOSITY[0m] [[36m--log_file [33mLOG_FILE[0m]
                                    [[36m--log_level [33mLOG_LEVEL[0m] [[36m--zip [33mZIP[0m |
                                    [36m--ql_path [33mQL_PATH[0m] [[36m--odir [33mODIR[0m] [[36m--rmzip[0m]
    
    Script to install PypeIt QL calibration files
    
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
      [1;36m--zip[0m [1;33mZIP[0m             Zip file of the full QL_CALIB directory downloaded from
                            the PypeIt Google Drive (default: None)
      [1;36m--ql_path[0m [1;33mQL_PATH[0m     An existing directory to symlink as the QL_CALIB
                            directory. (default: None)
      [1;36m--odir[0m [1;33mODIR[0m           The directory in which to extract the zip file. Ignored
                            if a direct path is provided using --ql_path. (default:
                            current working directory)
      [1;36m--rmzip[0m               Remove the downloaded zip file (default: False)
    