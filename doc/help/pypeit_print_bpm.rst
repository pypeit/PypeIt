.. code-block:: console

    $ pypeit_print_bpm -h
    [1;34musage: [0m[1;35mpypeit_print_bpm[0m [[32m-h[0m] [[32m-v [33mVERBOSITY[0m] [[36m--log_file [33mLOG_FILE[0m]
                            [[36m--log_level [33mLOG_LEVEL[0m] [[36m--file [33mFILE[0m] [[36m--det [33mDET[0m]
                            [32mbit[0m
    
    Print out an informative description of a bad pixel masked value. Usually, you
    should run pypeit_show_2dspec --showmask first to see the bad pixel mask values.
    Then, call this script with the BPM value that you want to findmore information
    about.
    
    [1;34mpositional arguments:[0m
      [1;32mbit[0m                   Bad pixel mask value to describe in plain text
    
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
      [1;36m--file[0m [1;33mFILE[0m           PypeIt spec2d file to use for the description(optional).
                            If provided, the bitmask contained in the spec2d file
                            will be used to describe the bad pixel mask value. If
                            not provided, the default pypeit bad pixel mask will be
                            used. (default: None)
      [1;36m--det[0m [1;33mDET[0m             Detector name or number. If a number, the name is
                            constructed assuming the reduction is for a single
                            detector. If a string, it must match the name of the
                            detector object (e.g., DET01 for a detector, MSC01 for a
                            mosaic). This is not required, and the value is
                            acceptable. Default is 1. (default: 1)
    