.. code-block:: console

    $ pypeit_chk_alignments -h
    [1;34musage: [0m[1;35mpypeit_chk_alignments[0m [[32m-h[0m] [[32m-v [33mVERBOSITY[0m] [[36m--log_file [33mLOG_FILE[0m]
                                 [[36m--log_level [33mLOG_LEVEL[0m] [[36m--chname [33mCHNAME[0m]
                                 [[36m--try_old[0m]
                                 [32mfile[0m
    
    Display Alignment image and the trace data
    
    [1;34mpositional arguments:[0m
      [1;32mfile[0m                  PypeIt Alignment file [e.g. Alignment_A_1_DET01.fits]
    
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
      [1;36m--chname[0m [1;33mCHNAME[0m       Channel name for image in Ginga (default: Alignments)
      [1;36m--try_old[0m             Attempt to load old datamodel versions. A crash may
                            ensue.. (default: False)
    