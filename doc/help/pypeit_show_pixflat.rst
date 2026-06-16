.. code-block:: console

    $ pypeit_show_pixflat -h
    [1;34musage: [0m[1;35mpypeit_show_pixflat[0m [[32m-h[0m] [[32m-v [33mVERBOSITY[0m] [[36m--log_file [33mLOG_FILE[0m]
                               [[36m--log_level [33mLOG_LEVEL[0m] [[36m--det [33mDET [DET ...][0m]
                               [32mfile[0m
    
    Show an archived Pixel Flat image in a ginga window.
    
    [1;34mpositional arguments:[0m
      [1;32mfile[0m                  Pixel Flat filename, e.g.
                            pixelflat_keck_lris_blue.fits.gz
    
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
      [1;36m--det[0m [1;33mDET [DET ...][0m   Detector(s) to show. If more than one, list the
                            detectors as, e.g. --det 1 2 to show detectors 1 and 2.
                            If not provided, all detectors will be shown. (default:
                            None)
    