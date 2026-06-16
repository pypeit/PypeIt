.. code-block:: console

    $ pypeit_chk_tilts -h
    [1;34musage: [0m[1;35mpypeit_chk_tilts[0m [[32m-h[0m] [[32m-v [33mVERBOSITY[0m] [[36m--log_file [33mLOG_FILE[0m]
                            [[36m--log_level [33mLOG_LEVEL[0m] [[36m--mpl[0m] [[36m--show_traces[0m]
                            [[36m--try_old[0m]
                            [32mfile[0m
    
    Display Tiltimg image and 2D fitted tilts in Ginga viewer or Matplotlib window.
    Tiltimg file must be in the same directory as Tilts.
    
    [1;34mpositional arguments:[0m
      [1;32mfile[0m                  PypeIt Tilts file [e.g. Tilt_A_1_01.fits]
    
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
      [1;36m--mpl[0m                 Use a matplotlib window instead of ginga to show the
                            tilts. Faster plotting. (default: False)
      [1;36m--show_traces[0m         Show the traced tilts. This slows down the plotting
                            (mostly in Ginga). If not set, only the fitted, masked
                            and rejected in the fit tilts are shown. (default:
                            False)
      [1;36m--try_old[0m             Attempt to load old datamodel versions. A crash may
                            ensue.. (default: False)
    