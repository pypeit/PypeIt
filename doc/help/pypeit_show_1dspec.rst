.. code-block:: console

    $ pypeit_show_1dspec -h
    [1;34musage: [0m[1;35mpypeit_show_1dspec[0m [[32m-h[0m] [[32m-v [33mVERBOSITY[0m] [[36m--log_file [33mLOG_FILE[0m]
                              [[36m--log_level [33mLOG_LEVEL[0m] [[36m--list[0m] [[36m--exten [33mEXTEN[0m |
                              [36m--obj [33mOBJ[0m] [[36m--extract [33m{BOX,OPT}[0m] [[36m--flux[0m] [[32m-m[0m]
                              [32mfile[0m
    
    Show a 1D spectrum
    
    [1;34mpositional arguments:[0m
      [1;32mfile[0m                  PypeIt spec1d file (this script does not work with
                            coadd_1dspec output spectra).
    
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
      [1;36m--list[0m                Instead of plotting any spectra, simply list the
                            extensions with spectra (default: False)
      [1;36m--exten[0m [1;33mEXTEN[0m         Number of the extension to plot (default: 1)
      [1;36m--obj[0m [1;33mOBJ[0m             Extension (object) name to plot, e.g.
                            SPAT0424-SLIT0000-DET01 (default: None)
      [1;36m--extract[0m [1;33m{BOX,OPT}[0m   Method used to extract the spectrum (default: OPT)
      [1;36m--flux[0m                Show the flux-calibrated spectrum (if available)
                            (default: False)
      [1;32m-m[0m, [1;36m--unmasked[0m        Only show the unmasked data. (default: True)
    