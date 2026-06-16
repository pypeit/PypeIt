.. code-block:: console

    $ pypeit_run_to_calibstep -h
    [1;34musage: [0m[1;35mpypeit_run_to_calibstep[0m [[32m-h[0m] [[32m-v [33mVERBOSITY[0m] [[36m--log_file [33mLOG_FILE[0m]
                                   [[36m--log_level [33mLOG_LEVEL[0m]
                                   [[36m--science_frame [33mSCIENCE_FRAME[0m]
                                   [[36m--calib_group [33mCALIB_GROUP[0m] [[36m--det [33mDET[0m]
                                   [[32m-r [33mREDUX_PATH[0m] [[32m-s[0m]
                                   [32mpypeit_file[0m [32mstep[0m
    
    Run PypeIt to a single calibration step for an input frame
    
    [1;34mpositional arguments:[0m
      [1;32mpypeit_file[0m           PypeIt reduction file (must have .pypeit extension)
      [1;32mstep[0m                  Calibration step to perform. Valid steps are: align,
                            arc, bias, bpm, dark, flats, scattlight, slits, tiltimg,
                            tilts, wv_calib
    
    [1;34moptions:[0m
      [1;32m-h[0m, [1;36m--help[0m            show this help message and exit
      [1;32m-v[0m, [1;36m--verbosity[0m [1;33mVERBOSITY[0m
                            Verbosity level, which must be 0, 1, or 2. Level 0
                            includes warning and error messages, level 1 adds
                            informational messages, and level 2 adds debugging
                            messages and the calling sequence.
      [1;36m--log_file[0m [1;33mLOG_FILE[0m   Name for the log file. If set to "default", a default
                            name is used. If None, a log file is not produced.
      [1;36m--log_level[0m [1;33mLOG_LEVEL[0m
                            Verbosity level for the log file. If a log file is
                            produce and this is None, the file log will match the
                            console stream log.
      [1;36m--science_frame[0m [1;33mSCIENCE_FRAME[0m
                            Raw science frame to reduce as listed in your PypeIt
                            file, e.g. b28.fits.gz. Either this or the calib_group
                            must be provided
      [1;36m--calib_group[0m [1;33mCALIB_GROUP[0m
                            Calibration group ID to reduce. Either this or the frame
                            must be provided
      [1;36m--det[0m [1;33mDET[0m             Detector to reduce
      [1;32m-r[0m, [1;36m--redux_path[0m [1;33mREDUX_PATH[0m
                            Path to directory for the reduction. Only advised for
                            testing
      [1;32m-s[0m, [1;36m--show[0m            Show reduction steps via plots (which will block further
                            execution until clicked on) and outputs to ginga.
                            Requires remote control ginga session via "ginga
                            --modules=RC,SlitWavelength &"
    