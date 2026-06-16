.. code-block:: console

    $ pypeit_chk_scattlight -h
    [1;34musage: [0m[1;35mpypeit_chk_scattlight[0m [[32m-h[0m] [[32m-v [33mVERBOSITY[0m] [[36m--log_file [33mLOG_FILE[0m]
                                 [[36m--log_level [33mLOG_LEVEL[0m] [[36m--spec2d [33mSPEC2D[0m]
                                 [[36m--det [33mDET[0m] [[36m--mask [33mMASK[0m] [[36m--try_old[0m]
                                 [32mfile[0m [32mslits[0m
    
    Display the scattered light image in a Ginga viewer
    
    [1;34mpositional arguments:[0m
      [1;32mfile[0m                  PypeIt Scattered Light file [e.g.
                            ScatteredLight_A_0_DET01.fits.gz]
      [1;32mslits[0m                 Slits calibration file [e.g. Slits_A_0_DET01.fits.gz]
    
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
      [1;36m--spec2d[0m [1;33mSPEC2D[0m       PypeIt science spec2d file (default: None)
      [1;36m--det[0m [1;33mDET[0m             Detector name or number. If a number, the name is
                            constructed assuming the reduction is for a single
                            detector. If a string, it must match the name of the
                            detector object (e.g., DET01 for a detector, MSC01 for a
                            mosaic). (default: 1)
      [1;36m--mask[0m [1;33mMASK[0m           If True, the detector pixels that are considered on the
                            slit will be masked to highlight the scattered light
                            regions (default: False)
      [1;36m--try_old[0m             Attempt to load old datamodel versions. A crash may
                            ensue.. (default: False)
    