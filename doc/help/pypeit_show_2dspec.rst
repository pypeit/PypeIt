.. code-block:: console

    $ pypeit_show_2dspec -h
    [1;34musage: [0m[1;35mpypeit_show_2dspec[0m [[32m-h[0m] [[32m-v [33mVERBOSITY[0m] [[36m--log_file [33mLOG_FILE[0m]
                              [[36m--log_level [33mLOG_LEVEL[0m] [[36m--list[0m] [[36m--det [33mDET[0m]
                              [[36m--spat_id [33mSPAT_ID[0m] [[36m--maskID [33mMASKID[0m]
                              [[36m--showmask [33m[SHOWMASK ...][0m] [[36m--removetrace[0m] [[36m--embed[0m]
                              [[36m--ignore_extract_mask[0m] [[36m--channels [33mCHANNELS[0m]
                              [[36m--prefix [33mPREFIX[0m] [[36m--no_clear[0m] [[36m--try_old[0m]
                              [32mfile[0m
    
    Display sky subtracted, spec2d image in a ginga viewer.
    
    [1;34mpositional arguments:[0m
      [1;32mfile[0m                  Path to a PypeIt spec2d file
    
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
      [1;36m--list[0m                List the extensions only? (default: False)
      [1;36m--det[0m [1;33mDET[0m             Detector name or number. If a number, the name is
                            constructed assuming the reduction is for a single
                            detector. If a string, it must match the name of the
                            detector object (e.g., DET01 for a detector, MSC01 for a
                            mosaic). If not set, the first available detectorin the
                            spec2d file will be shown (default: None)
      [1;36m--spat_id[0m [1;33mSPAT_ID[0m     Restrict plotting to this slit (PypeIt ID notation)
                            (default: None)
      [1;36m--maskID[0m [1;33mMASKID[0m       Restrict plotting to this maskID (default: None)
      [1;36m--showmask[0m [1;33m[SHOWMASK ...][0m
                            Include a channel showing the mask. If no arguments are
                            provided, the mask bit values are provided directly. You
                            can also specify one or more mask flags used to
                            construct an image identifying which pixels are flagged
                            by any of these issues. E.g., to show pixels flagged by
                            the instrument specific bad-pixel mask or cosmic arrays,
                            use --showmask BPM CR . See
                            https://pypeit.readthedocs.io/en/release/out_masks.html
                            for the list of flags. (default: None)
      [1;36m--removetrace[0m         Do not overplot traces in the skysub, sky_resid, and
                            resid channels (default: False)
      [1;36m--embed[0m               Upon completion embed in ipython shell (default: False)
      [1;36m--ignore_extract_mask[0m
                            Ignore the extraction mask (default: False)
      [1;36m--channels[0m [1;33mCHANNELS[0m   Only show a subset of the channels (0-indexed), e.g. 1,3
                            (default: None)
      [1;36m--prefix[0m [1;33mPREFIX[0m       Channel name prefix [lets you display more than one set]
                            (default: )
      [1;36m--no_clear[0m            Do *not* clear all existing tabs (default: True)
      [1;36m--try_old[0m             Attempt to load old datamodel versions. A crash may
                            ensue.. (default: False)
    