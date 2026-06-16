.. code-block:: console

    $ pypeit_reduce_by_step -h
    [1;34musage: [0m[1;35mpypeit_reduce_by_step[0m [[32m-h[0m] [[32m-v [33mVERBOSITY[0m] [[36m--log_file [33mLOG_FILE[0m]
                                 [[36m--log_level [33mLOG_LEVEL[0m] [[36m--det [33mDET[0m] [[36m--show[0m]
                                 [32mpypeit_file[0m [32mframe[0m [32m{process,findobj,extract}[0m
    
    Run one of the PypeIt reduction steps on a single frame (and detector)
    
    [1;34mpositional arguments:[0m
      [1;32mpypeit_file[0m           PypeIt reduction file (must have .pypeit extension)
      [1;32mframe[0m                 Raw science/standard frame to reduce as listed in your
                            PypeIt file, e.g. b28.fits.gz.
      [1;32m{process,findobj,extract}[0m
                            Reduction step to perform. Must be "process" to perform
                            basic image processing (bias subtraction, field
                            flattening, etc), "findobj" to perform object detection
                            and sky subtraction, or "extract" to extract 1D spectra.
    
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
      [1;36m--det[0m [1;33mDET[0m             Single detector number or Mosaic tuple. The Mosaic tuple
                            must include the parentheses and be provided as a
                            string, e.g. "(1,2)". Required, but the list of options
                            is provided if nothing is provided.
      [1;36m--show[0m                Show reduction steps via plots (which will block further
                            execution until clicked on) and outputs to ginga.
                            Requires remote control ginga session via "ginga
                            --modules=RC,SlitWavelength &"
    