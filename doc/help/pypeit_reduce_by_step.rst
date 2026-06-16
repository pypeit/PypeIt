.. code-block:: console

    $ pypeit_reduce_by_step -h
    usage: pypeit_reduce_by_step [-h] [-v VERBOSITY] [--log_file LOG_FILE]
                                 [--log_level LOG_LEVEL] [--det DET] [--show]
                                 pypeit_file frame {process,findobj,extract}
    
    Run one of the PypeIt reduction steps on a single frame (and detector)
    
    positional arguments:
      pypeit_file           PypeIt reduction file (must have .pypeit extension)
      frame                 Raw science/standard frame to reduce as listed in your
                            PypeIt file, e.g. b28.fits.gz.
      {process,findobj,extract}
                            Reduction step to perform. Must be "process" to perform
                            basic image processing (bias subtraction, field
                            flattening, etc), "findobj" to perform object detection
                            and sky subtraction, or "extract" to extract 1D spectra.
    
    options:
      -h, --help            show this help message and exit
      -v, --verbosity VERBOSITY
                            Verbosity level, which must be 0, 1, or 2. Level 0
                            includes warning and error messages, level 1 adds
                            informational messages, and level 2 adds debugging
                            messages and the calling sequence.
      --log_file LOG_FILE   Name for the log file. If set to "default", a default
                            name is used. If None, a log file is not produced.
      --log_level LOG_LEVEL
                            Verbosity level for the log file. If a log file is
                            produce and this is None, the file log will match the
                            console stream log.
      --det DET             Single detector number or Mosaic tuple. The Mosaic tuple
                            must include the parentheses and be provided as a
                            string, e.g. "(1,2)". Required, but the list of options
                            is provided if nothing is provided.
      --show                Show reduction steps via plots (which will block further
                            execution until clicked on) and outputs to ginga.
                            Requires remote control ginga session via "ginga
                            --modules=RC,SlitWavelength &"
    