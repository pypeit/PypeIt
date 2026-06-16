.. code-block:: console

    $ pypeit_chk_alignments -h
    usage: pypeit_chk_alignments [-h] [-v VERBOSITY] [--log_file LOG_FILE]
                                 [--log_level LOG_LEVEL] [--chname CHNAME]
                                 [--try_old]
                                 file
    
    Display Alignment image and the trace data
    
    positional arguments:
      file                  PypeIt Alignment file [e.g. Alignment_A_1_DET01.fits]
    
    options:
      -h, --help            show this help message and exit
      -v, --verbosity VERBOSITY
                            Verbosity level, which must be 0, 1, or 2. Level 0
                            includes warning and error messages, level 1 adds
                            informational messages, and level 2 adds debugging
                            messages and the calling sequence. (default: 2)
      --log_file LOG_FILE   Name for the log file. If set to "default", a default
                            name is used. If None, a log file is not produced.
                            (default: None)
      --log_level LOG_LEVEL
                            Verbosity level for the log file. If a log file is
                            produce and this is None, the file log will match the
                            console stream log. (default: None)
      --chname CHNAME       Channel name for image in Ginga (default: Alignments)
      --try_old             Attempt to load old datamodel versions. A crash may
                            ensue.. (default: False)
    