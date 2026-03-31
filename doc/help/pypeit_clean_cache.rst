.. code-block:: console

    $ pypeit_clean_cache -h
    usage: pypeit_clean_cache [-h] [-v VERBOSITY] [--log_file LOG_FILE]
                              [--log_level LOG_LEVEL] [-p PATTERN [PATTERN ...]]
                              [--all] [--clear] [-l]
    
    View/Remove fils in the PypeIt data cache
    
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
      -p, --pattern PATTERN [PATTERN ...]
                            Remove any files matching the provided pattern. If
                            combined with --version, this selects only files
                            downloaded from the identified GitHub versoin. If the
                            version is not specified, any file matching the provided
                            pattern(s) are removed. (default: None)
      --all                 By default, the presence of any of the listed patterns
                            yields a match. This flag requires all patterns to be
                            present for a match. (default: False)
      --clear               BEWARE: Removes all data from the pypeit cache. Use of
                            this option ignores the --pattern options. (default:
                            False)
      -l, --list            Only list the contents of the cache. (default: False)
    