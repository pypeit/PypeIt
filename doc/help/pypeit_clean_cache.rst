.. code-block:: console

    $ pypeit_clean_cache -h
    [1;34musage: [0m[1;35mpypeit_clean_cache[0m [[32m-h[0m] [[32m-v [33mVERBOSITY[0m] [[36m--log_file [33mLOG_FILE[0m]
                              [[36m--log_level [33mLOG_LEVEL[0m] [[32m-p [33mPATTERN [PATTERN ...][0m]
                              [[36m--all[0m] [[36m--clear[0m] [[32m-l[0m]
    
    View/Remove fils in the PypeIt data cache
    
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
      [1;32m-p[0m, [1;36m--pattern[0m [1;33mPATTERN [PATTERN ...][0m
                            Remove any files matching the provided pattern. If
                            combined with --version, this selects only files
                            downloaded from the identified GitHub versoin. If the
                            version is not specified, any file matching the provided
                            pattern(s) are removed. (default: None)
      [1;36m--all[0m                 By default, the presence of any of the listed patterns
                            yields a match. This flag requires all patterns to be
                            present for a match. (default: False)
      [1;36m--clear[0m               BEWARE: Removes all data from the pypeit cache. Use of
                            this option ignores the --pattern options. (default:
                            False)
      [1;32m-l[0m, [1;36m--list[0m            Only list the contents of the cache. (default: False)
    