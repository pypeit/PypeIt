.. code-block:: console

    $ pypeit_clean_cache -h
    usage: pypeit_clean_cache [-h] [-p PATTERN [PATTERN ...]] [--all] [--clear] [-l]
    
    View/Remove fils in the PypeIt data cache
    
    options:
      -h, --help            show this help message and exit
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
    