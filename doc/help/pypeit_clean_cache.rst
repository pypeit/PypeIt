.. code-block:: console

    $ pypeit_clean_cache -h
    usage: pypeit_clean_cache [-h] [-p PATTERN [PATTERN ...]]
                              [-v VERSION [VERSION ...]] [--remove_all] [-l]
    
    View/Remove fils in the PypeIt data cache
    
    options:
      -h, --help            show this help message and exit
      -p PATTERN [PATTERN ...], --pattern PATTERN [PATTERN ...]
                            Remove any files matching the provided pattern. If
                            combined with --version, this selects only files
                            downloaded from the identified GitHub versoin. If the
                            version is not specified, any file matching the provided
                            pattern(s) are removed. (default: None)
      -v VERSION [VERSION ...], --version VERSION [VERSION ...]
                            Remove files associated one or more provided tags,
                            branches, or commit references on GitHub. These must be
                            an exact match to the relevant GitHub reference. If
                            combined with --pattern, this selects the GitHub
                            reference for the files found. If no files are
                            specified, all files associated with the given reference
                            are removed. Note this is only relevant for the files on
                            GitHub, not s3. For files on s3, do not specify the
                            version. (default: None)
      --remove_all          BEWARE: Removes all data from the pypeit cache. Use of
                            this option ignores the --pattern and --version options.
                            (default: False)
      -l, --list            Only list the contents of the cache. (default: False)
    