.. code-block:: console

    $ pypeit_install_ql_calibs -h
    usage: pypeit_install_ql_calibs [-h] [-v VERBOSITY] [--log_file LOG_FILE]
                                    [--log_level LOG_LEVEL] [--zip ZIP |
                                    --ql_path QL_PATH] [--odir ODIR] [--rmzip]
    
    Script to install PypeIt QL calibration files
    
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
      --zip ZIP             Zip file of the full QL_CALIB directory downloaded from
                            the PypeIt Google Drive (default: None)
      --ql_path QL_PATH     An existing directory to symlink as the QL_CALIB
                            directory. (default: None)
      --odir ODIR           The directory in which to extract the zip file. Ignored
                            if a direct path is provided using --ql_path. (default:
                            current working directory)
      --rmzip               Remove the downloaded zip file (default: False)
    