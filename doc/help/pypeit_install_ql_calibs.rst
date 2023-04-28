.. code-block:: console

    $ pypeit_install_ql_calibs -h
    usage: pypeit_install_ql_calibs [-h] [--zip ZIP | --ql_path QL_PATH]
                                    [--odir ODIR] [--rmzip]
    
    Script to install PypeIt QL calibration files
    
    options:
      -h, --help         show this help message and exit
      --zip ZIP          Zip file of the full QL_CALIB directory downloaded from the
                         PypeIt Google Drive (default: None)
      --ql_path QL_PATH  An existing directory to symlink as the QL_CALIB directory.
                         (default: None)
      --odir ODIR        The directory in which to extract the zip file. Ignored if
                         a direct path is provided using --ql_path. (default:
                         current working directory)
      --rmzip            Remove the downloaded zip file (default: False)
    