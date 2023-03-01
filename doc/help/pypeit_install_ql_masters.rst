.. code-block:: console

    $ pypeit_install_ql_masters -h
    usage: pypeit_install_ql_masters [-h] [--zip ZIP | --ql_path QL_PATH]
                                     [--odir ODIR] [--rmzip]
    
    Script to install PypeIt QL Master files
    
    optional arguments:
      -h, --help         show this help message and exit
      --zip ZIP          Zip file of the full QL_MASTERS directory downloaded from
                         the PypeIt Google Drive (default: None)
      --ql_path QL_PATH  An existing directory to symlink as the QL_MASTERS
                         directory. (default: None)
      --odir ODIR        The directory in which to extract the zip file. Ignored if
                         a direct path is provided using --ql_path. (default:
                         /Users/westfall/Work/packages/pypeit/doc)
      --rmzip            Remove the downloaded zip file (default: False)
    