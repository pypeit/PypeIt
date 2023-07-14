.. code-block:: console

    $ pypeit_install_telluric -h
    usage: pypeit_install_telluric [-h] [--force_update] [--local_file]
                                   files [files ...]
    
    Script to download/install PypeIt telluric files
    
    positional arguments:
      files           Filename(s) of the TelFits files to be downloaded from the
                      Cloud and installed in the PypeIt cache
    
    options:
      -h, --help      show this help message and exit
      --force_update  Force download of latest version of the telluric grid
                      (default: False)
      --local_file    This is a local file (downloaded or created) to be installed
                      in the cache (default: False)
    