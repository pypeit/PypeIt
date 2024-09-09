.. code-block:: console

    $ pypeit_chk_flexure -h
    usage: pypeit_chk_flexure [-h] (--spec | --spat) [--try_old]
                              input_file [input_file ...]
    
    Print QA on flexure to the screen
    
    positional arguments:
      input_file  One or more PypeIt spec2d or spec1d file
    
    options:
      -h, --help  show this help message and exit
      --spec      Check the spectral flexure (default: False)
      --spat      Check the spatial flexure (default: False)
      --try_old   Attempt to load old datamodel versions. A crash may ensue..
                  (default: False)
    