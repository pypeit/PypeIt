.. code-block:: console

    $ pypeit_compile_wvarxiv -h
    usage: pypeit_compile_wvarxiv [-h] [--append] wvarxiv_folder instrument grating
    
    Read in a set of wxarxiv solutions from Identify and compile them into a single
    fits file to be used with the reidentify method.
    
    positional arguments:
      wvarxiv_folder  Location of the WVarxiv files
      instrument      Name of instrument. e.g. keck_lris_blue, keck_deimos,
                      gemini_gmos_south_ham
      grating         Instrument grating name. E.g. b600, r400, 600_10000.
    
    options:
      -h, --help      show this help message and exit
      --append        Append to an existing file for this instrument. (default:
                      False)
    