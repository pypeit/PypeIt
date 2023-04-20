.. code-block:: console

    $ pypeit_multislit_flexure -h
    usage: pypeit_multislit_flexure [-h] [--clobber] [--debug] flex_file outroot
    
    Calculate and apply flexure corrections for 1D spectra produced by PypeIt.
    
    positional arguments:
      flex_file   File to guide flexure corrections for this multi-slit mode.  This
                  file must have the following format:
                   
                  flexure read
                    filename
                    spec1dfile1
                    spec1dfile2
                       ...    
                  flexure end
                   
                   
      outroot     Output fileroot for the flexure fits saved as FITS.
    
    options:
      -h, --help  show this help message and exit
      --clobber   Clobber output files
      --debug     show debug plots?
    