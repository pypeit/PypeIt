.. code-block:: console

    $ pypeit_flux_calib -h
    usage: pypeit_flux_calib [-h] [--debug] [--par_outfile] flux_file
    
    Flux calibrate spectra
    
    positional arguments:
      flux_file      File to guide fluxing process.
                     This file must have the following format: 
                     
                     flux read
                       spec1dfile1 sensfile
                       spec1dfile2
                          ...    
                          ...    
                     flux end
                     
                         OR   
                     
                     flux read
                       spec1dfile1 sensfile1
                       spec1dfile2 sensfile2
                       spec1dfile3 sensfile3
                          ...    
                     flux end
                     
                     That is, you must specify either a sensfile for all spec1dfiles on the first line, or 
                     create a two column list of spec1dfiles and corresponding sensfiles
                     
    
    optional arguments:
      -h, --help     show this help message and exit
      --debug        show debug plots?
      --par_outfile  Output to save the parameters
    