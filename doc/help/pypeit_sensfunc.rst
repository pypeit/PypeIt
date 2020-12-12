.. code-block:: console

    $ pypeit_sensfunc -h
    usage: pypeit_sensfunc [-h] [--algorithm {UVIS,IR}] [--multi MULTI]
                           [-o OUTFILE] [-s SENS_FILE] [--debug] [--debug_init]
                           [--par_outfile PAR_OUTFILE]
                           spec1dfile
    
    Compute a sensitivity function
    
    positional arguments:
      spec1dfile            spec1d file for the standard that will be used to
                            compute sensitivity function
    
    optional arguments:
      -h, --help            show this help message and exit
      --algorithm {UVIS,IR}
                            Override the default algorithm for computing the sensitivity function. 
                            Note that it is not possible to set --algorithm and simultaneously use a .sens file with
                            the --sens_file option. If you are using a .sens file set the algorithm there via:
                            
                                [sensfunc]
                                     algorithm = IR
                            
                            The algorithm options are:
                            
                                UVIS = Should be used for data with lambda < 7000A.
                                No detailed model of telluric absorption but corrects for atmospheric extinction.
                            
                                IR   = Should be used for data with lambbda > 7000A.
                                Peforms joint fit for sensitivity function and telluric absorption using HITRAN models.
                            
      --multi MULTI         List of detector numbers to splice together for instruments with multiple detectors
                            arranged in the spectral direction, e.g. --multi = '3,7'
                            Note that it is not possible to set --multi and 
                            simultaneously use a .sens file with the --sens_file option.
                            If you are using a .sens file set the multi_spec_det param there via:
                            
                                     [sensfunc]
                                          multi_spec_det = 3,7
                            
      -o OUTFILE, --outfile OUTFILE
                            Ouput file for sensitivity function. If not specified,
                            the sensitivity function will be written out to a
                            standard filename in the current working directory,
                            i.e. if the standard spec1d file is named
                            spec1d_b24-Feige66_KASTb_foo.fits the sensfunc will be
                            written to sens_b24-Feige66_KASTb_foo.fits
      -s SENS_FILE, --sens_file SENS_FILE
                            Configuration file to change default sensivity
                            function parameters
      --debug               show debug plots?
      --debug_init          debug the initilization of the sensfunc + telluric fit
                            for the IR algorithm
      --par_outfile PAR_OUTFILE
                            Name of outut file to save the parameters used by the
                            fit
    