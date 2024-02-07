.. code-block:: console

    $ pypeit_flux_calib -h
    usage: pypeit_flux_calib [-h] [--par_outfile] [-v VERBOSITY] [--try_old]
                             flux_file
    
    Flux calibrate 1D spectra produced by PypeIt
    
    positional arguments:
      flux_file             File to guide fluxing process.  This file must have the
                            following format:
                             
                            flux read
                                 filename | sensfile
                              spec1dfile1 | sensfile1
                              spec1dfile2 | 
                                 ...    
                            flux end
                             
                            OR
                             
                            flux read
                                 filename | sensfile
                              spec1dfile1 | sensfile1
                              spec1dfile2 | sensfile2
                              spec1dfile3 | sensfile3
                                 ...    
                            flux end
                             
                            OR
                             
                            [fluxcalib]
                              use_archived_sens = True
                            flux read
                                 filename
                              spec1dfile1
                              spec1dfile2
                              spec1dfile3
                                 ...    
                            flux end
                             
                            That is, you must specify either a sensfile for all
                            spec1dfiles on the first line, specify one sensfile for
                            each spec1dfile, or specify no sensfiles and use an
                            archived one.
                            Archived sensfiles are available for the following
                            spectrographs: keck_deimos
                             
    
    options:
      -h, --help            show this help message and exit
      --par_outfile         Output to save the parameters
      -v VERBOSITY, --verbosity VERBOSITY
                            Verbosity level between 0 [none] and 2 [all]. Default:
                            1. Level 2 writes a log with filename
                            flux_calib_YYYYMMDD-HHMM.log
      --try_old             Attempt to load old datamodel versions. A crash may
                            ensue..
    