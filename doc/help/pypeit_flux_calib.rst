.. code-block:: console

    $ pypeit_flux_calib -h
    usage: pypeit_flux_calib [-h] [-v VERBOSITY] [--log_file LOG_FILE]
                             [--log_level LOG_LEVEL] [--par_outfile] [--try_old]
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
      -v, --verbosity VERBOSITY
                            Verbosity level, which must be 0, 1, or 2. Level 0
                            includes warning and error messages, level 1 adds
                            informational messages, and level 2 adds debugging
                            messages and the calling sequence.
      --log_file LOG_FILE   Name for the log file. If set to "default", a default
                            name is used. If None, a log file is not produced.
      --log_level LOG_LEVEL
                            Verbosity level for the log file. If a log file is
                            produce and this is None, the file log will match the
                            console stream log.
      --par_outfile         Output to save the parameters
      --try_old             Attempt to load old datamodel versions. A crash may
                            ensue..
    