.. code-block:: console

    $ pypeit_sensfunc -h
    usage: pypeit_sensfunc [-h] [--extr {OPT,BOX}] [--algorithm {UVIS,IR}]
                           [--multi MULTI] [-o OUTFILE] [-s SENS_FILE] [-f]
                           [--debug] [--par_outfile PAR_OUTFILE] [-v VERBOSITY]
                           spec1dfile
    
    Compute a sensitivity function
    
    positional arguments:
      spec1dfile            spec1d file for the standard that will be used to
                            compute the sensitivity function
    
    options:
      -h, --help            show this help message and exit
      --extr {OPT,BOX}      Override the default extraction method used for
                            computing the sensitivity function.  Note that it is not
                            possible to set --extr and simultaneously use a .sens
                            file with the --sens_file option. If you are using a
                            .sens file, set the algorithm there via:
                             
                                [sensfunc]
                                     extr = BOX
                             
                            The extraction options are: OPT or BOX
      --algorithm {UVIS,IR}
                            Override the default algorithm for computing the
                            sensitivity function.  Note that it is not possible to
                            set --algorithm and simultaneously use a .sens file with
                            the --sens_file option. If you are using a .sens file,
                            set the algorithm there via:
                             
                                [sensfunc]
                                     algorithm = IR
                             
                            The algorithm options are:
                             
                            UVIS = Should be used for data with lambda < 7000A.  No
                            detailed model of telluric absorption but corrects for
                            atmospheric extinction.
                             
                            IR = Should be used for data with lambbda > 7000A.
                            Performs joint fit for sensitivity function and telluric
                            absorption using HITRAN models.
                             
      --multi MULTI         List of detector numbers to splice together for
                            instruments with multiple detectors arranged in the
                            spectral direction, e.g. --multi = '3,7'.  Note that it
                            is not possible to set --multi and simultaneously use a
                            .sens file with the --sens_file option.  If you are
                            using a .sens file, set the multi_spec_det param there
                            via:
                             
                                [sensfunc]
                                    multi_spec_det = 3,7
                             
      -o OUTFILE, --outfile OUTFILE
                            Output file for sensitivity function. If not specified,
                            the sensitivity function will be written out to a
                            standard filename in the current working directory, i.e.
                            if the standard spec1d file is named
                            spec1d_b24-Feige66_KASTb_foo.fits the sensfunc will be
                            written to sens_b24-Feige66_KASTb_foo.fits. A QA file
                            will also be written as
                            sens_spec1d_b24-Feige66_KASTb_foo_QA.pdf and a file
                            showing throughput plots to
                            sens_spec1d_b24-Feige66_KASTb_foo_throughput.pdf. The
                            same extensions for QA and throughput will be used if
                            outfile is provided but with .fits trimmed off if it is
                            in the filename.
      -s SENS_FILE, --sens_file SENS_FILE
                            Configuration file with sensitivity function parameters
      -f, --use_flat        Use the extracted spectrum of the flatfield calibration
                            to estimate the blaze function when generating the
                            sensitivity function. This is helpful to account for
                            small scale undulations in the sensitivity function. The
                            spec1dfile must contain the extracted flatfield response
                            in order to use this option. This spectrum is extracted
                            by default, unless you did not compute a pixelflat
                            frame. Note that it is not possible to set --use_flat
                            and simultaneously use a .sens file with the --sens_file
                            option. If you are using a .sens file, set the use_flat
                            flag with the argument:
                             
                                [sensfunc]
                                     use_flat = True
      --debug               show debug plots?
      --par_outfile PAR_OUTFILE
                            Name of output file to save the parameters used by the
                            fit
      -v VERBOSITY, --verbosity VERBOSITY
                            Verbosity level between 0 [none] and 2 [all]. Default:
                            1. Level 2 writes a log with filename sensfunc_YYYYMMDD-
                            HHMM.log
    