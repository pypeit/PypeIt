.. code-block:: console

    $ pypeit_collate_1d -h
    usage: pypeit_collate_1d [-h] [--spec1d_files [SPEC1D_FILES ...]]
                             [--par_outfile PAR_OUTFILE] [--outdir OUTDIR]
                             [--spec1d_outdir SPEC1D_OUTDIR] [--tolerance TOLERANCE]
                             [--match_using MATCH_USING] [--dry_run] [--ignore_flux]
                             [--flux] [--exclude_slit_bm EXCLUDE_SLIT_BM]
                             [--exclude_serendip] [--wv_rms_thresh WV_RMS_THRESH]
                             [--refframe {observed,heliocentric,barycentric}]
                             [--chk_version] [-v VERBOSITY]
                             [input_file]
    
    Flux/Coadd multiple 1d spectra from multiple nights and prepare a directory for
    the KOA.
    
    positional arguments:
      input_file            (Optional) File for guiding the collate process.
                            Parameters in this file are overidden by the command
                            line. The file must have the following format:
                             
                            [collate1d]
                              tolerance             <tolerance>
                              outdir                <directory to place output files>
                              spec1d_outdir         <directory to place modified spec1ds, if any>
                              exclude_slit_trace_bm <slit types to exclude>
                              exclude_serendip      If set serendipitous objects are skipped.
                              match_using           Whether to match using "pixel" or
                                                    "ra/dec"
                              dry_run               If set the matches are displayed
                                                    without any processing
                              flux                  Flux calibrate using archived sensfuncs.
                              ignore_flux           Ignore any flux calibration information in
                                                    spec1d files.
                              wv_rms_thresh         If set, any objects with a wavelength rms > than the input
                                                    value are skipped, else all wavelength rms values are accepted.
                              refframe              Perform reference frame correction prior to coadding.
                                                    Options are ['observed', 'heliocentric', 'barycentric']. Defaults to None.
                             
                            spec1d read
                            <path to spec1d files, wildcards allowed>
                            ...
                            end
    
    options:
      -h, --help            show this help message and exit
      --spec1d_files [SPEC1D_FILES ...]
                            One or more spec1d files to flux/coadd/archive. Can
                            contain wildcards
      --par_outfile PAR_OUTFILE
                            Output to save the parameters
      --outdir OUTDIR       The path where all coadded output files and report files
                            will be placed. Defaults to the current directory.
      --spec1d_outdir SPEC1D_OUTDIR
                            The path where all modified spec1d files are placed.
                            These are only created if flux calibration or refframe
                            correction are asked for. Defaults to overwriting
                            existing spec1ds.
      --tolerance TOLERANCE
                            The tolerance used when comparing the coordinates of
                            objects. If two objects are within this distance from
                            each other, they are considered the same object. If
                            match_using is 'ra/dec' (the default) this is an angular
                            distance. The defaults units are arcseconds but other
                            units supported by astropy.coordinates.Angle can be used
                            (`e.g.`, '0.003d' or '0h1m30s'). If match_using is
                            'pixel' this is a float.
      --match_using MATCH_USING
                            Determines how 1D spectra are matched as being the same
                            object. Must be either 'pixel' or 'ra/dec'.
      --dry_run             If set, the script will display the matching File and
                            Object Ids but will not flux, coadd or archive.
      --ignore_flux         If set, the script will only coadd non-fluxed spectra
                            even if flux data is present. Otherwise fluxed spectra
                            are coadded if all spec1ds have been fluxed calibrated.
      --flux                If set, the script will flux calibrate using archived
                            sensfuncs before coadding.
      --exclude_slit_bm EXCLUDE_SLIT_BM
                            A list of slit trace bitmask bits that should be
                            excluded. Comma separated.
      --exclude_serendip    Whether to exclude SERENDIP objects from collating.
      --wv_rms_thresh WV_RMS_THRESH
                            If set, any objects with a wavelength RMS > this value
                            are skipped, else all wavelength RMS values are
                            accepted.
      --refframe {observed,heliocentric,barycentric}
                            Perform reference frame correction prior to coadding.
                            Options are: observed, heliocentric, barycentric
      --chk_version         If True enforce strict PypeIt version checking to ensure
                            that all files were created with the current version of
                            PypeIt. If set to False, the code will attempt to read
                            out-of-date files and keep going. Beware (!!) that this
                            can lead to unforeseen bugs that either cause the code
                            to crash or lead to erroneous results. I.e., you really
                            need to know what you are doing if you set this to
                            False!
      -v VERBOSITY, --verbosity VERBOSITY
                            Verbosity level between 0 [none] and 2 [all]. Default:
                            1. Level 2 writes a log with filename
                            collate_1d_YYYYMMDD-HHMM.log
    