.. code-block:: console

    $ pypeit_collate_1d -h
    usage: pypeit_collate_1d [-h] [--spec1d_files [SPEC1D_FILES ...]]
                             [--par_outfile PAR_OUTFILE] [--tolerance TOLERANCE]
                             [--match MATCH] [--dry_run] [--archive_dir ARCHIVE_DIR]
                             [--exclude_slit_bm [EXCLUDE_SLIT_BM ...]]
                             [--exclude_serendip]
                             [input_file]
    
    Flux/Coadd multiple 1d spectra from multiple nights and prepare a directory for
    the KOA.
    
    positional arguments:
      input_file            (Optional) File for guiding the collate process.
                            Parameters in this file are overidden by the command
                            line. The file must have the following format:
                             
                            [collate1d]
                              tolerance             <tolerance>
                              archive_root          <directory for archive files>
                              exclude_slit_trace_bm <slit types to exclude>
                              exclude_serendip      If set serendipitous objects are skipped.
                              match_using           Whether to match using "pixel" or
                                                    "ra/dec"
                              dry_run               If set the matches are displayed
                                                    without any processing
                             
                            spec1d read
                            <path to spec1d files, wildcards allowed>
                            ...
                            end
    
    optional arguments:
      -h, --help            show this help message and exit
      --spec1d_files [SPEC1D_FILES ...]
                            One or more spec1d files to flux/coadd/archive. Can
                            contain wildcards
      --par_outfile PAR_OUTFILE
                            Output to save the parameters
      --tolerance TOLERANCE
                            The tolerance used when comparing the coordinates of
                            objects. If two objects are within this distance from
                            each other, they are considered the same object. If
                            match_using is 'ra/dec' (the default) this is an angular
                            distance. The defaults units are arcseconds but other
                            units supported by astropy.coordinates.Angle can be
                            used(e.g. '0.003d' or '0h1m30s'). If match_using is
                            'pixel' this is a float.
      --match MATCH         Determines how 1D spectra are matched as being the same
                            object. Must be either 'pixel' or 'ra/dec'.
      --dry_run             If set, the script will display the matching File and
                            Object Ids but will not flux, coadd or archive.
      --archive_dir ARCHIVE_DIR
                            The path where files and metadata will be archived.
      --exclude_slit_bm [EXCLUDE_SLIT_BM ...]
                            A list of slit trace bitmask bits that should be
                            excluded.
      --exclude_serendip    Whether to exclude SERENDIP objects from collating.
    