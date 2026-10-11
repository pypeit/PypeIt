.. code-block:: console

    $ pypeit_osiris_rectify -h
    usage: pypeit_osiris_rectify [-h] [-v VERBOSITY] [--log_file LOG_FILE]
                                 [--log_level LOG_LEVEL] [--bkg BKG] -m MATRIX
                                 [--method {drp,lstsq}] [--niter NITER]
                                 [--relax RELAX] [--use_var] [--no_refrows]
                                 [--no_bpm] [-o OUTFILE] [--overwrite]
                                 raw_file
    
    Rectify a raw Keck/OSIRIS frame (optionally minus a background frame) with a
    Keck rectification matrix, producing one row per lenslet slice. Output: a FITS
    file with the rectified image (ADU), its inverse variance and mask (1 = bad),
    with shape (2048 detector columns, nslice) in the raw column order, and the raw
    header.
    
    positional arguments:
      raw_file              Raw OSIRIS frame (H2RG era, >= 2016)
    
    options:
      -h, --help            show this help message and exit
      -v, --verbosity VERBOSITY
                            Verbosity level, which must be 0, 1, or 2. Level 0
                            includes warning and error messages, level 1 adds
                            informational messages, and level 2 adds debugging
                            messages and the calling sequence. (default: 2)
      --log_file LOG_FILE   Name for the log file. If set to "default", a default
                            name is used. If None, a log file is not produced.
                            (default: None)
      --log_level LOG_LEVEL
                            Verbosity level for the log file. If a log file is
                            produce and this is None, the file log will match the
                            console stream log. (default: None)
      --bkg BKG             Raw background frame subtracted from raw_file before
                            rectification (A-B) (default: None)
      -m, --matrix MATRIX   Rectification matrix file, or a directory of matrices
                            from which the matrix is selected by filter, scale and
                            date
      --method {drp,lstsq}  Solver; default from the osiris_rect_method parameter
                            (drp) (default: None)
      --niter NITER         Iterations of the drp solver; default from
                            osiris_rect_niter (40) (default: None)
      --relax RELAX         Relaxation of the drp solver; default from
                            osiris_rect_relax (1.0) (default: None)
      --use_var             Use the VAR extension(s) for the pixel noise instead of
                            the PypeIt noise model (gain, read noise and shot noise)
                            (default: False)
      --no_refrows          Do not subtract the per-column median of the
                            unilluminated rows 2032-2044 (DRP "Adjust Channel
                            Levels" for the H2RG) (default: False)
      --no_bpm              Do not mask the Keck hot/dead pixels and the DQ != 9
                            pixels (default: False)
      -o, --outfile OUTFILE
                            Output file; default <raw_file stem>_rect.fits in the
                            current directory (default: None)
      --overwrite           Overwrite an existing output file (default: False)
    