.. code-block:: console

    $ pypeit_setup_datacube -h
    usage: pypeit_setup_datacube [-h] [-v VERBOSITY] [--log_file LOG_FILE]
                                 [--log_level LOG_LEVEL] [--sensfile SENSFILE]
                                 [--whitelight_range WHITELIGHT_RANGE]
                                 [--manual MANUAL] [--fwhm FWHM]
                                 [--snr_thresh SNR_THRESH]
                                 [--spatial_delta SPATIAL_DELTA] [--det DET]
                                 [-o] [--append]
                                 pypeit_file target

    Prepare .coadd3d and .extract files for point-source datacube work.

    positional arguments:
      pypeit_file           PypeIt reduction file.
      target                Target name, e.g. J0750+6927.

    options:
      -h, --help            show this help message and exit
      -v VERBOSITY, --verbosity VERBOSITY
                            Verbosity level, which must be 0, 1, or 2. Level 0
                            includes warning and error messages, level 1 adds
                            informational messages, and level 2 adds debugging
                            messages and the calling sequence. (default: 2)
      --log_file LOG_FILE   Name for the log file. If set to "default", a default
                            name is used. If None, a log file is not produced.
                            (default: default)
      --log_level LOG_LEVEL
                            Verbosity level for the log file. If a log file is
                            produce and this is None, the file log will match the
                            console stream log. (default: None)
      --sensfile SENSFILE   Optional sensitivity function file. If omitted, the
                            generated .coadd3d file will produce unfluxed cubes.
                            (default: None)
      --whitelight_range WHITELIGHT_RANGE, --wl_range WHITELIGHT_RANGE
                            White-light wavelength range, e.g. 9400,10000.
                            Defaults to None,None. (default: None,None)
      --manual MANUAL       Manual extraction position in x:y coordinates. If
                            provided, the extract file uses
                            opt_prof_method=user_gauss. (default: None)
      --fwhm FWHM           Extraction FWHM in arcsec. (default: 1.1)
      --snr_thresh SNR_THRESH
                            S/N threshold for automatic source finding in the
                            extraction file. (default: 4.0)
      --spatial_delta SPATIAL_DELTA
                            Output cube spatial sampling in arcsec. Defaults to
                            0.678924 for keck_kcrm; otherwise omitted. (default:
                            None)
      --det DET             Detector number. (default: 1)
      -o, --overwrite       Overwrite an existing .extract file. The .coadd3d file
                            is always refreshed. (default: False)
      --append              Append newly reduced spec2d files to an existing
                            .coadd3d file without changing any other lines. The
                            .extract file is left unchanged. (default: False)
