.. code-block:: console

    $ pypeit_identify -h
    usage: pypeit_identify [-h] [--lamps LAMPS] [-s] [--wmin WMIN] [--wmax WMAX]
                           [--slits SLITS] [-m] [-n] [--det DET] [--rmstol RMSTOL]
                           [--fwhm FWHM] [--sigdetect SIGDETECT] [--pixtol PIXTOL]
                           [--linear] [--force_save] [--rescale_resid]
                           [-v VERBOSITY] [--try_old]
                           arc_file slits_file
    
    Launch PypeIt pypeit_identify tool, display extracted Arc, and load linelist.
    
    positional arguments:
      arc_file              PypeIt Arc file
      slits_file            PypeIt Slits file
    
    options:
      -h, --help            show this help message and exit
      --lamps LAMPS         Comma separated list of calibration lamps (no spaces)
                            (default: None)
      -s, --solution        Load a wavelength solution from the arc_file (if it
                            exists) (default: False)
      --wmin WMIN           Minimum wavelength range (default: 3000.0)
      --wmax WMAX           Maximum wavelength range (default: 50000.0)
      --slits SLITS         Which slit to load for wavelength calibration. Format
                            should be [0,1,...] for multiple slits, 0 for only one
                            slit. If creating a new WaveCalib with the -n flag, this
                            is not necessary. (default: 0)
      -m, --multi           Set this flag to create wavelength solutions for muliple
                            slits (default: False)
      -n, --new_sol         Set this flag to construct a new WaveCalib file, rather
                            than using the exising one (default: False)
      --det DET             Detector index (default: 1)
      --rmstol RMSTOL       RMS tolerance (default: 0.1)
      --fwhm FWHM           FWHM for line finding (default: 4.0)
      --sigdetect SIGDETECT
                            sigma detection for line finding (default: None)
      --pixtol PIXTOL       Pixel tolerance for Auto IDs (default: 0.1)
      --linear              Show the spectrum in linear (rather than log) scale
                            (default: False)
      --force_save          Save the solutions, despite the RMS (default: False)
      --rescale_resid       Rescale the residual plot to include all points?
                            (default: False)
      -v VERBOSITY, --verbosity VERBOSITY
                            Verbosity level between 0 [none] and 2 [all]. Default:
                            1. Level 2 writes a log with filename identify_YYYYMMDD-
                            HHMM.log (default: 1)
      --try_old             Attempt to load old datamodel versions. A crash may
                            ensue.. (default: False)
    