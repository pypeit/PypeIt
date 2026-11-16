.. code-block:: console

    $ pypeit_chk_noise_1dspec -h
    usage: pypeit_chk_noise_1dspec [-h] [-v VERBOSITY] [--log_file LOG_FILE]
                                   [--log_level LOG_LEVEL] [--fileformat FILEFORMAT]
                                   [--extraction EXTRACTION] [--ploterr] [--step]
                                   [--z [Z ...]] [--maskdef_objname MASKDEF_OBJNAME]
                                   [--pypeit_name PYPEIT_NAME] [--wavemin WAVEMIN]
                                   [--wavemax WAVEMAX] [--plot_or_save PLOT_OR_SAVE]
                                   [--try_old]
                                   [files ...]
    
    Examine the noise in a PypeIt spectrum
    
    positional arguments:
      files                 PypeIt spec1d file(s) (default: None)
    
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
      --fileformat FILEFORMAT
                            Is this coadd1d or spec1d? (default: spec1d)
      --extraction EXTRACTION
                            If spec1d, which extraction? opt or box (default: opt)
      --ploterr             Plot noise spectrum (default: False)
      --step                Use `steps-mid` as linestyle (default: False)
      --z [Z ...]           Object redshift (default: None)
      --maskdef_objname MASKDEF_OBJNAME
                            MASKDEF_OBJNAME of the target that you want to plot. If
                            maskdef_objname is not provided, nor a pypeit_name, all
                            the 1D spectra in the file(s) will be plotted. (default:
                            None)
      --pypeit_name PYPEIT_NAME
                            PypeIt name of the target that you want to plot. If
                            pypeit_name is not provided, nor a maskdef_objname, all
                            the 1D spectra in the file(s) will be plotted. (default:
                            None)
      --wavemin WAVEMIN     Wavelength min. This is for selecting a region of the
                            spectrum to analyze. (default: None)
      --wavemax WAVEMAX     Wavelength max.This is for selecting a region of the
                            spectrum to analyze. (default: None)
      --plot_or_save PLOT_OR_SAVE
                            Do you want to save to disk or open a plot in a mpl
                            window. If you choose save, a folder called
                            spec1d*_noisecheck will be created and all the relevant
                            plot will be placed there. (default: plot)
      --try_old             Attempt to load old datamodel versions. A crash may
                            ensue.. (default: False)
    