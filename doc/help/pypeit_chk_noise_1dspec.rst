.. code-block:: console

    $ pypeit_chk_noise_1dspec -h
    [1;34musage: [0m[1;35mpypeit_chk_noise_1dspec[0m [[32m-h[0m] [[32m-v [33mVERBOSITY[0m] [[36m--log_file [33mLOG_FILE[0m]
                                   [[36m--log_level [33mLOG_LEVEL[0m] [[36m--fileformat [33mFILEFORMAT[0m]
                                   [[36m--extraction [33mEXTRACTION[0m] [[36m--ploterr[0m] [[36m--step[0m]
                                   [[36m--z [33m[Z ...][0m] [[36m--maskdef_objname [33mMASKDEF_OBJNAME[0m]
                                   [[36m--pypeit_name [33mPYPEIT_NAME[0m] [[36m--wavemin [33mWAVEMIN[0m]
                                   [[36m--wavemax [33mWAVEMAX[0m] [[36m--plot_or_save [33mPLOT_OR_SAVE[0m]
                                   [[36m--try_old[0m]
                                   [32m[files ...][0m
    
    Examine the noise in a PypeIt spectrum
    
    [1;34mpositional arguments:[0m
      [1;32mfiles[0m                 PypeIt spec1d file(s) (default: None)
    
    [1;34moptions:[0m
      [1;32m-h[0m, [1;36m--help[0m            show this help message and exit
      [1;32m-v[0m, [1;36m--verbosity[0m [1;33mVERBOSITY[0m
                            Verbosity level, which must be 0, 1, or 2. Level 0
                            includes warning and error messages, level 1 adds
                            informational messages, and level 2 adds debugging
                            messages and the calling sequence. (default: 2)
      [1;36m--log_file[0m [1;33mLOG_FILE[0m   Name for the log file. If set to "default", a default
                            name is used. If None, a log file is not produced.
                            (default: None)
      [1;36m--log_level[0m [1;33mLOG_LEVEL[0m
                            Verbosity level for the log file. If a log file is
                            produce and this is None, the file log will match the
                            console stream log. (default: None)
      [1;36m--fileformat[0m [1;33mFILEFORMAT[0m
                            Is this coadd1d or spec1d? (default: spec1d)
      [1;36m--extraction[0m [1;33mEXTRACTION[0m
                            If spec1d, which extraction? opt or box (default: opt)
      [1;36m--ploterr[0m             Plot noise spectrum (default: False)
      [1;36m--step[0m                Use `steps-mid` as linestyle (default: False)
      [1;36m--z[0m [1;33m[Z ...][0m           Object redshift (default: None)
      [1;36m--maskdef_objname[0m [1;33mMASKDEF_OBJNAME[0m
                            MASKDEF_OBJNAME of the target that you want to plot. If
                            maskdef_objname is not provided, nor a pypeit_name, all
                            the 1D spectra in the file(s) will be plotted. (default:
                            None)
      [1;36m--pypeit_name[0m [1;33mPYPEIT_NAME[0m
                            PypeIt name of the target that you want to plot. If
                            pypeit_name is not provided, nor a maskdef_objname, all
                            the 1D spectra in the file(s) will be plotted. (default:
                            None)
      [1;36m--wavemin[0m [1;33mWAVEMIN[0m     Wavelength min. This is for selecting a region of the
                            spectrum to analyze. (default: None)
      [1;36m--wavemax[0m [1;33mWAVEMAX[0m     Wavelength max.This is for selecting a region of the
                            spectrum to analyze. (default: None)
      [1;36m--plot_or_save[0m [1;33mPLOT_OR_SAVE[0m
                            Do you want to save to disk or open a plot in a mpl
                            window. If you choose save, a folder called
                            spec1d*_noisecheck will be created and all the relevant
                            plot will be placed there. (default: plot)
      [1;36m--try_old[0m             Attempt to load old datamodel versions. A crash may
                            ensue.. (default: False)
    