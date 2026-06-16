.. code-block:: console

    $ pypeit_identify -h
    [1;34musage: [0m[1;35mpypeit_identify[0m [[32m-h[0m] [[32m-v [33mVERBOSITY[0m] [[36m--log_file [33mLOG_FILE[0m]
                           [[36m--log_level [33mLOG_LEVEL[0m] [[36m--lamps [33mLAMPS[0m] [[32m-s[0m]
                           [[36m--wmin [33mWMIN[0m] [[36m--wmax [33mWMAX[0m] [[36m--slits [33mSLITS[0m] [[32m-m[0m] [[32m-n[0m]
                           [[36m--det [33mDET[0m] [[36m--rmstol [33mRMSTOL[0m] [[36m--fwhm [33mFWHM[0m]
                           [[36m--sigdetect [33mSIGDETECT[0m] [[36m--pixtol [33mPIXTOL[0m] [[36m--linear[0m]
                           [[36m--force_save[0m] [[36m--rescale_resid[0m] [[36m--try_old[0m]
                           [32marc_file[0m [32mslits_file[0m
    
    Launch PypeIt pypeit_identify tool, display extracted Arc, and load linelist.
    
    [1;34mpositional arguments:[0m
      [1;32marc_file[0m              PypeIt Arc file
      [1;32mslits_file[0m            PypeIt Slits file
    
    [1;34moptions:[0m
      [1;32m-h[0m, [1;36m--help[0m            show this help message and exit
      [1;32m-v[0m, [1;36m--verbosity[0m [1;33mVERBOSITY[0m
                            Verbosity level, which must be 0, 1, or 2. Level 0
                            includes warning and error messages, level 1 adds
                            informational messages, and level 2 adds debugging
                            messages and the calling sequence. (default: 2)
      [1;36m--log_file[0m [1;33mLOG_FILE[0m   Name for the log file. If set to "default", a default
                            name is used. If None, a log file is not produced.
                            (default: default)
      [1;36m--log_level[0m [1;33mLOG_LEVEL[0m
                            Verbosity level for the log file. If a log file is
                            produce and this is None, the file log will match the
                            console stream log. (default: None)
      [1;36m--lamps[0m [1;33mLAMPS[0m         Comma separated list of calibration lamps (no spaces)
                            (default: None)
      [1;32m-s[0m, [1;36m--solution[0m        Load a wavelength solution from the arc_file (if it
                            exists) (default: False)
      [1;36m--wmin[0m [1;33mWMIN[0m           Minimum wavelength range (default: 3000.0)
      [1;36m--wmax[0m [1;33mWMAX[0m           Maximum wavelength range (default: 50000.0)
      [1;36m--slits[0m [1;33mSLITS[0m         Which slit to load for wavelength calibration. Format
                            should be [0,1,...] for multiple slits, 0 for only one
                            slit. If creating a new WaveCalib with the -n flag, this
                            is not necessary. (default: 0)
      [1;32m-m[0m, [1;36m--multi[0m           Set this flag to create wavelength solutions for muliple
                            slits (default: False)
      [1;32m-n[0m, [1;36m--new_sol[0m         Set this flag to construct a new WaveCalib file, rather
                            than using the exising one (default: False)
      [1;36m--det[0m [1;33mDET[0m             Detector index (default: 1)
      [1;36m--rmstol[0m [1;33mRMSTOL[0m       RMS tolerance (default: 0.1)
      [1;36m--fwhm[0m [1;33mFWHM[0m           FWHM for line finding (default: 4.0)
      [1;36m--sigdetect[0m [1;33mSIGDETECT[0m
                            sigma detection for line finding (default: None)
      [1;36m--pixtol[0m [1;33mPIXTOL[0m       Pixel tolerance for Auto IDs (default: 0.1)
      [1;36m--linear[0m              Show the spectrum in linear (rather than log) scale
                            (default: False)
      [1;36m--force_save[0m          Save the solutions, despite the RMS (default: False)
      [1;36m--rescale_resid[0m       Rescale the residual plot to include all points?
                            (default: False)
      [1;36m--try_old[0m             Attempt to load old datamodel versions. A crash may
                            ensue.. (default: False)
    