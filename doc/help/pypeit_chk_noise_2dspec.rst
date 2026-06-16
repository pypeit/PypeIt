.. code-block:: console

    $ pypeit_chk_noise_2dspec -h
    [1;34musage: [0m[1;35mpypeit_chk_noise_2dspec[0m [[32m-h[0m] [[32m-v [33mVERBOSITY[0m] [[36m--log_file [33mLOG_FILE[0m]
                                   [[36m--log_level [33mLOG_LEVEL[0m] [[36m--det [33mDET[0m] [[36m--z [33m[Z ...][0m]
                                   [[36m--maskdef_id [33mMASKDEF_ID[0m] [[36m--pypeit_id [33mPYPEIT_ID[0m]
                                   [[36m--pad [33mPAD[0m] [[36m--aspect_ratio [33mASPECT_RATIO[0m]
                                   [[36m--wavemin [33mWAVEMIN[0m] [[36m--wavemax [33mWAVEMAX[0m]
                                   [[36m--mode [33mMODE[0m] [[36m--list[0m] [[36m--try_old[0m]
                                   [32m[files ...][0m
    
    Examine the noise in a PypeIt slit/order
    
    [1;34mpositional arguments:[0m
      [1;32mfiles[0m                 PypeIt spec2d file(s) (default: None)
    
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
      [1;36m--det[0m [1;33mDET[0m             Detector name or number. If a number, the name is
                            constructed assuming the reduction is for a single
                            detector. If a string, it must match the name of the
                            detector object (e.g., DET01 for a detector, MSC01 for a
                            mosaic). (default: 1)
      [1;36m--z[0m [1;33m[Z ...][0m           Object redshift (default: None)
      [1;36m--maskdef_id[0m [1;33mMASKDEF_ID[0m
                            MASKDEF_ID of the slit that you want to plot. If
                            maskdef_id is not provided, nor a pypeit_id, all the 2D
                            spectra in the file(s) will be plotted. (default: None)
      [1;36m--pypeit_id[0m [1;33mPYPEIT_ID[0m
                            PypeIt ID of the slit that you want to plot. If
                            pypeit_id is not provided, nor a maskdef_id, all the 2D
                            spectra in the file(s) will be plotted. (default: None)
      [1;36m--pad[0m [1;33mPAD[0m             Padding for the selected slit. Negative value will trim.
                            (default: -5)
      [1;36m--aspect_ratio[0m [1;33mASPECT_RATIO[0m
                            Aspect ratio when plotting the spec2d (default: 3)
      [1;36m--wavemin[0m [1;33mWAVEMIN[0m     Wavelength min. This is for selecting a region of the
                            spectrum to analyze. (default: None)
      [1;36m--wavemax[0m [1;33mWAVEMAX[0m     Wavelength max. This is for selecting a region of the
                            spectrum to analyze. (default: None)
      [1;36m--mode[0m [1;33mMODE[0m           Options are: [plot, save, print]. "plot" will open a
                            plot in a mpl window. "save" will create a folder called
                            spec2d*_noisecheck where all the relevant plots will be
                            placed. "print" will cause the check noise values to be
                            printed in the terminal. (default: plot)
      [1;36m--list[0m                List the extensions only? (default: False)
      [1;36m--try_old[0m             Attempt to load old datamodel versions. A crash may
                            ensue.. (default: False)
    