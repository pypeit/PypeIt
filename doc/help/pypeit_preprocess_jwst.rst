.. code-block:: console

    $ pypeit_preprocess_jwst -h
    usage: pypeit_preprocess_jwst [-h] [-v VERBOSITY] [--log_file LOG_FILE]
                                  [--log_level LOG_LEVEL] [-o OUTPUT_DIR]
                                  [--crds_path CRDS_PATH] [--pid PID [PID ...]]
                                  [--obs OBS [OBS ...]] [--vis VIS [VIS ...]]
                                  [--vg VG [VG ...]] [--s S [S ...]] [--a A [A ...]]
                                  [--exp EXP [EXP ...]] [-l]
                                  [--source_type {POINT,EXTENDED}]
                                  [--det1_config DET1_CONFIG]
                                  [--spec2_config SPEC2_CONFIG] [--save_config]
                                  [--overwrite_stage1] [--overwrite_stage2]
                                  raw_dir
    
    Preprocess JWST data through the JWST Detector1Pipeline (Stage 1) and
    Spec2Pipeline (Stage 2) prior to ingestion by PypeIt.
    
    positional arguments:
      raw_dir               Directory containing the raw uncalibrated JWST files
                            (*_uncal.fits)
    
    options:
      -h, --help            show this help message and exit
      -v, --verbosity VERBOSITY
                            Verbosity level, which must be 0, 1, or 2. Level 0
                            includes warning and error messages, level 1 adds
                            informational messages, and level 2 adds debugging
                            messages and the calling sequence.
      --log_file LOG_FILE   Name for the log file. If set to "default", a default
                            name is used. If None, a log file is not produced.
      --log_level LOG_LEVEL
                            Verbosity level for the log file. If a log file is
                            produce and this is None, the file log will match the
                            console stream log.
      -o, --output_dir OUTPUT_DIR
                            Directory for the JWST pipeline output files. Default is
                            jwst_output in the current working directory.
      --crds_path CRDS_PATH
                            Path to the CRDS cache directory. If not provided, the
                            CRDS_PATH environment variable must be set.
      --pid PID [PID ...]   Select only files matching these JWST program ID
                            number(s) (e.g., 01234 01235). Based on the JWST file
                            naming convention: jw<ppppp><ooo><vvv>_... If not
                            provided, files are not filtered by program ID.
      --obs OBS [OBS ...]   Select only files matching these JWST observation
                            number(s) (e.g., 003 004). Based on the JWST file naming
                            convention: jw<ppppp><ooo><vvv>_... If not provided,
                            files are not filtered by observation number.
      --vis VIS [VIS ...]   Select only files matching these JWST visit number(s)
                            (e.g., 001 002). Based on the JWST file naming
                            convention: jw<ppppp><ooo><vvv>_... If not provided,
                            files are not filtered by visit number.
      --vg VG [VG ...]      Select only files matching these JWST visit group
                            number(s) (e.g., 01 02). Based on the JWST file naming
                            convention: jw<ppppp><ooo><vvv>_<gg><s><aa>_... If not
                            provided, files are not filtered by visit group.
      --s S [S ...]         Select only files matching these JWST parallel sequence
                            ID(s) (e.g., 1 2). 1 = prime, 2-5 = parallel. Based on
                            the JWST file naming convention:
                            jw<ppppp><ooo><vvv>_<gg><s><aa>_... If not provided,
                            files are not filtered by parallel sequence.
      --a A [A ...]         Select only files matching these JWST activity number(s)
                            (base-36, e.g., 01 02). Based on the JWST file naming
                            convention: jw<ppppp><ooo><vvv>_<gg><s><aa>_... If not
                            provided, files are not filtered by activity number.
      --exp EXP [EXP ...]   Select only files matching these JWST exposure number(s)
                            (e.g., 00001 00002). Based on the JWST file naming
                            convention: jw<ppppp><ooo><vvv>_<gg><s><aa>_<eeeee>_...
                            If not provided, files are not filtered by exposure
                            number.
      -l, --list            List the files that would be processed and exit. No
                            preprocessing is performed.
      --source_type {POINT,EXTENDED}
                            Source type for the JWST Spec2Pipeline srctype step.
                            Default is EXTENDED.
      --det1_config DET1_CONFIG
                            JSON file with parameters for the JWST Detector1Pipeline
                            (Stage 1).  Parameters provided in this file are merged
                            with the defaults, overriding any matching keys.  If not
                            provided, the built-in defaults are used.  Use
                            --save_config to generate template files with the
                            default parameters.  Example JSON format:
                             
                                {
                                    "clean_flicker_noise": {
                                        "skip": false,
                                        "fit_method": "fft",
                                        "n_sigma": 2
                                    }
                                }
      --spec2_config SPEC2_CONFIG
                            JSON file with parameters for the JWST Spec2Pipeline
                            (Stage 2).  Parameters provided in this file are merged
                            with the defaults, overriding any matching keys.  If not
                            provided, the built-in defaults are used.  Use
                            --save_config to generate template files with the
                            default parameters.  Example JSON format:
                             
                                {
                                    "extract_1d": {"skip": false},
                                    "bkg_subtract": {"skip": false}
                                }
      --save_config         Save the default Detector1Pipeline and Spec2Pipeline
                            parameter dictionaries to JSON files (det1_config.json,
                            spec2_config.json) in the output directory and exit. No
                            preprocessing is performed when this flag is set. These
                            files can then be edited and passed back via
                            --det1_config and --spec2_config. Note that even without
                            this flag, the parameters used for each run are always
                            saved to the output directory.
      --overwrite_stage1    Overwrite existing Stage 1 output files
      --overwrite_stage2    Overwrite existing Stage 2 output files
    