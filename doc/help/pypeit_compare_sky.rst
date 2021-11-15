.. code-block:: console

    $ pypeit_compare_sky -h
    usage: pypeit_compare_sky [-h] [--exten EXTEN] [--optimal]
                              [--scale_user SCALE_USER] [--test]
                              file skyfile
    
    Compare the extracted sky spectrum against an archived sky model maintained by
    PypeIt.
    
    positional arguments:
      file                  spec1d Spectral file
      skyfile               Archived PypeIt sky file (e.g. paranal_sky.fits)
    
    optional arguments:
      -h, --help            show this help message and exit
      --exten EXTEN         FITS extension (default: None)
      --optimal             Show Optimal? Default is boxcar (default: False)
      --scale_user SCALE_USER
                            Scale user spectrum by a factor (default: 1.0)
      --test                Load files but do not show plot (default: False)
    