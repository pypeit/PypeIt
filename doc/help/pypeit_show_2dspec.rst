.. code-block:: console

    $ pypeit_show_2dspec -h
    usage: pypeit_show_2dspec [-h] [--list] [--det DET] [--showmask]
                              [--removetrace] [--embed] [--ignore_extract_mask]
                              file
    
    Display sky subtracted, spec2d image in a Ginga viewer. Run above the Science/
    folder
    
    positional arguments:
      file                  PYPIT spec2d file
    
    optional arguments:
      -h, --help            show this help message and exit
      --list                List the extensions only? (default: False)
      --det DET             Detector number (default: 1)
      --showmask            Overplot masked pixels (default: False)
      --removetrace         Do not overplot traces in the skysub, sky_resid and
                            resid channels (default: False)
      --embed               Upon completion embed in ipython shell (default:
                            False)
      --ignore_extract_mask
                            Ignore the extraction mask (default: False)
    