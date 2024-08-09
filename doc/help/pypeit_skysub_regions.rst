.. code-block:: console

    $ pypeit_skysub_regions -h
    usage: pypeit_skysub_regions [-h] [--det DET] [-o] [-i] [-f] [-s] [-v VERBOSITY]
                                 [--try_old]
                                 file
    
    Display a spec2d frame and interactively define the sky regions using a GUI. Run
    in the same folder as your .pypeit file
    
    positional arguments:
      file                  spec2d file
    
    options:
      -h, --help            show this help message and exit
      --det DET             Detector (default: 1)
      -o, --overwrite       Overwrite any existing files/directories (default:
                            False)
      -i, --initial         Use initial slit edges? (default: False)
      -f, --flexure         Use flexure corrected slit edges? (default: False)
      -s, --standard        List standard stars as well? (default: False)
      -v VERBOSITY, --verbosity VERBOSITY
                            Verbosity level between 0 [none] and 2 [all]. Default:
                            1. Level 2 writes a log with filename
                            skysub_regions_YYYYMMDD-HHMM.log (default: 1)
      --try_old             Attempt to load old datamodel versions. A crash may
                            ensue.. (default: False)
    