.. code-block:: console

    $ pypeit_find_objects -h
    usage: pypeit_find_objects [-h] [--list] [--det DET] [--old] file
    
    Display sky subtracted, spec2d image in theinteractive object finding GUI. Run
    abovethe Science/ folder
    
    positional arguments:
      file        PYPEIT spec2d file
    
    optional arguments:
      -h, --help  show this help message and exit
      --list      List the extensions only? (default: False)
      --det DET   Detector (default: 1)
      --old       Used old slit tracing (default: False)
    