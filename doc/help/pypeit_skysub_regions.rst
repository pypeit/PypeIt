.. code-block:: console

    $ pypeit_skysub_regions -h
    usage: pypeit_skysub_regions [-h] [--det DET] [-o] [-i] [-f] [-s] file
    
    Display a Raw science image and interactively definethe sky regions using a
    GUI. Run in the same folderas your .pypeit file
    
    positional arguments:
      file             PypeIt file
    
    optional arguments:
      -h, --help       show this help message and exit
      --det DET        Detector (default: 1)
      -o, --overwrite  Overwrite any existing files/directories (default: False)
      -i, --initial    Use initial slit edges? (default: False)
      -f, --flexure    Use flexure corrected slit edges? (default: False)
      -s, --standard   List standard stars as well? (default: False)
    