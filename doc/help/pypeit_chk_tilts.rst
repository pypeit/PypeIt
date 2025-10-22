.. code-block:: console

    $ pypeit_chk_tilts -h
    usage: pypeit_chk_tilts [-h] [--mpl] [--show_traces] [--try_old] file
    
    Display Tiltimg image and 2D fitted tilts in Ginga viewer or Matplotlib window.
    Tiltimg file must be in the same directory as Tilts.
    
    positional arguments:
      file           PypeIt Tilts file [e.g. Tilt_A_1_01.fits]
    
    options:
      -h, --help     show this help message and exit
      --mpl          Use a matplotlib window instead of ginga to show the tilts.
                     Faster plotting. (default: False)
      --show_traces  Show the traced tilts. This slows down the plotting (mostly in
                     Ginga). If not set, only the fitted, masked and rejected in the
                     fit tilts are shown. (default: False)
      --try_old      Attempt to load old datamodel versions. A crash may ensue..
                     (default: False)
    