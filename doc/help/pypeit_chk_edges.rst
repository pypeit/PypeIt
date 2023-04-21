.. code-block:: console

    $ pypeit_chk_edges -h
    usage: pypeit_chk_edges [-h] [--mpl] [--try_old] trace_file
    
    Display trace image and edge traces
    
    positional arguments:
      trace_file  PypeIt Edges file [e.g. Edges_A_0_DET01.fits.gz]
    
    options:
      -h, --help  show this help message and exit
      --mpl       Use a matplotlib window instead of ginga to show the trace
                  (default: False)
      --try_old   Attempt to load old datamodel versions. A crash may ensue..
                  (default: False)
    