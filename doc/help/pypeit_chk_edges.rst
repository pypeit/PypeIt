.. code-block:: console

    $ pypeit_chk_edges -h
    usage: pypeit_chk_edges [-h] [--chname CHNAME] [--mpl] [--try_old] trace_file
    
    Display MasterEdges image and trace data
    
    positional arguments:
      trace_file       PypeIt Master Trace file [e.g. MasterEdges_A_01_aa.fits.gz]
    
    optional arguments:
      -h, --help       show this help message and exit
      --chname CHNAME  Channel name for image in Ginga (default: MTrace)
      --mpl            Use a matplotlib window instead of ginga to show the trace
                       (default: False)
      --try_old        Attempt to load old datamodel versions. A crash may ensue..
                       (default: False)
    