.. code-block:: console

    $ pypeit_chk_alignments -h
    usage: pypeit_chk_alignments [-h] [--chname CHNAME] [--try_old] file
    
    Display Alignment image and the trace data
    
    positional arguments:
      file             PypeIt Alignment file [e.g. Alignment_A_1_DET01.fits]
    
    options:
      -h, --help       show this help message and exit
      --chname CHNAME  Channel name for image in Ginga (default: Alignments)
      --try_old        Attempt to load old datamodel versions. A crash may ensue..
                       (default: False)
    