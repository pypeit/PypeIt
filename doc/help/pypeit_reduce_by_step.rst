.. code-block:: console

    $ pypeit_reduce_by_step -h
    usage: pypeit_reduce_by_step [-h] [--det DET] [--show] pypeit_file frame step
    
    Run one of the PypeIt reduction steps on a single frame (and detector)
    
    positional arguments:
      pypeit_file  PypeIt reduction file (must have .pypeit extension)
      frame        Raw science/standard frame to reduce as listed in your PypeIt
                   file, e.g. b28.fits.gz.
      step         Reduction step to perform. Must be "process" to perform basic
                   image processing (bias subtraction, field flattening, etc),
                   "findobj" to perform object detection and initial sky
                   subtraction, or "extract" to extract 1D spectra.
    
    options:
      -h, --help   show this help message and exit
      --det DET    Detector number or Mosaic tuple. Required, but the list of
                   options is provided if nothing is provided.
      --show       Show reduction steps via plots (which will block further
                   execution until clicked on) and outputs to ginga. Requires remote
                   control ginga session via "ginga --modules=RC,SlitWavelength &"
    