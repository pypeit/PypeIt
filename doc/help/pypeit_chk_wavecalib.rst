.. code-block:: console

    $ pypeit_chk_wavecalib -h
    usage: pypeit_chk_wavecalib [-h] [--try_old] input_file [input_file ...]
    
    Print QA on Wavelength Calib to the screen
    
    positional arguments:
      input_file  One or more PypeIt WaveCalib file [e.g. WaveCalib_A_1_DET01.fits]
                  or spec2d file
    
    options:
      -h, --help  show this help message and exit
      --try_old   Attempt to load old datamodel versions. A crash may ensue..
                  (default: False)
    