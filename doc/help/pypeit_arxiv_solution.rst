.. code-block:: console

    $ pypeit_arxiv_solution -h
    usage: pypeit_arxiv_solution [-h] [-s SLIT] [-v VERBOSITY] [--try_old]
                                 file binning
    
    Read in a MasterWaveCalib solution and convert it into the format required for
    the PypeIt full template archive
    
    positional arguments:
      file                  MasterWaveCalib file
      binning               Spectral binning
    
    options:
      -h, --help            show this help message and exit
      -s, --slit SLIT       Slit number to use (default: 0)
      -v, --verbosity VERBOSITY
                            Verbosity level between 0 [none] and 2 [all]. Default:
                            1. Level 2 writes a log with filename
                            make_arxiv_solution_YYYYMMDD-HHMM.log (default: 1)
      --try_old             Attempt to load old datamodel versions. A crash may
                            ensue.. (default: False)
    