.. code-block:: console

    $ pypeit_coadd_1dspec -h
    usage: pypeit_coadd_1dspec [-h] [--debug] [--show] [--par_outfile PAR_OUTFILE]
                               [-v VERBOSITY]
                               coadd1d_file
    
    Coadd 1D spectra produced by PypeIt
    
    positional arguments:
      coadd1d_file          File to guide coadding process. This file must have the
                            following format (see docs for further details including
                            the use of paths):
                             
                            [coadd1d]
                               coaddfile='output_filename.fits' # Optional
                               sensfuncfile = 'sensfunc.fits' # Required only for Echelle
                             
                               coadd1d read
                                    filename | obj_id
                                 spec1dfile1 | objid1
                                 spec1dfile2 | objid2
                                 spec1dfile3 | objid3
                                    ...    
                               coadd1d end
                             
                             OR the coadd1d read/end block can look like
                             
                              coadd1d read
                                    filename | obj_id
                                 spec1dfile1 | objid 
                                 spec1dfile2 | 
                                 spec1dfile3 | 
                                 ...    
                              coadd1d end
                             
                            That is the coadd1d block must be a two column list of
                            spec1dfiles and objids, but you can specify only a
                            single objid for all spec1dfiles on the first line
                             
                            Where:
                             
                            spec1dfile: full path to a PypeIt spec1dfile
                             
                            objid: the object identifier. To determine the objids
                            inspect the spec1d_*.txt files or run pypeit_show_1dspec
                            spec1dfile --list
                             
                            If the coaddfile is not given the output file will be
                            placed in the same directory as the first spec1d file.
                             
    
    optional arguments:
      -h, --help            show this help message and exit
      --debug               show debug plots?
      --show                show QA during coadding process
      --par_outfile PAR_OUTFILE
                            Output to save the parameters
      -v VERBOSITY, --verbosity VERBOSITY
                            Verbosity level between 0 [none] and 2 [all]. Default:
                            1. Level 2 writes a log with filename
                            coadd_1dspec_YYYYMMDD-HHMM.log
    