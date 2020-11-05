.. code-block:: console

    $ pypeit_coadd_1dspec -h
    usage: pypeit_coadd_1dspec [-h] [--debug] [--show] [--par_outfile]
                               [--test_spec_path TEST_SPEC_PATH]
                               coadd1d_file
    
    Parse
    
    positional arguments:
      coadd1d_file          File to guide coadding process. This file must have the following format: 
                            
                            [coadd1d]
                               coaddfile='output_filename.fits'
                               sensfuncfile = 'sensfunc.fits' # Required only for Echelle
                            
                               coadd1d read
                                 spec1dfile1 objid1
                                 spec1dfile2 objid2
                                 spec1dfile3 objid3
                                    ...    
                               coadd1d end
                            
                                     OR the coadd1d read/end block can look like 
                            
                              coadd1d read
                                 spec1dfile1 objid 
                                 spec1dfile2 
                                 spec1dfile3 
                                 ...    
                              coadd1d end
                            
                            That is the coadd1d block must either be a two column list of spec1dfiles and objids,
                            or you can specify a single objid for all spec1dfiles on the first line
                            
                            Where: 
                            
                               spec1dfile -- full path to a PypeIt spec1dfile
                               objid      -- is the object identifier. To determine the objids inspect the 
                                             spec1d_*.txt files or run pypeit_show_1dspec spec1dfile --list
                            
    
    optional arguments:
      -h, --help            show this help message and exit
      --debug               show debug plots?
      --show                show QA during coadding process
      --par_outfile         Output to save the parameters
      --test_spec_path TEST_SPEC_PATH
                            Path for testing
    