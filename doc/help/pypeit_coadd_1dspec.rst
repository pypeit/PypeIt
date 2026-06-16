.. code-block:: console

    $ pypeit_coadd_1dspec -h
    [1;34musage: [0m[1;35mpypeit_coadd_1dspec[0m [[32m-h[0m] [[32m-v [33mVERBOSITY[0m] [[36m--log_file [33mLOG_FILE[0m]
                               [[36m--log_level [33mLOG_LEVEL[0m] [[36m--debug[0m] [[36m--show[0m]
                               [[36m--par_outfile [33mPAR_OUTFILE[0m]
                               [32mcoadd1d_file[0m
    
    Coadd 1D spectra produced by PypeIt
    
    [1;34mpositional arguments:[0m
      [1;32mcoadd1d_file[0m          File to guide coadding process.
                             
                            ------------------------
                                   MultiSlit
                            ------------------------
                             
                            For coadding Multislit spectra the file must have the
                            following format (see docs for further details including
                            the use of paths):
                             
                            [coadd1d]
                               coaddfile='output_filename.fits' # Optional
                             
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
                             
                            ------------------------
                                     Echelle
                            ------------------------
                             
                            For coadding Echelle spectra the file must have the
                            following format (see docs for further details):
                             
                            [coadd1d]
                               coaddfile='output_filename.fits' # Optional
                             
                               coadd1d read
                                    filename | obj_id | sensfile  | setup_id 
                                 spec1dfile1 | objid1 | sensfile1 | setup_id1
                                 spec1dfile2 | objid2 | sensfile2 | setup_id2
                                 spec1dfile3 | objid3 | sensfile3 | setup_id3
                                    ...    
                               coadd1d end
                             
                             OR the coadd1d read/end block can look like
                             
                              coadd1d read
                                    filename | obj_id | sensfile  | setup_id
                                 spec1dfile1 | objid1 | sensfile  | setup_id
                                 spec1dfile2 |        |           |         
                                 spec1dfile3 |        |           |         
                                    ...    
                               coadd1d end
                             
                            That is the coadd1d block is a four column list of
                            spec1dfiles, objids, sensitivity function files, and
                            setup_ids, but you can specify only a single objid,
                            sensfile, and setup_id for all spec1dfiles on the first
                            line
                             
                            Here:
                             
                            spec1dfile: full path to a PypeIt spec1dfile
                             
                            objid: the object identifier (see details above)
                             
                            sensfile: full path to a PypeIt sensitivity function
                            file for the echelle setup in question
                             
                            setup_id: string identifier for the echelle setup in
                            question, i.e. 'VIS', 'NIR', or 'UVB'
                             
                            If the coaddfile is not given the output file will be
                            placed in the same directory as the first spec1d file.
                             
    
    [1;34moptions:[0m
      [1;32m-h[0m, [1;36m--help[0m            show this help message and exit
      [1;32m-v[0m, [1;36m--verbosity[0m [1;33mVERBOSITY[0m
                            Verbosity level, which must be 0, 1, or 2. Level 0
                            includes warning and error messages, level 1 adds
                            informational messages, and level 2 adds debugging
                            messages and the calling sequence.
      [1;36m--log_file[0m [1;33mLOG_FILE[0m   Name for the log file. If set to "default", a default
                            name is used. If None, a log file is not produced.
      [1;36m--log_level[0m [1;33mLOG_LEVEL[0m
                            Verbosity level for the log file. If a log file is
                            produce and this is None, the file log will match the
                            console stream log.
      [1;36m--debug[0m               show debug plots?
      [1;36m--show[0m                show QA during coadding process
      [1;36m--par_outfile[0m [1;33mPAR_OUTFILE[0m
                            Output to save the parameters
    