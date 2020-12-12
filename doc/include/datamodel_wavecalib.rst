
Version 1.0.0

===============  =========  ==========  ============================================================
HDU Name         Obj Type   Array Type  Description                                                 
===============  =========  ==========  ============================================================
``PYP_SPEC``     str                    PypeIt spectrograph name                                    
``ARC_SPECTRA``  ndarray    floating    2D array: 1D extracted spectra, slit by slit (nspec, nslits)
``NSLITS``       int                    Total number of slits.  This can include masked slits       
``SPAT_IDS``     ndarray    integer     Slit spat_ids. Named distinctly from that in WaveFit        
``STRPAR``       str                    Parameters as a string                                      
``WV_FIT2D``     PypeItFit              2D wavelength solution (echelle)                            
``WV_FITS``      ndarray    WaveFit     WaveFit to each 1D wavelength solution                      
===============  =========  ==========  ============================================================
