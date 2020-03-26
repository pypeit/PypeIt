
Version 1.0.0

==============  ========  ==========  =============================================================
HDU Name        Obj Type  Array Type  Description                                                  
==============  ========  ==========  =============================================================
``PYP_SPEC``    str                   PypeIt spectrograph name                                     
``COEFFS``      ndarray   floating    2D coefficents for the fit.One set per slit/order (3D array).
``FUNC2D``      str                   Function used for the 2D fit                                 
``NSLIT``       int                   Number of slits                                              
``SLITCEN``     ndarray   floating    Location of the slit center.  2D array (nspec, nslit)        
``SPAT_ORDER``  ndarray   integer     Order for spatial fit                                        
``SPEC_ORDER``  ndarray   integer     Order for spectral fit                                       
``TILTS``       ndarray   floating    Image of the tilts (nspec, nspat)                            
==============  ========  ==========  =============================================================
