
Version 1.1.0

================  ========  ==========  =============================================================
HDU Name          Obj Type  Array Type  Description                                                  
================  ========  ==========  =============================================================
``PYP_SPEC``      str                   PypeIt spectrograph name                                     
``COEFFS``        ndarray   floating    2D coefficents for the fit.One set per slit/order (3D array).
``FUNC2D``        str                   Function used for the 2D fit                                 
``NSLIT``         int                   Total number of slits.  This can include masked slits        
``SPAT_FLEXURE``  float                 Flexure shift from the input TiltImage                       
``SPAT_ID``       ndarray   integer     Slit spat_id                                                 
``SPAT_ORDER``    ndarray   integer     Order for spatial fit (nslit)                                
``SPEC_ORDER``    ndarray   integer     Order for spectral fit (nslit)                               
================  ========  ==========  =============================================================
