
Version 1.1.0

===========================  ========  ==========  =================================================
HDU Name                     Obj Type  Array Type  Description                                      
===========================  ========  ==========  =================================================
``PYP_SPEC``                 str                   PypeIt spectrograph name                         
``ILLUMFLAT_BPM``            ndarray   integer     Mirrors SlitTraceSet mask for Flat-specific flags
``ILLUMFLAT_RAW``            ndarray   floating    Processed, combined illum flats                  
``ILLUMFLAT_SPAT_BSPLINES``  ndarray   bspline     B-spline models for illum flat                   
``PIXELFLAT_BPM``            ndarray   integer     Mirrors SlitTraceSet mask for Flat-specific flags
``PIXELFLAT_MODEL``          ndarray   floating    Model flat                                       
``PIXELFLAT_NORM``           ndarray   floating    Normalized pixel flat                            
``PIXELFLAT_RAW``            ndarray   floating    Processed, combined pixel flats                  
``PIXELFLAT_SPAT_BSPLINES``  ndarray   bspline     B-spline models for pixel flat                   
``PIXELFLAT_SPEC_ILLUM``     ndarray   floating    Relative spectral illumination                   
``SPAT_ID``                  ndarray   integer     Slit spat_id                                     
===========================  ========  ==========  =================================================
