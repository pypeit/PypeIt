
Version 1.1.2

===========================  ========  ==========  ====================================================================================
HDU Name                     Obj Type  Array Type  Description                                                                         
===========================  ========  ==========  ====================================================================================
``PYP_SPEC``                 str                   PypeIt spectrograph name                                                            
``ILLUMFLAT_BPM``            ndarray   integer     Mirrors SlitTraceSet mask for flat-specific flags                                   
``ILLUMFLAT_FINECORR``       ndarray   PypeItFit   PypeIt 2D polynomial fits to the fine correction of the spatial illumination profile
``ILLUMFLAT_RAW``            ndarray   floating    Processed, combined illum flats                                                     
``ILLUMFLAT_SPAT_BSPLINES``  ndarray   bspline     B-spline models for illum flat                                                      
``PIXELFLAT_BPM``            ndarray   integer     Mirrors SlitTraceSet mask for flat-specific flags                                   
``PIXELFLAT_FINECORR``       ndarray   PypeItFit   PypeIt 2D polynomial fits to the fine correction of the spatial illumination profile
``PIXELFLAT_MODEL``          ndarray   floating    Model flat                                                                          
``PIXELFLAT_NORM``           ndarray   floating    Normalized pixel flat                                                               
``PIXELFLAT_RAW``            ndarray   floating    Processed, combined pixel flats                                                     
``PIXELFLAT_SPAT_BSPLINES``  ndarray   bspline     B-spline models for pixel flat                                                      
``PIXELFLAT_SPEC_ILLUM``     ndarray   floating    Relative spectral illumination                                                      
``PIXELFLAT_WAVEIMG``        ndarray   floating    Waveimage for pixel flat                                                            
``SPAT_ID``                  ndarray   integer     Slit spat_id                                                                        
===========================  ========  ==========  ====================================================================================
