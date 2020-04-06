
Version 1.1.0

=================  ========  ==========  =====================================================
HDU Name           Obj Type  Array Type  Description                                          
=================  ========  ==========  =====================================================
``PYP_SPEC``       str                   PypeIt spectrograph name                             
``BPMFLATS``       ndarray   integer     Mirrors SlitTraceSet mask for the Flat-specific flags
``FLAT_MODEL``     ndarray   floating    Model flat                                           
``PIXELFLAT``      ndarray   floating    Pixel normalized flat                                
``PROCFLAT``       ndarray   floating    Processed, combined flats                            
``SPAT_BSPLINES``  ndarray   bspline     B-spline models for Illumination flat                
``SPAT_ID``        ndarray   integer     Slit spat_id                                         
=================  ========  ==========  =====================================================
