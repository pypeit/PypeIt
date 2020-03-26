

Version: 1.0.0

=============  ========  ==========  ===================================
Obj            Obj Type  Array Type  Description                        
=============  ========  ==========  ===================================
``PYP_SPEC``   str                   PypeIt: Spectrograph name          
``ext_mode``   str                   Extraction mode (options: BOX, OPT)
``flux``       ndarray   floating    Flux/counts array                  
``fluxed``     bool                  Fluxed?                            
``ivar``       ndarray   floating    Inverse variance array             
``mask``       ndarray   integer     Mask array (0=Good???)             
``obj_model``  ndarray   floating    Object model for tellurics         
``telluric``   ndarray   floating    Telluric model                     
``wave``       ndarray   floating    Wavelength array                   
=============  ========  ==========  ===================================
