

Version: 1.0.0

==============  ========  ==========  ========================================================
Obj Key         Obj Type  Array Type  Description                                             
==============  ========  ==========  ========================================================
``PYP_SPEC``    str                   PypeIt: Spectrograph name                               
``ext_mode``    str                   Extraction mode (options: BOX, OPT)                     
``flux``        ndarray   floating    Flux array in units of counts/s or 10^-17 erg/s/cm^2/Ang
``fluxed``      bool                  Boolean indicating if the spectrum is fluxed.           
``ivar``        ndarray   floating    Inverse variance array (matches units of flux)          
``mask``        ndarray   integer     Mask array (1=Good,0=Bad)                               
``obj_model``   ndarray   floating    Object model for tellurics                              
``spect_meta``  dict                  header dict                                             
``telluric``    ndarray   floating    Telluric model                                          
``wave``        ndarray   floating    Wavelength array (Ang)                                  
==============  ========  ==========  ========================================================
