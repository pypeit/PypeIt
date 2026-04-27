

Version: 1.0.2

=================  ========  ==========  ==========================================================================================================================================
Obj Key            Obj Type  Array Type  Description                                                                                                                               
=================  ========  ==========  ==========================================================================================================================================
``PYP_SPEC``       str                   ``PypeIt`` spectrograph designation                                                                                                       
``ext_mode``       str                   Extraction mode (options: BOX, OPT)                                                                                                       
``flux``           ndarray   floating    Flux array in units of counts/s or 10^-17 erg/s/cm^2/Ang; see ``fluxed``                                                                  
``fluxed``         bool                  Boolean indicating if the spectrum is fluxed.                                                                                             
``ivar``           ndarray   floating    Inverse variance array (matches units of flux)                                                                                            
``mask``           ndarray   integer     Mask array (1=Good,0=Bad)                                                                                                                 
``obj_model``      ndarray   floating    Object model for tellurics                                                                                                                
``sigma``          ndarray   floating    One sigma noise array, equivalent to 1/sqrt(ivar) (matches units of flux)                                                                 
``spect_meta``     dict                  header dict                                                                                                                               
``telluric``       ndarray   floating    Telluric model                                                                                                                            
``wave``           ndarray   floating    Wavelength array (angstroms in vacuum), weighted by pixel contributions                                                                   
``wave_grid_mid``  ndarray   floating    Wavelength (angstroms in vacuum) evaluated at the bin centers of a grid that is uniformly spaced in either lambda or log10-lambda/velocity
=================  ========  ==========  ==========================================================================================================================================
