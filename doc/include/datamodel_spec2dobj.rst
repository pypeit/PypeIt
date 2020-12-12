

Version: 1.0.3

====================  =================  ==========  ================================================================================================================================================================================
Obj Key               Obj Type           Array Type  Description                                                                                                                                                                     
====================  =================  ==========  ================================================================================================================================================================================
``bpmmask``           ndarray            integer     2D bad-pixel mask for the image                                                                                                                                                 
``det``               int                            Detector index                                                                                                                                                                  
``detector``          DetectorContainer              Detector DataContainer                                                                                                                                                          
``imgbitm``           str                            List of BITMASK keys from ImageBitMask                                                                                                                                          
``ivarmodel``         ndarray            floating    2D ivar model image (float32)                                                                                                                                                   
``ivarraw``           ndarray            floating    2D processed inverse variance image (float32)                                                                                                                                   
``objmodel``          ndarray            floating    2D object model image (float32)                                                                                                                                                 
``scaleimg``          ndarray            floating    2D multiplicative scale image that has been applied to the science image (float32)                                                                                              
``sci_spat_flexure``  float                          Shift, in spatial pixels, between this image and SlitTrace                                                                                                                      
``sci_spec_flexure``  Table                          Global shift of the spectrum to correct for spectralflexure (pixels). This is based on the sky spectrum atthe center of each slit                                               
``sciimg``            ndarray            floating    2D processed science image (float32)                                                                                                                                            
``skymodel``          ndarray            floating    2D sky model image (float32)                                                                                                                                                    
``slits``             SlitTraceSet                   SlitTraceSet defining the slits                                                                                                                                                 
``tilts``             ndarray            floating    2D tilts image (float64)                                                                                                                                                        
``vel_corr``          float                          Relativistic velocity correction for wavelengths                                                                                                                                
``vel_type``          str                            Type of reference frame correction (if any). Options are listed in the routine: WavelengthSolutionPar.valid_reference_frames() Current list: observed, heliocentric, barycentric
``waveimg``           ndarray            floating    2D wavelength image in vacuum (float64)                                                                                                                                         
====================  =================  ==========  ================================================================================================================================================================================
