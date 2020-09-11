

Version: 1.0.2

====================  =================  ==========  ==================================================================================
Obj Key               Obj Type           Array Type  Description                                                                       
====================  =================  ==========  ==================================================================================
``bpmmask``           ndarray            integer     2D bad-pixel mask for the image                                                   
``det``               int                            Detector index                                                                    
``detector``          DetectorContainer              Detector DataContainer                                                            
``imgbitm``           str                            List of BITMASK keys from ImageBitMask                                            
``ivarmodel``         ndarray            floating    2D ivar model image (float32)                                                     
``ivarraw``           ndarray            floating    2D processed inverse variance image (float32)                                     
``objmodel``          ndarray            floating    2D object model image (float32)                                                   
``scaleimg``          ndarray            floating    2D multiplicative scale image that has been applied to the science image (float32)
``sci_spat_flexure``  float                          Shift, in spatial pixels, between this image and SlitTrace                        
``sciimg``            ndarray            floating    2D processed science image (float32)                                              
``skymodel``          ndarray            floating    2D sky model image (float32)                                                      
``slits``             SlitTraceSet                   SlitTraceSet defining the slits                                                   
``tilts``             ndarray            floating    2D tilts image (float64)                                                          
``waveimg``           ndarray            floating    2D wavelength image (float64)                                                     
====================  =================  ==========  ==================================================================================
