

Version: 1.0.0

====================  =================  ==========  ==========================================================
Obj                   Obj Type           Array Type  Description                                               
====================  =================  ==========  ==========================================================
``bpmmask``           ndarray            integer     2D bad-pixel mask for the image                           
``det``               int                            Detector index                                            
``detector``          DetectorContainer              Detector DataContainer                                    
``imgbitm``           str                            List of BITMASK keys from ImageBitMask                    
``ivarmodel``         ndarray            floating    2D ivar model image                                       
``ivarraw``           ndarray            floating    2D processed inverse variance image                       
``objmodel``          ndarray            floating    2D object model image                                     
``sci_spat_flexure``  float                          Shift, in spatial pixels, between this image and SlitTrace
``sciimg``            ndarray            floating    2D processed science image                                
``skymodel``          ndarray            floating    2D sky model image                                        
``slits``             SlitTraceSet                   SlitTraceSet defining the slits                           
``tilts``             ndarray            floating    2D tilts image                                            
``waveimg``           ndarray            floating    2D wavelength image                                       
====================  =================  ==========  ==========================================================
