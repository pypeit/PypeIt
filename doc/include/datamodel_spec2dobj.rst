

Version: 1.0.0

================  =================  ==========  ==========================================================
Obj               Obj Type           Array Type  Description                                               
================  =================  ==========  ==========================================================
``bpmmask``       ndarray            integer     2D bad-pixel mask for the image                           
``det``           int                            Detector index                                            
``detector``      DetectorContainer              Detector DataContainer                                    
``ivarmodel``     ndarray            floating    2D ivar model image                                       
``ivarraw``       ndarray            floating    2D processed inverse variance image                       
``objmodel``      ndarray            floating    2D object model image                                     
``sciimg``        ndarray            floating    2D processed science image                                
``skymodel``      ndarray            floating    2D sky model image                                        
``slits``         SlitTraceSet                   SlitTraceSet defining the slits                           
``spat_flexure``  float                          Shift, in spatial pixels, between this image and SlitTrace
``waveimg``       ndarray            floating    2D wavelength image                                       
================  =================  ==========  ==========================================================
