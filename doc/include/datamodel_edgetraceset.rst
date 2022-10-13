
Version 1.0.1

==============  ==========  ===============  ===================================================================================================================================================================
HDU Name        Obj Type    Array Type       Description                                                                                                                                                        
==============  ==========  ===============  ===================================================================================================================================================================
``PYP_SPEC``    str                          PypeIt spectrograph name                                                                                                                                           
``DISPNAME``    str                          Spectrograph disperser name.  Primarily needed for reloading an existing MasterEdge file.                                                                          
``EDGE_CEN``    ndarray     float, floating  (Floating-point) Measured spatial coordinate of the edge traces for each spectral pixel.  Shape is (Nspec,Ntrace).                                                 
``EDGE_ERR``    ndarray     float, floating  Error in the measured spatial coordinate edge traces.                                                                                                              
``EDGE_FIT``    ndarray     float, floating  The best-fit model result for the trace edge.                                                                                                                      
``EDGE_MSK``    ndarray     int32            Bitmask for the edge trace positions.                                                                                                                              
``FITTYPE``     str                          An informational string identifying the type of model used to fit the trace data.  Either ``pca`` for a PCA decomposition or the polynomial function type and order
``LEFT_PCA``    TracePCA                     The PCA decomposition of the left-edge traces.  Ignored if the left_right_pca parameter is False.                                                                  
``MASKDEF_ID``  ndarray     int, integer     slitmask ID number for the edge traces. IDs are for, respectively, left and right edges.                                                                           
``NSPAT``       int                          Image pixels in the spatial direction.                                                                                                                             
``NSPEC``       int                          Image pixels in the spectral direction.                                                                                                                            
``ORDERID``     ndarray     int, integer     For echelle spectrographs, this is the order ID number for the edge traces.  Negative and positive IDs are for, respectively, left and right edges.                
``PCA``         TracePCA                     The PCA decomposition of all edge traces.  Ignored if the left_right_pca parameter is True.                                                                        
``PCATYPE``     str                          String identifier for the measurements used to construct the PCA (center or fit)                                                                                   
``RIGHT_PCA``   TracePCA                     The PCA decomposition of the right-edge traces.  Ignored if the left_right_pca parameter is False.                                                                 
``SOBELSIG``    ndarray     float, floating  Sobel-filtered image used to detect edges                                                                                                                          
``TRACEBPM``    ndarray     bool, np.bool    Bad-pixel mask for trace image                                                                                                                                     
``TRACEID``     ndarray     int, integer     ID number for the edge traces.  Negative and positive IDs are for, respectively, left and right edges.                                                             
``TRACEIMG``    TraceImage                   Image used to construct the edge traces.                                                                                                                           
==============  ==========  ===============  ===================================================================================================================================================================
