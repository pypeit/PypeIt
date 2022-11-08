
**Version**: 1.1.0

==============  ================  =================  =======================================
Attribute       Type              Array Type         Description                            
==============  ================  =================  =======================================
``PYP_SPEC``    str                                  PypeIt spectrograph name               
``alignframe``  `numpy.ndarray`_  `numpy.floating`_  Processed, combined alignment frames   
``nalign``      int                                  Number of alignment traces in each slit
``nslits``      int                                  The number of slits                    
``nspec``       int                                  The number of spectral elements        
``spat_id``     `numpy.ndarray`_  `numpy.integer`_   Slit spat_id                           
``traces``      `numpy.ndarray`_  `numpy.floating`_  Traces of the alignment frame          
==============  ================  =================  =======================================
