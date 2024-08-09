
Version 1.1.5

=====================  ==============================  =========  =============================================================================================================
HDU Name               HDU Type                        Data Type  Description                                                                                                  
=====================  ==============================  =========  =============================================================================================================
``PRIMARY``            `astropy.io.fits.PrimaryHDU`_   ...        Empty data HDU.  Contains basic header information.                                                          
``SLITS``              `astropy.io.fits.BinTableHDU`_  ...        All data from the :class:`~pypeit.slittrace.SlitTraceSet` datamodel, except ``maskdef_designtab`` if present.
``MASKDEF_DESIGNTAB``  `astropy.io.fits.BinTableHDU`_  ...        Table with slitmask design and object info                                                                   
=====================  ==============================  =========  =============================================================================================================


MASKDEF_DESIGNTAB content

================  =========  ====================================================================================================
Column            Data Type  Description                                                                                         
================  =========  ====================================================================================================
``TRACEID``       int64      Trace ID Number                                                                                     
``TRACESROW``     int64      Spectral row for provided left and right edges.                                                     
``TRACELPIX``     float64    Spatial pixel coordinate for left edge                                                              
``TRACERPIX``     float64    Spatial pixel coordinate for right edge                                                             
``SPAT_ID``       int64      ID Number assigned by the pypeline to each slit                                                     
``MASKDEF_ID``    int64      Slit ID Number from slit-mask design                                                                
``SLITLMASKDEF``  float64    Left edge of the slit in pixel from slit-mask design before x-correlation                           
``SLITRMASKDEF``  float64    Right edge of the slit in pixel from slit-mask design before x-correlation                          
``SLITRA``        float64    Right ascension of the slit center (deg)                                                            
``SLITDEC``       float64    Declination of the slit center (deg)                                                                
``SLITLEN``       float64    Slit length (arcsec)                                                                                
``SLITWID``       float64    Slit width (arcsec)                                                                                 
``SLITPA``        float64    Slit position angle onsky (deg from N through E)                                                    
``ALIGN``         int16      Slit used for alignment (1-yes; 0-no), not target observations.                                     
``OBJID``         int64      Object ID Number                                                                                    
``OBJRA``         float64    Right ascension of the object (deg)                                                                 
``OBJDEC``        float64    Declination of the object (deg)                                                                     
``OBJNAME``       str        Object name assigned by the observer                                                                
``OBJMAG``        float64    Object magnitude provided by the observer                                                           
``OBJMAG_BAND``   str        Band of the magnitude provided by the observer                                                      
``OBJ_TOPDIST``   float64    Projected distance (in arcsec) of the object from the left edge of the slit (in PypeIt orientation).
``OBJ_BOTDIST``   float64    Projected distance (in arcsec) of the object from the right edge of the slit (in PypeIt orientation)
================  =========  ====================================================================================================
