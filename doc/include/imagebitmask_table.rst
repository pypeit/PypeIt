==========  ==========  =============  ===========================================================
Bit Name    Bit Number  Decimal Value  Description                                                
==========  ==========  =============  ===========================================================
BPM         0           1              Component of the instrument-specific bad pixel mask        
CR          1           2              Cosmic ray detected                                        
SATURATION  2           4              Saturated pixel                                            
MINCOUNTS   3           8              Pixel below the instrument-specific minimum counts         
OFFSLITS    4           16             Pixel does not belong to any slit                          
IS_NAN      5           32             Pixel value is undefined                                   
IVAR0       6           64             Inverse variance is undefined                              
IVAR_NAN    7           128            Inverse variance is NaN                                    
EXTRACT     8           256            Pixel masked during local skysub and extraction            
BADSCALE    9           512            Bad image rescaling operation (e.g., flat-field value <= 0)
STCKMASK    10          1024           All pixels masked in image stack                           
USER        11          2048           Pixel masked by user                                       
==========  ==========  =============  ===========================================================
