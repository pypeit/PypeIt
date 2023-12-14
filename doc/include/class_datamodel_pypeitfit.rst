
**Version**: 1.0.0

===========  ================  =================  =====================================================================================
Attribute    Type              Array Type         Description                                                                          
===========  ================  =================  =====================================================================================
``fitc``     `numpy.ndarray`_  `numpy.floating`_  Fit coefficients                                                                     
``fitcov``   `numpy.ndarray`_  `numpy.floating`_  Covariance of the coefficients                                                       
``func``     str                                  Fit function (polynomial, legendre, chebyshev, polynomial2d, legendre2d)             
``gpm``      `numpy.ndarray`_  `numpy.integer`_   Mask (1=good)                                                                        
``maxx``     float                                maximum value in the array (or the right limit for a legendre / chebyshev polynomial)
``maxx2``    float                                Same as maxx for the second independent variable x2                                  
``minx``     float                                minimum value in the array (or the left limit for a legendre / chebyshev polynomial) 
``minx2``    float                                Same as minx for the second independent variable x2                                  
``order``    `numpy.ndarray`_  `numpy.integer`_   The order of the polynomial to be used in the fitting. This is a 2d array for 2d fits
``success``  int                                  Flag indicating whether fit was successful (success=1) or if it failed (success=0)   
``weights``  `numpy.ndarray`_  `numpy.floating`_  Weights.  Often the same as invvar                                                   
``x2``       `numpy.ndarray`_  `numpy.floating`_  x2 inputs, second independent variable                                               
``xval``     `numpy.ndarray`_  `numpy.floating`_  x inputs                                                                             
``yval``     `numpy.ndarray`_  `numpy.floating`_  y inputs                                                                             
===========  ================  =================  =====================================================================================
