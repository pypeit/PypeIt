
**Version**: 1.0.0

===============  ================  =================  ===================================================================
Attribute        Type              Array Type         Description                                                        
===============  ================  =================  ===================================================================
``breakpoints``  `numpy.ndarray`_  `numpy.floating`_  Breakpoint locations                                               
``coeff``        `numpy.ndarray`_  `numpy.floating`_  Output fit coefficients                                            
``funcname``     str                                  Function type for the 2nd variable (when x2 is specified)          
``icoeff``       `numpy.ndarray`_  `numpy.floating`_  Cholesky band matrix used to solve for the bspline coefficients    
``mask``         `numpy.ndarray`_  `numpy.bool`_      Output mask                                                        
``nord``         int                                  Order of the bspline fit                                           
``npoly``        int                                  Order of polynomial to fit over 2nd variable (when x2 is specified)
``xmax``         float                                Normalization maximum for x2                                       
``xmin``         float                                Normalization minimum for x2                                       
===============  ================  =================  ===================================================================
